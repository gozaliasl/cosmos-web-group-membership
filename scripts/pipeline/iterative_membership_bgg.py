#!/usr/bin/env python
"""
Iterative membership + BGG selection with a flexible aperture (up to
--radius, default 750 kpc), to avoid missing a more massive galaxy near a
group whose catalog center is offset from the true dynamical/light center
(miscentering) or that simply falls just outside a smaller hard-coded radius.

Per group, outer loop:
  1. Run the existing probabilistic membership + redshift refinement
     (determine_membership_dztier.refine_group_redshift) centered on the
     current center (initially the catalog Ra/Dec), with the *full* search
     radius (default 750 kpc, not the previous 500 kpc default) -- so the
     mass-floor BGG step (select_bgg, no spatial restriction) already sees
     any massive galaxy within that wider net.
  2. Select the BGG from the resulting members (mass-floored hybrid
     mass+R-luminosity rank score, soft 300 kpc centrality preference).
  3. If the BGG position differs from the current center by more than
     --recenter-tol-kpc, recenter on the BGG and repeat (up to
     --max-recenter-iter times) -- this specifically catches a massive
     neighbor that would fall outside the aperture from the *original*
     (possibly miscentered) catalog position but inside it once centered on
     the true light/mass concentration.
  4. Stop when the BGG position stabilizes or the iteration cap is reached.

Mass/magnitude (mass_final/MR_final, Webb-native LP_mass_med_PDF/LP_MR_phys
with a COSMOS2020 fallback for true-margin members) is precomputed once for
the whole galaxy catalog up front, so each outer iteration is cheap.

Usage:
    python iterative_membership_bgg.py --catalog both --radius 750
"""

import argparse
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from determine_membership_dztier import (  # noqa: E402
    C_KMS, SIGMA_V_CEIL_KMS, SIGMA_V_FLOOR_KMS,
    angular_radius_arcsec, apply_membership_cut, find_candidates_in_aperture,
    load_galaxy_catalog, load_group_catalog,
)
from merge_photz_properties import PROPERTY_COLUMNS, load_photz_catalog  # noqa: E402
from select_bgg import select_bgg_for_group  # noqa: E402

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data"
SPECZ_FITS = DATA_DIR / "specz" / "Webb_Specz_Feb2026_plus_COSMOS_field.fits"
CW_ALL_FITS = DATA_DIR / "group-catalog" / "cosmos_web_groups_catalog_refined_z.fits"
CW_HCG_FITS = DATA_DIR / "group-catalog" / "Py18_Groups.fits"
COSMOS2020_FITS = DATA_DIR / "galaxy_catalog_cosmos2020" / "COSMOS2020_CLASSIC_R1_v2.2_p3.fits"
OUTPUT_DIR = BASE_DIR / "outputs" / "results" / "membership_dztier"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

M_SUN_R_AB = 4.42


def attach_mass_properties(galaxies: pd.DataFrame, match_radius_arcsec: float = 0.5) -> pd.DataFrame:
    """Precompute mass_final/MR_final/mass_source for the whole galaxy catalog once."""
    galaxies = galaxies.copy()
    galaxies["mass_final"] = np.nan
    galaxies["MR_final"] = np.nan
    galaxies["mass_source"] = "none"

    print("Attaching stellar mass / R-band magnitude (Webb-native LePhare)...")
    photz = load_photz_catalog()
    valid = photz.dropna(subset=["RA_MODEL", "DEC_MODEL"])
    photz_coords = SkyCoord(ra=valid["RA_MODEL"].values * u.deg, dec=valid["DEC_MODEL"].values * u.deg)
    gal_coords = SkyCoord(ra=galaxies["RA"].values * u.deg, dec=galaxies["DEC"].values * u.deg)
    idx, sep2d, _ = gal_coords.match_to_catalog_sky(photz_coords)
    matched = np.asarray(sep2d < match_radius_arcsec * u.arcsec)
    matched_photz = valid.iloc[idx].reset_index(drop=True)
    galaxies.loc[matched, "mass_final"] = matched_photz.loc[matched, "LP_mass_med_PDF"].values
    galaxies.loc[matched, "MR_final"] = matched_photz.loc[matched, "LP_MR_phys"].values
    galaxies.loc[matched, "mass_source"] = "Webb_LePhare"
    print(f"  Matched to Webb photz catalog: {matched.sum()}/{len(galaxies)}")

    needs_fill = galaxies["mass_final"].isna() & (galaxies["source_catalog"] == "COSMOS_full_field")
    n_needs_fill = int(needs_fill.sum())
    print(f"Filling remaining COSMOS_full_field galaxies from COSMOS2020: {n_needs_fill}")
    if n_needs_fill > 0:
        t = Table.read(COSMOS2020_FITS)
        c20 = t["ALPHA_J2000", "DELTA_J2000", "lp_mass_med", "lp_MR"].to_pandas()
        c20 = c20.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])
        c20_coords = SkyCoord(ra=c20["ALPHA_J2000"].values * u.deg, dec=c20["DELTA_J2000"].values * u.deg)
        fill_gal = galaxies[needs_fill]
        fill_coords = SkyCoord(ra=fill_gal["RA"].values * u.deg, dec=fill_gal["DEC"].values * u.deg)
        idx2, sep2d2, _ = fill_coords.match_to_catalog_sky(c20_coords)
        matched2 = np.asarray(sep2d2 < match_radius_arcsec * u.arcsec)
        matched_c20 = c20.iloc[idx2].reset_index(drop=True)

        fill_idx = fill_gal.index
        galaxies.loc[fill_idx[matched2], "mass_final"] = matched_c20.loc[matched2, "lp_mass_med"].values
        galaxies.loc[fill_idx[matched2], "MR_final"] = matched_c20.loc[matched2, "lp_MR"].values
        galaxies.loc[fill_idx[matched2], "mass_source"] = "COSMOS2020_LePhare"
        print(f"  Filled from COSMOS2020: {int(matched2.sum())}/{n_needs_fill}")

    print(f"mass_source breakdown:\n{galaxies['mass_source'].value_counts()}")
    return galaxies


def separation_kpc(ra1, dec1, ra2, dec2, z):
    sep_arcsec = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg).separation(
        SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)).arcsec
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
    return sep_arcsec * kpc_per_arcsec


def refine_redshift_at_center(center_ra, center_dec, z0, galaxies, radius_kpc,
                               prob_threshold, default_sigma_v_kms, max_iter=10, tol_kms=25.0):
    candidates = find_candidates_in_aperture(center_ra, center_dec, z0, galaxies, radius_kpc)
    if len(candidates) == 0:
        return candidates, z0, 0, False

    z_current = z0
    sigma_v_current = default_sigma_v_kms
    converged = False
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        scored = apply_membership_cut(candidates, z_current, sigma_v_current, prob_threshold)
        members = scored[scored["is_member"]]
        if len(members) < 3:
            return scored, z_current, n_iter, len(members) >= 1
        z_new = np.median(members["zfin"].values)
        dv_shift = C_KMS * (z_new - z_current) / (1 + z_current)
        member_v = C_KMS * (members["zfin"].values - z_new) / (1 + z_new)
        sigma_v_current = float(np.clip(np.std(member_v), SIGMA_V_FLOOR_KMS, SIGMA_V_CEIL_KMS))
        z_current = z_new
        if abs(dv_shift) < tol_kms:
            converged = True
            break
    scored = apply_membership_cut(candidates, z_current, sigma_v_current, prob_threshold)
    return scored, z_current, n_iter, converged


def iterative_group(group_id, ra0, dec0, z0, galaxies, radius_kpc, prob_threshold, default_sigma_v_kms,
                     mass_floor_n, centrality_kpc, centrality_tol,
                     max_recenter_iter, recenter_tol_kpc):
    center_ra, center_dec = ra0, dec0
    z_current = z0
    bgg_row = None
    n_recenter = 0
    center_converged = False

    for outer in range(1, max_recenter_iter + 1):
        n_recenter = outer
        scored, z_refined, n_z_iter, z_converged = refine_redshift_at_center(
            center_ra, center_dec, z_current, galaxies, radius_kpc, prob_threshold, default_sigma_v_kms,
        )
        if len(scored) == 0:
            return None, z_current, center_ra, center_dec, outer, False, False

        members = scored[scored["is_member"]].copy()
        if len(members) == 0:
            return None, z_refined, center_ra, center_dec, outer, z_converged, False

        members["Group_ID"] = group_id
        scored_bgg = select_bgg_for_group(members, mass_floor_n, centrality_kpc, centrality_tol)
        bgg_candidates = scored_bgg[scored_bgg["is_bgg"]]

        if len(bgg_candidates) == 0:
            # No mass/magnitude info at all -- keep current members, stop recentering
            return scored_bgg, z_refined, center_ra, center_dec, outer, z_converged, False

        bgg_row = bgg_candidates.iloc[0]
        offset_kpc = separation_kpc(center_ra, center_dec, bgg_row["RA"], bgg_row["DEC"], z_refined)

        if offset_kpc < recenter_tol_kpc:
            center_converged = True
            return scored_bgg, z_refined, bgg_row["RA"], bgg_row["DEC"], outer, z_converged, center_converged

        center_ra, center_dec = bgg_row["RA"], bgg_row["DEC"]
        z_current = z_refined

    return scored_bgg, z_refined, center_ra, center_dec, n_recenter, z_converged, center_converged


def process_catalog(groups, galaxies, catalog_name, id_col, ra_col, dec_col, z_col,
                     radius_kpc, prob_threshold, default_sigma_v_kms,
                     mass_floor_n, centrality_kpc, centrality_tol,
                     max_recenter_iter, recenter_tol_kpc, sample_size=None,
                     radius_col=None, radius_factor=3.0, radius_min=100.0, radius_max=250.0):
    """
    radius_col: if given (e.g. 'Rad' for CW-HCG, arcsec), the per-group search
    radius is clip(radius_factor * Rad_kpc, radius_min, radius_max) instead of
    the fixed `radius_kpc`. Used to keep compact-group apertures physically
    sized to the group itself, rather than the general-group 750 kpc default
    (which would be 5-20x the actual size of a Hickson-type compact group and
    risks pulling in / recentering onto unrelated massive interlopers).
    """
    if sample_size is not None:
        groups = groups.head(sample_size)

    members_list = []
    summary_list = []

    for _, row in tqdm(groups.iterrows(), total=len(groups), desc=f"Processing {catalog_name}"):
        group_id = row[id_col]
        group_ra, group_dec, group_z0 = row[ra_col], row[dec_col], row[z_col]
        if group_z0 > 4.0 or pd.isna(group_z0):
            continue

        if radius_col is not None and radius_col in row.index and pd.notna(row[radius_col]):
            kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z0).value / 60.0
            rad_kpc = row[radius_col] * kpc_per_arcsec
            this_radius_kpc = float(np.clip(radius_factor * rad_kpc, radius_min, radius_max))
        else:
            this_radius_kpc = radius_kpc

        scored_bgg, z_refined, final_ra, final_dec, n_recenter, z_converged, center_converged = iterative_group(
            group_id, group_ra, group_dec, group_z0, galaxies, this_radius_kpc, prob_threshold, default_sigma_v_kms,
            mass_floor_n, centrality_kpc, centrality_tol, max_recenter_iter, recenter_tol_kpc,
        )

        if scored_bgg is None:
            members = pd.DataFrame()
        else:
            members = scored_bgg[scored_bgg["is_member"]].copy()

        if len(members) > 0:
            members["Group_Ra_catalog"] = group_ra
            members["Group_Dec_catalog"] = group_dec
            members["Group_Ra_final"] = final_ra
            members["Group_Dec_final"] = final_dec
            members["Group_z_catalog"] = group_z0
            members["Group_z_refined"] = z_refined
            members["catalog"] = catalog_name
            members_list.append(members)

        n_specz = int((members["redshift_type"] == "spec-z").sum()) if len(members) else 0
        n_photz = int((members["redshift_type"] == "photo-z").sum()) if len(members) else 0
        recenter_offset_kpc = (separation_kpc(group_ra, group_dec, final_ra, final_dec, z_refined)
                                if len(members) else np.nan)

        summary_list.append({
            "Group_ID": group_id,
            "search_radius_kpc": this_radius_kpc,
            "Group_Ra_catalog": group_ra,
            "Group_Dec_catalog": group_dec,
            "Group_Ra_final": final_ra,
            "Group_Dec_final": final_dec,
            "recenter_offset_kpc": recenter_offset_kpc,
            "Group_z_catalog": group_z0,
            "z_refined": z_refined,
            "n_recenter_iter": n_recenter,
            "z_converged": z_converged,
            "center_converged": center_converged,
            "n_members_total": len(members),
            "n_members_specz": n_specz,
            "n_members_photoz": n_photz,
            "catalog": catalog_name,
        })

    all_members = pd.concat(members_list, ignore_index=True) if members_list else pd.DataFrame()
    summary = pd.DataFrame(summary_list)
    return all_members, summary


def save(df, path_stub):
    if len(df) == 0:
        print(f"  (empty, skipping) {path_stub}")
        return
    csv_path = path_stub.with_suffix(".csv")
    fits_path = path_stub.with_suffix(".fits")
    df.to_csv(csv_path, index=False)
    Table.from_pandas(df).write(fits_path, format="fits", overwrite=True)
    print(f"  Saved: {csv_path.name}, {fits_path.name} (N={len(df)})")


def main():
    parser = argparse.ArgumentParser(description="Iterative membership + BGG selection with flexible (up to 750 kpc) aperture and BGG recentering")
    parser.add_argument("--catalog", choices=["all", "hcg", "both"], default="both")
    parser.add_argument("--radius", type=float, default=750.0, help="Search radius in kpc (flexible aperture)")
    parser.add_argument("--prob-threshold", type=float, default=0.05)
    parser.add_argument("--default-sigma-v-kms", type=float, default=500.0)
    parser.add_argument("--mass-floor-n", type=int, default=2)
    parser.add_argument("--centrality-kpc", type=float, default=300.0)
    parser.add_argument("--centrality-tol", type=float, default=0.05)
    parser.add_argument("--max-recenter-iter", type=int, default=3)
    parser.add_argument("--recenter-tol-kpc", type=float, default=20.0)
    # CW-HCG (Hickson-type compact groups) are physically tens-to-~100 kpc across --
    # NOT general/extended groups. Using the CW-All 750 kpc flexible aperture for them
    # would be 5-20x their actual size and risks pulling in / recentering onto unrelated
    # massive interlopers. Instead, scale the HCG search radius off the group's own
    # Rad column (arcsec -> kpc), clamped to a compact-group-appropriate range.
    parser.add_argument("--hcg-radius-col", type=str, default="Rad",
                        help="Py18_Groups.fits column (arcsec) defining each compact group's own scale")
    parser.add_argument("--hcg-radius-factor", type=float, default=3.0,
                        help="HCG search radius = clip(this * Rad_kpc, hcg-radius-min, hcg-radius-max)")
    parser.add_argument("--hcg-radius-min-kpc", type=float, default=100.0)
    parser.add_argument("--hcg-radius-max-kpc", type=float, default=250.0)
    parser.add_argument("--hcg-centrality-kpc", type=float, default=100.0,
                        help="Tighter centrality preference radius for compact groups (vs 300 kpc for CW-All)")
    parser.add_argument("--hcg-recenter-tol-kpc", type=float, default=10.0,
                        help="Tighter recenter convergence tolerance for compact groups (vs 20 kpc for CW-All)")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--group-catalog-all", type=str, default=None)
    parser.add_argument("--group-catalog-hcg", type=str, default=None)
    parser.add_argument("--xray-detected-only-all", type=str, default=None)
    parser.add_argument("--xray-detected-only-hcg", type=str, default=None)
    parser.add_argument("--output-prefix", type=str, default="iterative_")
    args = parser.parse_args()

    print("Loading galaxy catalog (spec-z + photo-z tiers)...")
    galaxies = load_galaxy_catalog(SPECZ_FITS)
    galaxies = attach_mass_properties(galaxies)

    def restrict_to_detected(groups, xray_path, id_col):
        p = Path(xray_path)
        xray_df = Table.read(p).to_pandas() if p.suffix.lower() == ".fits" else pd.read_csv(p)
        detected_ids = set(xray_df["Group_ID"].astype(str))
        restricted = groups[groups[id_col].astype(str).isin(detected_ids)].copy()
        print(f"  Restricted to X-ray-detected groups: {len(restricted)}/{len(groups)}")
        return restricted

    sample_size = 10 if args.test else None
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    common_kwargs = dict(
        radius_kpc=args.radius, prob_threshold=args.prob_threshold, default_sigma_v_kms=args.default_sigma_v_kms,
        mass_floor_n=args.mass_floor_n, centrality_kpc=args.centrality_kpc, centrality_tol=args.centrality_tol,
        max_recenter_iter=args.max_recenter_iter, recenter_tol_kpc=args.recenter_tol_kpc, sample_size=sample_size,
    )

    if args.catalog in ["all", "both"]:
        cw_all = load_group_catalog(Path(args.group_catalog_all) if args.group_catalog_all else CW_ALL_FITS)
        if args.xray_detected_only_all:
            cw_all = restrict_to_detected(cw_all, args.xray_detected_only_all, "Group_ID")
        members_all, summary_all = process_catalog(
            cw_all, galaxies, "CW-All", id_col="Group_ID", ra_col="Ra", dec_col="Dec", z_col="z", **common_kwargs,
        )
        print(f"  CW-All: {len(summary_all)} groups, {len(members_all)} members, "
              f"recentered>0: {(summary_all['n_recenter_iter']>1).sum()}, "
              f"center_converged: {summary_all['center_converged'].sum()}")
        save(members_all, OUTPUT_DIR / f"{args.output_prefix}cw_all_members_{timestamp}")
        save(summary_all, OUTPUT_DIR / f"{args.output_prefix}cw_all_summary_{timestamp}")

    if args.catalog in ["hcg", "both"]:
        cw_hcg = load_group_catalog(Path(args.group_catalog_hcg) if args.group_catalog_hcg else CW_HCG_FITS)
        if args.xray_detected_only_hcg:
            cw_hcg = restrict_to_detected(cw_hcg, args.xray_detected_only_hcg, "Group_ID")
        hcg_kwargs = dict(common_kwargs)
        hcg_kwargs.update(
            centrality_kpc=args.hcg_centrality_kpc,
            recenter_tol_kpc=args.hcg_recenter_tol_kpc,
            radius_col=args.hcg_radius_col,
            radius_factor=args.hcg_radius_factor,
            radius_min=args.hcg_radius_min_kpc,
            radius_max=args.hcg_radius_max_kpc,
        )
        print(f"  CW-HCG uses compact-group-scaled apertures (from '{args.hcg_radius_col}'), "
              f"clipped to [{args.hcg_radius_min_kpc}, {args.hcg_radius_max_kpc}] kpc -- "
              f"NOT the CW-All {args.radius} kpc flexible aperture")
        members_hcg, summary_hcg = process_catalog(
            cw_hcg, galaxies, "CW-HCG", id_col="Group_ID", ra_col="Ra", dec_col="Dec", z_col="z", **hcg_kwargs,
        )
        print(f"  CW-HCG: {len(summary_hcg)} groups, {len(members_hcg)} members, "
              f"recentered>0: {(summary_hcg['n_recenter_iter']>1).sum()}, "
              f"center_converged: {summary_hcg['center_converged'].sum()}")
        save(members_hcg, OUTPUT_DIR / f"{args.output_prefix}cw_hcg_members_{timestamp}")
        save(summary_hcg, OUTPUT_DIR / f"{args.output_prefix}cw_hcg_summary_{timestamp}")

    print("\nDone.")


if __name__ == "__main__":
    main()
