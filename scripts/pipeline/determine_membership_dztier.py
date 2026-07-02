#!/usr/bin/env python
"""
Group membership + redshift refinement using the merged Webb_Specz_Feb2026 catalog.

Data model (data/specz/Webb_Specz_Feb2026.fits, columns RA, DEC, zfin, dz):
Each galaxy has one final redshift (zfin) and one redshift error (dz), which
already encodes the quality tier assigned by the spec-z release:
    dz = 0.001  -> spec-z, quality 3 & 4 (best)
    dz = 0.002  -> spec-z, quality 2
    dz = 0.003  -> spec-z, quality 1, consistent with photo-z
    dz > 0.003  -> photometric redshift (dz = photo-z uncertainty)
There is no separate photo-z catalog to join against: rows with dz > 0.003
*are* the photo-z members, living in the same table as the spec-z rows.

Membership criterion (probabilistic, replaces an earlier hard n_sigma + fixed
1000 km/s cut): a fixed velocity window is a reasonable dynamical cut for real
spec-z (group sigma_v here runs ~300-800 km/s), but is badly mismatched to
photo-z, whose per-galaxy uncertainty (median 1-sigma ~3000+ km/s) dwarfs any
group's actual dynamics -- a flat 1000 km/s cut let only ~23% of true photo-z
members through by chance, and a flat looser cut would swamp low-z groups
with interlopers. Instead each candidate gets a continuous probability:

    P_z = exp(-0.5 * (dz / dz_gal)^2)          # uses the galaxy's OWN redshift error
    P_v = exp(-0.5 * (dv / sigma_v_group)^2)   # dynamical likelihood using the group's
                                                # own measured sigma_v (or a default
                                                # when not yet measured)
    membership_prob = P_z * P_v
    is_member = membership_prob >= prob_threshold

This degrades gracefully: a precise spec-z galaxy at the group velocity gets
prob~1; an imprecise photo-z galaxy is downweighted smoothly by its own
uncertainty rather than being coin-flipped by an arbitrary shared threshold.

Redshift refinement:
    Iteratively recentre the group redshift (and sigma_v_group) on the
    current probabilistic members and re-score, until the shift in velocity
    is below `tol_kms` or `max_iter` is reached (interloper rejection akin to
    the VRF/gapper approach used elsewhere in this repo, but keyed on each
    galaxy's own dz-tier error and a probabilistic membership score rather
    than a fixed cut).

Usage:
    python determine_membership_dztier.py --catalog both --radius 500
    python determine_membership_dztier.py --catalog all --test
"""

import argparse
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

warnings.filterwarnings("ignore")

C_KMS = 299792.458

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / "data"
SPECZ_FITS = DATA_DIR / "specz" / "Webb_Specz_Feb2026.fits"
CW_ALL_FITS = DATA_DIR / "group-catalog" / "cosmos_web_groups_catalog_refined_z.fits"
CW_HCG_FITS = DATA_DIR / "group-catalog" / "Py18_Groups.fits"
OUTPUT_DIR = BASE_DIR / "outputs" / "results" / "membership_dztier"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SPECZ_QUALITY_MAX_DZ = 0.003  # dz <= this -> spec-z; dz > this -> photo-z


def load_galaxy_catalog(path: Path = SPECZ_FITS) -> pd.DataFrame:
    df = Table.read(path).to_pandas()
    for col in df.select_dtypes(include=["object"]).columns:
        if df[col].map(lambda v: isinstance(v, bytes)).any():
            df[col] = df[col].str.decode("utf-8")
    df["redshift_type"] = np.where(df["dz"] <= SPECZ_QUALITY_MAX_DZ, "spec-z", "photo-z")
    return df


def load_group_catalog(path: Path) -> pd.DataFrame:
    return Table.read(path).to_pandas()


def angular_radius_arcsec(z: float, radius_kpc: float) -> float:
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
    return radius_kpc / kpc_per_arcsec


def find_candidates_in_aperture(group_ra, group_dec, group_z, galaxies, radius_kpc):
    """Spatial cut only (redshift cut applied afterwards per iteration)."""
    radius_arcsec = angular_radius_arcsec(group_z, radius_kpc)
    group_coord = SkyCoord(ra=group_ra * u.deg, dec=group_dec * u.deg)
    gal_coords = SkyCoord(ra=galaxies["RA"].values * u.deg, dec=galaxies["DEC"].values * u.deg)
    sep = group_coord.separation(gal_coords)
    within = sep < radius_arcsec * u.arcsec
    if within.sum() == 0:
        return pd.DataFrame()
    cand = galaxies[within].copy()
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    cand["sep_arcsec"] = sep[within].arcsec
    cand["sep_kpc"] = cand["sep_arcsec"] * kpc_per_arcsec
    return cand


SIGMA_V_FLOOR_KMS = 150.0   # avoid runaway collapse when few members coincide closely
SIGMA_V_CEIL_KMS = 2500.0   # avoid runaway growth from early-iteration outliers


def apply_membership_cut(candidates, group_z, sigma_v_kms, prob_threshold):
    dz = candidates["zfin"].values - group_z
    dv = C_KMS * dz / (1 + group_z)
    dz_norm = dz / (1 + group_z)
    dz_gal = candidates["dz"].values

    p_z = np.exp(-0.5 * (dz / dz_gal) ** 2)
    p_v = np.exp(-0.5 * (dv / sigma_v_kms) ** 2)
    membership_prob = p_z * p_v

    out = candidates.copy()
    out["dv"] = dv
    out["dz_norm"] = dz_norm
    out["n_sigma_offset"] = np.abs(dz) / dz_gal
    out["P_z"] = p_z
    out["P_v"] = p_v
    out["membership_prob"] = membership_prob
    out["is_member"] = membership_prob >= prob_threshold
    return out


def refine_group_redshift(group_ra, group_dec, group_z0, galaxies, radius_kpc,
                           prob_threshold=0.05, default_sigma_v_kms=500.0,
                           max_iter=10, tol_kms=25.0):
    """
    Iteratively refine the group center redshift (and its velocity dispersion)
    using probabilistic membership.

    Returns
    -------
    members : DataFrame
        Final candidate table (all galaxies in aperture, `is_member` flag,
        `membership_prob`, dv/dz_norm computed relative to the refined redshift).
    z_refined : float
    n_iter : int
    converged : bool
    """
    candidates = find_candidates_in_aperture(group_ra, group_dec, group_z0, galaxies, radius_kpc)
    if len(candidates) == 0:
        return candidates, group_z0, 0, False

    z_current = group_z0
    sigma_v_current = default_sigma_v_kms
    converged = False
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        scored = apply_membership_cut(candidates, z_current, sigma_v_current, prob_threshold)
        members = scored[scored["is_member"]]
        if len(members) < 3:
            # Not enough members to refine; keep last valid scoring
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


def process_catalog(groups, galaxies, catalog_name, id_col, ra_col, dec_col, z_col,
                     radius_kpc, prob_threshold, default_sigma_v_kms, sample_size=None):
    if sample_size is not None:
        groups = groups.head(sample_size)

    members_list = []
    summary_list = []

    for _, row in tqdm(groups.iterrows(), total=len(groups), desc=f"Processing {catalog_name}"):
        group_id = row[id_col]
        group_ra, group_dec, group_z0 = row[ra_col], row[dec_col], row[z_col]

        if group_z0 > 4.0 or pd.isna(group_z0):
            continue

        scored, z_refined, n_iter, converged = refine_group_redshift(
            group_ra, group_dec, group_z0, galaxies, radius_kpc,
            prob_threshold=prob_threshold, default_sigma_v_kms=default_sigma_v_kms,
        )

        if len(scored) > 0:
            members = scored[scored["is_member"]].copy()
        else:
            members = pd.DataFrame()

        if len(members) > 0:
            members["Group_ID"] = group_id
            members["Group_Ra"] = group_ra
            members["Group_Dec"] = group_dec
            members["Group_z_catalog"] = group_z0
            members["Group_z_refined"] = z_refined
            members["catalog"] = catalog_name
            members_list.append(members)

        n_specz = int((members["redshift_type"] == "spec-z").sum()) if len(members) else 0
        n_photz = int((members["redshift_type"] == "photo-z").sum()) if len(members) else 0

        sigma_v = np.nan
        if n_specz + n_photz >= 3:
            member_z = members["zfin"].values
            member_v = C_KMS * (member_z - z_refined) / (1 + z_refined)
            sigma_v = np.std(member_v)

        summary_list.append({
            "Group_ID": group_id,
            "Group_Ra": group_ra,
            "Group_Dec": group_dec,
            "Group_z_catalog": group_z0,
            "z_refined": z_refined,
            "dv_refinement_kms": C_KMS * (z_refined - group_z0) / (1 + group_z0),
            "n_refine_iter": n_iter,
            "refinement_converged": converged,
            "n_members_total": len(members),
            "n_members_specz": n_specz,
            "n_members_photoz": n_photz,
            "sigma_v_kms": sigma_v,
            "catalog": catalog_name,
        })

    all_members = pd.concat(members_list, ignore_index=True) if members_list else pd.DataFrame()
    summary = pd.DataFrame(summary_list)
    return all_members, summary


# X-ray columns worth carrying into the membership summary catalog
XRAY_MERGE_COLUMNS = [
    "Is_Detected", "SNR", "Significance_Sigma",
    "Flux_erg_cm2_s", "Flux_Error",
    "Luminosity_erg_s", "Luminosity_Error",
    "Temperature_keV", "Temperature_Error",
    "R200_kpc", "R500_kpc",
    "Log10_M200_Temp", "Log10_M200_Temp_Error",
    "Log10_M200_Luminosity", "Log10_M200_Luminosity_Error",
    "Is_Projected_Contaminated", "Contamination_Severity",
    "Is_Suspected_False_Positive",
]


def merge_xray_properties(summary: pd.DataFrame, xray_path: str) -> pd.DataFrame:
    """Left-merge X-ray properties (flux/L_X/T/mass/flags) onto the membership summary by Group_ID."""
    p = Path(xray_path)
    xray_df = Table.read(p).to_pandas() if p.suffix.lower() == ".fits" else pd.read_csv(p)
    cols = ["Group_ID"] + [c for c in XRAY_MERGE_COLUMNS if c in xray_df.columns]
    xray_df = xray_df[cols].copy()
    xray_df["Group_ID"] = xray_df["Group_ID"].astype(str)
    summary = summary.copy()
    summary["Group_ID"] = summary["Group_ID"].astype(str)
    merged = summary.merge(xray_df, on="Group_ID", how="left")
    print(f"  Merged X-ray properties from {p.name}: "
          f"{merged['Is_Detected'].notna().sum() if 'Is_Detected' in merged.columns else 0}/{len(merged)} groups matched")
    return merged


def add_quality_flag(summary: pd.DataFrame) -> pd.DataFrame:
    """
    Add an overall `quality_flag` (1=good, 2=marginal, 3=poor) and `quality_reason`.

    Membership-only criteria:
      - poor (3):    n_members_total < 3 (no redshift refinement ran) or refinement did not converge
      - marginal (2): 3 <= n_members_total < 8
      - good (1):     n_members_total >= 8 and refinement_converged

    If X-ray columns are present (from merge_xray_properties), a group is downgraded
    (capped at marginal) when Is_Suspected_False_Positive or Is_Projected_Contaminated
    is True, or Significance_Sigma < 3.
    """
    summary = summary.copy()
    flag = np.full(len(summary), 1, dtype=int)
    reason = np.full(len(summary), "", dtype=object)

    poor_mask = (summary["n_members_total"] < 3) | (~summary["refinement_converged"].astype(bool))
    flag[poor_mask.values] = 3
    reason[poor_mask.values] = "too_few_members_or_not_converged"

    marginal_mask = (~poor_mask) & (summary["n_members_total"] < 8)
    flag[marginal_mask.values] = 2
    reason[marginal_mask.values] = "few_members"

    good_mask = (~poor_mask) & (~marginal_mask)
    reason[good_mask.values] = "ok"

    if "Is_Suspected_False_Positive" in summary.columns:
        downgrade = summary["Is_Suspected_False_Positive"].fillna(False).astype(bool).values
        flag[downgrade & (flag == 1)] = 2
        reason[downgrade] = np.where(reason[downgrade] == "", "suspected_false_positive",
                                      reason[downgrade] + "+suspected_false_positive")

    if "Is_Projected_Contaminated" in summary.columns:
        downgrade = summary["Is_Projected_Contaminated"].fillna(False).astype(bool).values
        flag[downgrade & (flag == 1)] = 2
        reason[downgrade] = np.where(reason[downgrade] == "", "projected_contamination",
                                      reason[downgrade] + "+projected_contamination")

    if "Significance_Sigma" in summary.columns:
        low_sig = (summary["Significance_Sigma"] < 3).fillna(False).values
        flag[low_sig & (flag == 1)] = 2
        reason[low_sig] = np.where(reason[low_sig] == "", "low_xray_significance",
                                    reason[low_sig] + "+low_xray_significance")

    summary["quality_flag"] = flag
    summary["quality_reason"] = reason
    return summary


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
    parser = argparse.ArgumentParser(
        description="Membership determination + redshift refinement on Webb_Specz_Feb2026 (spec-z + photo-z tiers combined)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--catalog", choices=["all", "hcg", "both"], default="both",
                        help="all=CW-All, hcg=CW-HCG (a.k.a. CW-CHG), both=run both")
    parser.add_argument("--radius", type=float, default=500.0, help="Search radius in kpc")
    parser.add_argument("--prob-threshold", type=float, default=0.05,
                        help="Minimum membership_prob = P_z * P_v to flag is_member")
    parser.add_argument("--default-sigma-v-kms", type=float, default=500.0,
                        help="Initial/fallback group velocity dispersion (km/s) for P_v before enough members are found")
    parser.add_argument("--test", action="store_true", help="Process only 10 groups per catalog")
    parser.add_argument("--group-catalog-all", type=str, default=None,
                        help="Override CW-All group catalog FITS path")
    parser.add_argument("--group-catalog-hcg", type=str, default=None,
                        help="Override CW-HCG group catalog FITS path")
    parser.add_argument("--xray-detected-only-all", type=str, default=None,
                        help="Path to an X-ray catalog FITS/CSV for CW-All (e.g. xray_detected_groups.fits); "
                             "restricts processing to the Group_IDs present in it")
    parser.add_argument("--xray-detected-only-hcg", type=str, default=None,
                        help="Path to an X-ray catalog FITS/CSV for CW-HCG; restricts to its Group_IDs")
    parser.add_argument("--output-prefix", type=str, default="",
                        help="Optional prefix for output filenames (e.g. 'xray_detected_')")
    parser.add_argument("--merge-xray-all", type=str, default=None,
                        help="Path to xray_catalog/xray_detected_groups FITS/CSV for CW-All; "
                             "merges X-ray properties (flux, L_X, T, mass, flags) into the summary")
    parser.add_argument("--merge-xray-hcg", type=str, default=None,
                        help="Path to xray_catalog/xray_detected_groups FITS/CSV for CW-HCG; "
                             "merges X-ray properties into the summary")
    parser.add_argument("--galaxy-catalog", type=str, default=None,
                        help="Override the galaxy catalog (default: data/specz/Webb_Specz_Feb2026.fits). "
                             "Pass Webb_Specz_Feb2026_plus_COSMOS_field.fits to include full-COSMOS-field "
                             "margin coverage (has a `source_catalog` column: 'Webb' or 'COSMOS_full_field').")
    args = parser.parse_args()

    print("Loading galaxy catalog (spec-z + photo-z tiers)...")
    galaxies = load_galaxy_catalog(Path(args.galaxy_catalog) if args.galaxy_catalog else SPECZ_FITS)
    print(f"  Total galaxies: {len(galaxies)} "
          f"(spec-z: {(galaxies['redshift_type']=='spec-z').sum()}, "
          f"photo-z: {(galaxies['redshift_type']=='photo-z').sum()})")
    if "source_catalog" in galaxies.columns:
        print(f"  By source: {galaxies['source_catalog'].value_counts().to_dict()}")

    sample_size = 10 if args.test else None
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    def restrict_to_detected(groups: pd.DataFrame, xray_path: str, id_col: str) -> pd.DataFrame:
        p = Path(xray_path)
        xray_df = Table.read(p).to_pandas() if p.suffix.lower() == ".fits" else pd.read_csv(p)
        detected_ids = set(xray_df["Group_ID"].astype(str))
        restricted = groups[groups[id_col].astype(str).isin(detected_ids)].copy()
        print(f"  Restricted to X-ray-detected groups: {len(restricted)}/{len(groups)} "
              f"(from {p.name})")
        return restricted

    if args.catalog in ["all", "both"]:
        print("\nLoading CW-All group catalog...")
        cw_all = load_group_catalog(Path(args.group_catalog_all) if args.group_catalog_all else CW_ALL_FITS)
        print(f"  Groups: {len(cw_all)}")
        if args.xray_detected_only_all:
            cw_all = restrict_to_detected(cw_all, args.xray_detected_only_all, "Group_ID")
        members_all, summary_all = process_catalog(
            cw_all, galaxies, "CW-All",
            id_col="Group_ID", ra_col="Ra", dec_col="Dec", z_col="z",
            radius_kpc=args.radius, prob_threshold=args.prob_threshold, default_sigma_v_kms=args.default_sigma_v_kms,
            sample_size=sample_size,
        )
        print(f"  CW-All: {len(summary_all)} groups processed, {len(members_all)} total members")
        if args.merge_xray_all:
            summary_all = merge_xray_properties(summary_all, args.merge_xray_all)
        summary_all = add_quality_flag(summary_all)
        save(members_all, OUTPUT_DIR / f"{args.output_prefix}cw_all_members_{timestamp}")
        save(summary_all, OUTPUT_DIR / f"{args.output_prefix}cw_all_summary_{timestamp}")

    if args.catalog in ["hcg", "both"]:
        print("\nLoading CW-HCG (CW-CHG) group catalog...")
        cw_hcg = load_group_catalog(Path(args.group_catalog_hcg) if args.group_catalog_hcg else CW_HCG_FITS)
        print(f"  Groups: {len(cw_hcg)}")
        if args.xray_detected_only_hcg:
            cw_hcg = restrict_to_detected(cw_hcg, args.xray_detected_only_hcg, "Group_ID")
        members_hcg, summary_hcg = process_catalog(
            cw_hcg, galaxies, "CW-HCG",
            id_col="Group_ID", ra_col="Ra", dec_col="Dec", z_col="z",
            radius_kpc=args.radius, prob_threshold=args.prob_threshold, default_sigma_v_kms=args.default_sigma_v_kms,
            sample_size=sample_size,
        )
        print(f"  CW-HCG: {len(summary_hcg)} groups processed, {len(members_hcg)} total members")
        if args.merge_xray_hcg:
            summary_hcg = merge_xray_properties(summary_hcg, args.merge_xray_hcg)
        summary_hcg = add_quality_flag(summary_hcg)
        save(members_hcg, OUTPUT_DIR / f"{args.output_prefix}cw_hcg_members_{timestamp}")
        save(summary_hcg, OUTPUT_DIR / f"{args.output_prefix}cw_hcg_summary_{timestamp}")

    print("\nDone.")


if __name__ == "__main__":
    main()
