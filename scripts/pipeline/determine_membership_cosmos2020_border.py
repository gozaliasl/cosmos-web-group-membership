#!/usr/bin/env python
"""
Supplement group membership near the COSMOS-Web footprint edge using COSMOS2020.

The Webb_Specz_Feb2026 catalog only covers the COSMOS-Web footprint, so groups
whose 500 kpc search aperture crosses the field edge can lose true members that
fall just outside it. This script:

  1. Builds a footprint proxy as the RA/Dec convex hull of the CW-All + CW-HCG
     group positions.
  2. Flags "border groups": those whose search aperture (radius_kpc, converted
     to arcsec at the group's redshift) comes within `--margin-arcmin` of the
     hull boundary.
  3. For border groups only, additionally searches COSMOS2020
     (data/galaxy_catalog_cosmos2020/COSMOS2020_CLASSIC_R1_v2.2_p3.fits) for
     candidate members in the same aperture, using COSMOS2020's own photo-z
     (lp_zBEST) and per-galaxy error (half the lp_zPDF 68% interval), applying
     the same 3-sigma + |dv|<max-velocity criterion as the Webb-based pipeline.
  4. Merges these into the existing member catalog, tagged
     `member_source` = 'Webb' or 'COSMOS2020', deduplicating any COSMOS2020
     candidate that positionally matches an existing Webb member (kept as
     'Webb').

Usage:
    python determine_membership_cosmos2020_border.py \
        --summary-file outputs/results/membership_dztier/xray_detected_cw_all_summary_*.csv \
        --members-file outputs/results/membership_dztier/xray_detected_cw_all_members_*_with_photz.csv \
        --catalog-name CW-All --margin-arcmin 1.0
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from tqdm import tqdm

warnings.filterwarnings("ignore")

C_KMS = 299792.458
BASE_DIR = Path(__file__).resolve().parent.parent.parent
COSMOS2020_FITS = BASE_DIR / "data" / "galaxy_catalog_cosmos2020" / "COSMOS2020_CLASSIC_R1_v2.2_p3.fits"


def build_footprint_hull(all_group_positions: np.ndarray) -> ConvexHull:
    """RA/Dec convex hull of group positions, as a proxy for the COSMOS-Web footprint."""
    return ConvexHull(all_group_positions)


def point_to_hull_boundary_distance_deg(point: np.ndarray, hull: ConvexHull) -> float:
    """Minimum distance (deg) from a point to the nearest hull edge segment."""
    pts = hull.points
    min_dist = np.inf
    for simplex in hull.simplices:
        a, b = pts[simplex[0]], pts[simplex[1]]
        ab = b - a
        ab_len2 = np.dot(ab, ab)
        if ab_len2 == 0:
            d = np.linalg.norm(point - a)
        else:
            t = np.clip(np.dot(point - a, ab) / ab_len2, 0.0, 1.0)
            proj = a + t * ab
            d = np.linalg.norm(point - proj)
        min_dist = min(min_dist, d)
    return min_dist


def angular_radius_arcsec(z: float, radius_kpc: float) -> float:
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
    return radius_kpc / kpc_per_arcsec


def flag_border_groups(summary: pd.DataFrame, hull: ConvexHull, radius_kpc: float, margin_arcmin: float) -> pd.DataFrame:
    summary = summary.copy()
    margin_arcsec = margin_arcmin * 60.0
    is_border = np.zeros(len(summary), dtype=bool)
    dist_to_edge_arcsec = np.full(len(summary), np.nan)

    for i, row in enumerate(summary.itertuples()):
        z = row.z_refined if np.isfinite(row.z_refined) else row.Group_z_catalog
        aperture_arcsec = angular_radius_arcsec(z, radius_kpc)
        point = np.array([row.Group_Ra, row.Group_Dec])
        d_deg = point_to_hull_boundary_distance_deg(point, hull)
        d_arcsec = d_deg * 3600.0
        dist_to_edge_arcsec[i] = d_arcsec
        if d_arcsec < aperture_arcsec + margin_arcsec:
            is_border[i] = True

    summary["dist_to_footprint_edge_arcsec"] = dist_to_edge_arcsec
    summary["is_border_group"] = is_border
    return summary


def load_cosmos2020() -> pd.DataFrame:
    print(f"Loading COSMOS2020: {COSMOS2020_FITS}")
    cols = ["ALPHA_J2000", "DELTA_J2000", "lp_zBEST", "lp_zPDF_l68", "lp_zPDF_u68", "lp_mass_med"]
    t = Table.read(COSMOS2020_FITS)
    available = [c for c in cols if c in t.colnames]
    df = t[available].to_pandas()
    df = df.rename(columns={"ALPHA_J2000": "RA", "DELTA_J2000": "DEC", "lp_zBEST": "zfin"})
    df["ez"] = (df["lp_zPDF_u68"] - df["lp_zPDF_l68"]) / 2.0
    df = df[np.isfinite(df["zfin"]) & (df["zfin"] > 0) & np.isfinite(df["ez"]) & (df["ez"] > 0)]
    print(f"  Loaded {len(df)} valid COSMOS2020 galaxies (of {len(t)})")
    return df


def find_cosmos2020_candidates(group_ra, group_dec, group_z, cosmos2020, radius_kpc, n_sigma, max_velocity):
    radius_arcsec = angular_radius_arcsec(group_z, radius_kpc)
    group_coord = SkyCoord(ra=group_ra * u.deg, dec=group_dec * u.deg)
    gal_coords = SkyCoord(ra=cosmos2020["RA"].values * u.deg, dec=cosmos2020["DEC"].values * u.deg)
    sep = group_coord.separation(gal_coords)
    within = sep < radius_arcsec * u.arcsec
    if within.sum() == 0:
        return pd.DataFrame()
    cand = cosmos2020[within].copy()
    dz = cand["zfin"].values - group_z
    dv = C_KMS * dz / (1 + group_z)
    member = (np.abs(dz) < n_sigma * cand["ez"].values) & (np.abs(dv) < max_velocity)
    cand["dv"] = dv
    cand["dz_norm"] = dz / (1 + group_z)
    cand["is_member"] = member
    return cand[cand["is_member"]]


def main():
    parser = argparse.ArgumentParser(description="Supplement border-group membership using COSMOS2020")
    parser.add_argument("--summary-file", type=str, required=True)
    parser.add_argument("--members-file", type=str, required=True)
    parser.add_argument("--catalog-name", type=str, required=True, choices=["CW-All", "CW-HCG"])
    parser.add_argument("--radius", type=float, default=500.0, help="Search radius in kpc (must match membership run)")
    parser.add_argument("--margin-arcmin", type=float, default=1.0, help="Footprint-edge margin in arcmin")
    parser.add_argument("--n-sigma", type=float, default=3.0)
    parser.add_argument("--max-velocity", type=float, default=1000.0)
    parser.add_argument("--dedup-radius-arcsec", type=float, default=1.0,
                        help="Positional match radius to dedupe COSMOS2020 candidates against existing Webb members")
    args = parser.parse_args()

    summary_path = Path(args.summary_file)
    members_path = Path(args.members_file)
    summary = pd.read_csv(summary_path) if summary_path.suffix == ".csv" else Table.read(summary_path).to_pandas()
    members = pd.read_csv(members_path) if members_path.suffix == ".csv" else Table.read(members_path).to_pandas()
    print(f"Loaded {len(summary)} groups, {len(members)} existing members ({args.catalog_name})")

    # Footprint hull from both catalogs' group positions if available, else this catalog only
    hull_points = summary[["Group_Ra", "Group_Dec"]].dropna().values
    hull = build_footprint_hull(hull_points)
    print(f"Footprint hull built from {len(hull_points)} group positions ({len(hull.simplices)} edges)")

    summary = flag_border_groups(summary, hull, args.radius, args.margin_arcmin)
    n_border = summary["is_border_group"].sum()
    print(f"Border groups (aperture within {args.margin_arcmin} arcmin of footprint edge): {n_border}/{len(summary)}")

    if "member_source" not in members.columns:
        members["member_source"] = "Webb"

    if n_border == 0:
        print("No border groups found; nothing to supplement.")
    else:
        cosmos2020 = load_cosmos2020()
        webb_coords_by_group = {}
        new_members_list = []

        for row in tqdm(summary[summary["is_border_group"]].itertuples(), total=n_border, desc="Border groups"):
            group_id = row.Group_ID
            z = row.z_refined if np.isfinite(row.z_refined) else row.Group_z_catalog

            cand = find_cosmos2020_candidates(
                row.Group_Ra, row.Group_Dec, z, cosmos2020,
                args.radius, args.n_sigma, args.max_velocity,
            )
            if len(cand) == 0:
                continue

            existing = members[members["Group_ID"].astype(str) == str(group_id)]
            if len(existing) > 0:
                existing_coords = SkyCoord(ra=existing["RA"].values * u.deg, dec=existing["DEC"].values * u.deg)
                cand_coords = SkyCoord(ra=cand["RA"].values * u.deg, dec=cand["DEC"].values * u.deg)
                _, sep2d, _ = cand_coords.match_to_catalog_sky(existing_coords)
                is_dup = sep2d < args.dedup_radius_arcsec * u.arcsec
                cand = cand[~is_dup]

            if len(cand) == 0:
                continue

            cand = cand.copy()
            cand["Group_ID"] = group_id
            cand["Group_Ra"] = row.Group_Ra
            cand["Group_Dec"] = row.Group_Dec
            cand["Group_z_catalog"] = row.Group_z_catalog
            cand["Group_z_refined"] = z
            cand["catalog"] = args.catalog_name
            cand["member_source"] = "COSMOS2020"
            cand["redshift_type"] = "photo-z"
            new_members_list.append(cand)

        if new_members_list:
            new_members = pd.concat(new_members_list, ignore_index=True)
            print(f"New COSMOS2020 members added: {len(new_members)} across "
                  f"{new_members['Group_ID'].nunique()} border groups")
            members = pd.concat([members, new_members], ignore_index=True, sort=False)
        else:
            print("No new non-duplicate COSMOS2020 members found in border groups.")

    # Normalize mixed dtypes introduced by concatenating Webb + COSMOS2020 rows
    # (bool columns become object/NaN-mixed after concat with rows lacking them)
    for bool_col in ["photz_matched", "is_member"]:
        if bool_col in members.columns:
            members[bool_col] = members[bool_col].fillna(False).astype(bool)

    out_stub = members_path.parent / f"{members_path.stem.replace('_with_photz', '')}_with_cosmos2020"
    members.to_csv(out_stub.with_suffix(".csv"), index=False)
    Table.from_pandas(members).write(out_stub.with_suffix(".fits"), format="fits", overwrite=True)
    print(f"Saved: {out_stub.with_suffix('.csv').name}, {out_stub.with_suffix('.fits').name}")

    summary_out_stub = summary_path.parent / f"{summary_path.stem}_border_flagged"
    summary.to_csv(summary_out_stub.with_suffix(".csv"), index=False)
    print(f"Saved: {summary_out_stub.with_suffix('.csv').name}")


if __name__ == "__main__":
    main()
