#!/usr/bin/env python
"""
Re-derive total group stellar mass and re-select the BGG using the physically
motivated R200 (from cosmos-web-xray-igm/scripts/compute_group_r200.py:
X-ray-detected -> X-ray R200; spec-z-rich non-detected -> dynamical R200
via the HZ-AGN-calibrated sigma_v-M200 relation; else -> stacked-X-ray R200)
instead of the flat 750 kpc / compact-group aperture used for the initial
membership search.

For each group:
  1. Restrict members to sep_kpc <= R200_kpc (the group's own physically
     motivated boundary), from the existing iterative membership selection --
     no RVF/probabilistic membership rerun, just a spatial re-filter of
     members already found within the (larger) 750 kpc / compact-group net.
  2. Total stellar mass = sum(10**mass_final) over those R200 members with a
     valid mass_final (Webb LePhare LP_mass_med_PDF, or COSMOS2020 lp_mass_med
     fallback for true-margin members).
  3. BGG re-selected within the R200-restricted member set using the same
     mass-floored hybrid rank-score method (select_bgg.select_bgg_for_group),
     with the centrality-preference radius scaled to the group's own R200
     (min(300 kpc, 0.5*R200)) rather than a fixed 300/100 kpc.

Usage:
    python recompute_mass_bgg_with_r200.py --catalog-name CW-All \
        --members-file outputs/results/membership_dztier/iterative_cw_all_members_*.csv \
        --r200-file ../cosmos-web-xray-igm/outputs/results/r200_catalog/r200_catalog_cw_all.csv
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).parent))
from select_bgg import select_bgg_for_group  # noqa: E402

warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser(description="Re-derive total stellar mass + BGG using physically motivated R200")
    parser.add_argument("--catalog-name", type=str, required=True)
    parser.add_argument("--members-file", type=str, required=True)
    parser.add_argument("--r200-file", type=str, required=True)
    parser.add_argument("--mass-floor-n", type=int, default=2)
    parser.add_argument("--centrality-tol", type=float, default=0.05)
    parser.add_argument("--max-centrality-kpc", type=float, default=300.0,
                        help="Cap on the R200-scaled centrality-preference radius")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    members_path = Path(args.members_file)
    members = pd.read_csv(members_path) if members_path.suffix == ".csv" else Table.read(members_path).to_pandas()
    r200 = pd.read_csv(args.r200_file)
    members["Group_ID"] = members["Group_ID"].astype(str)
    r200["Group_ID"] = r200["Group_ID"].astype(str)
    print(f"Loaded {len(members)} members, {r200['Group_ID'].nunique()} groups with R200")

    r200_idx = r200.set_index("Group_ID")

    group_rows = []
    r200_members_list = []
    for group_id, group_df in members.groupby("Group_ID", sort=False):
        if group_id not in r200_idx.index or not np.isfinite(r200_idx.loc[group_id, "R200_kpc"]):
            continue
        r200_kpc = r200_idx.loc[group_id, "R200_kpc"]
        r200_method = r200_idx.loc[group_id, "r200_method"]

        r200_members = group_df[group_df["sep_kpc"] <= r200_kpc].copy()
        n_total = len(r200_members)
        n_specz = int((r200_members["redshift_type"] == "spec-z").sum())

        valid_mass = r200_members["mass_final"].dropna()
        total_mstar = float(np.sum(10 ** valid_mass)) if len(valid_mass) > 0 else np.nan
        log_total_mstar = np.log10(total_mstar) if total_mstar > 0 else np.nan

        centrality_kpc = min(args.max_centrality_kpc, 0.5 * r200_kpc)
        bgg_mstar, bgg_ra, bgg_dec, bgg_sep_kpc = np.nan, np.nan, np.nan, np.nan
        if n_total > 0:
            scored = select_bgg_for_group(r200_members, args.mass_floor_n, centrality_kpc, args.centrality_tol)
            bgg_row = scored[scored["is_bgg"]]
            if len(bgg_row) > 0:
                bgg_row = bgg_row.iloc[0]
                bgg_mstar = bgg_row["mass_final"]
                bgg_ra, bgg_dec = bgg_row["RA"], bgg_row["DEC"]
                bgg_sep_kpc = bgg_row["sep_kpc"]
            r200_members_list.append(scored)

        group_rows.append(dict(
            Group_ID=group_id, catalog=args.catalog_name, R200_kpc=r200_kpc, r200_method=r200_method,
            n_members_within_r200=n_total, n_specz_within_r200=n_specz,
            log_total_Mstar_Msun=log_total_mstar, total_Mstar_Msun=total_mstar,
            log_BGG_Mstar_Msun=bgg_mstar, BGG_RA=bgg_ra, BGG_DEC=bgg_dec, BGG_sep_kpc=bgg_sep_kpc,
            centrality_kpc_used=centrality_kpc,
        ))

    summary = pd.DataFrame(group_rows)
    r200_members_all = pd.concat(r200_members_list, ignore_index=True) if r200_members_list else pd.DataFrame()

    print(f"\nGroups processed: {len(summary)}")
    print(f"  Median total Mstar: 10^{summary['log_total_Mstar_Msun'].median():.2f} Msun")
    print(f"  Median BGG Mstar: 10^{summary['log_BGG_Mstar_Msun'].median():.2f} Msun")
    print(f"  Median members within R200: {summary['n_members_within_r200'].median():.1f}")
    print(f"  Groups with 0 members within R200: {(summary['n_members_within_r200']==0).sum()}")

    out_stub = Path(args.output) if args.output else members_path.parent / f"{members_path.stem}_r200recomputed"
    summary_csv = out_stub.parent / f"{out_stub.name}_summary.csv"
    summary_fits = out_stub.parent / f"{out_stub.name}_summary.fits"
    summary.to_csv(summary_csv, index=False)
    Table.from_pandas(summary).write(summary_fits, format="fits", overwrite=True)
    if len(r200_members_all) > 0:
        members_csv = out_stub.parent / f"{out_stub.name}_members.csv"
        members_fits = out_stub.parent / f"{out_stub.name}_members.fits"
        r200_members_all.to_csv(members_csv, index=False)
        Table.from_pandas(r200_members_all).write(members_fits, format="fits", overwrite=True)
    print(f"\nSaved: {summary_csv.name}, {summary_fits.name}, and members files")


if __name__ == "__main__":
    main()
