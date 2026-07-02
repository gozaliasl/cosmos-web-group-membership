#!/usr/bin/env python
"""
Select the Brightest Group Galaxy (BGG) for each group.

Algorithm (mass-floored hybrid rank score with soft centrality preference):

    1. Mass floor: restrict candidates to the `--mass-floor-n` (default 2)
       most massive members by mass_final (LP_mass_med_PDF for Webb-native
       members; COSMOS2020 lp_mass_med fallback for footprint-margin members
       -- see fill_cosmos2020_properties.py). This prevents a low-mass but
       UV/blue-bright galaxy from being picked purely on luminosity.
    2. Hybrid score per candidate (within the mass-floor pool only):
           score = 0.5 * mass_rank(percentile) + 0.5 * (-MR_final)_rank(percentile)
       i.e. equal weight between stellar-mass rank and R-band luminosity rank
       (brighter = higher rank, since MR_final is an absolute magnitude and
       more negative = brighter).
    3. Centrality preference: if the best-scoring candidate within
       `--centrality-kpc` (default 300 kpc) of the group center scores within
       `--centrality-tol` (default 0.05) of the global best (unrestricted)
       score, the inner galaxy is preferred; otherwise the global best wins.

Net effect: the BGG is the more R-band-luminous of the ~2 most massive
members, with a preference for the central one when scores are close.

Requires `mass_final`, `MR_final`, `sep_kpc` columns (mass_final/MR_final
from fill_cosmos2020_properties.py; sep_kpc from determine_membership_dztier.py).

Usage:
    python select_bgg.py --members-file outputs/results/membership_dztier/fullfield_prob_cw_all_members_*_massfilled.csv
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

warnings.filterwarnings("ignore")


def select_bgg_for_group(group_df: pd.DataFrame, mass_floor_n: int,
                          centrality_kpc: float, centrality_tol: float) -> pd.DataFrame:
    group_df = group_df.copy()
    group_df["is_bgg"] = False
    group_df["bgg_score"] = np.nan

    candidates = group_df.dropna(subset=["mass_final", "MR_final"])
    if len(candidates) == 0:
        return group_df

    # Step 1: mass floor
    pool = candidates.nlargest(mass_floor_n, "mass_final").copy()

    # Step 2: hybrid percentile-rank score within the mass-floor pool
    mass_pct = pool["mass_final"].rank(pct=True, ascending=True, method="average")
    mr_pct = (-pool["MR_final"]).rank(pct=True, ascending=True, method="average")
    pool["bgg_score"] = 0.5 * mass_pct + 0.5 * mr_pct
    group_df.loc[pool.index, "bgg_score"] = pool["bgg_score"]

    # Step 3: global best within the mass-floor pool
    global_best_idx = pool["bgg_score"].idxmax()
    global_best_score = pool.loc[global_best_idx, "bgg_score"]

    # Step 4: soft centrality preference
    bgg_idx = global_best_idx
    used_centrality_preference = False
    inner = pool[pool["sep_kpc"] <= centrality_kpc]
    if len(inner) > 0:
        inner_best_idx = inner["bgg_score"].idxmax()
        inner_best_score = inner.loc[inner_best_idx, "bgg_score"]
        if (global_best_score - inner_best_score) <= centrality_tol and inner_best_idx != global_best_idx:
            bgg_idx = inner_best_idx
            used_centrality_preference = True

    group_df.loc[bgg_idx, "is_bgg"] = True
    group_df.attrs["bgg_mass_floor_pool_size"] = len(pool)
    group_df.attrs["bgg_used_centrality_preference"] = used_centrality_preference
    return group_df


def main():
    parser = argparse.ArgumentParser(description="Select BGG per group: mass-floored hybrid mass+R-luminosity rank score, soft centrality preference")
    parser.add_argument("--members-file", type=str, required=True)
    parser.add_argument("--mass-floor-n", type=int, default=2,
                        help="Restrict BGG candidates to the N most massive members before scoring")
    parser.add_argument("--centrality-kpc", type=float, default=300.0,
                        help="Radius (kpc) defining the 'central' candidate pool for the soft preference")
    parser.add_argument("--centrality-tol", type=float, default=0.05,
                        help="Prefer the best central candidate if its score is within this of the global best")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    members_path = Path(args.members_file)
    members = (Table.read(members_path).to_pandas() if members_path.suffix.lower() == ".fits"
               else pd.read_csv(members_path))
    print(f"Loaded {len(members)} members from {members_path.name}")

    for col in ["mass_final", "MR_final", "sep_kpc", "Group_ID"]:
        if col not in members.columns:
            raise SystemExit(f"Required column '{col}' not found in members file")

    groups_out = []
    n_centrality_pref = 0
    n_no_candidates = 0
    for group_id, group_df in members.groupby("Group_ID", sort=False):
        scored = select_bgg_for_group(group_df, args.mass_floor_n, args.centrality_kpc, args.centrality_tol)
        if scored["is_bgg"].sum() == 0:
            n_no_candidates += 1
        elif scored.attrs.get("bgg_used_centrality_preference"):
            n_centrality_pref += 1
        groups_out.append(scored)

    out = pd.concat(groups_out, ignore_index=True)
    n_groups = out["Group_ID"].nunique()
    n_bgg = int(out["is_bgg"].sum())
    print(f"\nBGG selected for {n_bgg}/{n_groups} groups")
    print(f"  Groups where centrality preference flipped the choice: {n_centrality_pref}")
    print(f"  Groups with no valid mass+MR candidates (no BGG assigned): {n_no_candidates}")

    out_stub = Path(args.output) if args.output else members_path.parent / f"{members_path.stem}_bgg"
    out.to_csv(out_stub.with_suffix(".csv"), index=False)
    Table.from_pandas(out).write(out_stub.with_suffix(".fits"), format="fits", overwrite=True)
    print(f"\nSaved: {out_stub.with_suffix('.csv').name}, {out_stub.with_suffix('.fits').name}")


if __name__ == "__main__":
    main()
