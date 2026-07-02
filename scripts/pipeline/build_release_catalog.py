#!/usr/bin/env python
"""
membership_release_v1.0 -- build the definitive production membership
catalogs (CW-All, CW-HCG).

Production methodology (membership_pipeline_v1.0) is UNCHANGED:
  - membership: probabilistic P_z*P_v model, fixed 750 kpc search aperture
    (CW-All) / compact-group-scaled aperture (CW-HCG)
  - BGG: mass-floored hybrid rank score (select_bgg.py)
  - dynamical mass: gapper sigma_v estimator, Munari et al. (2013) interim
    calibration (compute_group_r200.py, cosmos-web-xray-igm repo)
  - X-ray properties: xray_pipeline_v1.1_production
    (cosmos-web-xray-igm/outputs/release_v1.1)

This script only MERGES already-produced v1.0 outputs (membership, BGG,
R200/M200/sigma_v/method) into a single release table, VERIFIES catalog
integrity, and ADDS informational-only quality columns from the v2
investigation (membership-v2 branch, Project B, D18) -- it does not alter
any production value and does not re-run or modify the membership,
BGG, or dynamical-mass algorithms. The v2 informational columns are
computed via membership_v2_dynamics.py's `fixed750` configuration only
(the config that exactly reproduces v1.0 by construction).

Usage:
    python build_release_catalog.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
MEMBERSHIP_DIR = BASE_DIR / "outputs" / "results" / "membership_dztier"
R200_DIR = XRAY_REPO / "outputs" / "results" / "r200_catalog"
XRAY_RELEASE_DIR = XRAY_REPO / "outputs" / "release_v1.1"
V2_DIR = BASE_DIR / "outputs" / "results" / "membership_v2"

RELEASE_DIR = BASE_DIR / "outputs" / "release" / "membership_release_v1.0"
RELEASE_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(Path(__file__).parent))
from select_bgg import select_bgg_for_group  # noqa: E402

NSPEC_DEGRADATION_LO, NSPEC_DEGRADATION_HI = 8, 14
MIN_SPECZ_FOR_DYNAMICS = 5


def latest_wide_aperture(pattern: str) -> Path:
    candidates = [p for p in MEMBERSHIP_DIR.glob(pattern) if "_r200recomputed_" not in p.name]
    if not candidates:
        raise FileNotFoundError(f"No file matching {pattern} in {MEMBERSHIP_DIR}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def nspec_regime(n_specz):
    if pd.isna(n_specz):
        return "n/a"
    return "degradation-regime" if NSPEC_DEGRADATION_LO <= n_specz <= NSPEC_DEGRADATION_HI else "standard"


def dynamical_mass_confidence(row):
    """Per Phase 7 (membership-v2, D18): default LOW for SPECZ-only method
    groups; XRAY-anchored methods are not dynamical-mass-dependent so are
    marked n/a (X-ray mass, not dynamical mass, is the science value)."""
    method = row["r200_method"]
    if method in ("XRAY", "STACKED_XRAY"):
        return "n/a (X-ray-derived mass, not dynamical)"
    if method == "XRAY+SPECZ" and not row.get("r200_disagreement", False):
        return "medium (X-ray-anchored, dynamical agreement confirmed)"
    if method == "XRAY+SPECZ" and row.get("r200_disagreement", False):
        return "low (X-ray/dynamical disagreement -- X-ray value adopted)"
    if method == "SPECZ":
        if row.get("nspec_regime") == "degradation-regime":
            return "low (spec-z-only, degradation-regime richness)"
        return "low (spec-z-only)"
    return "n/a"


def sigma_v_quality(row):
    n = row.get("n_specz_members", np.nan)
    if pd.isna(n) or n < MIN_SPECZ_FOR_DYNAMICS:
        return "n/a (below n_specz>=5 floor)"
    if row.get("hit_iteration_cap", False) or row.get("r200_exceeds_search_aperture", False):
        return "low (aperture/iteration flag set)"
    if n >= 15:
        return "medium (n_specz>=15, but see nspec_regime)"
    if nspec_regime(n) == "degradation-regime":
        return "low (degradation-regime richness)"
    return "medium"


def bgg_stability_for_group(members_all_group: pd.DataFrame, apertures_kpc: list) -> str:
    """Re-scores select_bgg_for_group at each of 3 apertures (matching
    Phase 3's independent, re-scored BGG-stability test) and reports
    whether the identity is stable. Informational only."""
    v1_bgg = members_all_group[members_all_group["is_bgg"] == True]  # noqa: E712
    if len(v1_bgg) == 0:
        return "n/a (no BGG assigned)"
    v1_ra, v1_dec = v1_bgg.iloc[0]["RA"], v1_bgg.iloc[0]["DEC"]

    agreements = []
    for ap in apertures_kpc:
        trimmed = members_all_group[members_all_group["sep_kpc"] <= ap]
        if len(trimmed) == 0:
            continue
        scored = select_bgg_for_group(trimmed, mass_floor_n=2, centrality_kpc=300.0, centrality_tol=0.05)
        if scored["is_bgg"].sum() == 0:
            agreements.append(False)
            continue
        new_row = scored[scored.is_bgg].iloc[0]
        agreements.append(bool(np.isclose(new_row.RA, v1_ra) and np.isclose(new_row.DEC, v1_dec)))
    if not agreements:
        return "n/a"
    frac = np.mean(agreements)
    return "stable" if frac == 1.0 else ("mostly_stable" if frac >= 0.5 else "unstable")


def build_catalog(catalog_name: str, xray_slug: str, membership_prefix: str, apertures_for_bgg_test: list):
    print(f"\n{'='*70}\n{catalog_name}\n{'='*70}")

    members_path = latest_wide_aperture(f"{membership_prefix}members_*.csv")
    summary_path = latest_wide_aperture(f"{membership_prefix}summary_*.csv")
    members = pd.read_csv(members_path)
    summary = pd.read_csv(summary_path)
    members["Group_ID"] = members["Group_ID"].astype(str)
    summary["Group_ID"] = summary["Group_ID"].astype(str)
    print(f"  Membership source: {members_path.name} ({len(members)} members, {len(summary)} groups)")

    r200 = pd.read_csv(R200_DIR / f"r200_catalog_{xray_slug}.csv")
    r200["Group_ID"] = r200["Group_ID"].astype(str)
    print(f"  R200/M200/sigma_v source: r200_catalog_{xray_slug}.csv ({len(r200)} groups)")

    v2_path = V2_DIR / f"membership_v2_{xray_slug}_fixed750.csv"
    v2 = pd.read_csv(v2_path) if v2_path.exists() else None
    v2["Group_ID"] = v2["Group_ID"].astype(str) if v2 is not None else None
    print(f"  v2 informational-quality source: {v2_path.name if v2 is not None else 'MISSING'}")

    # ---------------- group-level release table ----------------
    group_tbl = summary.merge(r200, on="Group_ID", how="outer", suffixes=("", "_r200"))
    group_tbl["nspec_regime"] = group_tbl["n_specz_members"].apply(nspec_regime)
    group_tbl["dynamical_mass_confidence"] = group_tbl.apply(dynamical_mass_confidence, axis=1)
    group_tbl["sigma_v_quality"] = group_tbl.apply(sigma_v_quality, axis=1)

    if v2 is not None:
        v2_cols = v2.set_index("Group_ID")[["contamination_score", "phase_space_outlier_frac",
                                             "dynamical_reliability_score", "membership_confidence"]]
        v2_cols = v2_cols.rename(columns={"membership_confidence": "membership_quality"})
        group_tbl = group_tbl.merge(v2_cols, on="Group_ID", how="left")
    else:
        for c in ["contamination_score", "phase_space_outlier_frac", "dynamical_reliability_score", "membership_quality"]:
            group_tbl[c] = np.nan

    bgg_stability_map = {}
    for gid, gdf in members.groupby("Group_ID"):
        bgg_stability_map[gid] = bgg_stability_for_group(gdf, apertures_for_bgg_test)
    group_tbl["BGG_stability"] = group_tbl["Group_ID"].map(bgg_stability_map)

    # ---------------- verification ----------------
    report = verify_catalog(catalog_name, members, group_tbl)

    # ---------------- write outputs ----------------
    members_stub = f"{catalog_name.replace('-', '_')}_members_release_v1.0"
    groups_stub = f"{catalog_name.replace('-', '_')}_groups_release_v1.0"
    members.to_csv(RELEASE_DIR / f"{members_stub}.csv", index=False)
    Table.from_pandas(members).write(RELEASE_DIR / f"{members_stub}.fits", format="fits", overwrite=True)
    group_tbl.to_csv(RELEASE_DIR / f"{groups_stub}.csv", index=False)
    Table.from_pandas(group_tbl).write(RELEASE_DIR / f"{groups_stub}.fits", format="fits", overwrite=True)
    print(f"  Saved: {members_stub}.[csv,fits], {groups_stub}.[csv,fits]")

    return group_tbl, members, report


def verify_catalog(catalog_name: str, members: pd.DataFrame, group_tbl: pd.DataFrame) -> dict:
    print(f"\n  --- verification: {catalog_name} ---")
    report = {}

    report["n_groups"] = group_tbl["Group_ID"].nunique()
    report["group_id_duplicates"] = int(group_tbl["Group_ID"].duplicated().sum())

    members_only = members[members["is_member"] == True]  # noqa: E712
    dup_key = members_only[["Group_ID", "RA", "DEC"]].duplicated().sum()
    report["duplicate_member_rows_within_group"] = int(dup_key)

    # cross-group multiplicity: a galaxy assigned as a member of >1 group.
    # NOT a bug -- v1.0's probabilistic P_z*P_v membership model is
    # non-exclusive by design (nearby groups' apertures/velocity windows can
    # overlap), but this must be reported transparently, not silently
    # treated as "no duplicates."
    per_galaxy_group_count = members_only.groupby(["RA", "DEC"])["Group_ID"].nunique()
    report["unique_galaxies_total"] = int(len(per_galaxy_group_count))
    report["galaxies_member_of_multiple_groups"] = int((per_galaxy_group_count > 1).sum())
    report["galaxies_member_of_multiple_groups_frac"] = float((per_galaxy_group_count > 1).mean())

    bgg_counts = members_only.groupby("Group_ID")["is_bgg"].sum()
    report["groups_with_zero_bgg"] = int((bgg_counts == 0).sum())
    report["groups_with_multiple_bgg"] = int((bgg_counts > 1).sum())
    report["groups_with_exactly_one_bgg"] = int((bgg_counts == 1).sum())

    bgg_rows = members_only[members_only["is_bgg"] == True]  # noqa: E712
    bgg_multiplicity = bgg_rows.groupby(["RA", "DEC"])["Group_ID"].nunique()
    report["bgg_galaxies_used_as_bgg_for_multiple_groups"] = int((bgg_multiplicity > 1).sum())
    report["bgg_galaxies_total"] = int(len(bgg_multiplicity))

    n_specz_check = members_only[members_only["redshift_type"] == "spec-z"].groupby("Group_ID").size()
    summary_specz = group_tbl.set_index("Group_ID")["n_members_specz"] if "n_members_specz" in group_tbl.columns else None
    if summary_specz is not None:
        merged_check = pd.DataFrame({"recount": n_specz_check}).join(summary_specz, how="outer").fillna(0)
        mismatch = (merged_check["recount"] != merged_check["n_members_specz"]).sum()
        report["n_specz_recount_mismatches"] = int(mismatch)
    else:
        report["n_specz_recount_mismatches"] = "n/a (column missing)"

    # sigma_v is legitimately N/A for XRAY-only and STACKED_XRAY method groups
    # (no dynamical estimate computed for them by design) -- only flag it missing
    # where a method that SHOULD carry sigma_v (SPECZ, XRAY+SPECZ) lacks it.
    dyn_methods = group_tbl["r200_method"].isin(["SPECZ", "XRAY+SPECZ"])
    report["sigma_v_kms_missing_for_dynamical_methods"] = int(group_tbl.loc[dyn_methods, "sigma_v_kms"].isna().sum())
    report["n_groups_dynamical_methods"] = int(dyn_methods.sum())

    for col in ["R200_kpc", "M200_Msun"]:
        if col in group_tbl.columns:
            has_method = group_tbl["r200_method"].notna()
            missing = group_tbl.loc[has_method, col].isna().sum()
            report[f"{col}_missing_where_method_assigned"] = int(missing)
            if missing > 0:
                bad_ids = group_tbl.loc[has_method & group_tbl[col].isna(), "Group_ID"].tolist()
                report[f"{col}_missing_group_ids"] = bad_ids

    method_null = int(group_tbl["r200_method"].isna().sum())
    report["groups_with_null_r200_method"] = method_null

    print(f"    n_groups: {report['n_groups']}")
    print(f"    Group_ID duplicates: {report['group_id_duplicates']}")
    print(f"    duplicate member (Group_ID,RA,DEC) rows: {report['duplicate_member_rows_within_group']}")
    print(f"    galaxies that are members of >1 group (non-exclusive membership, by design): "
          f"{report['galaxies_member_of_multiple_groups']} / {report['unique_galaxies_total']} "
          f"({100*report['galaxies_member_of_multiple_groups_frac']:.1f}%)")
    print(f"    groups with 0 BGG / >1 BGG / exactly 1 BGG: "
          f"{report['groups_with_zero_bgg']} / {report['groups_with_multiple_bgg']} / {report['groups_with_exactly_one_bgg']}")
    print(f"    BGG galaxies that serve as BGG for >1 group (inherits cross-group multiplicity above): "
          f"{report['bgg_galaxies_used_as_bgg_for_multiple_groups']} / {report['bgg_galaxies_total']}")
    print(f"    n_specz recount mismatches vs. summary: {report['n_specz_recount_mismatches']}")
    print(f"    groups with null r200_method: {report['groups_with_null_r200_method']}")
    print(f"    sigma_v missing among SPECZ/XRAY+SPECZ (dynamical) methods: "
          f"{report['sigma_v_kms_missing_for_dynamical_methods']} / {report['n_groups_dynamical_methods']}")
    for col in ["R200_kpc", "M200_Msun"]:
        k = f"{col}_missing_where_method_assigned"
        if k in report:
            print(f"    {col} missing where method assigned: {report[k]}"
                  + (f" (Group_IDs: {report[col + '_missing_group_ids']})" if report[k] > 0 else ""))

    return report


def main():
    all_group_tbl, all_members, all_report = build_catalog(
        "CW-All", "cw_all", "iterative_cw_all_", apertures_for_bgg_test=[750.0, 700.0, 650.0])
    hcg_group_tbl, hcg_members, hcg_report = build_catalog(
        "CW-HCG", "cw_hcg", "iterative_cw_hcg_", apertures_for_bgg_test=[250.0, 200.0, 150.0])

    import json
    with open(RELEASE_DIR / "verification_report.json", "w") as f:
        json.dump(dict(CW_All=all_report, CW_HCG=hcg_report), f, indent=2, default=str)
    print(f"\nVerification report written to {RELEASE_DIR / 'verification_report.json'}")
    print(f"\nAll release outputs in {RELEASE_DIR}")


if __name__ == "__main__":
    main()
