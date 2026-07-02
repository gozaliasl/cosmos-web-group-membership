#!/usr/bin/env python
"""
membership-v2 Phase 3: scientific validation of the configurable
aperture/estimator module against independent quantities.

Reuses membership_v2_dynamics.py output (run that first). Produces, per
config, for CW-All:
  - sigma_v vs. independent X-ray M200 correlation (log-log Pearson/Spearman)
  - richness (AMICO LAMBDA_STAR, independent of this pipeline) vs. M200,dyn
    correlation, plus the same richness vs. X-ray M200 correlation on the
    identical subsample as a reference point
  - BGG identity stability under each config's final (trimmed) aperture,
    by actually re-scoring select_bgg_for_group on the trimmed member set
    and comparing identity (RA/DEC) to the v1.0 production BGG -- not a
    self-referential subset comparison.

Explicitly does NOT attempt a true completeness/purity estimate: no
ground-truth membership exists in the real data (this gap was flagged in
REFEREE_REVIEW.md and remains open; see Phase 5 for the planned resolution
via realistic mock lightcones). `projected_member_fraction` (in the
membership_v2_dynamics.py output) is reported here as an explicitly labeled
proxy only -- the fraction of the search-aperture candidate pool retained
after aperture trimming, not a purity or completeness measure.

Usage:
    python membership_v2_phase3_validation.py --config fixed750,r200x_trim,r200dyn_adaptive,r200dyn_adaptive_biweight,r200dyn_adaptive_vclip3
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
MEMBERSHIP_DIR = BASE_DIR / "outputs" / "results" / "membership_dztier"
V2_DIR = BASE_DIR / "outputs" / "results" / "membership_v2"

sys.path.insert(0, str(BASE_DIR / "scripts" / "pipeline"))
from select_bgg import select_bgg_for_group  # noqa: E402


def latest_wide_aperture(pattern: str) -> Path:
    candidates = [p for p in MEMBERSHIP_DIR.glob(pattern) if "_r200recomputed_" not in p.name]
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def load_richness() -> pd.Series:
    rc = Table.read(XRAY_REPO / "data" / "group-catalog" / "cosmos_web_groups_catalog_refined_z_dztier.fits").to_pandas()
    rc["Group_ID"] = rc["Group_ID"].astype(str)
    return rc.set_index("Group_ID")["LAMBDA_STAR"]


def load_xray_m200() -> pd.Series:
    xray = Table.read(XRAY_REPO / "outputs" / "results" / "cw_all" / "xray_catalog.fits").to_pandas()
    xray["Group_ID"] = xray["Group_ID"].astype(str)
    xray = xray[xray["Is_Detected"].astype(bool)]
    return xray.set_index("Group_ID")["M200_Temp_Msun"]


def sigma_v_vs_xray_m200(config: str, xray_m200: pd.Series) -> dict:
    df = pd.read_csv(V2_DIR / f"membership_v2_cw_all_{config}.csv")
    df["Group_ID"] = df["Group_ID"].astype(str)
    df = df.join(xray_m200, on="Group_ID")
    m = df[(df.n_specz >= 5) & np.isfinite(df.sigma_v_kms) & (df.sigma_v_kms > 0) &
           np.isfinite(df.M200_Temp_Msun) & (df.M200_Temp_Msun > 0)]
    if len(m) < 4:
        return dict(config=config, n=len(m), pearson_r=np.nan, pearson_p=np.nan, spearman_r=np.nan)
    r, p = pearsonr(np.log10(m.sigma_v_kms), np.log10(m.M200_Temp_Msun))
    rs, _ = spearmanr(m.sigma_v_kms, m.M200_Temp_Msun)
    return dict(config=config, n=len(m), pearson_r=r, pearson_p=p, spearman_r=rs)


def richness_vs_mass(config: str, richness: pd.Series, xray_m200: pd.Series) -> dict:
    df = pd.read_csv(V2_DIR / f"membership_v2_cw_all_{config}.csv")
    df["Group_ID"] = df["Group_ID"].astype(str)
    df = df.join(richness, on="Group_ID").join(xray_m200, on="Group_ID")

    m_dyn = df[(df.n_specz >= 5) & np.isfinite(df.M200_Msun) & (df.M200_Msun > 0) &
               np.isfinite(df.LAMBDA_STAR) & (df.LAMBDA_STAR > 0)]
    r_dyn, p_dyn = (pearsonr(np.log10(m_dyn.LAMBDA_STAR), np.log10(m_dyn.M200_Msun))
                    if len(m_dyn) >= 4 else (np.nan, np.nan))

    m_xray = df[(df.n_specz >= 5) & np.isfinite(df.M200_Temp_Msun) & (df.M200_Temp_Msun > 0) &
                np.isfinite(df.LAMBDA_STAR) & (df.LAMBDA_STAR > 0)]
    r_xray, p_xray = (pearsonr(np.log10(m_xray.LAMBDA_STAR), np.log10(m_xray.M200_Temp_Msun))
                       if len(m_xray) >= 4 else (np.nan, np.nan))

    return dict(config=config, n_richness_M200dyn=len(m_dyn), r_richness_M200dyn=r_dyn, p_richness_M200dyn=p_dyn,
                n_richness_M200xray=len(m_xray), r_richness_M200xray=r_xray, p_richness_M200xray=p_xray)


def bgg_stability(config: str, members_all: pd.DataFrame, v1_bgg: pd.DataFrame) -> dict:
    v2 = pd.read_csv(V2_DIR / f"membership_v2_cw_all_{config}.csv")
    v2["Group_ID"] = v2["Group_ID"].astype(str)
    v2 = v2[v2.n_specz >= 5].set_index("Group_ID")["aperture_used_kpc"]

    agree, total = 0, 0
    for gid, aperture in v2.items():
        if gid not in v1_bgg.index or not np.isfinite(aperture):
            continue
        gdf = members_all[members_all.Group_ID == gid]
        trimmed = gdf[gdf.sep_kpc <= aperture]
        if len(trimmed) == 0:
            continue
        scored = select_bgg_for_group(trimmed, mass_floor_n=2, centrality_kpc=300.0, centrality_tol=0.05)
        if scored["is_bgg"].sum() == 0:
            continue
        new_row = scored[scored.is_bgg].iloc[0]
        old_row = v1_bgg.loc[gid]
        same = bool(np.isclose(new_row.RA, old_row.RA) and np.isclose(new_row.DEC, old_row.DEC))
        agree += int(same)
        total += 1
    return dict(config=config, bgg_agreement_n=agree, bgg_agreement_total=total,
                bgg_agreement_frac=(agree / total if total else np.nan))


def projected_member_fraction_summary(config: str) -> dict:
    df = pd.read_csv(V2_DIR / f"membership_v2_cw_all_{config}.csv")
    m = df[df.n_specz >= 5]
    return dict(config=config,
                median_projected_member_fraction=float(np.nanmedian(m.projected_member_fraction)),
                median_contamination_score=float(np.nanmedian(m.contamination_score)),
                median_dynamical_reliability_score=float(np.nanmedian(m.dynamical_reliability_score)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str,
                         default="fixed750,r200x_trim,r200dyn_adaptive,r200dyn_adaptive_biweight,r200dyn_adaptive_vclip3")
    args = parser.parse_args()
    configs = [c.strip() for c in args.config.split(",")]

    richness = load_richness()
    xray_m200 = load_xray_m200()

    members_path = latest_wide_aperture("iterative_cw_all_members_*.csv")
    members_all = pd.read_csv(members_path)
    members_all["Group_ID"] = members_all["Group_ID"].astype(str)
    v1_bgg = members_all[members_all["is_bgg"] == True].set_index("Group_ID")  # noqa: E712

    sigma_rows, rich_rows, bgg_rows, proxy_rows = [], [], [], []
    for c in configs:
        sigma_rows.append(sigma_v_vs_xray_m200(c, xray_m200))
        rich_rows.append(richness_vs_mass(c, richness, xray_m200))
        bgg_rows.append(bgg_stability(c, members_all, v1_bgg))
        proxy_rows.append(projected_member_fraction_summary(c))

    sigma_df = pd.DataFrame(sigma_rows)
    rich_df = pd.DataFrame(rich_rows)
    bgg_df = pd.DataFrame(bgg_rows)
    proxy_df = pd.DataFrame(proxy_rows)

    print("\n=== sigma_v vs independent X-ray M200 ===")
    print(sigma_df.to_string(index=False))
    print("\n=== richness (AMICO LAMBDA_STAR) vs M200,dyn and vs X-ray M200 (reference) ===")
    print(rich_df.to_string(index=False))
    print("\n=== BGG stability under v2 apertures (re-scored, not subset comparison) ===")
    print(bgg_df.to_string(index=False))
    print("\n=== projected-member-fraction / contamination / reliability proxies (NOT purity/completeness) ===")
    print(proxy_df.to_string(index=False))

    sigma_df.to_csv(V2_DIR / "phase3_sigmav_xray_m200.csv", index=False)
    rich_df.to_csv(V2_DIR / "phase3_richness_mass.csv", index=False)
    bgg_df.to_csv(V2_DIR / "phase3_bgg_stability.csv", index=False)
    proxy_df.to_csv(V2_DIR / "phase3_quality_proxies.csv", index=False)
    print(f"\nSaved 4 summary tables to {V2_DIR}")


if __name__ == "__main__":
    main()
