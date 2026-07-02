#!/usr/bin/env python
"""
membership-v2 Phase 4: does the v2 quality-metric toolkit explain, or at
least further localize, the Nspec-dependent sigma_v-M200 degradation found
in the v1.0 audit (technical_validation_report.md Section 8-11;
REFEREE_REVIEW.md)?

Three tests, all on real CW-All data (no simulation involved -- this script
is entirely observational; TNG-side follow-up remains Phase 5's job):

  T1. Do the new per-group quality metrics (contamination_score,
      phase_space_outlier_frac, members_outside_r200_frac) scale with
      n_specz? (fixed750 config, search-aperture pool)
  T2. Is contamination systematically different between X-ray-detected and
      non-X-ray-detected groups at fixed richness (selection-bias probe)?
  T3. Does adaptive R200,dyn trimming (which measurably lowers contamination
      relative to fixed750) also repair the sigma_v-M200 correlation in the
      Nspec bins where v1.0 found it degrade? Re-run the same Nspec-binned
      correlation test v1.0 used (diagnose_nspec_contradiction.py) under
      the fixed750 and r200dyn_adaptive configs side-by-side.

Caution: several Nspec bins have very small n (down to n=4 at Nspec>=15).
Bin-level r/p values here are suggestive, not individually significant
findings on their own -- see the written interpretation in
docs/membership_v2_phase4_report.md for the epistemic framing.

Usage:
    python membership_v2_phase4_causal.py
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.stats import pearsonr, spearmanr, mannwhitneyu

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
V2_DIR = BASE_DIR / "outputs" / "results" / "membership_v2"

NSPEC_BINS = [(5, 7), (8, 9), (10, 14), (15, 999)]


def load_xray():
    xray = Table.read(XRAY_REPO / "outputs" / "results" / "cw_all" / "xray_catalog.fits").to_pandas()
    xray["Group_ID"] = xray["Group_ID"].astype(str)
    return xray


def t1_metrics_vs_nspec(df: pd.DataFrame) -> pd.DataFrame:
    m = df[df.n_specz >= 5]
    rows = []
    for metric in ["contamination_score", "phase_space_outlier_frac", "members_outside_r200_frac"]:
        r, p = spearmanr(m.n_specz, m[metric])
        rows.append(dict(metric=metric, spearman_r=r, p=p, n=len(m)))
    return pd.DataFrame(rows)


def t2_xray_selection_probe(df: pd.DataFrame, xray: pd.DataFrame) -> dict:
    xray_det = xray.set_index("Group_ID")["Is_Detected"].astype(bool)
    m = df[df.n_specz >= 5].join(xray_det.rename("is_xray_detected"), on="Group_ID")
    m["is_xray_detected"] = m["is_xray_detected"].fillna(False)
    xd = m[m.is_xray_detected]
    nxd = m[~m.is_xray_detected]
    if len(xd) < 3 or len(nxd) < 3:
        return dict(n_xray=len(xd), n_non_xray=len(nxd), median_contam_xray=np.nan,
                     median_contam_non_xray=np.nan, mannwhitney_p=np.nan)
    u, p = mannwhitneyu(xd.contamination_score.dropna(), nxd.contamination_score.dropna())
    return dict(n_xray=len(xd), n_non_xray=len(nxd),
                median_contam_xray=float(xd.contamination_score.median()),
                median_contam_non_xray=float(nxd.contamination_score.median()),
                median_nspec_xray=float(xd.n_specz.median()),
                median_nspec_non_xray=float(nxd.n_specz.median()),
                mannwhitney_p=p)


def t3_nspec_binned_correlation(config: str, xray_m200: pd.Series) -> pd.DataFrame:
    df = pd.read_csv(V2_DIR / f"membership_v2_cw_all_{config}.csv")
    df["Group_ID"] = df["Group_ID"].astype(str)
    df = df.join(xray_m200, on="Group_ID")
    rows = []
    for lo, hi in NSPEC_BINS:
        m = df[(df.n_specz >= lo) & (df.n_specz <= hi) & np.isfinite(df.sigma_v_kms) & (df.sigma_v_kms > 0) &
               np.isfinite(df.M200_Temp_Msun) & (df.M200_Temp_Msun > 0)]
        if len(m) >= 4:
            r, p = pearsonr(np.log10(m.sigma_v_kms), np.log10(m.M200_Temp_Msun))
        else:
            r, p = np.nan, np.nan
        rows.append(dict(config=config, nspec_bin=f"{lo}-{hi}", n=len(m), pearson_r=r, p=p,
                          median_contamination=float(m.contamination_score.median()) if len(m) else np.nan))
    return pd.DataFrame(rows)


def main():
    xray = load_xray()
    xray_m200 = xray[xray.Is_Detected.astype(bool)].set_index("Group_ID")["M200_Temp_Msun"]

    df_fixed = pd.read_csv(V2_DIR / "membership_v2_cw_all_fixed750.csv")
    df_fixed["Group_ID"] = df_fixed["Group_ID"].astype(str)

    print("=== T1: quality metrics vs n_specz (fixed750) ===")
    t1 = t1_metrics_vs_nspec(df_fixed)
    print(t1.to_string(index=False))
    t1.to_csv(V2_DIR / "phase4_t1_metrics_vs_nspec.csv", index=False)

    print("\n=== T2: contamination_score, X-ray-detected vs non-detected ===")
    t2 = t2_xray_selection_probe(df_fixed, xray)
    print(pd.Series(t2))
    pd.DataFrame([t2]).to_csv(V2_DIR / "phase4_t2_xray_selection_probe.csv", index=False)

    print("\n=== T3: Nspec-binned sigma_v-M200 correlation, fixed750 vs r200dyn_adaptive ===")
    t3 = pd.concat([t3_nspec_binned_correlation("fixed750", xray_m200),
                     t3_nspec_binned_correlation("r200dyn_adaptive", xray_m200)], ignore_index=True)
    print(t3.to_string(index=False))
    t3.to_csv(V2_DIR / "phase4_t3_nspec_binned_correlation.csv", index=False)


if __name__ == "__main__":
    main()
