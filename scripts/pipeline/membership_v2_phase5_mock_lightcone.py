#!/usr/bin/env python
"""
membership-v2 Phase 5: realistic mock lightcone + sigma_v-M200 calibration
diagnostic test.

Per D14 (DECISION_LOG_V2.md): this is a *calibration-diagnostic* test, not
a generic realism upgrade. Ground truth (ACTUAL halo M200, ACTUAL member
peculiar velocities) comes directly from the TNG50-1 + TNG100-1 "enriched
local catalog" (external, not part of this repo -- see --tng-dir), so we
are testing recovery against real simulation ground truth, not a re-derived
one.

Pipeline per mock group (built from one TNG halo with >= min_richness
member subhalos):
  1. True members: TNG's own RA_mock/DEC_mock (pre-computed mock sky
     projection) + vel_z_kms (real simulation peculiar velocity) ->
     z_obs_true = z_snap + vel_z_kms/c*(1+z_snap).
  2. Spectroscopic incompleteness: randomly retain a `completeness_frac`
     subset of true members (COSMOS spec-z completeness is far from 100%;
     tested at 0.3 / 0.5 / 0.7).
  3. Redshift uncertainty: add Gaussian velocity-space noise
     (z_err_kms, default 150 km/s, matching the spec-z precision floor used
     throughout the production pipeline, SIGMA_V_FLOOR_KMS) to each
     retained true member's observed velocity.
  4. Field contamination: inject interloper subhalos drawn from OTHER,
     unrelated TNG halos at the same snapshot. Interlopers are given: (a) a
     projected sky position drawn uniformly within the group's search
     aperture (representing genuine chance projection), and (b) an observed
     velocity offset drawn from Uniform(-2500, +2500) km/s around the group
     redshift -- i.e. uniformly filling the production pipeline's actual
     velocity acceptance window (SIGMA_V_CEIL_KMS=2500), representing a
     photo-z-consistent field population that was never excluded at the
     photo-z stage (matches the real-data finding that photo-z 1-sigma
     uncertainty is itself ~3000+ km/s -- D1 in the v1.0 Decision Log).
     `n_field_per_true_member` controls the contamination rate (default 0.5,
     i.e. one interloper injected for every two true members on average).
  5. Membership selection: reproduces v1.0's actual default (fixed
     750 kpc aperture + P_v probabilistic velocity gate,
     P_v = exp(-0.5*(dv/sigma_prior)^2), sigma_prior=500 km/s,
     prob_threshold=0.05 -- same functional form and same default constants
     as determine_membership_dztier.py's apply_membership_cut(); this is
     the actual production algorithm's velocity-gating logic, applied here
     to the mock catalog).
  6. Recovered sigma_v: gapper estimator (v1.0's choice) on the final
     selected member set.
  7. Mass recovery: apply EACH of 4 candidate sigma_v-M200 calibrations to
     the recovered sigma_v, compare to the TRUE logM200_msun (ground truth,
     not derived).

Calibrations tested
--------------------
  munari2013_galaxies : Munari et al. (2013), galaxy tracers. Currently the
                         production interim relation (D7).
                         sigma_15=1177.0 km/s, alpha=0.364.
  evrard2008_dm        : Evrard et al. (2008), dark-matter-particle
                         calibration (ApJ 672, 122). sigma_DM,15=1082.9
                         km/s, alpha=0.3361. External, independent of this
                         project's simulation suite.
  tng_selfcal_all      : This-work self-calibration, power-law fit of
                         log10(sigma_v) vs. log10(M200) on ALL TNG halos
                         with >= min_richness members (the same sample used
                         for the mock lightcone itself) -- a genuinely
                         simulation-derived, in-sample calibration.
  tng_selfcal_groups   : Same fit restricted to the group-scale mass range
                         only (log M200 < 13.5, matching this project's
                         actual regime) -- the "group-specific, lower-mass"
                         calibration requested; fit directly from TNG rather
                         than an external literature value, to avoid citing
                         an external group-scale calibration without full
                         confidence in its exact published parameters.

Usage:
    python membership_v2_phase5_mock_lightcone.py --tng-dir /Volumes/extHD/tng_local_catalog
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as cosmo
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
OUTPUT_DIR = BASE_DIR / "outputs" / "results" / "membership_v2_phase5"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(XRAY_REPO / "src"))
from xray_analysis.mass_estimation import calculate_radius_from_mass  # noqa: E402

C_KMS = 299792.458
SEARCH_RADIUS_KPC = 750.0        # matches v1.0 CW-All default (fixed750)
SIGMA_PRIOR_KMS = 500.0          # matches v1.0's default_sigma_v_kms
PROB_THRESHOLD = 0.05            # matches v1.0's prob_threshold
Z_ERR_KMS_DEFAULT = 150.0        # matches SIGMA_V_FLOOR_KMS (spec-z precision floor)
VELOCITY_CEIL_KMS = 2500.0       # matches v1.0's SIGMA_V_CEIL_KMS
MIN_RICHNESS = 5                 # matches MIN_SPECZ_FOR_DYNAMICS


def gapper_sigma_v(v: np.ndarray) -> float:
    v = np.sort(v)
    n = len(v)
    if n < 2:
        return np.nan
    gaps = np.diff(v)
    w = np.arange(1, n) * (n - np.arange(1, n))
    return float(np.sqrt(np.pi) / (n * (n - 1)) * np.sum(w * gaps))


# --------------------------------------------------------------------------
# candidate calibrations: all take (sigma_v_kms, z) -> M200_Msun
# --------------------------------------------------------------------------

def _munari2013(sigma_v_kms, z):
    hz = (cosmo.H(z) / cosmo.H0).value
    return (1e15 / hz) * (sigma_v_kms / 1177.0) ** (1.0 / 0.364)


def _evrard2008(sigma_v_kms, z):
    hz = (cosmo.H(z) / cosmo.H0).value
    return (1e15 / hz) * (sigma_v_kms / 1082.9) ** (1.0 / 0.3361)


def _make_selfcal(sigma_15: float, alpha: float):
    def _fn(sigma_v_kms, z):
        hz = (cosmo.H(z) / cosmo.H0).value
        return (1e15 / hz) * (sigma_v_kms / sigma_15) ** (1.0 / alpha)
    return _fn


def fit_selfcal(df: pd.DataFrame) -> tuple:
    """Fit sigma(M,z) = sigma_15 * [h(z)*M/1e15]^alpha via linear regression
    in log space: log10(sigma) = log10(sigma_15) + alpha*log10(h(z)*M/1e15)."""
    hz = (cosmo.H(df["snapshot_redshift"].values) / cosmo.H0).value
    x = np.log10(hz * 10 ** (df["logM200_msun"].values - 15.0))
    y = np.log10(df["veldisp_halo_kms"].values)
    valid = np.isfinite(x) & np.isfinite(y)
    alpha, log_sigma15 = np.polyfit(x[valid], y[valid], 1)
    sigma_15 = 10 ** log_sigma15
    return sigma_15, alpha


# --------------------------------------------------------------------------
# mock lightcone construction
# --------------------------------------------------------------------------

def load_tng(tng_dir: Path) -> pd.DataFrame:
    df100 = pd.read_parquet(tng_dir / "tng100-1_local_catalog_enriched.parquet")
    df50 = pd.read_parquet(tng_dir / "tng50-1_local_catalog_enriched.parquet")
    df = pd.concat([df100, df50], ignore_index=True)
    df["halo_key"] = df["snapshot"].astype(str) + "_" + df["halo_id"].astype(str)
    return df


def kpc_per_arcsec(z: float) -> float:
    return cosmo.kpc_proper_per_arcmin(z).value / 60.0


def build_mock_group(true_members: pd.DataFrame, field_pool: pd.DataFrame, z_group: float,
                      completeness_frac: float, z_err_kms: float, n_field_per_true_member: float,
                      rng: np.random.Generator) -> pd.DataFrame:
    n_true = len(true_members)
    n_keep = max(1, int(round(n_true * completeness_frac)))
    kept = true_members.sample(n=min(n_keep, n_true), random_state=rng.integers(1e9))

    v_true = kept["vel_z_kms"].values + rng.normal(0, z_err_kms, size=len(kept))
    kpc_per_as = kpc_per_arcsec(z_group)
    ra0, dec0 = kept["RA_mock"].mean(), kept["DEC_mock"].mean()
    sep_kpc_true = np.hypot((kept["RA_mock"].values - ra0) * 3600 * kpc_per_as,
                             (kept["DEC_mock"].values - dec0) * 3600 * kpc_per_as)

    true_rows = pd.DataFrame(dict(dv_kms=v_true, sep_kpc=sep_kpc_true, is_true_member=True))

    n_field = rng.poisson(n_field_per_true_member * n_true)
    if n_field > 0 and len(field_pool) > 0:
        field = field_pool.sample(n=n_field, replace=True, random_state=rng.integers(1e9))
        sep_kpc_field = rng.uniform(0, SEARCH_RADIUS_KPC, size=n_field)
        dv_field = rng.uniform(-VELOCITY_CEIL_KMS, VELOCITY_CEIL_KMS, size=n_field)
        field_rows = pd.DataFrame(dict(dv_kms=dv_field, sep_kpc=sep_kpc_field, is_true_member=False))
        candidates = pd.concat([true_rows, field_rows], ignore_index=True)
    else:
        candidates = true_rows

    return candidates[candidates["sep_kpc"] <= SEARCH_RADIUS_KPC].copy()


def select_members(candidates: pd.DataFrame) -> pd.DataFrame:
    """Reproduces v1.0's P_v probabilistic velocity gate (same functional
    form/defaults as determine_membership_dztier.py's apply_membership_cut)."""
    p_v = np.exp(-0.5 * (candidates["dv_kms"].values / SIGMA_PRIOR_KMS) ** 2)
    return candidates[p_v > PROB_THRESHOLD].copy()


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def run_experiment(tng: pd.DataFrame, completeness_fracs: list, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    counts = tng.groupby("halo_key").size()
    rich_keys = counts[counts >= MIN_RICHNESS].index
    groups = tng[tng["halo_key"].isin(rich_keys)]

    sigma15_all, alpha_all = fit_selfcal(groups)
    low_mass = groups[groups["logM200_msun"] < 13.5]
    sigma15_grp, alpha_grp = fit_selfcal(low_mass) if len(low_mass) >= 20 else (sigma15_all, alpha_all)
    print(f"tng_selfcal_all:    sigma_15={sigma15_all:.1f} km/s, alpha={alpha_all:.4f} (n={len(groups)} halos' members)")
    print(f"tng_selfcal_groups: sigma_15={sigma15_grp:.1f} km/s, alpha={alpha_grp:.4f} "
          f"(n={len(low_mass)} halos' members, log M200 < 13.5)")

    calibrations = dict(
        munari2013_galaxies=_munari2013,
        evrard2008_dm=_evrard2008,
        tng_selfcal_all=_make_selfcal(sigma15_all, alpha_all),
        tng_selfcal_groups=_make_selfcal(sigma15_grp, alpha_grp),
    )

    rows = []
    for halo_key, true_members in groups.groupby("halo_key"):
        z_group = true_members["snapshot_redshift"].iloc[0]
        logM200_true = true_members["logM200_msun"].iloc[0]
        n_true = len(true_members)
        field_pool = tng[(tng["snapshot"] == true_members["snapshot"].iloc[0]) &
                          (tng["halo_key"] != halo_key)]

        for completeness in completeness_fracs:
            candidates = build_mock_group(true_members, field_pool, z_group, completeness,
                                           Z_ERR_KMS_DEFAULT, 0.5, rng)
            selected = select_members(candidates)
            n_specz_final = len(selected)
            if n_specz_final < MIN_RICHNESS:
                continue

            sigma_v_obs = gapper_sigma_v(selected["dv_kms"].values)
            purity = float(selected["is_true_member"].mean()) if n_specz_final > 0 else np.nan
            recall = float(selected["is_true_member"].sum() / n_true)

            row = dict(halo_key=halo_key, z_group=z_group, logM200_true=logM200_true,
                       n_true_members=n_true, completeness_frac=completeness,
                       n_specz_final=n_specz_final, sigma_v_obs_kms=sigma_v_obs,
                       purity=purity, recall=recall)
            for name, fn in calibrations.items():
                m_pred = fn(sigma_v_obs, z_group)
                row[f"logM200_pred_{name}"] = np.log10(m_pred) if m_pred > 0 else np.nan
            rows.append(row)

    out = pd.DataFrame(rows)
    print(f"\nBuilt {len(out)} mock (halo x completeness) observations from {len(rich_keys)} eligible TNG halos")
    return out, calibrations


def summarize_bias_scatter(df: pd.DataFrame, calibrations: dict, bin_col: str, bins: list) -> pd.DataFrame:
    rows = []
    for name in calibrations:
        pred_col = f"logM200_pred_{name}"
        resid = df[pred_col] - df["logM200_true"]
        for lo, hi in bins:
            m = df[(df[bin_col] >= lo) & (df[bin_col] < hi)]
            r = m[pred_col] - m["logM200_true"]
            r = r.replace([np.inf, -np.inf], np.nan).dropna()
            if len(r) < 3:
                rows.append(dict(calibration=name, bin_col=bin_col, bin=f"{lo}-{hi}", n=len(r),
                                  bias_dex=np.nan, scatter_dex=np.nan, catastrophic_frac=np.nan))
                continue
            rows.append(dict(calibration=name, bin_col=bin_col, bin=f"{lo}-{hi}", n=len(r),
                              bias_dex=float(np.median(r)), scatter_dex=float(np.std(r)),
                              catastrophic_frac=float((np.abs(r) > 0.5).mean())))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tng-dir", type=str, default="/Volumes/extHD/tng_local_catalog")
    args = parser.parse_args()

    tng = load_tng(Path(args.tng_dir))
    completeness_fracs = [0.3, 0.5, 0.7]
    df, calibrations = run_experiment(tng, completeness_fracs)
    df.to_csv(OUTPUT_DIR / "phase5_mock_recovery_pergroup.csv", index=False)

    print("\n=== overall bias/scatter/catastrophic-failure-rate per calibration ===")
    overall_rows = []
    for name in calibrations:
        pred_col = f"logM200_pred_{name}"
        r = (df[pred_col] - df["logM200_true"]).replace([np.inf, -np.inf], np.nan).dropna()
        overall_rows.append(dict(calibration=name, n=len(r), bias_dex=float(np.median(r)),
                                  scatter_dex=float(np.std(r)), catastrophic_frac=float((np.abs(r) > 0.5).mean())))
        r_pearson, p_pearson = pearsonr(df[pred_col].replace([np.inf, -np.inf], np.nan).dropna(),
                                         df.loc[df[pred_col].replace([np.inf, -np.inf], np.nan).dropna().index, "logM200_true"])
        overall_rows[-1]["pearson_r_pred_vs_true"] = r_pearson
        overall_rows[-1]["pearson_p"] = p_pearson
    overall = pd.DataFrame(overall_rows)
    print(overall.to_string(index=False))
    overall.to_csv(OUTPUT_DIR / "phase5_overall_bias_scatter.csv", index=False)

    print("\n=== by true halo mass (log M200) ===")
    mass_bins = [(11, 12.5), (12.5, 13.0), (13.0, 13.5), (13.5, 15.5)]
    by_mass = summarize_bias_scatter(df, calibrations, "logM200_true", mass_bins)
    print(by_mass.to_string(index=False))
    by_mass.to_csv(OUTPUT_DIR / "phase5_by_mass.csv", index=False)

    print("\n=== by richness (n_true_members) ===")
    rich_bins = [(5, 7), (8, 9), (10, 14), (15, 999)]
    by_rich = summarize_bias_scatter(df, calibrations, "n_true_members", rich_bins)
    print(by_rich.to_string(index=False))
    by_rich.to_csv(OUTPUT_DIR / "phase5_by_richness.csv", index=False)

    print("\n=== by redshift ===")
    z_bins = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.7), (1.7, 4.0)]
    by_z = summarize_bias_scatter(df, calibrations, "z_group", z_bins)
    print(by_z.to_string(index=False))
    by_z.to_csv(OUTPUT_DIR / "phase5_by_redshift.csv", index=False)

    print("\n=== by spectroscopic completeness ===")
    comp_bins = [(0.25, 0.35), (0.45, 0.55), (0.65, 0.75)]
    by_comp = summarize_bias_scatter(df, calibrations, "completeness_frac", comp_bins)
    print(by_comp.to_string(index=False))
    by_comp.to_csv(OUTPUT_DIR / "phase5_by_completeness.csv", index=False)

    print(f"\nAll outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
