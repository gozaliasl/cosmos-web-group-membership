#!/usr/bin/env python
"""
membership-v2 Phase 6: definitive benchmark + error-budget decomposition.

Builds on Phase 5's mock lightcone machinery (same TNG ground truth, same
observational-effects model, same production membership gate, same 4
calibrations). Two products:

  A. Full configuration benchmark: for every combination of {aperture mode,
     sigma_v estimator, calibration} x {completeness fraction}, report mass
     bias, scatter, catastrophic-failure rate, completeness (recall),
     purity, sigma_v precision, robustness (cross-seed scatter), and
     computational cost.

  B. Error-budget decomposition: a SEQUENTIAL ABLATION (not an orthogonal
     variance decomposition -- residuals are correlated across stages since
     the same halos are reused throughout, so this reports MARGINAL scatter
     added at each stage, clearly labeled as such, not a strict additive
     ANOVA-style split):

       S0 intrinsic       : calibration applied directly to TNG's own
                             veldisp_halo_kms (ground-truth, full-halo
                             dispersion, no member-sampling at all) -->
                             isolates factor 7, intrinsic halo dynamics /
                             calibration's own intrinsic scatter (these two
                             cannot be fully separated from each other with
                             this test alone -- see report).
       S1 + estimator      : switch to gapper(true members, full richness,
                             no incompleteness/z-error/aperture-cut/
                             contamination) --> marginal addition = factor 3,
                             sigma_v estimator sampling noise.
       S2 + incompleteness : subsample to completeness_frac=0.5 --> marginal
                             addition = factor 5.
       S3 + aperture       : apply the 750 kpc spatial cut (vs. an
                             effectively unbounded 5000 kpc cut) --> marginal
                             addition = factor 2.
       S4 + membership gate: apply the production P_v probabilistic velocity
                             gate (still no interlopers present) --> marginal
                             addition = factor 1 (membership selection logic
                             itself, distinct from aperture).
       S5 + contamination  : add field interlopers --> marginal addition =
                             factor 6, projection/interlopers.

     Factor 4 (calibration choice) is quantified separately as the
     spread in bias across the 4 calibrations at the final (S5) stage --
     it is a normalization choice, not a variance source in the same sense
     as S0-S5, so it is reported alongside rather than inside the sequential
     ladder.

Usage:
    python membership_v2_phase6_benchmark.py --tng-dir /Volumes/extHD/tng_local_catalog
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).parent))
import membership_v2_phase5_mock_lightcone as p5  # noqa: E402

OUTPUT_DIR = BASE_DIR / "outputs" / "results" / "membership_v2_phase6"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def resid_stats(pred, true):
    r = (pred - true).replace([np.inf, -np.inf], np.nan).dropna() if hasattr(pred, "replace") else pd.Series(pred - true).replace([np.inf, -np.inf], np.nan).dropna()
    if len(r) < 3:
        return dict(n=len(r), bias_dex=np.nan, scatter_dex=np.nan, catastrophic_frac=np.nan)
    return dict(n=len(r), bias_dex=float(np.median(r)), scatter_dex=float(np.std(r)),
                catastrophic_frac=float((np.abs(r) > 0.5).mean()))


# --------------------------------------------------------------------------
# A. full configuration benchmark
# --------------------------------------------------------------------------

def full_benchmark(tng, groups, calibrations, estimators, apertures, completeness_fracs, n_seeds=3):
    rows = []
    for aperture_kpc in apertures:
        for est_name, est_fn in estimators.items():
            for completeness in completeness_fracs:
                seed_biases = {name: [] for name in calibrations}
                seed_scatters = {name: [] for name in calibrations}
                purities, recalls, sigma_precisions, n_final_list = [], [], [], []
                t0 = time.time()
                for seed in range(n_seeds):
                    rng = np.random.default_rng(seed)
                    rows_this_seed = []
                    for halo_key, true_members in groups.groupby("halo_key"):
                        z_group = true_members["snapshot_redshift"].iloc[0]
                        logM200_true = true_members["logM200_msun"].iloc[0]
                        veldisp_true = true_members["veldisp_halo_kms"].iloc[0]
                        n_true = len(true_members)
                        field_pool = tng[(tng["snapshot"] == true_members["snapshot"].iloc[0]) &
                                          (tng["halo_key"] != halo_key)]

                        cand = _build_group_custom(true_members, field_pool, z_group, completeness,
                                                    p5.Z_ERR_KMS_DEFAULT, 0.5, aperture_kpc, rng)
                        sel = p5.select_members(cand)
                        if len(sel) < p5.MIN_RICHNESS:
                            continue
                        sigma_v_obs = est_fn(sel["dv_kms"].values)
                        if not np.isfinite(sigma_v_obs) or sigma_v_obs <= 0:
                            continue
                        purity = float(sel["is_true_member"].mean())
                        recall = float(sel["is_true_member"].sum() / n_true)
                        sigma_prec = (sigma_v_obs - veldisp_true) / veldisp_true if np.isfinite(veldisp_true) and veldisp_true > 0 else np.nan
                        row = dict(logM200_true=logM200_true, purity=purity, recall=recall, sigma_prec=sigma_prec, n_final=len(sel))
                        for name, fn in calibrations.items():
                            mp = fn(sigma_v_obs, z_group)
                            row[f"pred_{name}"] = np.log10(mp) if mp > 0 else np.nan
                        rows_this_seed.append(row)
                    seed_df = pd.DataFrame(rows_this_seed)
                    if len(seed_df) == 0:
                        continue
                    purities.append(seed_df["purity"].mean())
                    recalls.append(seed_df["recall"].mean())
                    sigma_precisions.append(seed_df["sigma_prec"].std())
                    n_final_list.append(seed_df["n_final"].mean())
                    for name in calibrations:
                        r = (seed_df[f"pred_{name}"] - seed_df["logM200_true"]).replace([np.inf, -np.inf], np.nan).dropna()
                        if len(r) >= 3:
                            seed_biases[name].append(float(np.median(r)))
                            seed_scatters[name].append(float(np.std(r)))
                elapsed = time.time() - t0

                for name in calibrations:
                    if not seed_biases[name]:
                        continue
                    rows.append(dict(
                        aperture_kpc=aperture_kpc, estimator=est_name, completeness_frac=completeness,
                        calibration=name,
                        bias_dex=float(np.mean(seed_biases[name])),
                        scatter_dex=float(np.mean(seed_scatters[name])),
                        robustness_bias_std_across_seeds=float(np.std(seed_biases[name])),
                        robustness_scatter_std_across_seeds=float(np.std(seed_scatters[name])),
                        purity=float(np.mean(purities)) if purities else np.nan,
                        completeness_recall=float(np.mean(recalls)) if recalls else np.nan,
                        sigma_v_precision_frac_scatter=float(np.mean(sigma_precisions)) if sigma_precisions else np.nan,
                        n_final_mean=float(np.mean(n_final_list)) if n_final_list else np.nan,
                        runtime_sec_per_seed=elapsed / n_seeds,
                    ))
    return pd.DataFrame(rows)


def _build_group_custom(true_members, field_pool, z_group, completeness_frac, z_err_kms,
                         n_field_per_true_member, aperture_kpc, rng):
    """Same as p5.build_mock_group but with a configurable aperture (p5's
    version hardcodes SEARCH_RADIUS_KPC)."""
    n_true = len(true_members)
    n_keep = max(1, int(round(n_true * completeness_frac)))
    kept = true_members.sample(n=min(n_keep, n_true), random_state=rng.integers(1e9))
    v_true = kept["vel_z_kms"].values + rng.normal(0, z_err_kms, size=len(kept))
    kpc_per_as = p5.kpc_per_arcsec(z_group)
    ra0, dec0 = kept["RA_mock"].mean(), kept["DEC_mock"].mean()
    sep_kpc_true = np.hypot((kept["RA_mock"].values - ra0) * 3600 * kpc_per_as,
                             (kept["DEC_mock"].values - dec0) * 3600 * kpc_per_as)
    true_rows = pd.DataFrame(dict(dv_kms=v_true, sep_kpc=sep_kpc_true, is_true_member=True))

    n_field = rng.poisson(n_field_per_true_member * n_true) if n_field_per_true_member > 0 else 0
    if n_field > 0 and len(field_pool) > 0:
        field = field_pool.sample(n=n_field, replace=True, random_state=rng.integers(1e9))
        sep_kpc_field = rng.uniform(0, aperture_kpc, size=n_field)
        dv_field = rng.uniform(-p5.VELOCITY_CEIL_KMS, p5.VELOCITY_CEIL_KMS, size=n_field)
        field_rows = pd.DataFrame(dict(dv_kms=dv_field, sep_kpc=sep_kpc_field, is_true_member=False))
        candidates = pd.concat([true_rows, field_rows], ignore_index=True)
    else:
        candidates = true_rows
    return candidates[candidates["sep_kpc"] <= aperture_kpc].copy()


# --------------------------------------------------------------------------
# B. sequential ablation / error-budget decomposition
# --------------------------------------------------------------------------

def ablation_ladder(tng, groups, calibrations, seed=42, completeness_frac=0.5):
    """
    Matched-sample ablation: the incompleteness retention mask is drawn ONCE
    per halo and reused unchanged at every downstream stage (S2-S6), and the
    per-stage results are additionally restricted to the intersection of
    halos that pass n>=5 at EVERY stage (`core_keys`), so the marginal
    variance added at each step is a genuine like-for-like comparison on the
    same halos/members, not a comparison across differently-resampled
    survivor sets.
    """
    rng = np.random.default_rng(seed)

    cal_name = "munari2013_galaxies"  # fixed reference calibration for the ladder itself
    cal_fn = calibrations[cal_name]

    # draw the completeness retention mask once per halo
    retained_idx = {}
    for halo_key, tm in groups.groupby("halo_key"):
        n_true = len(tm)
        n_keep = max(1, int(round(n_true * completeness_frac)))
        retained_idx[halo_key] = tm.sample(n=min(n_keep, n_true), random_state=rng.integers(1e9)).index

    def _stage_df(use_incompleteness, z_err_kms, aperture_kpc, use_gate, use_contamination):
        rows = []
        for halo_key, tm in groups.groupby("halo_key"):
            z, logM = tm["snapshot_redshift"].iloc[0], tm["logM200_msun"].iloc[0]
            base = tm.loc[retained_idx[halo_key]] if use_incompleteness else tm
            v = base["vel_z_kms"].values + rng.normal(0, z_err_kms, size=len(base)) if z_err_kms > 0 else base["vel_z_kms"].values
            kpc_per_as = p5.kpc_per_arcsec(z)
            ra0, dec0 = base["RA_mock"].mean(), base["DEC_mock"].mean()
            sep = np.hypot((base["RA_mock"].values - ra0) * 3600 * kpc_per_as,
                            (base["DEC_mock"].values - dec0) * 3600 * kpc_per_as)
            cand = pd.DataFrame(dict(dv_kms=v, sep_kpc=sep, is_true_member=True))
            cand = cand[cand["sep_kpc"] <= aperture_kpc]
            if use_contamination:
                field_pool = tng[(tng["snapshot"] == tm["snapshot"].iloc[0]) & (tng["halo_key"] != halo_key)]
                n_field = rng.poisson(0.5 * len(tm))
                if n_field > 0 and len(field_pool) > 0:
                    field = field_pool.sample(n=n_field, replace=True, random_state=rng.integers(1e9))
                    field_rows = pd.DataFrame(dict(
                        dv_kms=rng.uniform(-p5.VELOCITY_CEIL_KMS, p5.VELOCITY_CEIL_KMS, size=n_field),
                        sep_kpc=rng.uniform(0, aperture_kpc, size=n_field), is_true_member=False))
                    cand = pd.concat([cand, field_rows[field_rows["sep_kpc"] <= aperture_kpc]], ignore_index=True)
            sel = p5.select_members(cand) if use_gate else cand
            if len(sel) < p5.MIN_RICHNESS:
                continue
            sv = p5.gapper_sigma_v(sel["dv_kms"].values)
            if not np.isfinite(sv) or sv <= 0:
                continue
            mp = cal_fn(sv, z)
            rows.append(dict(halo_key=halo_key, logM200_true=logM, pred=np.log10(mp) if mp > 0 else np.nan))
        return pd.DataFrame(rows)

    stages = {}
    # S0: intrinsic -- calibration applied directly to TNG's own veldisp_halo_kms (no sampling at all)
    rows = []
    for halo_key, tm in groups.groupby("halo_key"):
        z, logM = tm["snapshot_redshift"].iloc[0], tm["logM200_msun"].iloc[0]
        vd = tm["veldisp_halo_kms"].iloc[0]
        if not np.isfinite(vd) or vd <= 0:
            continue
        mp = cal_fn(vd, z)
        rows.append(dict(halo_key=halo_key, logM200_true=logM, pred=np.log10(mp) if mp > 0 else np.nan))
    stages["S0_intrinsic"] = pd.DataFrame(rows)

    stages["S1_plus_estimator"] = _stage_df(False, 0.0, 5000.0, False, False)
    stages["S2_plus_incompleteness"] = _stage_df(True, 0.0, 5000.0, False, False)
    stages["S3_plus_zerror"] = _stage_df(True, p5.Z_ERR_KMS_DEFAULT, 5000.0, False, False)
    stages["S4_plus_aperture"] = _stage_df(True, p5.Z_ERR_KMS_DEFAULT, p5.SEARCH_RADIUS_KPC, False, False)
    stages["S5_plus_membership_gate"] = _stage_df(True, p5.Z_ERR_KMS_DEFAULT, p5.SEARCH_RADIUS_KPC, True, False)
    stages["S6_plus_contamination"] = _stage_df(True, p5.Z_ERR_KMS_DEFAULT, p5.SEARCH_RADIUS_KPC, True, True)

    # matched core sample: halos present in every stage (S0 has a different, superset-like
    # eligibility criterion -- veldisp_halo_kms coverage -- so the "core" is defined over S1-S6)
    core_keys = None
    for name, df in stages.items():
        if name == "S0_intrinsic":
            continue
        keys = set(df["halo_key"])
        core_keys = keys if core_keys is None else (core_keys & keys)

    ladder_rows = []
    prev_var = None
    for name, df in stages.items():
        df_matched = df[df["halo_key"].isin(core_keys)] if name != "S0_intrinsic" else df[df["halo_key"].isin(core_keys)]
        stats = resid_stats(df_matched["pred"], df_matched["logM200_true"])
        var = stats["scatter_dex"] ** 2 if np.isfinite(stats["scatter_dex"]) else np.nan
        marginal_var = (var - prev_var) if (prev_var is not None and np.isfinite(var) and np.isfinite(prev_var)) else np.nan
        ladder_rows.append(dict(stage=name, n_full_stage=len(df), n_matched_core=len(df_matched),
                                 **stats, variance_dex2=var, marginal_variance_added_dex2=marginal_var))
        prev_var = var
    return pd.DataFrame(ladder_rows), len(core_keys)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tng-dir", type=str, default="/Volumes/extHD/tng_local_catalog")
    args = parser.parse_args()

    tng = p5.load_tng(Path(args.tng_dir))
    counts = tng.groupby("halo_key").size()
    rich_keys = counts[counts >= p5.MIN_RICHNESS].index
    groups = tng[tng["halo_key"].isin(rich_keys)]

    sigma15_all, alpha_all = p5.fit_selfcal(groups)
    low_mass = groups[groups["logM200_msun"] < 13.5]
    sigma15_grp, alpha_grp = p5.fit_selfcal(low_mass)
    calibrations = dict(
        munari2013_galaxies=p5._munari2013,
        evrard2008_dm=p5._evrard2008,
        tng_selfcal_all=p5._make_selfcal(sigma15_all, alpha_all),
        tng_selfcal_groups=p5._make_selfcal(sigma15_grp, alpha_grp),
    )
    estimators = dict(gapper=p5.gapper_sigma_v)
    try:
        from astropy.stats import biweight_scale
        estimators["biweight"] = lambda v: float(biweight_scale(v)) if len(v) >= 2 else np.nan
    except ImportError:
        pass
    estimators["std"] = lambda v: float(np.std(v, ddof=1)) if len(v) >= 2 else np.nan

    print("=== B. Sequential ablation / error-budget decomposition ===")
    ladder, n_core = ablation_ladder(tng, groups, calibrations)
    print(f"(matched core sample: {n_core} halos present at every stage S1-S6)")
    print(ladder.to_string(index=False))
    ladder.to_csv(OUTPUT_DIR / "phase6_error_budget_ladder.csv", index=False)

    total_var = ladder["variance_dex2"].iloc[-1]
    budget = ladder.iloc[1:].copy()
    budget["pct_of_total_variance"] = 100 * budget["marginal_variance_added_dex2"] / total_var
    print("\n=== Error budget as % of total variance (S6 endpoint) ===")
    print(budget[["stage", "marginal_variance_added_dex2", "pct_of_total_variance"]].to_string(index=False))
    budget.to_csv(OUTPUT_DIR / "phase6_error_budget_pct.csv", index=False)

    # calibration-choice spread at the final (S6-equivalent) stage
    print("\n=== Calibration-choice spread (factor 4), at full-mock (S6-equivalent) conditions, 50% completeness ===")
    rng = np.random.default_rng(42)
    rows = []
    for halo_key, tm in groups.groupby("halo_key"):
        z, logM = tm["snapshot_redshift"].iloc[0], tm["logM200_msun"].iloc[0]
        field_pool = tng[(tng["snapshot"] == tm["snapshot"].iloc[0]) & (tng["halo_key"] != halo_key)]
        cand = _build_group_custom(tm, field_pool, z, 0.5, p5.Z_ERR_KMS_DEFAULT, 0.5, p5.SEARCH_RADIUS_KPC, rng)
        sel = p5.select_members(cand)
        if len(sel) < p5.MIN_RICHNESS:
            continue
        sv = p5.gapper_sigma_v(sel["dv_kms"].values)
        row = dict(logM200_true=logM)
        for name, fn in calibrations.items():
            mp = fn(sv, z)
            row[f"pred_{name}"] = np.log10(mp) if mp > 0 else np.nan
        rows.append(row)
    cal_df = pd.DataFrame(rows)
    cal_biases = {name: float(np.median((cal_df[f"pred_{name}"] - cal_df["logM200_true"]).replace([np.inf, -np.inf], np.nan).dropna()))
                  for name in calibrations}
    cal_spread = max(cal_biases.values()) - min(cal_biases.values())
    print(pd.Series(cal_biases))
    print(f"calibration-choice bias spread: {cal_spread:.3f} dex (this is factor 4's contribution -- a normalization spread, not a variance-ladder term)")
    pd.Series(cal_biases).to_csv(OUTPUT_DIR / "phase6_calibration_spread.csv")

    print("\n=== A. Full configuration benchmark (this is slower -- aperture x estimator x completeness x calibration) ===")
    apertures = [750.0]  # kept to the production default; aperture's effect is isolated in the ladder (S3->S4) instead
    completeness_fracs = [0.3, 0.5, 0.7]
    bench = full_benchmark(tng, groups, calibrations, estimators, apertures, completeness_fracs, n_seeds=3)
    print(bench.to_string(index=False))
    bench.to_csv(OUTPUT_DIR / "phase6_full_benchmark.csv", index=False)

    print(f"\nAll outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
