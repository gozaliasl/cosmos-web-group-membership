#!/usr/bin/env python
"""
membership-v2: configurable membership aperture + multi-estimator dynamics +
membership-quality metrics.

Phase 2 of the membership-v2 methodology program (see
cosmos-web-xray-igm/docs/membership_aperture_methodology.md and
REFEREE_REVIEW.md for the audit this builds on). This module does NOT modify
or replace v1.0 (iterative_membership_bgg.py, select_bgg.py,
determine_membership_dztier.py, or cosmos-web-xray-igm/scripts/compute_group_r200.py)
-- v1.0 remains the untouched production pipeline. This is a new,
side-by-side, fully configurable module operating on v1.0's already-selected
member pool (same input files), so results are directly comparable.

Everything here is driven by a MembershipV2Config; the CLI can run any number
of named presets in one invocation and writes one output file per preset.

Aperture modes
---------------
fixed            : use the original search aperture as-is (no trimming) --
                   reproduces v1.0 SPECZ/XRAY+SPECZ member pools exactly.
r200x_trim       : trim the member pool to sep_kpc <= R200,X (X-ray-derived
                   R200; requires an r200_xray_kpc column supplied via
                   --xray-r200-file). Groups without an X-ray R200 fall back
                   to 'fixed'.
r200dyn_adaptive : iterate the spatial cut to the *dynamically estimated*
                   R200 to convergence (this is v1.0's existing
                   dynamical_estimate() logic from compute_group_r200.py,
                   reimplemented here so it can be run with any of the four
                   sigma_v estimators below, not just gapper).

sigma_v estimators (all four computed and reported side-by-side every run)
---------------------------------------------------------------------------
gapper   : Beers, Flynn & Gebhardt (1990) gapper estimator (v1.0's choice).
biweight : Beers et al. (1990) biweight midvariance (astropy.stats.biweight_scale).
std      : ordinary sample standard deviation (ddof=1).
robust   : 1.4826 * median absolute deviation (MAD), a simple order-statistic
           robust estimator distinct from both gapper and biweight.

The aperture-mode iteration itself always uses the estimator named in
`sigma_v_estimator_for_iteration` (default gapper, matching v1.0); the other
three are computed once on the final member set for comparison, not used to
drive convergence, so the four columns are on an equal footing at the same
final aperture.

Usage:
    python membership_v2_dynamics.py --catalog both --config fixed750,r200x_trim,r200dyn_adaptive
"""

import argparse
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.stats import biweight_scale
from astropy.table import Table

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent  # cosmos-web-group-membership/
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
MEMBERSHIP_DIR = BASE_DIR / "outputs" / "results" / "membership_dztier"
OUTPUT_DIR = BASE_DIR / "outputs" / "results" / "membership_v2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Reuse v1.0's validated sigma_v-M200 relation and R200<->M200 conversion
# rather than re-deriving them -- this module changes *membership/aperture/
# estimator* methodology only, per the Project B scope (dynamical mass
# calibration itself, D7, is explicitly out of scope for this branch).
sys.path.insert(0, str(XRAY_REPO / "scripts"))
sys.path.insert(0, str(XRAY_REPO / "src"))
from compute_group_r200 import (  # noqa: E402
    r200_m200_from_sigma_v as _v1_r200_m200_from_sigma_v,
    gapper_sigma_v as _v1_gapper_sigma_v,
)

C_KMS = 299792.458
MIN_SPECZ_FOR_DYNAMICS = 5
DYNAMICS_MAX_ITER = 10
DYNAMICS_CONVERGE_TOL_FRAC = 0.10


# --------------------------------------------------------------------------
# sigma_v estimators
# --------------------------------------------------------------------------

def gapper_sigma_v(v: np.ndarray) -> float:
    return _v1_gapper_sigma_v(v)


def biweight_sigma_v(v: np.ndarray) -> float:
    if len(v) < 2:
        return np.nan
    try:
        return float(biweight_scale(v))
    except Exception:
        return np.nan


def std_sigma_v(v: np.ndarray) -> float:
    if len(v) < 2:
        return np.nan
    return float(np.std(v, ddof=1))


def robust_mad_sigma_v(v: np.ndarray) -> float:
    if len(v) < 2:
        return np.nan
    med = np.median(v)
    return float(1.4826 * np.median(np.abs(v - med)))


ESTIMATORS = {
    "gapper": gapper_sigma_v,
    "biweight": biweight_sigma_v,
    "std": std_sigma_v,
    "robust_mad": robust_mad_sigma_v,
}


def all_estimators(velocities_kms: np.ndarray) -> dict:
    return {name: fn(velocities_kms) for name, fn in ESTIMATORS.items()}


def r200_m200_from_sigma_v(sigma_v_kms, sigma_v_err_kms, z):
    """Thin wrapper on v1.0's validated (interim, Munari+2013) relation --
    not re-derived here; see module docstring."""
    return _v1_r200_m200_from_sigma_v(sigma_v_kms, sigma_v_err_kms, z)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

@dataclass
class MembershipV2Config:
    name: str
    aperture_mode: str = "fixed"          # fixed | r200x_trim | r200dyn_adaptive
    velocity_window_nsigma: float = None  # if set, clip members with |dv| > N*sigma_v after estimator convergence (None = no clip, matches v1.0)
    sigma_v_estimator_for_iteration: str = "gapper"
    min_specz: int = MIN_SPECZ_FOR_DYNAMICS


PRESETS = {
    "fixed750": MembershipV2Config(name="fixed750", aperture_mode="fixed"),
    "r200x_trim": MembershipV2Config(name="r200x_trim", aperture_mode="r200x_trim"),
    "r200dyn_adaptive": MembershipV2Config(name="r200dyn_adaptive", aperture_mode="r200dyn_adaptive"),
    "r200dyn_adaptive_biweight": MembershipV2Config(
        name="r200dyn_adaptive_biweight", aperture_mode="r200dyn_adaptive",
        sigma_v_estimator_for_iteration="biweight"),
    "r200dyn_adaptive_vclip3": MembershipV2Config(
        name="r200dyn_adaptive_vclip3", aperture_mode="r200dyn_adaptive",
        velocity_window_nsigma=3.0),
}


# --------------------------------------------------------------------------
# quality metrics
# --------------------------------------------------------------------------

def compute_quality_metrics(specz_in_search: pd.DataFrame, specz_final: pd.DataFrame,
                             sigma_v_final: float, R200_final: float, z_group: float) -> dict:
    """
    All metrics defined relative to the *search-aperture* candidate pool
    (specz_in_search) vs. the *final, converged* member set (specz_final),
    so they quantify what the aperture/estimator choice actually did.
    """
    n_search = len(specz_in_search)
    n_final = len(specz_final)

    if n_search == 0:
        return dict(projected_member_fraction=np.nan, members_outside_r200_frac=np.nan,
                     phase_space_outlier_frac=np.nan, velocity_tail_frac=np.nan,
                     contamination_score=np.nan, dynamical_reliability_score=np.nan,
                     membership_confidence=np.nan)

    projected_member_fraction = n_final / n_search

    if n_final == 0 or not np.isfinite(R200_final) or not np.isfinite(sigma_v_final) or sigma_v_final <= 0:
        return dict(projected_member_fraction=projected_member_fraction,
                     members_outside_r200_frac=np.nan, phase_space_outlier_frac=np.nan,
                     velocity_tail_frac=np.nan, contamination_score=np.nan,
                     dynamical_reliability_score=0.0, membership_confidence=0.0)

    members_outside_r200_frac = float((specz_final["sep_kpc"] > R200_final).mean())

    z_med = np.median(specz_final["zfin"].values)
    dv = C_KMS * (specz_final["zfin"].values - z_med) / (1 + z_med)
    r_norm = specz_final["sep_kpc"].values / R200_final
    v_norm = dv / sigma_v_final

    # phase-space outlier: escape-velocity-like envelope |v_norm| > sqrt(2)*(1 - r_norm/3)
    # for r_norm<=3, else always an outlier beyond r_norm=3 (caustic-motivated,
    # matches the corrected core-normalized phase-space diagnostic from the v1.0 audit)
    envelope = np.where(r_norm <= 3.0, np.sqrt(2.0) * (1.0 - r_norm / 3.0), 0.0)
    phase_space_outlier_frac = float((np.abs(v_norm) > np.maximum(envelope, 0.3)).mean())

    velocity_tail_frac = float((np.abs(v_norm) > 2.0).mean())

    contamination_score = float(np.clip(
        0.5 * members_outside_r200_frac + 0.3 * phase_space_outlier_frac + 0.2 * velocity_tail_frac,
        0.0, 1.0))

    # dynamical reliability: richness term (saturating at min_specz*3) x (1 - contamination)
    richness_term = min(1.0, n_final / (MIN_SPECZ_FOR_DYNAMICS * 3.0))
    dynamical_reliability_score = float(np.clip(richness_term * (1.0 - contamination_score), 0.0, 1.0))

    membership_confidence = float(np.clip(
        0.5 * projected_member_fraction + 0.5 * dynamical_reliability_score, 0.0, 1.0))

    return dict(projected_member_fraction=projected_member_fraction,
                members_outside_r200_frac=members_outside_r200_frac,
                phase_space_outlier_frac=phase_space_outlier_frac,
                velocity_tail_frac=velocity_tail_frac,
                contamination_score=contamination_score,
                dynamical_reliability_score=dynamical_reliability_score,
                membership_confidence=membership_confidence)


# --------------------------------------------------------------------------
# core per-group computation
# --------------------------------------------------------------------------

def dynamical_estimate_v2(members: pd.DataFrame, z_group: float, radius_used_kpc: float,
                           config: MembershipV2Config, r200_xray_kpc: float = np.nan) -> dict:
    empty = dict(sigma_v_kms=np.nan, sigma_v_err_kms=np.nan, R200_kpc=np.nan, R200_err_kpc=np.nan,
                 M200_Msun=np.nan, M200_err_Msun=np.nan, n_specz=0, n_dynamics_iter=0,
                 aperture_used_kpc=np.nan, sigma_v_gapper_kms=np.nan, sigma_v_biweight_kms=np.nan,
                 sigma_v_std_kms=np.nan, sigma_v_robust_mad_kms=np.nan)
    empty.update({k: np.nan for k in ("projected_member_fraction", "members_outside_r200_frac",
                                       "phase_space_outlier_frac", "velocity_tail_frac",
                                       "contamination_score", "dynamical_reliability_score",
                                       "membership_confidence")})

    specz_pool = members[members["redshift_type"] == "spec-z"]
    n_before = len(specz_pool)
    if n_before < config.min_specz:
        return dict(empty, n_specz=n_before)

    specz_search = specz_pool[specz_pool["sep_kpc"] <= radius_used_kpc] if "sep_kpc" in specz_pool.columns else specz_pool
    if len(specz_search) < config.min_specz:
        return dict(empty, n_specz=len(specz_search))

    est_fn = ESTIMATORS[config.sigma_v_estimator_for_iteration]

    def _compute(subset):
        z_med = np.median(subset["zfin"].values)
        v = C_KMS * (subset["zfin"].values - z_med) / (1 + z_med)
        sigma_v = est_fn(v)
        n = len(subset)
        sigma_v_err = sigma_v / np.sqrt(2 * (n - 1)) if n >= 3 else sigma_v / np.sqrt(max(n, 1))
        R200, R200_err, M200, M200_err = r200_m200_from_sigma_v(sigma_v, sigma_v_err, z_med)
        return sigma_v, sigma_v_err, R200, R200_err, M200, M200_err

    if config.aperture_mode == "fixed":
        specz_final = specz_search
        sigma_v, sigma_v_err, R200, R200_err, M200, M200_err = _compute(specz_final)
        aperture_final = radius_used_kpc
        n_iter = 1

    elif config.aperture_mode == "r200x_trim":
        if not np.isfinite(r200_xray_kpc) or r200_xray_kpc <= 0:
            specz_final = specz_search
            aperture_final = radius_used_kpc
        else:
            trimmed = specz_search[specz_search["sep_kpc"] <= r200_xray_kpc]
            if len(trimmed) < config.min_specz:
                specz_final = specz_search
                aperture_final = radius_used_kpc
            else:
                specz_final = trimmed
                aperture_final = r200_xray_kpc
        sigma_v, sigma_v_err, R200, R200_err, M200, M200_err = _compute(specz_final)
        n_iter = 1

    elif config.aperture_mode == "r200dyn_adaptive":
        specz_final = specz_search
        sigma_v, sigma_v_err, R200, R200_err, M200, M200_err = _compute(specz_final)
        current_aperture = radius_used_kpc
        n_iter = 1
        for _ in range(DYNAMICS_MAX_ITER - 1):
            if not np.isfinite(R200) or current_aperture <= 0:
                break
            R200_capped = min(R200, radius_used_kpc)
            if abs(R200_capped - current_aperture) / current_aperture <= DYNAMICS_CONVERGE_TOL_FRAC:
                break
            current_aperture = R200_capped
            new_pool = specz_search[specz_search["sep_kpc"] <= current_aperture]
            if len(new_pool) < config.min_specz:
                break
            specz_final = new_pool
            sigma_v, sigma_v_err, R200, R200_err, M200, M200_err = _compute(specz_final)
            n_iter += 1
        aperture_final = min(R200, radius_used_kpc) if np.isfinite(R200) else current_aperture

    else:
        raise ValueError(f"Unknown aperture_mode: {config.aperture_mode}")

    # optional velocity-window clip on the final member set (report-only unless configured)
    if config.velocity_window_nsigma is not None and np.isfinite(sigma_v) and sigma_v > 0 and len(specz_final) > config.min_specz:
        z_med = np.median(specz_final["zfin"].values)
        dv = C_KMS * (specz_final["zfin"].values - z_med) / (1 + z_med)
        keep = np.abs(dv) <= config.velocity_window_nsigma * sigma_v
        if keep.sum() >= config.min_specz:
            specz_final = specz_final[keep]
            sigma_v, sigma_v_err, R200, R200_err, M200, M200_err = _compute(specz_final)

    # all four estimators, evaluated on the SAME final member set for direct comparison
    z_med = np.median(specz_final["zfin"].values)
    v_final = C_KMS * (specz_final["zfin"].values - z_med) / (1 + z_med)
    all_est = all_estimators(v_final)

    quality = compute_quality_metrics(specz_search, specz_final, sigma_v, R200, z_group)

    return dict(sigma_v_kms=sigma_v, sigma_v_err_kms=sigma_v_err, R200_kpc=R200, R200_err_kpc=R200_err,
                M200_Msun=M200, M200_err_Msun=M200_err, n_specz=len(specz_final), n_dynamics_iter=n_iter,
                aperture_used_kpc=aperture_final,
                sigma_v_gapper_kms=all_est["gapper"], sigma_v_biweight_kms=all_est["biweight"],
                sigma_v_std_kms=all_est["std"], sigma_v_robust_mad_kms=all_est["robust_mad"],
                **quality)


# --------------------------------------------------------------------------
# catalog-level driver
# --------------------------------------------------------------------------

def latest_wide_aperture(pattern: str, directory: Path) -> Path:
    candidates = [p for p in directory.glob(pattern) if "_r200recomputed_" not in p.name]
    if not candidates:
        raise FileNotFoundError(f"No file matching {pattern} (excl. _r200recomputed_) in {directory}")
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def load_xray_r200(xray_slug: str) -> pd.DataFrame:
    xray_path = XRAY_REPO / "outputs" / "results" / xray_slug / "xray_catalog.fits"
    if not xray_path.exists():
        return pd.DataFrame(columns=["Group_ID", "r200_xray_kpc"])
    xray = Table.read(xray_path).to_pandas()
    xray["Group_ID"] = xray["Group_ID"].astype(str)
    out = xray[["Group_ID", "Is_Detected", "R200_kpc"]].copy()
    out = out[out["Is_Detected"].astype(bool)]
    out = out.rename(columns={"R200_kpc": "r200_xray_kpc"})
    return out[["Group_ID", "r200_xray_kpc"]]


def process_catalog(catalog_name: str, xray_slug: str, membership_prefix: str, config: MembershipV2Config) -> pd.DataFrame:
    members_path = latest_wide_aperture(f"{membership_prefix}members_*.csv", MEMBERSHIP_DIR)
    summary_path = latest_wide_aperture(f"{membership_prefix}summary_*.csv", MEMBERSHIP_DIR)
    members_all = pd.read_csv(members_path)
    summary = pd.read_csv(summary_path)
    members_all["Group_ID"] = members_all["Group_ID"].astype(str)
    summary["Group_ID"] = summary["Group_ID"].astype(str)
    summary_idx = summary.set_index("Group_ID")

    xray_r200 = load_xray_r200(xray_slug).set_index("Group_ID")

    rows = []
    for gid, group_members in members_all.groupby("Group_ID", sort=False):
        if gid not in summary_idx.index:
            continue
        z_group = summary_idx.loc[gid, "z_refined"]
        radius_used = summary_idx.loc[gid, "search_radius_kpc"] if "search_radius_kpc" in summary_idx.columns else 750.0
        r200x = xray_r200.loc[gid, "r200_xray_kpc"] if gid in xray_r200.index else np.nan

        dyn = dynamical_estimate_v2(group_members, z_group, radius_used, config, r200_xray_kpc=r200x)
        rows.append(dict(Group_ID=gid, Catalog=catalog_name, Redshift=z_group,
                          config=config.name, aperture_mode=config.aperture_mode,
                          n_specz_search=len(group_members[(group_members["redshift_type"] == "spec-z") &
                                                            (group_members["sep_kpc"] <= radius_used)]),
                          r200_xray_kpc=r200x, **dyn))

    out = pd.DataFrame(rows)
    n_valid = int(out["n_specz"].ge(config.min_specz).sum()) if len(out) else 0
    print(f"  [{config.name}] {catalog_name}: {len(out)} groups, {n_valid} with >= {config.min_specz} final spec-z members "
          f"(median sigma_v[{config.sigma_v_estimator_for_iteration}]="
          f"{np.nanmedian(out['sigma_v_kms']):.1f} km/s)" if n_valid else f"  [{config.name}] {catalog_name}: 0 valid groups")

    out_path = OUTPUT_DIR / f"membership_v2_{xray_slug}_{config.name}.csv"
    out.to_csv(out_path, index=False)
    return out


def main():
    parser = argparse.ArgumentParser(description="membership-v2: configurable aperture + multi-estimator dynamics + quality metrics")
    parser.add_argument("--catalog", choices=["all", "hcg", "both"], default="both")
    parser.add_argument("--config", type=str, default="fixed750,r200x_trim,r200dyn_adaptive",
                         help=f"Comma-separated preset names from: {list(PRESETS)}")
    args = parser.parse_args()

    configs = [PRESETS[c.strip()] for c in args.config.split(",")]
    catalogs = []
    if args.catalog in ("all", "both"):
        catalogs.append(("CW-All", "cw_all", "iterative_cw_all_"))
    if args.catalog in ("hcg", "both"):
        catalogs.append(("CW-HCG", "cw_hcg", "iterative_cw_hcg_"))

    for config in configs:
        print(f"\n{'='*70}\nConfig: {config.name} (aperture={config.aperture_mode}, "
              f"estimator={config.sigma_v_estimator_for_iteration}, "
              f"vclip={config.velocity_window_nsigma})\n{'='*70}")
        for catalog_name, xray_slug, prefix in catalogs:
            process_catalog(catalog_name, xray_slug, prefix, config)

    print(f"\nAll outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
