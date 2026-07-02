# membership-v2 — Phase 2 & 3 Interim Report

Branch: `membership-v2` (production `main` untouched; v1.0 tagged
`v1.0-audited-baseline` on both this repo and `cosmos-web-xray-igm`).

## What was built (Phase 2)

`scripts/pipeline/membership_v2_dynamics.py` — a new, side-by-side module,
not a modification of v1.0. It reuses v1.0's validated member pools and
sigma_v-M200 relation (Munari et al. 2013, still interim) rather than
re-deriving them, and adds:

- **Configurable aperture**: `fixed` (reproduces v1.0 exactly), `r200x_trim`
  (spatial cut to X-ray-derived R200), `r200dyn_adaptive` (iterate the
  spatial cut to the dynamically-estimated R200 to convergence, v1.0's own
  logic but usable with any estimator).
- **Four sigma_v estimators reported side-by-side** on the same final member
  set: gapper (v1.0's choice), biweight, standard deviation, robust-MAD.
- **Per-group quality metrics**: projected-member fraction, members-outside-
  R200 fraction, phase-space outlier fraction (escape-velocity-envelope
  test), velocity-tail fraction, a composite contamination score, a
  dynamical reliability score, and an overall membership confidence.
- **R200,dyn / M200,dyn with propagated uncertainties** (statistical +
  relation-scatter terms, same propagation approach as v1.0).

**Correctness check**: the `fixed750` config reproduces v1.0's own headline
number exactly (CW-All median sigma_v = 247.6 km/s), confirming the module
is built correctly on v1.0 rather than silently diverging from it.

Run on the full CW-All and CW-HCG catalogs (small-subset-first, then full,
per the agreed workflow): CW-All 1669 groups (112 with >=5 final spec-z
members across all configs), CW-HCG 788 groups (11 with >=5 members).

## Phase 3 findings (real data, CW-All)

### 1. sigma_v vs. independent X-ray M200

| config | n | Pearson r (log-log) | p | Spearman r |
|---|---|---|---|---|
| fixed750 | 63 | -0.006 | 0.965 | -0.135 |
| r200x_trim | 63 | -0.018 | 0.886 | -0.126 |
| r200dyn_adaptive | 63 | **0.156** | 0.222 | -0.069 |
| r200dyn_adaptive_biweight | 50 | 0.093 | 0.522 | -0.088 |
| r200dyn_adaptive_vclip3 | 62 | 0.151 | 0.240 | -0.082 |

**None reach statistical significance.** Adaptive R200,dyn trimming gives
the least-bad (but still non-significant) result. Neither the biweight
estimator nor an added 3-sigma velocity clip improves on plain adaptive
trimming with gapper.

### 2. Richness (AMICO LAMBDA_STAR, independent of this pipeline) vs. mass

| config | n | r(richness, M200,dyn) | p | — | r(richness, M200,X-ray), same subsample | p |
|---|---|---|---|---|---|---|
| fixed750 | 112 | -0.064 | 0.501 | | 0.641 | 1.6e-8 |
| r200x_trim | 112 | -0.079 | 0.408 | | 0.641 | 1.6e-8 |
| r200dyn_adaptive | 112 | 0.012 | 0.904 | | 0.641 | 1.6e-8 |
| r200dyn_adaptive_biweight | 90 | -0.167 | 0.117 | | 0.641 | 1.6e-8 |
| r200dyn_adaptive_vclip3 | 111 | 0.021 | 0.828 | | 0.641 | 1.6e-8 |

**This is the most informative result of Phase 3.** An independent richness
estimator (AMICO's LAMBDA_STAR, computed by a separate detection pipeline)
correlates strongly and significantly with X-ray-derived M200 (r=0.64,
p=1.6e-8) but shows **no significant correlation with M200,dyn under any of
the five membership/aperture/estimator configurations tested** (|r|<0.17 in
all cases). Since richness-mass and X-ray mass-mass both check out
independently, but the dynamical mass fails against both, this localizes
the problem to the **sigma_v -> M200 relation itself** (currently the
interim Munari et al. 2013 calibration, D7 in `technical_validation_report.md`),
not to membership aperture or sigma_v estimator choice. This is consistent
with, and sharpens, the earlier audit's conclusion that D7 is the highest-
priority open item — no combination of aperture/estimator tested here fixes
the dynamical-mass calibration problem.

### 3. BGG stability under v2 apertures

Unlike the earlier v1.0-cycle check (which compared BGG identity between
two overlapping member-set configurations — a self-referential test flagged
in `REFEREE_REVIEW.md`), this check **re-scores** `select_bgg_for_group` on
each config's actually-trimmed member set and compares the result to the
v1.0 production BGG by position:

| config | agreement |
|---|---|
| fixed750 | 112/112 = 100.0% |
| r200x_trim | 112/112 = 100.0% |
| r200dyn_adaptive (+biweight, +vclip3) | 109/112 = 97.3% |

BGG selection is robust to the aperture-definition choice — a genuinely
independent confirmation (not just a stronger restatement) of the earlier
"~95%" finding.

### 4. Quality-metric proxies — explicitly not completeness/purity

Median `projected_member_fraction` = 1.0 for every config (expected: for
`n_specz>=5` groups, the dynamics-trimming step rarely removes members below
the floor). Median contamination score drops modestly from 0.33 (fixed750)
to 0.26-0.30 (adaptive configs) — consistent with adaptive trimming removing
some large-radius/high-velocity outliers, but this is a proxy defined
relative to our own search aperture, **not a measurement against ground
truth**. No true completeness or purity estimate exists for the real-data
pipeline (this gap was flagged in `REFEREE_REVIEW.md` and remains open —
resolving it is deferred to Phase 5's mock lightcones, which will provide
actual ground-truth membership to test against).

## Bottom line so far

Aperture choice and sigma_v estimator choice, tested exhaustively across
five configurations, do **not** produce a significant sigma_v-M200
correlation in the real data, and BGG selection is robust regardless. The
new richness cross-check adds real evidence that the bottleneck is the
sigma_v-M200 calibration relation (D7), not membership methodology. This
argues for prioritizing Phase 4/8 (self-calibration or resolving D7) over
further membership-side tuning, but that determination should be confirmed
via Phase 4/5 before making a final Phase 7 recommendation.

## Reproducibility

```
python scripts/pipeline/membership_v2_dynamics.py --catalog both \
    --config fixed750,r200x_trim,r200dyn_adaptive,r200dyn_adaptive_biweight,r200dyn_adaptive_vclip3
python scripts/pipeline/membership_v2_phase3_validation.py
```
Outputs: `outputs/results/membership_v2/*.csv` (gitignored; scripts are
committed and fully reproduce these tables from the existing v1.0 member
catalogs and the X-ray/richness catalogs in `cosmos-web-xray-igm`).
