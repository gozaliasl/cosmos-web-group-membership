# membership-v2 — Phase 4 Report: does better membership methodology repair the Nspec-dependent degradation?

Branch: `membership-v2`. Builds directly on Phase 2/3
(`docs/membership_v2_phase2_3_report.md`) and the v1.0-cycle audit
(`cosmos-web-xray-igm/docs/technical_validation_report.md` Sections 8-11;
`REFEREE_REVIEW.md`).

## Question

The unresolved issue carried over from the v1.0 audit: the real
CW-All(spec-z, X-ray) overlap sample shows the sigma_v-M200 correlation
*degrade* with increasing spectroscopic richness (Nspec) beyond ~N=7-8,
opposite to what TNG mocks predict. The v1.0 audit's evidence-supported (not
proven) interpretation was residual interloper contamination in
larger-Nspec groups. Phase 4 asks: **does the v2 toolkit (configurable
aperture + explicit contamination/quality metrics) let us do better than
"evidence-supported interpretation"?**

This report covers the observational half of Phase 4 only (selection,
projection, X-ray-selection, contamination — all testable on real data).
Testing the "limitations of the current simulations" branch requires Phase
5's realistic mock lightcones and is deferred there.

## T1 — do the new quality metrics scale with Nspec?

| metric | Spearman r | p | n |
|---|---|---|---|
| contamination_score | 0.019 | 0.844 | 112 |
| phase_space_outlier_frac | **0.189** | **0.046** | 112 |
| members_outside_r200_frac | -0.033 | 0.730 | 112 |

Only `phase_space_outlier_frac` shows a significant (if weak) rise with
Nspec. The composite `contamination_score` — which blends this with
members-outside-R200 and velocity-tail fraction — washes this signal out,
because the other two components don't share the trend. This is a useful,
concrete result: **not all plausible "contamination" proxies actually track
Nspec** — only the specifically phase-space-based one does, mildly.

## T2 — is X-ray selection introducing the contamination?

One hypothesis carried from the v1.0 audit was that X-ray selection (flux-
limited, Malmquist/Eddington-biased) could preferentially select groups in
denser environments, compounding interloper contamination for the X-ray+
spec-z overlap subsample specifically.

| | X-ray-detected | non-X-ray-detected |
|---|---|---|
| n | 63 | 49 |
| median n_specz | 8.0 | 7.0 |
| median contamination_score | 0.320 | 0.410 |

Mann-Whitney p=0.14 (not significant), and the *direction* is opposite to
the hypothesis: X-ray-detected groups have **lower**, not higher, median
contamination at similar richness. This doesn't support the X-ray-selection-
compounding hypothesis as tested here; it neither confirms nor strongly
refutes it given the non-significance, but it removes it as a leading
explanation pending a larger sample.

## T3 — does lowering contamination via better trimming repair the correlation? (key result)

Re-running the same Nspec-binned sigma_v-M200 correlation test from the
v1.0 audit, this time also under the v2 `r200dyn_adaptive` config
(spatially trims to the dynamically-converged R200 rather than the fixed
750 kpc search net):

| config | Nspec bin | n | Pearson r | p | median contamination |
|---|---|---|---|---|---|
| fixed750 | 5-7 | 25 | +0.421 | 0.036 | 0.320 |
| fixed750 | 8-9 | 10 | -0.709 | 0.022 | 0.305 |
| fixed750 | 10-14 | 15 | -0.582 | 0.023 | 0.264 |
| fixed750 | 15+ | 13 | -0.049 | 0.874 | 0.339 |
| r200dyn_adaptive | 5-7 | 39 | +0.389 | 0.014 | 0.160 |
| r200dyn_adaptive | 8-9 | 12 | -0.587 | 0.045 | 0.206 |
| r200dyn_adaptive | 10-14 | 8 | **-0.814** | 0.014 | **0.075** |
| r200dyn_adaptive | 15+ | 4 | +0.366 | 0.634 | 0.339 |

The qualitative pattern (positive at low Nspec, strongly negative at
Nspec=8-14, uncertain/small-n at Nspec>=15) **persists under adaptive
trimming**. In the Nspec=10-14 bin specifically, adaptive trimming lowers
the median contamination score by a factor of ~3.5 (0.264 -> 0.075) — a
real, substantial reduction in the measured contamination proxy — yet the
correlation gets **more** negative (-0.582 -> -0.814), not less.

**Caveat**: bin sizes are small (n=4-15 per bin), especially the two most
diagnostic ones (n=8, n=4) — none of these per-bin results should be read
as independently conclusive, and the r200dyn_adaptive Nspec=15+ bin (n=4) is
too small to interpret at all.

## Interpretation

**Direct observation**: reducing the phase-space/spatial contamination
proxy via adaptive trimming does not measurably improve, and if anything
coincides with a worse, sigma_v-M200 correlation in the richest bin tested.

**Evidence-supported interpretation, not proven**: this is more consistent
with the sigma_v-M200 relation itself (the interim Munari et al. 2013
calibration, D7) being miscalibrated for this sample's mass/richness regime
than with residual, fixable interloper contamination being the dominant
driver. This sharpens (does not contradict) Phase 3's independent finding
that richness correlates with X-ray M200 but not with any tested
M200,dyn configuration.

**Remaining uncertainty**: with n=112 groups total and single-digit counts
in the most diagnostic Nspec bins, a genuine change in what physically
causes the degradation (versus small-number noise in the bin-by-bin
correlation itself) cannot be fully distinguished from real data alone.
Phase 5's realistic mock lightcone is still needed to settle whether this is
truly a calibration problem (testable by injecting a known M200-sigma_v
relation and checking recovery) or reflects a real physical effect this
Nspec-binned test cannot isolate with the current sample size.

## Recommendation carried into Phase 5/7

Given T3, Phase 5's mock-lightcone effort should be designed to explicitly
test **both** hypotheses side-by-side: (a) inject a known, correct
sigma_v-M200 relation and check whether the *current* membership/aperture
pipeline (any of the 5 configs) recovers it under realistic observational
selection, and (b) if it does not, that isolates the problem to
methodology/selection rather than to D7's calibration specifically.

## Reproducibility

```
python scripts/pipeline/membership_v2_phase4_causal.py
```
Requires Phase 2's `membership_v2_dynamics.py` outputs to already exist for
`fixed750` and `r200dyn_adaptive`.
