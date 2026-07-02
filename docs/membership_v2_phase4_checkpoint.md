# membership-v2 — Phase 4 Checkpoint

Branch: `membership-v2`. Consolidates `membership_v2_phase2_3_report.md` and
`membership_v2_phase4_report.md` before Phase 5 begins.

## 1. What was tested

- Five membership/aperture/estimator configurations (`fixed750`, `r200x_trim`,
  `r200dyn_adaptive`, `r200dyn_adaptive_biweight`, `r200dyn_adaptive_vclip3`)
  run side-by-side on the full CW-All and CW-HCG member catalogs, all built
  on v1.0's validated member pools (v1.0 itself untouched).
- sigma_v-M200 correlation against independent X-ray mass, for every config.
- Richness (AMICO LAMBDA_STAR, an independently-computed quantity) against
  both M200,dyn (every config) and X-ray M200 (reference), on matched
  subsamples.
- BGG identity stability under each config's actual re-trimmed member set
  (not a subset-overlap proxy).
- Per-group quality metrics (contamination score, phase-space outlier
  fraction, members-outside-R200 fraction) against n_specz, against X-ray
  detection status, and against the Nspec-binned sigma_v-M200 correlation
  pattern first identified in the v1.0 audit.

## 2. What improved

- **BGG selection** is confirmed robust to aperture-definition choice
  (100% agreement for fixed750/r200x_trim, 97.3% for the adaptive configs),
  now via a genuinely independent re-scoring test rather than the earlier
  self-referential comparison flagged in `REFEREE_REVIEW.md`.
- **Measured contamination** (phase-space outlier fraction and the
  composite contamination score) is reduced by adaptive R200,dyn trimming,
  substantially so in the Nspec=10-14 bin (factor ~3.5 lower).
- The v2 module itself is a genuine methodological improvement in
  *transparency*: every group now carries an explicit, reproducible
  contamination/reliability/confidence score, where v1.0 had none.

## 3. What did not improve

- **The sigma_v-M200 correlation itself.** No config among the five reaches
  statistical significance against independent X-ray mass (best case:
  r200dyn_adaptive, r=0.156, p=0.222; worst case: fixed750, r=-0.006,
  p=0.965).
- **Richness-M200,dyn.** No config shows a significant richness-mass
  correlation (|r|<0.17 for all five), despite the same richness showing a
  strong, highly significant correlation with X-ray M200 on the identical
  subsample (r=0.64, p=1.6e-8).
- **The Nspec-dependent degradation pattern itself.** Adaptive trimming
  measurably lowers contamination in the Nspec=10-14 bin yet the
  correlation there gets *more* negative (-0.582 -> -0.814), not less.

## 4. What was ruled out (as leading explanations, on current evidence)

- **Membership aperture definition** as the primary driver: five different
  aperture strategies produce materially the same non-significant
  sigma_v-M200 outcome.
- **sigma_v estimator choice** as the primary driver: gapper, biweight, and
  a 3-sigma velocity clip all give statistically indistinguishable, equally
  non-significant results.
- **Simple, aperture-fixable interloper contamination** as a *complete*
  explanation: reducing the measured contamination proxy via better
  trimming does not repair, and in the most diagnostic bin coincides with a
  worsened, correlation.
- **X-ray selection actively compounding contamination**, as tested:
  X-ray-detected groups show *lower*, not higher, median contamination than
  non-detected groups at similar richness (not significant, but wrong-signed
  relative to the hypothesis).

None of these are proven false in every possible sense -- they are ruled out
specifically as the *primary, dominant* explanation given the evidence
gathered so far, per the same evidentiary standard applied throughout the
v1.0 audit.

## 5. What remains unresolved

- Whether the Nspec-dependent degradation reflects (a) a genuinely
  miscalibrated sigma_v-M200 relation, (b) a real physical effect not
  captured by the current contamination proxies, or (c) small-number noise
  in the bin-by-bin correlation test itself (the two most diagnostic bins
  have n=8 and n=4 groups) cannot be fully distinguished using real data
  alone at the current sample size.
- True completeness/purity of membership selection remains unmeasured (no
  ground truth exists in the real data) -- unchanged from the v1.0 audit's
  own limitation.

## 6. Why the Munari et al. (2013) sigma_v-M200 calibration is now the leading suspect

Three independent lines of evidence converge on the dynamical-mass
calibration (D7 in the v1.0 Decision Log) rather than membership
methodology:

1. **Richness cross-check (Phase 3)**: an entirely independent mass proxy
   (AMICO richness) recovers a strong, expected correlation with X-ray mass
   (r=0.64, p=1.6e-8) but fails against *every* dynamical-mass configuration
   tested (|r|<0.17). If membership/aperture choice were the dominant
   problem, at least one of the five configurations tested should have
   produced a detectable richness-mass signal; none did.
2. **Contamination-reduction test (Phase 4, T3)**: substantially lowering
   the measured contamination proxy via adaptive trimming does not improve,
   and in the richest diagnostic bin coincides with a worse, sigma_v-M200
   correlation -- inconsistent with contamination being the dominant,
   fixable driver.
3. **The calibration was already flagged, independently, as broken and
   interim** before this branch existed (D7: the original HZ-AGN
   calibration reproduces masses ~4 dex too high on its own calibration
   sample; the current Munari et al. 2013 substitute is an admitted
   placeholder, calibrated on cluster-scale halos (sigma_v ~500-1500 km/s)
   and extrapolated to this sample's group-scale regime (sigma_v ~100-300
   km/s) -- a materially different regime from where it was derived.

Individually, none of these three proves the calibration is the cause. Taken
together, they are the most parsimonious explanation consistent with all
observations gathered in Phases 2-4: membership methodology has been
exhausted as an explanation across five materially different
configurations, while the calibration was already independently suspect for
unrelated reasons predating this investigation.

## 7. What Phase 5 needs to confirm or reject this

Phase 5's mock lightcone must be designed as a **calibration-diagnostic
test**, not just a realism upgrade:

1. Inject a **known, correct** synthetic sigma_v-M200 relation into the mock
   lightcone (i.e. simulate galaxies whose true halo M200 and true velocity
   dispersion follow a specified, ground-truth relation).
2. Apply realistic observational selection on top: field
   foreground/background contamination, spectroscopic incompleteness,
   COSMOS-like redshift errors, and an X-ray flux-limited selection
   function consistent with the real survey.
3. Run the *actual* v1.0 and v2 membership/aperture pipelines (all five
   configs) on the resulting mock catalog and check recovery:
   - If the pipeline recovers the injected ground-truth relation
     (statistically, within the expected scatter) under realistic
     selection, the real-data degradation is attributable to the
     **calibration** (D7), not the membership methodology -- confirming
     this checkpoint's leading hypothesis.
   - If the pipeline systematically fails to recover even the known-correct
     relation under realistic selection, the problem is (at least partly)
     **methodological/selection-driven** after all, and this checkpoint's
     conclusion should be revised.
   - A specific sub-test: reproduce the Nspec=10-14 "adaptive trimming makes
     it worse" result on the mock. If the mock shows the same sign flip
     under injected ground truth, that is strong evidence of a
     methodology-side effect Phase 4 could not otherwise detect real-data
     side; if the mock does *not* reproduce it, that strengthens the
     calibration-bottleneck interpretation.

## Decision log update

See `DECISION_LOG_V2.md` (new, this branch): **do not adopt adaptive
R200,dyn trimming as a pipeline default.** It measurably reduces the
contamination proxy but does not improve (and in the most diagnostic bin,
coincides with a worse) mass-recovery correlation -- adopting it now would
trade a clearer contamination metric for no demonstrated gain in accuracy.
It remains available as an experimental, opt-in configuration
(`--config r200dyn_adaptive`) pending Phase 5's calibration test.
