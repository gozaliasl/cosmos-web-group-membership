# membership-v2 — Phase 7: Final Recommendation

Branch: `membership-v2`. Consolidates Phases 2-6
(`membership_v2_phase2_3_report.md`, `membership_v2_phase4_report.md` +
`membership_v2_phase4_checkpoint.md`, `membership_v2_phase5_report.md`,
`membership_v2_phase6_report.md`, `DECISION_LOG_V2.md` D11-D17).

## Framing

This is not a "we found the best dynamical-mass method" result. Phase 6's
error-budget decomposition shows that **intrinsic halo dynamical scatter
(>=46% of total variance) and finite-N sigma_v estimation noise (a further
~56% relative increase) dominate the error budget, and neither is fixable
by membership methodology, aperture choice, estimator choice, or
calibration choice.** The correct, evidence-based conclusion is:

> At the current richness (n_specz typically 5-15) and spectroscopic
> completeness of this survey, individual-group dynamical masses are
> intrinsically low-confidence. The pipeline should prioritize quality
> flags, uncertainty reporting, and ensemble/statistical use of M200,dyn
> rather than over-interpreting individual group values.

Membership/aperture tuning (Phases 2-4) is **not** the dominant limitation.
It was tested exhaustively (5 aperture/estimator configurations against
real data, then again in the Phase 6 ablation) and found to be a real but
secondary contributor. The dominant limitations, confirmed independently by
both the real-data richness cross-check (Phase 3) and the mock-lightcone
error budget (Phase 6), are **intrinsic halo dynamics and finite-N sigma_v
estimation noise**.

---

## 1. Final recommended production methodology for Project B

**No change to the production default.** v1.0's existing methodology
(fixed 750 kpc search aperture, gapper sigma_v estimator, Munari et al.
2013 interim calibration, D1's probabilistic P_z*P_v membership gate)
remains the recommended production configuration. Nothing tested in
membership-v2 demonstrates a superior alternative worth adopting as
default (D13, D15). What changes is **how the output is used and reported**,
not the algorithm that produces it:

- Every group carrying a dynamical mass (SPECZ or XRAY+SPECZ method) should
  carry an explicit, propagated per-group uncertainty already
  incorporating the Phase 6 scatter floor (see Section 7), not just the
  formal statistical error on the gapper estimate.
- Individual M200,dyn values should be flagged as low-confidence by
  default (Section 6) and only used in aggregate/statistical contexts
  unless a specific group has corroborating evidence (e.g. X-ray agreement,
  XRAY+SPECZ method).

## 2. Which v2 features should be adopted now

- **The v2 quality-metric toolkit** (`membership_v2_dynamics.py`'s
  per-group contamination score, phase-space outlier fraction,
  members-outside-R200 fraction, dynamical reliability score, membership
  confidence) should be adopted as an **additional, informational output**
  alongside the existing v1.0 catalog columns -- not as a selection or
  trimming criterion, but as a transparency tool that lets downstream users
  see which groups are more or less reliable. This is a low-risk addition:
  it changes no existing numbers, only adds columns.
- **Multi-estimator sigma_v reporting** (gapper, biweight, std, robust-MAD
  computed side-by-side) should be adopted for the same reason --
  informational, not decision-driving. Phase 6 showed gapper/std are
  comparably robust and biweight is measurably worse (larger scatter) for
  this richness regime, which itself is a useful, publishable methodological
  finding worth keeping visible in the catalog rather than only in this
  report.
- **The re-scored BGG-stability check** (Phase 3) should be adopted as a
  standard release-checklist item (add to `VALIDATION_CHECKLIST.md`'s BGG
  section) since it is a stronger, independent test than the original
  subset-overlap check.

## 3. Which v2 features should remain experimental

- **R200,X-trim and adaptive R200,dyn aperture modes** remain
  experimental/opt-in only (D13, unchanged). They do not improve
  sigma_v-M200 recovery and should not be used to select or trim
  production catalog membership.
- **The 4 tested sigma_v-M200 calibrations beyond Munari** (Evrard 2008,
  both TNG self-calibrations) remain experimental/diagnostic only (D15,
  D17). None is a demonstrated improvement; do not substitute any of them
  for Munari in production.
- **The mock-lightcone and benchmark scripts themselves**
  (`membership_v2_phase5_mock_lightcone.py`,
  `membership_v2_phase6_benchmark.py`) remain research/diagnostic tools,
  not part of the production catalog pipeline -- they depend on an
  external TNG data product not shipped with this repo and are not
  intended to run per-release.

## 4. Which v1.0 choices are validated and should remain unchanged

- **Gapper as the sigma_v estimator** -- now validated on comparative,
  not just literature-precedent, grounds (Phase 6: gapper/std scatter
  ~0.55 dex vs. biweight ~0.76 dex under realistic mock conditions).
- **750 kpc fixed search aperture (CW-All) / compact-group-scaled aperture
  (CW-HCG)** -- tested against two alternatives (R200,X-trim, adaptive
  R200,dyn) and found not to be a significant limiting factor (Phases 3-4,
  6).
- **The probabilistic P_z*P_v membership gate (D1)** -- the Phase 6 ladder
  shows it contributes ~0% marginal variance beyond aperture + estimator +
  incompleteness effects; it is doing its job (removing genuine outliers)
  without adding noise.
- **BGG selection algorithm** -- re-confirmed robust (97-100% identity
  stability) under three independently re-scored aperture definitions
  (Phase 3), a stronger test than the original v1.0-cycle check.
- **MIN_SPECZ_FOR_DYNAMICS=5 floor and Nspec>=10 quality flag (D6)** --
  unchanged; still the correct interim mitigation given the degradation
  pattern is not resolved by any tested methodology (Phase 4).
- **Munari et al. (2013) as the documented interim calibration (D7)** --
  unchanged; not superseded by any alternative tested (D15).

## 5. Which quantities are safe for science use

- **X-ray-derived quantities** (flux, luminosity, temperature, M200_Temp,
  R200_kpc for XRAY-method groups) -- entirely out of scope for
  Project B/membership-v2; unaffected by any finding here. (Project A
  covers their own validation status separately.)
- **Group membership (is_member) and richness (n_specz, n_photz)** -- the
  probabilistic membership selection itself is well-validated (v1.0 audit
  + Phase 3/4 re-confirmation) and safe for science use, including as an
  input to richness-based analyses (confirmed by the independent AMICO
  richness cross-check, Phase 3).
- **BGG identity** -- safe for science use; robust across all tested
  aperture/methodology variants.
- **M200 for the XRAY and XRAY+SPECZ (agreement case) methods** -- safe,
  since these are anchored to the independently-validated X-ray mass, with
  the dynamical estimate only contributing when it agrees within
  uncertainty.
- **STATISTICAL / ENSEMBLE use of M200,dyn** (e.g. stacked or binned
  sigma_v-M200 relations, population-level richness-mass trends, sample
  medians) -- appropriate; the Phase 6 error budget describes per-group
  noise that averages down with sample size in the usual way, and ensemble
  trends are not subject to the same individual-group unreliability.

## 6. Which quantities should be treated as low-confidence

- **Individual-group M200,dyn, R200,dyn, and sigma_v for SPECZ-only method
  groups** -- explicitly low-confidence per group. Phase 6: even under
  idealized (contamination-free) conditions, per-group scatter is ~0.5 dex
  (>3x in mass) with 37-43% of groups more than 0.5 dex from the true
  value; under realistic conditions, scatter is similar (~0.5-0.6 dex) with
  bias additionally shifting by 0.3-0.9 dex depending on calibration and
  contamination level.
- **M200,dyn for the XRAY+SPECZ disagreement case** (`r200_disagreement=True`)
  -- already flagged in v1.0, and this finding reinforces that flag's
  importance; the X-ray value should be preferred in these cases (as v1.0
  already does).
- **Any group with Nspec in the range where the real-data correlation was
  shown to invert** (approximately Nspec=8-14, Phase 4) -- these should not
  be treated as *more* reliable than lower-Nspec groups despite higher
  member counts; richness alone does not track dynamical-mass reliability
  in this sample.
- **Groups at low true mass (log M200<12.5) or high redshift (z>1)** -- the
  mock-lightcone benchmark shows catastrophic-failure rates of 70-90%
  (low-mass) and 60-80% (high-z) for every calibration tested; these
  regimes should be treated as the least reliable part of the dynamical-mass
  catalog.

## 7. Required flags for the catalog

The following should be added to (or confirmed already present in) the
production catalog schema, per this investigation's findings:

- `dynamical_mass_confidence` (new, categorical: high/medium/low) --
  derived from n_specz, contamination_score (v2 toolkit, Section 2),
  and whether the group falls in the low-mass/high-z regime flagged in
  Section 6. Default to **low** for any SPECZ-only-method group; this is
  the single most important new flag motivated by this investigation.
- `sigma_v_scatter_floor_dex` (new, fixed reference value = 0.5, or the
  per-regime value from Phase 6's mass/z/completeness-binned tables) --
  an explicit, documented irreducible uncertainty to attach to any
  M200,dyn value, distinct from (and typically larger than) the formal
  gapper statistical error already reported.
- `n_specz_regime_flag` (new, categorical: "degradation-regime" for
  Nspec 8-14, else "standard") -- flags the specific richness range where
  Phase 4 found the real-data correlation to invert, so downstream users
  do not treat these groups as more reliable due to higher member counts.
- `r200_disagreement` (existing, v1.0) -- confirmed important; no change,
  but its role in down-weighting dynamical mass should be stated explicitly
  in any paper using this catalog.
- `membership_confidence`, `contamination_score`,
  `dynamical_reliability_score` (new, from the v2 toolkit, Section 2) --
  adopted as informational columns per Section 2.

## 8. Final decision log update

See `DECISION_LOG_V2.md` D18 (added below): membership-v2's overall
conclusion is recorded as **"validated the production default, quantified
its limitations, did not identify a superior alternative worth adopting."**
This is itself the correct outcome of a rigorous methodology investigation
-- absence of a better alternative, backed by an exhaustive, quantified
search, is a strong and useful result, not a null one.

## 9. Publication-ready summary paragraph for the methods section

> Group membership was assigned using a probabilistic model combining
> photometric-redshift consistency and radial-velocity offset (Section
> [X]), with brightest-group-galaxy selection via a mass-floored,
> centrality-weighted hybrid ranking (Section [X]). We assessed the
> robustness of this methodology, including alternative search apertures
> (fixed vs. dynamically-estimated R200), alternative velocity-dispersion
> estimators (gapper, biweight, standard deviation), and the adopted
> sigma_v-M200 calibration, against both the observed spectroscopic sample
> and a mock lightcone built from IllustrisTNG halos with injected
> observational selection (spectroscopic incompleteness, redshift
> uncertainty, and field contamination). This analysis shows that
> membership and aperture choices are not the dominant source of
> uncertainty in the dynamical mass estimates; rather, the error budget is
> dominated by intrinsic scatter in the halo dynamical state and by
> finite-sample noise in the velocity-dispersion estimator at the
> richnesses typical of this survey (n_spec ~ 5-15 members per group). We
> therefore report dynamical masses with an explicit, conservative
> uncertainty floor (Section [X]) and recommend their use in ensemble or
> statistical analyses; individual-group dynamical masses, particularly for
> systems without independent X-ray confirmation, should be interpreted
> with corresponding caution.

## 10. Remaining work before any standalone methodology paper

1. **Mock-observed recalibration** (D16/D17): fit a new sigma_v-M200
   relation directly against mock-*observed* (post-selection) sigma_v,
   with explicit mass- and redshift-dependent terms, rather than a global
   power law fit to ground truth. This is the highest-value remaining item
   and was explicitly deferred, not attempted, in Phase 5-6.
2. **Extend the mock lightcone's realism**: the current field-contamination
   model (Phase 5) is a stated, labeled modeling choice (interlopers
   uniformly filling the velocity-acceptance window), not a direct
   measurement of COSMOS's actual interloper population; validating or
   replacing this recipe with a measured/simulated LSS-consistent
   interloper rate would strengthen the Phase 5/6 conclusions.
3. **Verify the mock's richness/mass coverage matches the real survey**:
   the TNG sample used (868 halos, z<~1.7-2) should be checked for
   representativeness against CW-All/CW-HCG's actual n_specz and redshift
   distributions before quoting Phase 6's error-budget percentages as
   survey-wide numbers rather than regime-specific diagnostics.
4. **Re-run the completeness/purity/error-budget benchmark with the
   dynamical_mass_confidence flag's actual thresholds once defined**
   (Section 7), to confirm the flag correctly separates high- from
   low-reliability groups in the real catalog, not just in the mock.
5. **Resolve or explicitly carry forward Paper I/II's dependency on this
   finding**: any Paper I/II result that uses individual-group M200,dyn
   values (rather than ensemble/statistical use) should be revisited in
   light of Section 6's low-confidence flagging -- this connects directly
   to `FUTURE_WORK.md` item 1 (Blocking) in the v1.0 documentation set.
6. **Formal write-up of the biweight-vs-gapper comparison** (Phase 6) as a
   citable, standalone methodological result, since it is a genuine,
   comparative (not just precedent-based) justification for the estimator
   choice that would strengthen a methods paper.
