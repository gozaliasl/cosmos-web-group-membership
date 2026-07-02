# Decision Log — membership-v2

Continues the numbering of the v1.0 Decision Log
(`cosmos-web-xray-igm/docs/technical_validation_report.md`, D1-D10). This
log covers decisions made on the `membership-v2` branch only; v1.0's
decisions and production defaults are unchanged.

| # | Decision | Status | Evidence |
|---|---|---|---|
| D11 | Build `membership_v2_dynamics.py` as a new, side-by-side, fully configurable module rather than editing v1.0's scripts in place. | Accepted | v1.0 must remain completely reproducible (explicit branch requirement); verified `fixed750` config exactly reproduces v1.0's own numbers. |
| D12 | Support 3 aperture modes (fixed, R200,X-trim, adaptive dynamical R200) and 4 sigma_v estimators (gapper, biweight, std, robust-MAD), all reported side-by-side rather than picking one. | Accepted | Needed to test membership methodology as a possible explanation for the Nspec-dependent sigma_v-M200 degradation (v1.0 audit, open item). |
| D13 | **Do not adopt adaptive R200,dyn trimming (or any tested v2 config) as the pipeline default.** Keep `fixed750` (= v1.0's existing behavior) as the default; all v2 configs remain opt-in/experimental. | **Accepted** | Phase 3/4: adaptive trimming measurably lowers the contamination proxy (~3.5x in the Nspec=10-14 bin) but does not improve, and in that same bin coincides with a *worse*, sigma_v-M200 correlation (-0.582 -> -0.814). No config among the five tested reaches a statistically significant sigma_v-M200 or richness-M200,dyn correlation. Adopting any of them now would trade pipeline complexity for no demonstrated accuracy gain. |
| D14 | Reframe Phase 5 (mock lightcones) specifically as a **calibration-diagnostic test** (inject a known sigma_v-M200 relation, check recovery under realistic selection) rather than a general realism upgrade. | Accepted | Phase 3+4 evidence (richness cross-check; contamination-reduction test) converges on the sigma_v-M200 calibration (D7, Munari et al. 2013 interim) as the leading suspect over membership methodology -- Phase 5 needs to be able to confirm or reject this specifically, not just produce a more realistic mock in general. See `membership_v2_phase4_checkpoint.md` Section 7 for the exact test design required. |
| D15 | **Do not replace Munari et al. (2013) with Evrard et al. (2008) or either TNG self-calibration.** Keep it as the documented interim placeholder (D7 unchanged). | **Accepted** | Phase 5 mock-lightcone test (4 calibrations x idealized + contaminated cases, 731 mock observations from 868 real TNG halos): no calibration is a consistent improvement across conditions -- TNG self-calibrations perform best in the idealized (contamination-free) case (bias -0.11 to -0.12 dex) but worst under realistic contamination (bias +0.58 to +0.64 dex); Munari performs worst idealized (-0.23 dex) but least-badly contaminated (+0.44 dex). All four show ~0.5 dex scatter and 37-65% catastrophic-failure rates depending on conditions -- none is a clear winner. |
| D16 | The dynamical-mass branch (any tested calibration) carries irreducible ~0.5 dex scatter and 40-65% per-group catastrophic-failure risk in this survey's regime (n_specz~5-15, group/poor-cluster mass, z up to ~2), **not primarily attributable to calibration choice**. A real fix requires recalibrating against mock-*observed* (post-contamination, post-incompleteness) sigma_v with explicit mass- and redshift-dependent terms, not a global power-law swap. | Accepted (finding, not yet actioned) | Phase 5, Section 3-4 (`membership_v2_phase5_report.md`): idealized-case bias/scatter alone (zero injected contamination) already shows -0.1 to -0.3 dex bias and ~0.5 dex scatter for all 4 calibrations; bias and catastrophic-failure rate both degrade sharply at low mass (log M200<12.5) and high z (z>1) for every calibration tested. Flagged as required future work, not implemented on this branch (would require a full mock-observed recalibration campaign, out of Phase 5's diagnostic scope). |

| D17 | **Do not fit a new calibration yet.** Phase 6's error-budget decomposition shows the two largest error sources (intrinsic halo dynamics, ~46%+ of variance; sigma_v estimator sampling noise, ~26-56% relative) are not fixable by any calibration choice, so a new fit alone cannot resolve the dominant error terms. | Accepted | `membership_v2_phase6_report.md`: matched-sample ablation ladder (n=39 core TNG halos) shows the S0 (intrinsic, true veldisp_halo_kms) floor already accounts for 46% of final total variance; S0->S1 (adding gapper sampling noise at realistic richness) adds a further, robust +56% relative increase. Any future recalibration effort should target the bias contribution (calibration + contamination combined, D16) rather than attempt to eliminate the scatter floor, which is irreducible with this pipeline's inputs. |

| D18 | **membership-v2 conclusion: validated the v1.0 production default, quantified its limitations, did not identify a superior alternative worth adopting as default.** Adopt the v2 quality-metric/multi-estimator toolkit as informational catalog columns only (not selection/trimming criteria). Add `dynamical_mass_confidence`, `sigma_v_scatter_floor_dex`, and `n_specz_regime_flag` to the catalog schema; default individual SPECZ-method M200,dyn to low confidence. Do not fit or adopt a new sigma_v-M200 calibration at this time. | **Accepted (final, Phase 7)** | `membership_v2_phase7_recommendation.md`, consolidating D11-D17: membership/aperture tuning (Phases 2-4) and calibration substitution (Phase 5) were both tested exhaustively and found to be secondary; Phase 6's error-budget decomposition shows intrinsic halo dynamics (>=46% of variance) and finite-N sigma_v sampling noise (~56% relative increase) dominate, and neither is fixable by any methodology tested on this branch. This is the correct, conservative conclusion of a rigorous investigation, not a null result. |

## Relationship to v1.0 Decision Log D5-D7

- D5 (v1.0): R200,X-based membership trimming tested as a diagnostic-only
  alternative, not adopted as default -- **D13 extends this same
  conclusion to the full v2 aperture/estimator space**, now with
  quantitative contamination-vs-correlation evidence rather than only the
  original sample-size argument.
- D6 (v1.0): Nspec>=10 flagged as higher-confidence despite no adopted fix
  for the underlying degradation -- **unchanged**; D13/Phase 4 show the
  degradation is not resolved by the v2 methodology, so the D6 quality flag
  remains the correct interim mitigation.
- D7 (v1.0): Munari et al. (2013) interim sigma_v-M200 calibration -- **now
  the leading, evidence-supported (not proven) suspect** for the
  degradation itself; D14 commits Phase 5 to directly testing this.
