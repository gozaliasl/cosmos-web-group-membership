# Project B Executive Summary — Membership Selection, BGG Selection, Velocity Dispersion & Dynamical Mass (membership-v2)

*A 15-minute orientation for a new collaborator. For full detail see the
phase reports, `DECISION_LOG_V2.md`, and `REFEREE_REVIEW_membership_v2.md`
on the `membership-v2` branch of `cosmos-web-group-membership`.*

---

## Objectives

The COSMOS-Web group catalog (v1.0, production) assigns group membership
via a probabilistic model, selects a brightest group galaxy (BGG), and — for
groups with enough spectroscopic members — estimates a dynamical mass from
the velocity dispersion (sigma_v) via an interim literature calibration
(Munari et al. 2013). An earlier audit of v1.0 found that this dynamical
mass branch behaves strangely: its correlation with independent X-ray mass
*degrades* as spectroscopic richness increases, opposite to what
simulations predict. Project B's goal was to determine, rigorously and
without prior commitment to an answer, **whether a better membership/
dynamical-mass methodology exists, and if not, what is actually causing the
problem.** The investigation was explicitly instructed not to just tune
parameters, but to test whether a demonstrably superior methodology could
be built — and to falsify its own conclusions wherever possible.

## Methodology

Work proceeded in seven phases on an isolated `membership-v2` branch,
with the v1.0 production pipeline left completely untouched and fully
reproducible throughout (tagged `v1.0-audited-baseline`):

1. **Branch setup** — committed the previously-uncommitted v1.0 pipeline
   scripts as a true baseline, then branched.
2. **Configurable module** (`membership_v2_dynamics.py`) — built a
   side-by-side (not a replacement) tool supporting 3 membership aperture
   modes (fixed, X-ray-R200-trimmed, adaptive dynamical R200), 4 sigma_v
   estimators (gapper, biweight, standard deviation, robust-MAD), and a
   suite of per-group quality metrics (contamination score, phase-space
   outlier fraction, dynamical reliability score). Verified to exactly
   reproduce v1.0's own numbers under matching settings.
3. **Real-data validation (Phase 3)** — tested all 5 aperture/estimator
   configurations against independent X-ray mass and an independent
   richness estimator (AMICO's LAMBDA_STAR), plus a genuinely re-scored
   (non-self-referential) BGG stability check.
4. **Causal follow-up (Phase 4)** — tested whether reducing measured
   contamination via better trimming repairs the correlation.
5. **Mock lightcone + calibration test (Phase 5)** — built mock
   observations from 868 real IllustrisTNG halos (real ground-truth mass
   and velocities), injected realistic spectroscopic incompleteness,
   redshift error, and field contamination, then tested 4 candidate
   sigma_v-M200 calibrations (Munari 2013, Evrard 2008, and two in-sample
   TNG self-calibrations) against true mass.
6. **Full benchmark + error-budget decomposition (Phase 6)** — a 3
   estimator x 3 completeness x 4 calibration grid, plus a matched-sample
   sequential ablation isolating how much each factor (estimator noise,
   incompleteness, aperture, membership gate, contamination, calibration)
   contributes to total error.
7. **Final recommendation (Phase 7)** — conservative synthesis, explicitly
   framed to avoid over-claiming a "best method."
8. **Independent referee review** — a self-skeptical audit of Phases 2-7
   that found and corrected two overclaims and one scope gap before this
   summary was written (see "Downgraded findings" below).

## Major experiments and what they showed

- **5-configuration real-data test**: no aperture/estimator combination
  produces a statistically significant sigma_v-M200 correlation (best case
  r=0.16, p=0.22; n~63 X-ray-detected groups with n_specz>=5).
- **Richness cross-check**: an independent richness estimator correlates
  strongly with X-ray mass (r=0.64, p=1.6e-8) but with *no* tested
  dynamical-mass configuration (|r|<0.17 for all five) — localizing the
  problem away from membership methodology.
- **Contamination-reduction test**: adaptive trimming cuts the measured
  contamination proxy by ~3.5x in the most diagnostic richness bin, yet the
  correlation there gets *worse*, not better — evidence against "fixable
  residual contamination" as the full story.
- **BGG re-scored stability test**: 97-100% identity agreement across three
  independently re-trimmed aperture definitions.
- **Mock-lightcone calibration test** (868 real TNG halos, 731 mock
  observations): even with **zero** injected contamination, all four
  calibrations show ~0.5 dex scatter and 37-43% catastrophic-failure rates;
  under realistic contamination, bias *flips sign* and calibration rankings
  invert entirely — no calibration is a consistent winner.
- **Error-budget ablation** (matched n=39 TNG halos): confirmed that
  finite-N sigma_v sampling noise is a large, robust, well-measured
  contributor (+56% relative variance vs. using the true velocity
  dispersion directly); a further scatter/mismatch term exists but its
  clean size was not established (see Downgraded findings).
- **Estimator comparison**: gapper and standard deviation perform
  comparably; biweight shows markedly larger scatter (0.76 vs. 0.55 dex) —
  the first genuinely comparative (not just precedent-based) justification
  for v1.0's original estimator choice.

## Validated findings (high confidence)

- **Membership aperture and sigma_v estimator choice are not the dominant
  limitation** of the dynamical-mass pipeline — tested exhaustively (5
  real-data configurations, then again in the mock ablation) and
  consistently found secondary.
- **BGG selection is robust** to how the group aperture is defined.
- **Individual-group dynamical masses (M200,dyn, R200,dyn, sigma_v) should
  be treated as low-confidence** and used in ensemble/statistical
  analyses, not interpreted group-by-group — the single best-supported
  conclusion of the entire project, backed independently by real-data
  non-significance, mock-based scatter/failure rates, and calibration-
  ranking instability.
- **Gapper is a comparatively justified choice** of sigma_v estimator for
  this survey's richness regime.
- **Field contamination's dominant effect is a mass bias, not scatter.**

## Downgraded findings (corrected after independent review — do not cite without caveats)

An internal referee-style audit of Phases 5-7 (performed before closing the
project) found two specific numeric claims were overstated, and corrected
them across all affected documents:

- ~~"Intrinsic halo dynamics accounts for >=46% of total variance"~~ — this
  was computed using only one calibration (Munari) applied to TNG's true
  velocity dispersion, entangling **calibration mismatch** with genuine
  **intrinsic scatter**, on a sample of only 39 halos with no uncertainty
  interval. **Withdrawn as a headline number**; now stated only as an
  upper bound pending a re-analysis across all four calibrations with
  bootstrap intervals.
- ~~"The membership gate contributes ~0% of the error budget"~~ — the mock
  test used a simplified, single-pass velocity gate, not the real
  production algorithm, which iterates up to 10 times and re-estimates
  each group's own velocity dispersion each round (confirmed by direct
  code inspection). **Marked as requiring re-validation** with a faithful
  reproduction of the real gate before being cited.
- **Scope gap**: no genuine external, published, group-specific
  sigma_v-M200 calibration was ever tested, despite this being explicitly
  requested. Two in-sample TNG self-calibrations were substituted as
  internal alternatives — a real, useful comparison, but not equivalent to
  testing an independent literature relation.

Neither correction changes the project's actual recommendation — both
withdrawn numbers were supporting detail, not load-bearing for the
qualitative conclusion above.

## Final recommendations

- **No change to the v1.0 production pipeline default**: fixed 750 kpc
  aperture, gapper estimator, Munari et al. (2013) interim calibration, and
  the existing probabilistic membership gate are all retained — not
  because any has been proven optimal, but because nothing tested here is
  a demonstrated improvement.
- **Adopt, as informational-only additions**: the v2 quality-metric
  toolkit (contamination score, dynamical reliability score, etc.) and
  multi-estimator sigma_v reporting — useful transparency, not selection
  criteria.
- **Add new catalog flags**: `dynamical_mass_confidence` (defaulting to
  low for spec-z-only groups), `sigma_v_scatter_floor_dex`, and
  `n_specz_regime_flag`.
- **Do not fit a new sigma_v-M200 calibration yet** — the dominant error
  terms identified so far are not calibration-fixable on their own.

## Limitations

- The mock lightcone's field-contamination model (injection rate,
  velocity-offset distribution) is a stated, untested modeling choice, not
  a direct measurement of COSMOS's real interloper population.
- The TNG sample's mass/richness/redshift coverage was never formally
  checked against the real survey's actual distribution.
- `veldisp_halo_kms` (used as ground truth in the mock) was never
  independently verified for its exact definition.
- No completeness/purity ground truth exists for the *real* (non-mock)
  catalog — this remains an open gap inherited from the original v1.0
  audit.
- Statistical rigor (bootstrap/jackknife intervals) was applied
  inconsistently across phases; the most consequential Phase 6 numbers
  currently lack them.

## Future work (priority order)

1. Re-run the disputed S0 ablation stage across all four calibrations with
   bootstrap intervals, to responsibly establish (or bound) an
   intrinsic-scatter-floor figure.
2. Rebuild the mock's membership gate as a faithful, iterating
   reproduction of the real production algorithm and re-test.
3. Verify `veldisp_halo_kms`'s exact definition against TNG documentation.
4. Sensitivity-sweep the field-contamination injection rate.
5. Obtain and test a genuine external, published group-specific
   sigma_v-M200 calibration.
6. If items 1-2 confirm it is still warranted, undertake a mock-*observed*
   (not ground-truth) recalibration with explicit mass- and
   redshift-dependent terms.
7. Revisit any Paper I/II result that relies on individual-group M200,dyn
   values in light of the low-confidence flagging recommended here.
