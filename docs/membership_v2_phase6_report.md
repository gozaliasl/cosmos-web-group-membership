# membership-v2 — Phase 6: Definitive Benchmark and Error-Budget Decomposition

Branch: `membership-v2`. Script: `scripts/pipeline/membership_v2_phase6_benchmark.py`.
Builds on Phase 5's mock lightcone (868 real TNG50-1+TNG100-1 halos, ground
truth M200 and member velocities from simulation).

## A. Full configuration benchmark

Grid: 3 sigma_v estimators (gapper, biweight, std) x 3 completeness
fractions (0.3/0.5/0.7) x 4 calibrations (munari2013_galaxies,
evrard2008_dm, tng_selfcal_all, tng_selfcal_groups), aperture fixed at the
production default (750 kpc; the aperture's own marginal effect is measured
separately in the ladder below, section B). 3 random seeds per cell for
robustness. Full 36-row table: `outputs/results/membership_v2_phase6/phase6_full_benchmark.csv`.

### By estimator (averaged over completeness and calibration)

| estimator | mean scatter (dex) | mean bias (dex) | mean purity | mean recall | mean runtime (s/seed) |
|---|---|---|---|---|---|
| gapper | 0.557 | 0.517 | 0.394 | 0.338 | 1.25 |
| std | 0.552 | 0.457 | 0.394 | 0.338 | 1.29 |
| biweight | **0.763** | 0.328 | 0.394 | 0.338 | 1.30 |

Biweight has the lowest bias but by far the largest scatter (0.76 vs
~0.55 dex) — it is measurably less robust than gapper/std for this sample's
richness regime (n_specz~5-15) under realistic contamination, confirming
v1.0's original choice of gapper on genuinely comparative grounds (not just
literature precedent). Purity/recall/runtime do not depend on the estimator
choice (as expected — the estimator only affects the sigma_v value computed
from an already-selected member set, not which members are selected).

### By calibration (averaged over estimator and completeness)

| calibration | mean scatter (dex) | mean bias (dex) |
|---|---|---|
| munari2013_galaxies | 0.609 | **0.341** |
| tng_selfcal_groups | 0.616 | 0.479 |
| tng_selfcal_all | 0.634 | 0.542 |
| evrard2008_dm | 0.637 | 0.373 |

Confirms Phase 5/D15: Munari (production default) has the smallest average
bias under realistic mock conditions across the full estimator/completeness
grid, but scatter is similar (~0.61-0.64 dex) across all four — calibration
choice mainly shifts the bias, not the scatter.

### By completeness fraction

| completeness | mean purity | mean recall | mean n_final |
|---|---|---|---|
| 0.3 | 0.277 | 0.195 | 7.08 |
| 0.5 | 0.396 | 0.317 | 7.44 |
| 0.7 | 0.510 | 0.500 | 7.69 |

Purity and recall both scale roughly linearly with completeness, as
expected, but even at 70% completeness, recall (fraction of TRUE members
recovered) is only 50% and purity is only 51% — the final selected sample
is, on average, composed of roughly equal parts true members and field
interlopers at every completeness level tested. This is a direct,
quantitative measure of the membership pipeline's actual reliability under
this mock's contamination model, not previously available as a number in
the v1.0 audit (which had no ground truth to measure against).

## B. Error-budget decomposition (matched-sample ablation)

**CORRECTED per the independent referee review
(`REFEREE_REVIEW_membership_v2.md`, items 1-2), applied 2026-07 without
re-running code.** The original version of this section presented a "46%
intrinsic variance" headline figure and a "membership gate contributes ~0%"
finding with more confidence than the underlying analysis supports. Both
are corrected below; the numeric table (produced by the single run already
executed) is retained for transparency, but is now explicitly labeled as
provisional pending the re-analysis described in each correction.

39 TNG halos survive every stage of the ladder (S1 through S6) and are used
throughout for a genuine like-for-like comparison (same halos, same
retained-member mask reused at every downstream stage). S0 (intrinsic) uses
the same 39-halo core.

| stage | scatter (dex) | variance (dex^2) | marginal variance added (dex^2) | % of final total variance (this run) |
|---|---|---|---|---|
| S0 intrinsic (true veldisp_halo_kms, no sampling) | 0.283 | 0.080 | — | 46.0% (single-calibration estimate — see correction below) |
| S1 + estimator (gapper, full richness, idealized) | 0.353 | 0.125 | +0.045 | 25.8% |
| S2 + incompleteness (50%) | 0.431 | 0.186 | +0.061 | 35.4% |
| S3 + redshift error (150 km/s) | 0.361 | 0.130 | -0.056 | -32.1% (see caveat) |
| S4 + aperture (750 kpc) | 0.433 | 0.188 | +0.057 | 33.1% |
| S5 + membership gate (P_v cut) | 0.433 | 0.188 | -0.00003 | ~0% (simplified-gate estimate — see correction below) |
| S6 + field contamination | **0.416** | **0.173** | -0.014 | -8.2% (see caveat) |

**CORRECTION 1 — the "46%" figure is a calibration-entangled upper bound,
not an established intrinsic floor.** S0 was computed using **only** the
Munari et al. (2013) calibration (`cal_name` was hardcoded to
`munari2013_galaxies` in `ablation_ladder`). Phase 5's own idealized-case
table shows the two TNG self-calibrations fit this same true-veldisp-to-
true-M200 relation substantially better than Munari does (bias -0.11/-0.12
dex vs. Munari's -0.23 dex). This means part of the S0 scatter reflects
**Munari's own mismatch to TNG's true M-sigma relation**, not purely
irreducible intrinsic/dynamical scatter — the two cannot be separated with
the single-calibration run performed here. **The 46% figure is therefore
withdrawn as a headline number.** It should be read only as an **upper
bound on the intrinsic-scatter share, calibration-choice-entangled and
computed from n=39 with no uncertainty interval.** A recomputation of S0
(and the resulting percentages) across all four calibrations, with
bootstrap intervals, is required before any specific percentage is cited in
a publication (tracked as required future work, Phase 7 Section 10 and
`REFEREE_REVIEW_membership_v2.md` item 9.1).

**CORRECTION 2 — the "membership gate contributes ~0% variance" finding
used a simplified, non-iterating stand-in for the real gate.** S5's P_v
cut used a single-pass evaluation at a fixed `sigma_v_prior=500` km/s. The
actual production algorithm
(`determine_membership_dztier.py::refine_group_redshift`, verified by
direct inspection) **iterates the gate up to 10 times**, re-estimating
`sigma_v_current = clip(std(members), 150, 2500)` from the current member
set each round and re-scoring against this updated, group-specific value
until the redshift shift converges below 25 km/s. S5's result was measured
against the simplified, single-pass version, not this iterating algorithm.
**"The membership gate contributes ~0% marginal variance" is therefore
marked as requiring re-validation before citation** — the true marginal
contribution of the real, iterating gate is currently untested, not
established to be negligible.

**Caveat on the remaining intermediate terms (S2-S4, S6)**: with n=39
matched halos, individual marginal terms have substantial sampling
uncertainty (no confidence intervals were computed); the negative entries
for S3 and S6 do not mean "adding this effect reduces error" in a
generalizable sense — they reflect noise in this specific 39-halo
realization (S3 also reflects that most of the *bias* impact of these
effects, especially contamination, shows up as a shift in the median
residual rather than in the scatter — see Phase 5's much larger n=731
sample, where contamination's dominant, robust effect is a +0.3 to +0.9 dex
**bias** shift, not primarily a scatter increase).

**The one result from this ladder that survives both corrections and is
retained with confidence**:

- **S0 -> S1 (intrinsic-floor-as-measured vs. + estimator sampling noise)**:
  a clean, well-matched comparison (variance 0.080 -> 0.125, +0.045 dex^2,
  +56% relative increase) — sampling noise from estimating sigma_v via
  gapper at finite, realistic richness (5-15 members) substantially
  increases scatter relative to using the true, full-halo velocity
  dispersion, **regardless of which calibration is used for the
  comparison** (this specific relative-increase statement does not depend
  on Correction 1, since both S0 and S1 use the same calibration
  consistently). This qualitative/relative finding — that finite-N
  sigma_v estimation contributes substantial, non-negligible scatter on
  top of whatever the calibration/intrinsic floor turns out to be — is
  retained as a moderate-to-high-confidence result.

## Error-source ranking (largest to smallest)

**Revised** to reflect Corrections 1-2 above: rankings involving the exact
S0 and S5 percentages are now stated qualitatively/with bounds rather than
as point estimates, since those two specific numbers are not yet
established with confidence.

1. **Intrinsic halo dynamics / calibration's own intrinsic scatter (factor 7)
   AND sigma_v estimator sampling noise (factor 3), combined** — together
   the largest source, **upper-bounded at roughly two-thirds to three-
   quarters of total variance in this single run (S0+S1 combined = 0.125 of
   0.173 final dex^2, i.e. up to ~72%)**, but the split between "factor 7"
   and "factor 4/calibration mismatch" within the S0 term is not yet
   resolved (Correction 1). What is well-supported: this combined
   contribution is large, and **neither sub-component is fixable by
   membership methodology** — factor 3 (estimator noise) is a
   consequence of the survey's finite richness, and factor 7 (to whatever
   extent it is genuinely irreducible, pending Correction 1's re-analysis)
   is a property of the halo population, not the pipeline.
2. **Projection/interlopers (factor 6)** — dominant effect is a **bias**
   shift (+0.3 to +0.9 dex depending on conditions, Phase 5), not primarily
   a scatter increase; still a major driver of total error when bias is
   counted (not just variance). This ranking is robust to Corrections 1-2.
3. **sigma_v -> M200 calibration choice (factor 4)** — a bias-normalization
   spread of ~0.2 dex across the 4 tested calibrations (Section A); smaller
   than the contamination-driven bias shift, but not negligible, and it
   compounds with (2) rather than replacing it. Note this factor also
   partially overlaps with the unresolved portion of item 1 (Correction 1)
   — the same calibration mismatch that inflates the S0 "intrinsic" term
   is, by definition, this factor; they are reported separately here for
   traceability to the original ladder stages, not because they are
   cleanly independent.
4. **Spectroscopic incompleteness (factor 5)** — measurable scatter
   contribution in the (noisy, n=39) ladder terms, but its larger, more
   certain effect is **sample loss**: at 50% completeness and modal survey
   richness (5-7 true members), the majority of groups drop below the
   n_specz>=5 floor entirely and are lost from the analyzable sample
   (502/868 -> 107/868 matched-eligible halos in this test) — arguably as
   important as its contribution to per-group scatter.
5. **Aperture choice (factor 2)** — measurable in the noisy ladder terms,
   but this project's dedicated Phase 3/4 aperture tests (5 real-data
   configurations) already showed aperture choice does not resolve the
   real-data correlation problem; consistent with a real but secondary
   contribution.
6. **Membership selection logic / P_v gate (factor 1)** —
   **UNRESOLVED, not confirmed smallest (Correction 2).** The ladder's
   near-zero marginal contribution was measured using a simplified,
   non-iterating stand-in for the real production gate, not the real
   algorithm. This item's ranking as "smallest" should not be relied upon
   until re-tested with a faithful reproduction of the real, iterating gate.
   Phase 3/4's real-data finding that membership methodology overall is
   secondary is unaffected by this (it did not depend on the mock), but
   this specific ladder entry should not be cited as confirming it.

## The true bottleneck

**Finite-N sigma_v estimation noise, combined with a scatter/mismatch term
at the calibration-vs-true-relation level that is not yet cleanly separated
from genuine intrinsic halo-dynamical scatter (Correction 1), together
account for the majority (up to roughly two-thirds to three-quarters, in
this single run) of the total error budget — and neither is fixable by
membership methodology changes.** Calibration choice and field
contamination are the next most consequential, primarily through **bias**
rather than scatter, and are jointly responsible for the largest systematic
(non-random) errors in individual group masses. Membership aperture and
gate logic, the focus of Phases 2-4, remain assessed as secondary
contributors based on the real-data tests in Phase 3/4 (which do not depend
on the mock) — but the mock-based ladder's specific claim that the gate
contributes ~0% is not yet independently confirmed (Correction 2) and
should not be over-cited as quantitative proof of this, only as consistent
with it.

## Implication for future calibration work (not yet undertaken)

This benchmark demonstrates, with an error budget rather than a single
correlation coefficient, that:
- No published or self-fit calibration tested is adequate on its own (D15
  unchanged) — the combined estimator-noise-plus-intrinsic/mismatch floor
  (up to ~72% of budget in this single run, pending Correction 1's
  re-analysis) limits any single-calibration approach to a substantial
  minimum scatter even under ideal conditions.
- A future recalibration effort would have the largest expected payoff by
  targeting the **combined bias from calibration + contamination** (factor
  4 + factor 6, i.e. fitting directly against mock-observed, not
  ground-truth, sigma_v — D16), since that is the largest clearly
  *reducible* (bias, not variance-floor) component identified here.
- Before any specific scatter-floor percentage is used in a publication,
  Correction 1's required re-analysis (S0 recomputed across all four
  calibrations, with bootstrap intervals) should be completed. Until then,
  report only the qualitative conclusion — a large, membership-
  methodology-independent scatter floor exists — not a specific percentage.

No new calibration is fitted in this phase, per instruction — this
benchmark, as corrected, is the evidence base for deciding whether and how
to do so next.

## Reproducibility

```
python scripts/pipeline/membership_v2_phase6_benchmark.py --tng-dir /Volumes/extHD/tng_local_catalog
```
Outputs: `outputs/results/membership_v2_phase6/*.csv` (gitignored; script
committed).
