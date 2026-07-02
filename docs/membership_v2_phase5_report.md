# membership-v2 — Phase 5 Report: mock lightcone + sigma_v-M200 calibration diagnostic

Branch: `membership-v2`. Implements D14 (`DECISION_LOG_V2.md`): a
calibration-diagnostic mock lightcone test, not a general realism exercise.
Script: `scripts/pipeline/membership_v2_phase5_mock_lightcone.py`.

## 1. What was built

- **Ground truth**: 868 TNG50-1+TNG100-1 halos with >=5 member subhalos
  (real simulation data, external to this repo, `/Volumes/extHD/tng_local_catalog`),
  spanning z=0-1.67+ and log M200 = 11.0-14.5. True M200 and true member
  peculiar velocities (`vel_z_kms`) come directly from the simulation, not
  re-derived.
- **Mock observational effects**, applied per group: spectroscopic
  incompleteness (30/50/70% retained), redshift uncertainty (150 km/s
  Gaussian, matching the production spec-z precision floor), and field
  contamination (interlopers from unrelated TNG halos at the same
  snapshot, positioned uniformly within the 750 kpc search aperture, with
  velocity offsets drawn uniformly across the production pipeline's actual
  velocity acceptance window).
- **Actual production membership gate** applied to the mock catalog: the
  same P_v probabilistic velocity cut, same functional form and same
  default constants (`SIGMA_PRIOR_KMS=500`, `PROB_THRESHOLD=0.05`) as
  `determine_membership_dztier.py`'s `apply_membership_cut()`.
- **Four sigma_v-M200 calibrations tested** against true M200: Munari et
  al. (2013) galaxy-tracer (current production default, D7), Evrard et al.
  (2008) dark-matter-particle (external, independent literature source),
  and two in-sample TNG self-calibrations (all halos; group-scale only,
  log M200<13.5) fit directly from this same TNG sample rather than citing
  an external group-scale value without full confidence in its exact
  published parameters.

731 mock (halo x completeness) observations were built and passed through
the full recovery chain.

**Scope limitation, added post-hoc per the independent referee review
(`REFEREE_REVIEW_membership_v2.md`, item 6)**: the task requested an
external, published, group-specific sigma_v-M200 calibration for
lower-mass systems as one of four calibration categories to test. This was
**not delivered**. `tng_selfcal_all` and `tng_selfcal_groups` are
**internal, in-sample fits to this project's own TNG data**, substituted
for that requirement to avoid citing an external group-scale relation
without full confidence in its exact published parameters (stated at the
time, but not flagged prominently enough in the original report). They
should be read as **internal alternatives/sensitivity checks, not as
literature replacements** — no genuine external group-specific calibration
has yet been tested against this pipeline. This is an open gap, not a
result.

## 2. Mass recovery results

### 2a. Idealized control (no contamination, no incompleteness, no z-error)

This isolates the calibrations' own accuracy from the observational-effects
model, using only the true member velocities at each halo's natural
richness:

| calibration | bias (dex) | scatter (dex) | catastrophic frac (\|resid\|>0.5 dex) |
|---|---|---|---|
| munari2013_galaxies | -0.233 | 0.509 | 0.383 |
| evrard2008_dm | -0.282 | 0.542 | 0.431 |
| tng_selfcal_all | -0.106 | 0.539 | 0.399 |
| tng_selfcal_groups | -0.121 | 0.518 | 0.367 |

**Even with zero injected observational contamination**, every calibration
under-predicts true M200 by 0.1-0.3 dex on median, with ~0.5 dex scatter
and 37-43% of individual groups landing more than 0.5 dex (>3x in mass)
from the truth. The two TNG self-calibrations (fit to this exact sample)
perform best, as expected, but still show large scatter and a >35%
catastrophic-failure rate. This scatter is not a calibration-choice
problem — it reflects the fundamental statistical noise of estimating
sigma_v (via gapper) from 5-15 real, individual (not idealized-Gaussian)
halo members, several of which are dynamically young / non-virialized in
the simulation itself.

By richness, idealized case (median bias, dex):

| n_true members | munari2013 | evrard2008 | tng_selfcal_all | tng_selfcal_groups |
|---|---|---|---|---|
| 5-7 | -0.17 | -0.22 | -0.05 | -0.05 |
| 8-9 | -0.18 | -0.23 | -0.05 | -0.06 |
| 10-14 | -0.31 | -0.34 | -0.17 | -0.20 |
| 15+ | -0.52 | -0.56 | -0.39 | -0.41 |

Counter-intuitively, bias *grows* with richness rather than shrinking —
richer TNG halos in this sample are more massive and more often
dynamically complex, and neither external calibration (fit on much larger,
more relaxed cluster samples) nor the TNG self-fit (a single power law
across the full mass range) captures this population's mass-dependent
behavior well. This is worth flagging as a real, non-obvious result, not
an artifact of the recovery method.

### 2b. Full mock (incompleteness + z-error + field contamination)

| calibration | bias (dex) | scatter (dex) | catastrophic frac | Pearson r(pred,true) |
|---|---|---|---|---|
| munari2013_galaxies | +0.443 | 0.543 | 0.512 | 0.677 |
| evrard2008_dm | +0.471 | 0.566 | 0.527 | 0.663 |
| tng_selfcal_all | +0.639 | 0.564 | 0.655 | 0.664 |
| tng_selfcal_groups | +0.576 | 0.549 | 0.596 | 0.673 |

Under the full mock (contamination + incompleteness + z-error), the bias
**flips sign** relative to the idealized case: all four calibrations now
substantially *over*-predict mass (field contamination inflates the
observed gapper sigma_v). **Munari (the current production default) shows
the smallest bias and lowest catastrophic-failure rate of the four under
this specific contamination model** — but this should not be read as
"Munari is the most accurate calibration": it is fit with a shallower
power-law index that happens to compress the contamination-inflated sigma_v
back down more than the steeper TNG self-fits do in this regime. All four
Pearson correlations between predicted and true mass are similar (r~0.66-0.68,
all highly significant, n=731) — the *linear ranking* of groups by mass is
recovered reasonably well by every calibration; it is the *absolute
normalization* that is unreliable, and differently so per calibration.

Trends common to all four calibrations (by mass, richness, redshift,
completeness — full tables in `outputs/results/membership_v2_phase5/`):
bias and catastrophic-failure rate are worst at the lowest true mass
(log M200<12.5: catastrophic fraction 72-88%), improve at higher mass
(log M200>13.5: catastrophic fraction 26-46%), worsen at higher redshift
(z>1: catastrophic fraction 62-79% vs 40-55% at z<0.5), and show only a
weak, non-monotonic dependence on the specific completeness fraction
tested (30/50/70% give similar results — contamination rate matters more
than completeness fraction in this model).

## 3. Interpretation

**Direct observation**: no tested calibration recovers true M200 with
useful accuracy under either the idealized or the realistic mock — bias
ranges from -0.3 to +0.6 dex and scatter is consistently ~0.5 dex across
all four. Under contamination, Munari performs least-badly of the four;
under idealized conditions, the TNG self-calibrations perform least-badly.

**Evidence-supported interpretation, not proven**: the dominant source of
sigma_v-M200 recovery error in this sample's regime (n_specz~5-15, group
to poor-cluster mass, z up to ~2) is not primarily "which published
calibration is used" — it is (a) intrinsic gapper-estimator variance at
low N, (b) genuine dynamical non-equilibrium in this simulated population
at these masses/redshifts, and (c) sensitivity of the *observed* sigma_v
to field contamination, which shifts the effective bias by a comparable
amount to the calibration-to-calibration spread itself. This is consistent
with, and now gives a mechanistic explanation for, the Phase 3/4 findings
that no membership/aperture configuration repaired the real-data
correlation: **membership methodology cannot fix a problem whose dominant
term is estimator variance plus true simulation-scale scatter, not
residual contamination alone.**

**Caveat on the contamination model itself**: the specific field-injection
recipe (interloper rate, velocity-offset distribution) is a stated modeling
choice, not a direct COSMOS measurement — the *sign flip* between idealized
and contaminated cases confirms contamination matters a great deal, but the
exact magnitude of the contaminated-case bias should not be over-read as a
precise prediction for the real survey; it is a diagnostic stress test.

## 4. Recommendation

**Should Munari remain?** Not as-is, unconditionally. It is not
demonstrably better than the alternatives tested in the idealized
(contamination-free) case — the TNG self-calibrations do noticeably better
there — but it does perform least-badly under the realistic contaminated
mock, largely by chance of its particular power-law index rather than by
demonstrated superior physical accuracy.

**Should it be replaced?** Not with any of the three alternatives tested
here as a direct swap — none is a clear, consistent improvement across
both the idealized and contaminated cases; each performs best under a
different (and equally plausible) subset of conditions.

**Should a new calibration be fitted?** Yes, but the evidence here argues
against fitting a single new power-law replacement (which is what
`tng_selfcal_all`/`tng_selfcal_groups` already represent, and which did not
solve the problem). A genuinely improved calibration would need to be
fitted **on mock-observed sigma_v (post-contamination, post-incompleteness)
rather than on true simulation sigma_v** — i.e. calibrated to the actual
observable, not the ground truth, which is a materially different
regression problem than any of the four tested here.

**Is a mass-dependent or redshift-dependent calibration required?** The
evidence says yes: bias and catastrophic-failure rate both degrade
substantially at low mass (log M200<12.5) and high redshift (z>1) for
every calibration tested, with no single global power law fitting the full
range well. A calibration with explicit mass and redshift terms (or
separate low-mass/high-z corrections) is better supported by this test than
retaining, or simply swapping, a single global power law.

**Overall**: keep Munari et al. (2013) as the *documented interim*
placeholder (unchanged from D7) rather than replacing it with any of the
three alternatives tested here. This is a decision to **retain an
acceptable interim default because no tested alternative is a demonstrated
improvement — not a finding that Munari is uniquely correct or optimal**
(only two internal TNG self-calibrations and one external cluster-scale
relation were compared against it; no external group-specific calibration
was tested, see the scope limitation above). This mock test provides direct
evidence that **the dynamical-mass branch of the pipeline (any calibration
tested) carries ~0.5 dex scatter and 40-65% catastrophic-failure risk for
individual groups in this survey's richness/mass/redshift regime**, and
that fixing this requires a mock-observed (not ground-truth) recalibration
with explicit mass/redshift dependence — a larger undertaking than swapping
one power law for another, and outside the scope that can be completed by
tuning membership/aperture methodology alone (consistent with D13).

## Reproducibility

```
python scripts/pipeline/membership_v2_phase5_mock_lightcone.py --tng-dir /Volumes/extHD/tng_local_catalog
```
Outputs: `outputs/results/membership_v2_phase5/*.csv` (gitignored; script is
committed).
