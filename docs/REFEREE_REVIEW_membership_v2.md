# Independent Referee Review — membership-v2 (Project B)

Reviewer stance: external referee for A&A/MNRAS, not the developer. Assumes
no prior conclusion in `membership_v2_phase2_3_report.md` through
`membership_v2_phase7_recommendation.md` and `DECISION_LOG_V2.md` (D11-D18)
is correct until checked. Goal: falsify the methodology if possible. No
code was modified to produce this review; two claims were checked directly
against the actual production script
(`determine_membership_dztier.py`) to verify fidelity claims made in the
mock lightcone documentation.

---

## 1. Are any conclusions stronger than the evidence supports?

**Yes, one materially.** The headline number carried into D17/D18 —
"intrinsic halo dynamics accounts for >=46% of final total variance" — is
computed from a **single ablation run using only the Munari et al. (2013)
calibration** as the fixed reference (`ablation_ladder`'s `cal_name =
"munari2013_galaxies"`, hardcoded). The S0 ("intrinsic") stage applies this
one calibration directly to TNG's own `veldisp_halo_kms`. But Phase 5's own
idealized-case table shows the **TNG self-calibrations fit this exact
relation far better** (bias -0.11/-0.12 dex vs. Munari's -0.23 dex). This
means the "46% floor" is not a clean measurement of *irreducible* intrinsic
scatter — it is entangled with **how well Munari's specific power-law
happens to match TNG's own M-sigma relation**, which the project's own
data show is imperfect. Re-running S0 with `tng_selfcal_all` as the
reference would very likely shrink the floor, and the true irreducible
component is *bounded above* by 46%, not established to equal it. The
report acknowledges a related but narrower caveat ("cannot be fully
separated... with this test alone") but the 46% figure is nonetheless
carried forward unqualified into two decision-log entries (D17, D18) and
the Phase 7 recommendation's headline framing. **This is the single most
important overclaim in the package.**

**Secondary instance**: Phase 3's claim that richness "correlates
significantly with X-ray M200 (r=0.64, p=1.6e-8) but with none of the
dynamical-mass configs" is accurately reported, but the *interpretive*
leap in Phase 4's Section 6 ("three independent lines of evidence
converge...") slightly overstates independence — two of the three lines
(the richness cross-check and the contamination-reduction test) both derive
from the same 63-112 group real-data sample, so they are not fully
statistically independent lines of evidence, only different measurements on
substantially overlapping data.

## 2. Are there hidden circular arguments?

**One confirmed, previously undocumented.** The Phase 5/6 mock's
membership-selection stage (`select_members()` in
`membership_v2_phase5_mock_lightcone.py`) applies a **single-pass,
fixed-sigma_v_prior=500 km/s** P_v gate. But the actual production
algorithm it claims to reproduce
(`determine_membership_dztier.py::refine_group_redshift`, verified by
direct inspection, lines 134-176) **iterates the gate up to 10 times**,
re-estimating `sigma_v_current = clip(std(members), 150, 2500)` from the
current member set each round and re-scoring against this updated,
group-specific value, converging when the redshift shift is below 25 km/s.
The mock's documentation states it "reproduces v1.0's actual production
membership gate... same functional form and same default constants" —
this is **only true of the first iteration**; the mock never reproduces the
convergence behavior that is the actual point of D1's probabilistic model.
This is not circular in the classic sense, but it is a **fidelity gap
between the claimed and actual test article**: the Phase 6 finding that
"the membership gate contributes ~0% marginal variance" (S5 in the ladder)
was measured on a simplified, non-iterating stand-in for the real gate, not
the real gate itself. The true marginal contribution of the real,
iterating gate is untested.

## 3. Are there statistical weaknesses?

**Yes, several, concentrated in Phase 6.**

- **n=39 core sample for the headline error-budget percentages.** The
  matched-sample ladder restricts to 39 TNG halos that survive every
  stage. At n=39, the standard error on a *variance* estimate is
  order sqrt(2/38) ~ 23% relative — yet the report presents point
  percentages (25.8%, 35.4%, -32.1%, 33.1%, ~0%, -8.2%) to three
  significant figures with no confidence intervals, and the two negative
  "marginal variance" entries are explicitly acknowledged as noise but not
  quantified with an interval that would let a reader judge whether *any*
  of the individual stage contributions (beyond the S0-S1 comparison) are
  distinguishable from zero.
- **No bootstrap or jackknife intervals anywhere in Phases 5-6**, a
  regression from the v1.0-cycle audit's own standard (which used
  bootstrap/jackknife throughout, per `technical_validation_report.md`).
  Given this project's own prior audit explicitly modeled good practice
  here, its absence in the newer, more consequential Phase 6 numbers is a
  real inconsistency in rigor across the project's history.
- **Catastrophic-failure threshold (|residual|>0.5 dex) is an unjustified,
  unvaried constant.** No sensitivity check on this threshold is reported;
  all catastrophic-fraction numbers throughout Phases 5-7 depend on it.
- **The single most consequential free parameter in the entire mock —
  the field-contamination injection rate (`n_field_per_true_member=0.5`)
  — was never varied.** Every Phase 5/6 result involving contamination
  (which is most of them, including the sign-flip finding and the D15
  calibration-ranking reversal) is conditioned on this one, stated-but-
  untested value. A referee would very likely ask for at minimum a 3-point
  sensitivity sweep (e.g. 0.2/0.5/1.0) before accepting the specific
  magnitude of the contamination-driven bias as more than illustrative.

## 4. Are any simulation assumptions unrealistic?

- **Field-interloper velocity model**: interlopers are given `dv_kms`
  drawn from `Uniform(-2500, +2500)` km/s, i.e. uniformly filling the
  production pipeline's *velocity-acceptance ceiling* by construction.
  This is explicitly labeled as a modeling choice in the Phase 5 docstring,
  which is good practice, but it is also somewhat self-fulfilling: a
  uniform distribution across exactly the window the gate can accept
  maximizes the fraction of injected interlopers that pass P_v, which
  will tend to produce a *larger* contamination-driven bias than a more
  physically motivated distribution (e.g. one peaked near dv=0, since real
  field galaxies at a given photo-z are not uniformly likely at any
  velocity offset — a Gaussian or LSS-correlation-motivated distribution
  would be more defensible). This should be treated as a plausible
  upper-bound stress test, not a calibrated estimate of the real
  contamination rate or its bias, and Phase 5/6's own text mostly respects
  this framing, but the specific numbers (e.g. "+0.3 to +0.9 dex bias")
  are being used somewhat more quantitatively in Phase 6/7 than this
  caveat supports.
- **`veldisp_halo_kms` (TNG's ground-truth dispersion used at S0) was never
  independently verified for its exact definition** (1D line-of-sight vs.
  3D; all subhalos vs. central+satellites only; mass-weighted or not).
  If this is a 3D dispersion while gapper (from `vel_z_kms` only) is
  inherently 1D/line-of-sight, part of the apparent "intrinsic scatter"
  at S0 could be a **systematic normalization mismatch** (a
  sqrt(3)-type factor or similar) rather than genuine astrophysical
  scatter — this compounds concern #1 above and was not checked.
  Given this is the anchor value for the single most important number in
  the whole Phase 6-7 conclusion, this is a significant, checkable gap.
- **Redshift error fixed at 150 km/s for all mock members regardless of
  magnitude, instrument, or redshift** — a simplification; real COSMOS
  spec-z precision varies by target and survey (DEIMOS/MOSFIRE/etc. have
  different precisions). Not sensitivity-tested.
- **TNG halo mass/redshift coverage was never checked against the real
  survey's actual n_specz/redshift distribution** before the Phase 6
  percentages were presented as descriptive of "this survey's regime" —
  this gap is self-identified in Phase 7 item 3 ("remaining work"), which
  is good practice, but it means the headline numbers in D17/D18 are
  currently *unverified as representative* of the real CW-All/CW-HCG
  catalog, a caveat that should arguably appear in the Phase 6/7 headline
  text itself, not only in a forward-looking to-do item.

## 5. Are there untested observational biases?

- **X-ray selection was tested for a contamination-score association in
  Phase 4 (T2, real data) but never incorporated into the Phase 5/6 mock**
  — the mock has no X-ray-flux-limited selection function applied at all,
  despite Phase 5's own docstring listing "X-ray selection" as one of the
  realistic effects to include. It was not implemented; the XRAY+SPECZ
  method-combination logic and its associated Malmquist/Eddington-bias
  risk (flagged as a Critical gap in `REFEREE_REVIEW.md`, the earlier
  v1.0-cycle audit) is entirely untested by this branch's mock lightcone,
  despite being explicitly in scope per the Phase 5 task description ("COSMOS
  observing strategy... X-ray selection").
- **No test of whether the real survey's own photo-z-driven candidate
  selection (the "search aperture" candidate pool before the P_v gate)
  has a redshift- or magnitude-dependent bias** that the mock's uniform
  interloper injection does not capture.

## 6. Are any literature comparisons incomplete or misleading?

**Yes — a stated user requirement was not actually fulfilled, only
substituted.** The task explicitly asked for "any group-specific
calibration available for lower-mass systems" as one of (at least) four
calibrations to inject. Phase 5 substituted a **self-fit TNG calibration**
for this slot instead of citing a genuine external group-scale literature
relation (e.g. a group-calibrated M-sigma relation from the group/cluster
literature), with the stated reason being "avoid citing an external
group-scale calibration without full confidence in its exact published
parameters." This is a defensible caution against fabricating a citation,
and it is transparently disclosed — but it means the deliverable does
**not** actually contain an independent, externally-calibrated,
group-specific test, which was the explicit request. A referee would
likely either require this be added properly (with a verified citation) or
require the substitution be stated more prominently as a scope limitation
in the headline recommendation (Phase 7), not only in Phase 5's methods
section.

## 7. Are there remaining implementation assumptions that were never validated?

- The membership-gate fidelity gap (Section 2 above) — the single-pass
  vs. iterating P_v gate discrepancy.
- The `veldisp_halo_kms` definition (Section 4 above).
- The fixed 0.5 field-contamination rate (Section 3 above).
- **`SIGMA_PRIOR_KMS=500.0` and `PROB_THRESHOLD=0.05` were confirmed by
  direct inspection to match production's actual defaults** (verified:
  `determine_membership_dztier.py` lines 109-135) — this specific check
  passes and should be credited as a correctly-validated assumption, not
  flagged as a gap.
- The mock's aperture-based sky projection uses TNG's pre-computed
  `RA_mock`/`DEC_mock` fields without independent verification of how
  these were generated (e.g., whether they preserve realistic angular
  clustering/covariance for a lightcone, or are a simplified per-snapshot
  projection) — not documented or checked in this branch.

## 8. Confidence classification of conclusions

| Conclusion | Confidence | Basis |
|---|---|---|
| Aperture choice is not the dominant driver of the real-data sigma_v-M200 problem (Phase 3/4) | **High** | Tested across 5 real-data configurations with consistent, reproducible non-significance; independently corroborated by the richness cross-check on real (not mock) data. |
| BGG selection is robust to aperture redefinition | **High** | Re-scored (not self-referential) test, 97-100% agreement, reproducible. |
| Gapper is comparatively more robust than biweight at this richness | **Moderate-High** | Real, reproducible mock-benchmark difference (0.55 vs 0.76 dex) but no significance interval computed, and dependent on the (unvalidated) contamination model. |
| Munari et al. (2013) should not be replaced with the 3 alternatives tested | **Moderate** | Consistent across idealized and contaminated mock cases in relative ranking terms, but absolute numbers depend on the untested contamination-rate parameter, and the "group-specific" comparison arm was not a genuine external calibration (Section 6). |
| Intrinsic halo dynamics account for >=46% of total variance | **Low-Moderate** | The specific number is calibration-choice-dependent (Section 1) and based on n=39 with no uncertainty interval (Section 3); the *qualitative* conclusion that intrinsic scatter is large and calibration-choice-independent is better supported than the specific 46% figure. |
| The membership gate contributes ~0% marginal variance | **Low** | Measured using a simplified, non-iterating stand-in for the real production gate (Section 2); untested whether the real, iterating gate behaves differently. |
| Individual-group M200,dyn should be treated as low-confidence | **High** | Robustly supported by every test in Phases 3-6 independently (real-data non-significance, mock scatter/catastrophic-failure rates, calibration-ranking instability) — this qualitative, conservative conclusion is the best-supported result in the whole package, appropriately so given it is also the most cautious one. |

## 9. Figures and tables a referee would request

- A **re-run of the S0/S1 ablation stage using each of the 4 calibrations**
  (not just Munari), to directly test/bound concern #1 (Section 1) — this
  is the single most important missing analysis in the package.
- A **sensitivity plot of the headline error-budget percentages vs. the
  field-contamination rate parameter** (e.g. 0.2/0.5/1.0/2.0 interlopers
  per true member), given how much rides on this untested constant.
- **Bootstrap confidence intervals** on every Phase 6 bias/scatter/
  catastrophic-fraction number, at minimum for the n=39 ladder.
- A **direct comparison plot of the mock's single-pass P_v gate vs. the
  real, iterating production gate** on a shared test sample, to quantify
  the fidelity gap identified in Section 2.
- A plot or table showing **TNG's n_specz/mass/redshift distribution
  overlaid on the real CW-All/CW-HCG survey's actual distribution**, to
  support (or correct) the implicit representativeness claim behind
  quoting Phase 6 percentages as descriptive of "this survey's regime."
- A **scatter plot of predicted vs. true log M200** for at least one
  calibration/condition combination (currently only summary statistics are
  reported anywhere in Phases 5-7; no recovery plot exists).

## 10. Additional experiments that would most strengthen the work

1. Re-run the S0 ladder stage with all 4 calibrations (addresses #1, #9;
   highest priority — cheap to do, directly bears on the headline number).
2. A contamination-rate sensitivity sweep (addresses #3, #4, #6).
3. Re-implement the mock's membership gate as a faithful, iterating
   reproduction of `refine_group_redshift`, and re-run the S5/S6 ladder
   stages (addresses #2, #7 — the current "gate contributes ~0%" claim
   should not be trusted until this is done).
4. Verify `veldisp_halo_kms`'s exact definition against TNG documentation
   or the enrichment script that produced it (addresses #4, #7).
5. Bootstrap/jackknife intervals on all Phase 6 headline numbers (addresses
   #3).
6. A genuine external group-scale sigma_v-M200 literature calibration,
   properly cited and verified, to replace or supplement the TNG
   self-calibration substitution (addresses #6).

---

## A. Executive summary

The membership-v2 investigation (Phases 1-7) is a substantially more
rigorous and more honestly self-critical piece of work than most
methodology validation efforts of this scope, and its central, cautious
conclusion — that individual-group dynamical masses are intrinsically
low-confidence at this survey's richness and should be used in
ensemble/statistical contexts rather than individually — is well
supported by multiple, largely independent lines of evidence and should
survive scrutiny. However, the specific **quantitative** headline claim
that intrinsic halo dynamics account for ">=46% of variance" is less
secure than its prominence in the decision log (D17, D18) suggests, for a
concrete, checkable reason (calibration-choice entanglement in the S0
measurement) and a statistical one (n=39, no uncertainty interval). A
second, previously undocumented fidelity gap — the mock's membership gate
does not actually reproduce the real, iterating production algorithm — 
undermines confidence in the specific "membership gate contributes ~0%"
sub-claim, though not the broader conclusion that membership methodology
is secondary. Neither issue invalidates the project's overall conservative
recommendation, but both should be resolved before the specific numbers are
used in a methodology paper.

## B. Strengths

- Consistent, cumulative, falsification-oriented investigation across
  Phases 2-7, with real, reproducible scripts at every stage (not just
  narrative claims).
- The final recommendation is explicitly conservative and resists the
  temptation to claim a "winning" configuration — a genuinely uncommon and
  creditable property for a methodology validation of this kind.
- Real ground truth (TNG simulation halos) used throughout Phase 5-6,
  not a synthetic population with assumed properties.
- Transparent, explicit labeling of modeling choices as modeling choices
  (e.g. the contamination-model caveat already present in Phase 5's own
  text) — the self-critical instinct is present throughout, even though
  this review finds it was not applied evenly to every subsequent claim.
- Decision log (D11-D18) provides full traceability of every methodological
  choice and its evidentiary basis.

## C. Weaknesses

- The 46% intrinsic-variance figure is calibration-choice-entangled and
  not yet isolated (Section 1, 9.1).
- The mock's membership gate does not faithfully reproduce the real,
  iterating production algorithm (Section 2).
- No uncertainty intervals on any Phase 6 headline statistic; n=39 core
  sample for the most consequential numbers (Section 3).
- The single most consequential mock parameter (contamination rate) was
  never varied (Section 3, 4).
- One of the four requested calibration categories (an external,
  group-specific literature relation) was not actually delivered, only
  substituted (Section 6).
- `veldisp_halo_kms`'s exact definition was never verified (Section 4, 7).

## D. Publication blockers

1. Re-run the S0 ladder stage across all 4 calibrations and report whether
   the ~46% figure is robust to this choice, or revise the headline claim
   to a range/qualitative statement if not.
2. Verify (or correct) the mock membership gate's fidelity to the real,
   iterating production algorithm before citing the "gate contributes ~0%"
   finding in any publication.
3. Verify `veldisp_halo_kms`'s definition against TNG/enrichment-script
   documentation, since it anchors the single most important number in the
   package.

## E. Nice-to-have improvements

- Contamination-rate sensitivity sweep.
- Bootstrap/jackknife intervals throughout Phase 6.
- A genuine, verified external group-scale calibration in place of the TNG
  self-calibration substitution.
- Recovery scatter plots (predicted vs. true M200) for at least the
  production calibration.
- Representativeness check of the TNG mock sample against the real
  survey's n_specz/mass/z distribution.

## F. Overall recommendation

**Needs major revisions** for the *specific quantitative claims* (D17's
46% figure and the membership-gate-contributes-~0% sub-claim) before they
can be cited in a methodology paper as currently stated — both rest on
identifiable, checkable, and moderately costly-to-fix issues (Section 1, 2)
rather than being fundamentally unsound.

The **qualitative, conservative conclusion** — that membership/aperture
methodology is not the dominant limitation, and that individual-group
dynamical masses should be treated as low-confidence and used in
ensemble/statistical contexts — is independently well-supported by
multiple lines of evidence robust to the weaknesses identified here, and
this framing (already adopted in Phase 7 per instruction) is **ready to
carry forward** into production catalog documentation and paper text
without further revision.
