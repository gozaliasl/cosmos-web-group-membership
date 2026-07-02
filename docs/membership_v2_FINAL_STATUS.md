# membership-v2 (Project B) — Final Corrected Status

Closes out Phases 1-7 after the independent referee review
(`REFEREE_REVIEW_membership_v2.md`) and the resulting correction pass
(applied to Phase 5, Phase 6, Phase 7, D17, D18; see git history on
`membership-v2` for the exact diffs — no code was changed, documentation
only).

## High-confidence claims (safe to cite as-is)

- Membership aperture choice (fixed vs. R200,X-trim vs. adaptive R200,dyn)
  is not the dominant driver of the real-data sigma_v-M200 problem —
  tested across 5 configurations on real data, reproducible and
  consistent (Phase 3/4).
- BGG selection is robust to aperture redefinition (97-100% identity
  agreement under a genuinely re-scored, non-self-referential test)
  (Phase 3).
- Gapper is comparatively more robust than biweight at this survey's
  richness (0.55 vs. 0.76 dex scatter in the mock benchmark) (Phase 6A).
- Individual-group dynamical masses (M200,dyn, R200,dyn, sigma_v for
  SPECZ-only-method groups) should be treated as low-confidence and used
  in ensemble/statistical contexts rather than individually — supported
  independently by real-data non-significance (Phase 3/4), mock-based
  scatter and catastrophic-failure rates measured directly on n=731 mock
  observations (Phase 5), and calibration-ranking instability across
  conditions (Phase 5/6). This is the best-supported conclusion in the
  entire investigation.
- Finite-N sigma_v estimator sampling noise is a real, substantial,
  well-matched contributor to total error (+56% relative variance
  increase, S0->S1 in the Phase 6 ladder, using a consistent calibration
  choice across both stages).
- Field contamination's dominant effect is a bias shift (0.3-0.9 dex),
  not primarily a scatter increase (Phase 5, n=731).

## Downgraded claims (do not cite the specific numbers below without the stated caveat)

- ~~"Intrinsic halo dynamics accounts for >=46% of total variance"~~ —
  **withdrawn as a headline figure.** It was computed using only the
  Munari et al. (2013) calibration in the S0 ladder stage, entangling
  calibration-to-true-relation mismatch with genuine intrinsic scatter,
  and rests on n=39 with no uncertainty interval. Now stated only as an
  **upper bound** ("up to ~72% combined with estimator noise, in the
  single mock run performed to date"), pending re-analysis across all 4
  calibrations with bootstrap intervals.
- ~~"The membership gate contributes ~0% marginal variance"~~ —
  **marked as requiring re-validation.** Measured against a simplified,
  single-pass, fixed-sigma_v_prior=500 mock gate; the real production
  algorithm iterates up to 10 times, re-estimating each group's own
  sigma_v each round (confirmed by direct code inspection). The true
  marginal contribution of the real, iterating gate is untested.
- "Munari et al. (2013) is the best-performing calibration" — was never
  claimed this strongly, but any paraphrase in that direction should be
  corrected to: Munari is **retained as an acceptable interim default
  because no tested alternative demonstrates a consistent improvement**,
  not because it has been shown to be uniquely correct or optimal.

## Requires future validation (open items, not yet resolved)

1. Re-run the S0 ladder stage across all 4 calibrations with bootstrap
   intervals, to establish a defensible intrinsic-scatter-floor estimate
   (or range) before any specific percentage is used in a publication.
2. Re-implement the mock's membership gate as a faithful, iterating
   reproduction of the real production algorithm and re-run the affected
   ladder stages.
3. Verify `veldisp_halo_kms`'s exact definition (line-of-sight vs. 3D;
   which subhalos it includes) against TNG/enrichment-script
   documentation, since it anchors the disputed S0 stage.
4. Sensitivity sweep of the field-contamination injection rate
   (currently a fixed, untested constant).
5. Obtain and test a genuine, verified external group-specific
   sigma_v-M200 calibration from the literature — not yet done; the TNG
   self-calibrations used so far are internal alternatives only.
6. Mock-observed (not ground-truth) recalibration with explicit mass- and
   redshift-dependent terms (D16/D17), if items 1-2 confirm it is still
   warranted after the scatter floor is properly isolated.
7. Verify the TNG mock's richness/mass/redshift coverage against the real
   CW-All/CW-HCG survey's actual distributions before treating any Phase 6
   percentage as descriptive of the real catalog rather than a regime
   diagnostic.

## Net effect on the Project B recommendation

**No change.** The qualitative, conservative recommendation (Phase 7,
D18) — keep the v1.0 production defaults unchanged, adopt the v2
quality-metric toolkit as informational-only, add explicit low-confidence
flagging for individual-group dynamical mass, and do not fit a new
calibration yet — does not depend on either withdrawn number and is
retained in full. What has changed is that the two specific quantitative
figures behind part of the justification are now correctly labeled as
provisional rather than established, with a concrete path to resolving
them logged as future work.
