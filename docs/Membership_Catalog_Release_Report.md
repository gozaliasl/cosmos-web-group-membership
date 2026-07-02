# Membership Catalog Release Report — membership_release_v1.0

Pipeline: `membership_pipeline_v1.0` (unchanged; membership-v2/Project B
remains validation-only, per instruction).
Release tag: `membership_release_v1.0`
Build script: `scripts/pipeline/build_release_catalog.py`
Outputs: `outputs/release/membership_release_v1.0/`

## Contents

| File | Rows | Description |
|---|---|---|
| `CW_All_groups_release_v1.0.[csv,fits]` | 1678 | One row per CW-All group: membership summary, R200/M200/sigma_v/method, informational quality columns |
| `CW_All_members_release_v1.0.[csv,fits]` | 20856 | Per-galaxy candidate/member table for CW-All |
| `CW_HCG_groups_release_v1.0.[csv,fits]` | 912 | One row per CW-HCG group, same schema |
| `CW_HCG_members_release_v1.0.[csv,fits]` | 1879 | Per-galaxy candidate/member table for CW-HCG |
| `verification_report.json` | — | Machine-readable form of the validation summary below |

## Production methodology (unchanged)

- **Membership**: probabilistic P_z*P_v model (`determine_membership_dztier.py`
  / `iterative_membership_bgg.py`), fixed 750 kpc search aperture (CW-All) /
  compact-group-scaled aperture (CW-HCG), with iterative BGG-recentering.
- **BGG selection**: mass-floored hybrid rank score, soft centrality
  preference (`select_bgg.py`).
- **Dynamical mass**: **gapper** sigma_v estimator, **Munari et al. (2013)**
  interim sigma_v-M200 calibration (`compute_group_r200.py`,
  `cosmos-web-xray-igm`).
- **Search aperture**: **fixed 750 kpc** (CW-All), unchanged.
- **X-ray properties**: `xray_pipeline_v1.1_production`
  (`cosmos-web-xray-igm/outputs/release_v1.1`) — confirmed byte-identical
  in `M200_Temp_Msun` to the `xray_catalog.fits` input already used by
  `compute_group_r200.py`, so the method-combination catalog used here is
  consistent with the official v1.1 X-ray release.

No methodology was modified, recalibrated, or replaced to produce this
release. No adaptive-aperture or experimental configuration was used.

## Validation summary

### Group ID / member ID integrity

| Check | CW-All | CW-HCG |
|---|---|---|
| `Group_ID` duplicates | 0 | 0 |
| Duplicate `(Group_ID, RA, DEC)` member rows | 0 | 0 |
| Groups with null `r200_method` | 0 | 0 |
| n_specz recount vs. summary mismatches | 0 | 0 |

Group IDs and member-row identity are fully unique and internally
consistent in both catalogs.

### BGG assignments

| Check | CW-All | CW-HCG |
|---|---|---|
| Groups with exactly 1 BGG | 1669 | 788 |
| Groups with 0 BGG | 9 | 124 |
| Groups with >1 BGG (within a single group) | 0 | 0 |

No group has more than one BGG assigned. Groups with zero BGG lack any
member with both `mass_final` and `MR_final` populated (v1.0's documented
BGG-candidate requirement) — this is expected behavior, not an error, and
is more common in CW-HCG (13.6% of groups) than CW-All (0.5%) because
compact-group members more often lack full photometric coverage.

### Cross-group galaxy multiplicity (transparently reported, not "fixed")

| Check | CW-All | CW-HCG |
|---|---|---|
| Unique galaxies flagged as members | 20383 | 1870 |
| Galaxies that are members of >1 group | 473 (2.3%) | 9 (0.5%) |
| BGG galaxies that serve as BGG for >1 group | 10 | 1 |

**This is a known, inherited characteristic of v1.0's probabilistic
membership model, not a defect introduced by this release.** The P_z*P_v
gate is evaluated independently per group; a galaxy near two nearby groups
can satisfy both groups' probability thresholds simultaneously, and the
model does not enforce mutual exclusivity across groups. A small number of
these dual-membership galaxies are consequently selected as BGG in more
than one group. Per instruction, this release does not modify the
membership algorithm to force exclusivity — the finding is reported here so
downstream users are aware of it. See Limitations and Recommended Usage
below for how to treat this.

### Richness / spectroscopic counts

Recomputed directly from the member table and compared against the
production summary table's own counts: **0 mismatches in both catalogs.**
`n_members_specz` and `n_members_photoz` in the group table are confirmed
accurate.

### sigma_v, sigma_v uncertainty, R200_dyn, M200_dyn

| Check | CW-All | CW-HCG |
|---|---|---|
| Groups with SPECZ or XRAY+SPECZ method (dynamical mass expected) | 112 | 11 |
| ...missing sigma_v/sigma_v_err among those | 0 | 0 |
| R200_kpc missing where a method is assigned | 0 | 1 (Group_ID 171103) |
| M200_Msun missing where a method is assigned | 0 | 1 (Group_ID 171103) |

CW-HCG Group_ID 171103 is `STACKED_XRAY` method at z=3.821, outside the
redshift range covered by the stacked-flux bins
(`stacked_xray_lookup`/`assign_stacked_bin` in `compute_group_r200.py`).
This is a pre-existing edge case in the production pipeline (an
out-of-coverage redshift for the stacking bins), not introduced by this
release. It affects exactly one group in the combined 2590-group release.

### Method breakdown

| Method | CW-All | CW-HCG |
|---|---|---|
| STACKED_XRAY | 1133 | 604 |
| XRAY | 433 | 297 |
| XRAY+SPECZ | 63 | 7 |
| SPECZ | 49 | 4 |
| **Total** | **1678** | **912** |

## Informational quality columns (v2, metadata only)

Added per Project B's D18 recommendation. **These columns do not affect
any production membership, BGG, sigma_v, R200, or M200 value** — they are
computed via `membership_v2_dynamics.py`'s `fixed750` configuration, which
was verified (Project B, Phase 2) to exactly reproduce v1.0's own numbers
by construction.

| Column | Meaning |
|---|---|
| `dynamical_mass_confidence` | Categorical (n/a / low / medium), defaults to low for any spec-z-only dynamical mass, per Project B's central finding that individual-group dynamical masses are intrinsically low-confidence at this survey's richness. |
| `sigma_v_quality` | Categorical, reflects richness regime and aperture/iteration flags. |
| `nspec_regime` | "degradation-regime" (Nspec 8-14, where the real-data sigma_v-M200 correlation was found to invert) vs. "standard". |
| `phase_space_outlier_frac` | Fraction of a group's final members that are phase-space outliers (escape-velocity-envelope test). |
| `membership_quality` | Composite membership-confidence score (0-1) from the v2 toolkit. |
| `BGG_stability` | stable / mostly_stable / unstable — whether the BGG identity is unchanged when the group is independently re-trimmed to 3 alternative apertures (a genuine re-scored test, not a self-referential one). |

`dynamical_mass_confidence` distribution (CW-All): 1566 n/a (X-ray-derived,
not dynamical), 43 medium (X-ray-anchored, dynamical agreement confirmed),
69 low (spec-z-only, including the degradation-regime subset). `BGG_stability`
(CW-All): 1634 stable, 14 mostly_stable, 21 unstable (out of 1669 assigned
BGGs).

## Limitations

- Cross-group galaxy multiplicity (473 CW-All / 9 CW-HCG galaxies) means
  per-galaxy analyses that assume exclusive group membership should
  de-duplicate explicitly, choosing a resolution rule appropriate to their
  science case (e.g. highest `membership_prob`, or nearest `sep_kpc`) —
  this release does not impose one.
- One CW-HCG group (171103) lacks R200/M200 due to a stacking-bin
  redshift-coverage gap.
- Individual-group dynamical mass (`M200_Msun`/`R200_kpc` for SPECZ and
  XRAY+SPECZ-disagreement methods) carries substantial, quantified
  uncertainty beyond the formal statistical error reported — see Project
  B's `membership_v2_phase7_recommendation.md` and
  `Project_B_EXECUTIVE_SUMMARY.md` for the full evidence base. This is
  flagged per-group via `dynamical_mass_confidence`.
- The Munari et al. (2013) sigma_v-M200 calibration remains an interim
  placeholder (v1.0 Decision Log D7), unchanged by this release.
- `dynamical_mass_confidence`/`BGG_stability`/etc. are informational
  columns computed by a research toolkit (Project B), not part of the
  audited v1.0 production algorithm; treat them as guidance, not as a
  validated selection criterion.

## Confidence levels

| Quantity | Confidence |
|---|---|
| Group ID / member ID integrity, uniqueness | High (fully verified, 0 mismatches) |
| BGG assignment (per-group uniqueness) | High (0 groups with >1 BGG) |
| Richness / n_specz / n_photz counts | High (recomputed and matched exactly) |
| Membership (`is_member`) | High (v1.0-audited; see `technical_validation_report.md`) |
| M200/R200 for XRAY and STACKED_XRAY methods | High (X-ray-anchored, `xray_release_v1.1`) |
| M200/R200 for XRAY+SPECZ (agreement case) | Medium-High (X-ray-anchored, dynamical cross-check passed) |
| sigma_v, R200_dyn, M200_dyn for SPECZ-only groups | **Low** (per-group; see Project B) |
| Cross-group multiplicity handling | N/A — reported, not resolved; downstream responsibility |

## Recommended scientific usage

- **Safe for direct use**: group membership, richness, BGG identity (with
  awareness of the small cross-group-multiplicity subset), and M200/R200
  for XRAY and STACKED_XRAY method groups.
- **Use with the `dynamical_mass_confidence` flag**: any analysis touching
  M200_Msun/R200_kpc for SPECZ or XRAY+SPECZ-disagreement groups should
  filter or weight by this flag rather than treating all groups uniformly.
- **Prefer ensemble/statistical use over individual-group interpretation**
  for dynamical mass, sigma_v, and R200_dyn — this is the central,
  evidence-based recommendation from Project B and applies directly to
  this release's SPECZ-method subset (112/1678 CW-All, 11/912 CW-HCG
  groups).
- **De-duplicate cross-group-member galaxies explicitly** before any
  per-galaxy (not per-group) statistical analysis.
- Do not substitute any v2 experimental configuration (adaptive aperture,
  alternative calibration) for this release's values — none is a
  demonstrated improvement (Project B, D13/D15).

## Freeze

This release is frozen as `membership_release_v1.0`. No further
methodology modification, recalibration, or adaptive-aperture change is
included. Any future production update requires a new, explicitly
versioned release.
