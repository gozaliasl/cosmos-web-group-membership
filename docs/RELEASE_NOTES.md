# Release Notes — Master Catalog v1.0

**Date**: 2026-07-02
**Status**: Official release for Papers I and II.

## What's in this release

The first combined release merging Project A (X-ray, `xray_release_v1.1`)
and Project B (membership/BGG/dynamical mass, `membership_release_v1.0`)
production catalogs with COSMOSWEB2025 photometric/SED/morphology
properties, COSMOS2020 fallback properties, VLA/radio classification, and
multi-method AGN flags, into Master Group and Master Galaxy catalogs for
CW-All and CW-HCG.

## Pipeline versions used

| Component | Version | Notes |
|---|---|---|
| Membership / BGG / dynamical mass | `membership_pipeline_v1.0` | Unchanged production defaults (fixed750 aperture, gapper, Munari et al. 2013) |
| X-ray properties | `xray_pipeline_v1.1_production` | Includes Bug #7 (duplicate config), #2 (WMAP5→Planck18 rescale), #4 (temperature-dependent ECF), #6 (background quality flag) |
| Photometry/SED/morphology | COSMOSWEB2025 (`galaxy_catalog.fits`) | Primary source for all galaxy-level properties |
| Fallback photometry | COSMOS2020 (`COSMOS2020_CLASSIC_R1_v2.2_p3.fits`) | Footprint-margin galaxies only |
| Spec-z / photo-z | `Webb_Specz_Feb2026_plus_COSMOS_field.fits` | Combined dz-tier catalog |
| Radio | VLA-COSMOS 3GHz counterpart catalog | HLAGN/MLAGN/SFG classification |
| Methodology validation | `membership-v2` (Project B) | **Validation only** — contributes informational quality columns; no experimental configuration used in production values |

## Known issues / limitations carried into this release

1. **Cross-group galaxy multiplicity** (473 CW-All / 9 CW-HCG galaxies
   belong to more than one group) is inherited from v1.0's non-exclusive
   probabilistic membership model. The Master Galaxy Catalog resolves this
   to one row per galaxy via a documented rule (highest `membership_prob`
   wins as the primary group); `n_group_memberships` and
   `all_group_ids_this_catalog` preserve the full information.
2. **One CW-HCG group** (Group_ID 171103, STACKED_XRAY, z=3.821) lacks
   R200/M200 due to a stacking-bin redshift-coverage gap.
3. **`COSMOS2020_AGN` is always NaN** — no independent COSMOS2020 AGN
   classification exists in the data holdings used by this pipeline.
4. **Individual-group dynamical mass (SPECZ method) is intrinsically
   low-confidence** at this survey's richness (Project B's central
   finding). Flagged via `dynamical_mass_confidence`, not removed from the
   catalog.
5. **Morphology** is a real Sersic-index/effective-radius/axis-ratio fit
   from COSMOSWEB2025 for galaxies within JWST coverage; footprint-margin
   galaxies only have the coarser COSMOS2020 `ACS_MU_CLASS` compactness
   classifier as a fallback (`morphology_source` documents which applies).

## What changed from the individual Project A / Project B releases

Nothing in the underlying X-ray, membership, BGG, or dynamical-mass values
was recalculated or modified. This release is a merge and aggregation
layer:

- Group-level: X-ray properties joined by `Group_ID`, BGG summary computed
  from the galaxy catalog, AGN/radio statistics aggregated per group.
- Galaxy-level: COSMOSWEB2025 SED/morphology properties and AGN/radio
  flags attached by positional cross-match (0.5 arcsec); COSMOS2020 used
  only as an explicit fallback.

## Recommended citation / usage note

Users should cite both `xray_release_v1.1` (Project A) and
`membership_release_v1.0` (Project B) as the underlying production
releases, and reference the Master Catalog v1.0 merge for the specific
file/column layout used in Papers I and II. See `QUALITY_FLAGS.md` for
required filtering before using individual-group dynamical mass or AGN
demographics that touch `COSMOS2020_AGN`.

## Reproducibility

```
python scripts/pipeline/build_release_catalog.py       # membership_release_v1.0
python scripts/pipeline/merge_photz_properties.py ...  # per catalog, run twice (CW-All, CW-HCG)
python scripts/pipeline/add_agn_flags.py ...            # per catalog, run twice
python scripts/pipeline/build_master_catalog.py         # final merge
```

Outputs: `outputs/release/master_catalog_v1.0/Master_{Group,Galaxy}_Catalog_{CW_ALL,CW_HCG}.[csv,fits]`

## Future work carried forward

See `Membership_Catalog_Release_Report.md` and
`Project_B_EXECUTIVE_SUMMARY.md` for the full prioritized list (mock-lightcone
recalibration diagnostics, external group-specific calibration testing,
etc.) — none of this is required to use the current release, but is
relevant context for any future version.
