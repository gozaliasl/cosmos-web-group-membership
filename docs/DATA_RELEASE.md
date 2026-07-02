# COSMOS-Web Group Catalog — Master Data Release v1.0

This is the official catalog release supporting Papers I and II. It merges:

- **Project A** (`cosmos-web-xray-igm`): X-ray group properties,
  `xray_pipeline_v1.1_production` (`outputs/release_v1.1`).
- **Project B** (`cosmos-web-group-membership`): membership, BGG selection,
  dynamical mass, `membership_pipeline_v1.0`
  (`outputs/release/membership_release_v1.0`), plus informational-only
  quality metadata from the Project B validation program (`membership-v2`).
- **COSMOSWEB2025** (`data/galaxy_catalog_photz/galaxy_catalog.fits`): the
  **primary** source for every galaxy-level photometric, SED, and
  morphological property (JWST-based 2025 COSMOS-Web catalog: LePhare
  SED fits, stellar mass/SFR/age, multi-band photometry, and real
  Sersic-index/effective-radius/axis-ratio structural fits).
- **COSMOS2020** (`data/galaxy_catalog_cosmos2020`): used **only** as the
  fallback for galaxies outside COSMOSWEB2025/JWST coverage
  (footprint-margin members), for stellar mass/absolute magnitude and a
  coarse morphology classifier.
- **Spec-z**: `data/specz/Webb_Specz_Feb2026_plus_COSMOS_field.fits` (the
  combined spec-z + photo-z dz-tier catalog membership was built from).
- **VLA** (VLA-COSMOS 3GHz counterpart catalog, Smolčić et al.) and
  **radio-AGN classification** (HLAGN/MLAGN/SFG).
- **AGN catalogs**: SED-AGN (fit against COSMOSWEB2025 photometry),
  X-ray-AGN (Chandra COSMOS-Legacy nuclear point-source match), VLA
  radio-AGN.

## Release contents

| Product | Rows (CW-All) | Rows (CW-HCG) |
|---|---|---|
| Master Group Catalog | 1678 | 912 |
| Master Galaxy Catalog | ~20,383 unique galaxies | ~1870 unique galaxies |

Files: `outputs/release/master_catalog_v1.0/Master_{Group,Galaxy}_Catalog_{CW_ALL,CW_HCG}.[csv,fits]`

CW-All and CW-HCG remain two separate catalogs (different group-finding
apertures, different target populations by design — see
`membership_aperture_methodology.md`), each internally self-consistent with
its own `Group_ID` namespace. They are not merged into a single
cross-catalog table.

## What every galaxy row contains

Group ID (primary, plus a list of all group memberships in this catalog if
more than one — see Limitations), BGG flag, X-ray group flag, radio galaxy
flag, VLA AGN flag, COSMOSWEB2025 (SED) AGN flag, COSMOS2020 AGN flag
(documented as unavailable — see `QUALITY_FLAGS.md`), X-ray AGN flag,
morphology (COSMOSWEB2025 Sersic index/effective radius/axis ratio,
primary; COSMOS2020 compactness class, fallback only), stellar mass
(COSMOSWEB2025/LePhare, with COSMOS2020 fallback where applicable),
redshift, and full spectroscopic/photo-z information (redshift type,
uncertainty).

## What every group row contains

X-ray properties (flux, luminosity, temperature, detection significance,
contamination flags — from `xray_release_v1.1`), dynamical properties
(sigma_v, R200_dyn, M200_dyn, method — from `membership_release_v1.0`),
richness (n_specz, n_photz), BGG summary (position, mass, magnitude,
redshift), AGN statistics (counts and fractions by AGN type), radio galaxy
statistics, and confidence flags (`dynamical_mass_confidence`,
`sigma_v_quality`, `BGG_stability`, `nspec_regime`).

## Version and provenance

- `membership_pipeline_v1.0` — unchanged production defaults: fixed 750 kpc
  search aperture, gapper sigma_v estimator, Munari et al. (2013) interim
  calibration.
- `xray_pipeline_v1.1_production` — includes the confirmed v1.1 bug fixes
  (duplicate config block, WMAP5→Planck18 rescaling, temperature-dependent
  ECF, background-quality flag).
- Project B's methodology validation (`membership-v2` branch) is
  **validation-only** — none of its experimental configurations (adaptive
  apertures, alternative calibrations) are used in this release.

See `COLUMN_DESCRIPTION.md` for the full column dictionary,
`QUALITY_FLAGS.md` for flag definitions and what is/isn't available, and
`RELEASE_NOTES.md` for version history and known issues.
