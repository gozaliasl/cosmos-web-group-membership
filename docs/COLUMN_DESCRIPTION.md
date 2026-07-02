# Column Dictionary — Master Group and Galaxy Catalogs (v1.0)

## Master Galaxy Catalog

| Column | Source | Description |
|---|---|---|
| `Catalog` | Project B | "CW-All" or "CW-HCG" |
| `Group_ID` | Project B | Primary group membership (see `n_group_memberships` for cross-group cases) |
| `RA`, `DEC` | COSMOSWEB2025-matched membership position | Sky position (deg) |
| `redshift` | Webb_Specz_Feb2026_plus_COSMOS_field.fits (`zfin`) | Adopted redshift |
| `redshift_uncertainty_dz` | same (`dz`) | Redshift uncertainty |
| `redshift_type` | same | `spec-z` or `photo-z` (dz-tier classification) |
| `is_bgg` | Project B (`select_bgg.py`) | BGG flag for the primary group |
| `is_bgg_any_group` | derived here | True if this galaxy is BGG in *any* of its group memberships (relevant only for the small cross-group-multiplicity subset) |
| `n_group_memberships` | derived here | Number of distinct groups (within this catalog) this galaxy is a member of |
| `all_group_ids_this_catalog` | derived here | Comma-separated list of all `Group_ID`s this galaxy belongs to |
| `membership_prob` | Project B (P_z*P_v model) | Membership probability used to select the primary group |
| `sep_kpc`, `dv` | Project B | Projected separation and velocity offset from the primary group center |
| `stellar_mass_final` | COSMOSWEB2025 (primary) / COSMOS2020 (fallback, footprint-margin) | Adopted stellar mass (log Msun); see `mass_source` |
| `stellar_mass_LePhare_med_PDF` | COSMOSWEB2025 LePhare fit | LePhare median-PDF stellar mass, direct (not fallback-merged) |
| `mass_source` | Project B | Which catalog `stellar_mass_final` came from |
| `MR_final` | COSMOSWEB2025 (primary) / COSMOS2020 (fallback) | Absolute R-band magnitude |
| `morphology_source` | derived here | "COSMOSWEB2025", "COSMOS2020_fallback", or "none" |
| `morphology_sersic_index`, `_err` | COSMOSWEB2025 (`SERSIC`) | Single-Sersic profile index, SE++ fit |
| `morphology_radius_arcsec`, `_err` | COSMOSWEB2025 (`RADIUS`) | Effective radius (arcsec) |
| `morphology_axis_ratio`, `_err` | COSMOSWEB2025 (`AXRATIO`) | Axis ratio (b/a) |
| `morphology_class_ACS_MU_CLASS` | COSMOS2020, fallback only | Coarse compactness-based star/extended-source class; **not** a Sersic morphology; populated only when COSMOSWEB2025 morphology is unmatched |
| `X_ray_AGN` | Chandra COSMOS-Legacy match | Nuclear X-ray point-source AGN |
| `COSMOSWEB2025_AGN` | SED fit against COSMOSWEB2025 photometry | SED-selected AGN |
| `VLA_AGN` | VLA 3GHz + radio classification | Radio-excess or (H/M)LAGN classified AGN |
| `radio_galaxy_flag` | VLA 3GHz counterpart match | Has any VLA 3GHz counterpart (AGN or star-forming) |
| `COSMOS2020_AGN` | **not available** | See `QUALITY_FLAGS.md` — no independent COSMOS2020 AGN classification exists in the data holdings used by this pipeline; always NaN |
| `any_AGN_flag` | derived here | True if any of the above AGN flags is True |
| `agn_flag` | Project B (`add_agn_flags.py`) | Combined categorical AGN classification string |
| `radio_type` | VLA classification | `SFG` / `MLAGN` / `HLAGN` / `unclassified` / `no_match` |
| `X_ray_group_flag` | derived here | True if this galaxy's primary group has an X-ray detection (method = XRAY or XRAY+SPECZ) |

## Master Group Catalog

| Column | Source | Description |
|---|---|---|
| `Group_ID`, `Catalog` | Project B | Identifiers |
| `search_radius_kpc` | Project B | Fixed 750 kpc (CW-All) / compact-group-scaled (CW-HCG) |
| `n_members_total`, `n_members_specz`, `n_members_photoz` | Project B | Richness |
| `z_refined`, `Group_Ra_final`, `Group_Dec_final` | Project B | Refined group center |
| `Is_Detected`, `Flux_erg_cm2_s`, `Luminosity_erg_s`, `Temperature_keV`, `Significance_Sigma`, `Background_Quality`, `Is_Projected_Contaminated`, `Contamination_Severity`, `Is_Suspected_False_Positive` | `xray_release_v1.1` | X-ray properties |
| `LAMBDA_STAR` | `xray_release_v1.1` (AMICO) | Independent richness estimate |
| `sigma_v_kms`, `sigma_v_err_kms` | Project B (`compute_group_r200.py`) | Gapper velocity dispersion, production default |
| `R200_kpc`, `R200_err_kpc`, `M200_Msun`, `M200_err_Msun` | Project B | Final adopted R200/M200 (method-combined) |
| `r200_method` | Project B | XRAY / SPECZ / STACKED_XRAY / XRAY+SPECZ |
| `r200_disagreement` | Project B | X-ray/dynamical mass disagreement flag |
| `BGG_RA`, `BGG_DEC`, `BGG_stellar_mass`, `BGG_MR`, `BGG_redshift`, `BGG_redshift_type` | derived here | BGG summary |
| `n_X_ray_AGN`, `n_COSMOSWEB2025_AGN`, `n_VLA_AGN`, `n_radio_galaxies`, `n_any_AGN`, `AGN_fraction`, `radio_galaxy_fraction` | derived here | AGN/radio statistics per group |
| `dynamical_mass_confidence`, `sigma_v_quality`, `nspec_regime`, `contamination_score`, `phase_space_outlier_frac`, `dynamical_reliability_score`, `membership_quality`, `BGG_stability` | Project B / `membership-v2` (informational metadata only) | Confidence flags — see `QUALITY_FLAGS.md` |

Units: masses in log10(Msun) unless noted `_Msun` (linear); fluxes in
erg/cm^2/s; luminosities in erg/s; temperatures in keV; radii in kpc unless
`_arcsec`.
