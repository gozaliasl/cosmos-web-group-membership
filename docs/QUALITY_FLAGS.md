# Quality Flags — Master Group and Galaxy Catalogs (v1.0)

## Galaxy-level flags

| Flag | Meaning | Reliability |
|---|---|---|
| `is_bgg` | Brightest group galaxy for its primary group | High — mass-floored hybrid rank algorithm, re-scored stability test passed (97-100% agreement under 3 independent aperture redefinitions) |
| `n_group_memberships > 1` | Galaxy is a member of more than one group | **Known, inherited, non-exclusive-by-design property of the probabilistic membership model** (473/20383 CW-All, 9/1870 CW-HCG galaxies). Not resolved automatically — analyses assuming exclusive membership must de-duplicate explicitly. |
| `morphology_source = "none"` | Neither COSMOSWEB2025 nor COSMOS2020 morphology matched | Treat morphology columns as missing, not zero |
| `X_ray_AGN`, `COSMOSWEB2025_AGN`, `VLA_AGN` | Independently derived AGN classifications | High — each from a distinct, established method (Chandra nuclear point-source match; SED fit; VLA 3GHz radio-excess/HLAGN/MLAGN classification) |
| `COSMOS2020_AGN` | **Always NaN** | **Not available.** No independent AGN classification column exists in `COSMOS2020_CLASSIC_R1_v2.2_p3.fits` in the data holdings used by this pipeline. This is reported as unavailable rather than approximated or fabricated. If a COSMOS2020-based AGN classification becomes available, this column should be populated in a future release, not silently left as a zero/False. |
| `radio_galaxy_flag` | Any VLA 3GHz counterpart (AGN or star-forming) | Medium — depends on VLA-COSMOS depth/completeness at the relevant flux limit; not independently re-validated here |

## Group-level flags

| Flag | Meaning | Reliability |
|---|---|---|
| `dynamical_mass_confidence` | n/a (X-ray-derived) / medium (X-ray+dynamical agreement) / low (spec-z-only) | **Central finding of Project B**: individual-group dynamical mass (SPECZ method) is intrinsically low-confidence at this survey's richness. Defaults to low for any spec-z-only group. See `membership_v2_phase7_recommendation.md` / `Project_B_EXECUTIVE_SUMMARY.md`. |
| `sigma_v_quality` | n/a (below floor) / low / medium | Reflects richness regime and aperture/iteration convergence flags; informational only, not a validated selection criterion |
| `nspec_regime` | standard / degradation-regime (Nspec 8-14) | The real-data richness range where the sigma_v-M200 correlation was found to invert (Project B, Phase 4); a HIGH member count in this range should NOT be read as higher reliability |
| `BGG_stability` | stable / mostly_stable / unstable | Whether BGG identity is unchanged when the group is independently re-trimmed to 3 alternative apertures; a genuinely re-scored test |
| `r200_disagreement` | X-ray and dynamical mass estimates disagree beyond combined uncertainty | When True, the X-ray value is adopted as the final `M200_Msun`/`R200_kpc` (already reflected in those columns; this flag documents why) |
| `Is_Projected_Contaminated`, `Contamination_Severity` | X-ray-side contamination flags (`xray_release_v1.1`) | See Project A's `RELEASE_VALIDATION_v1.1.md` for their derivation |
| `Is_Suspected_False_Positive` | X-ray-side flag | Groups with this True should be treated with caution for X-ray-derived quantities |

## Explicitly unavailable / not fabricated

- **`COSMOS2020_AGN`**: no source column exists; always NaN (see above).
- **True Sersic-index morphology for footprint-margin galaxies** (those
  outside COSMOSWEB2025/JWST coverage): only the coarse COSMOS2020
  `ACS_MU_CLASS` compactness classifier is available as a fallback; this is
  NOT equivalent to a Sersic fit and should not be used as one.
- **Completeness/purity of the real (non-mock) membership catalog**: no
  ground truth exists for the actual survey; only mock-lightcone-based
  estimates exist (Project B, Phase 5-6), and those apply to a simulated
  population, not directly to this catalog's real completeness/purity.

## How to use these flags

- For any analysis using **individual-group dynamical mass**: filter or
  weight by `dynamical_mass_confidence`; do not treat SPECZ-method M200 as
  equivalent in reliability to XRAY-method M200.
- For **per-galaxy statistics**: check `n_group_memberships`; de-duplicate
  using a rule appropriate to the science case (e.g. highest
  `membership_prob`) if exclusive membership is required.
- For **AGN demographics**: report `COSMOS2020_AGN` as "not available" in
  any AGN-fraction calculation that spans the full sample, rather than
  implicitly treating footprint-margin galaxies as non-AGN.
- For **morphology-dependent analyses**: filter on `morphology_source ==
  "COSMOSWEB2025"` if a genuine Sersic-based morphology is required; do not
  mix `morphology_class_ACS_MU_CLASS` (fallback) values into a
  Sersic-index-based analysis.
