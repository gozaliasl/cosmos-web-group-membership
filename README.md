# COSMOS-Web Group Membership, Redshift Refinement & Velocity Dispersion

**Membership selection, spectroscopic redshift refinement, and velocity dispersion measurements for COSMOS-Web galaxy groups**

## Overview

- Galaxy membership probabilities using spec-z + photo-z combination
- Redshift refinement (3σ + Δv < 1000 km/s criterion)
- Velocity dispersion (VD) estimation via VRF estimator
- Comparison of photo-z vs spec-z membership quality tiers (dz = 0.001/0.002/0.003)

## Structure

```
cosmos-web-group-membership/
├── src/
│   ├── determine_membership.py    # Main membership pipeline
│   ├── membership_functions.py    # Core membership algorithms
│   ├── radius_optimizer.py        # Adaptive radius optimisation
│   ├── VRFEstimator.py            # Velocity dispersion estimator
│   └── membership_funcs_improved.py
├── scripts/pipeline/
│   ├── run_membership_pipeline.py
│   ├── determine_photoz_membership.py
│   ├── determine_specz_membership.py
│   ├── interactive_specz_dashboard.py
│   └── create_publication_plots.py
└── outputs/
```

## Data inputs (not in git — see ../data/)

- `../data/specz/Webb_Specz_Feb2026.fits` — spec-z catalog with quality flags (dz tiers)
- `../data/group-catalog/cosmos_web_groups_catalog_refined_z.fits`
- `../data/galaxy_catalog_photz/galaxy_catalog.fits`

## Spec-z quality tiers

| dz   | Quality | N sources |
|------|---------|-----------|
| 0.001 | Spec quality 3 & 4 (best) | 13,370 |
| 0.002 | Spec quality 2 | 8,132 |
| 0.003 | Spec quality 1 + photo-z consistent | 4,484 |
| >0.003 | Photometric | 206,142 |
