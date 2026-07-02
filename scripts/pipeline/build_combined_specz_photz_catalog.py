#!/usr/bin/env python
"""
Build the combined galaxy catalog used for membership determination:
  - Prefer Webb_Specz_Feb2026.fits wherever a galaxy exists there (native
    COSMOS-Web catalog, highest priority).
  - For galaxies not in Webb_Specz_Feb2026 AND outside the COSMOS-Web
    footprint (true margin, not an internal gap/mask/chip-gap within the
    footprint), fall back to Full_COSMOS_Field_Specz_Photz_DR1.1.fits (built
    by build_full_cosmos_field_specz_photz.py from the COSMOS spec-z
    compilation DR1.1 + COSMOS2020, covering the full COSMOS field that
    COSMOS-Web sits within).

    Footprint test: a Delaunay triangulation over the full set of Webb galaxy
    positions. A full-field galaxy with no Webb positional match is only kept
    as supplement if it falls OUTSIDE this triangulated hull -- i.e. genuinely
    beyond the imaged COSMOS-Web area. Galaxies that are unmatched but still
    INSIDE the hull are internal gaps in Webb's own catalog (masking, chip
    gaps, quality cuts) and are deliberately NOT filled from COSMOS2020,
    since COSMOS-Web imaging *does* cover that position -- mixing in a
    different photometric system there would be inconsistent, not a genuine
    margin extension.

Every row is tagged `source_catalog` = 'Webb' or 'COSMOS_full_field' so it's
always possible to tell which galaxies came from which catalog downstream.

Usage:
    python build_combined_specz_photz_catalog.py
    python build_combined_specz_photz_catalog.py --match-radius-arcsec 1.0
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.spatial import Delaunay

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
WEBB_FITS = BASE_DIR / "data" / "specz" / "Webb_Specz_Feb2026.fits"
FULL_FIELD_FITS = BASE_DIR / "data" / "specz" / "Full_COSMOS_Field_Specz_Photz_DR1.1.fits"
OUTPUT_FITS = BASE_DIR / "data" / "specz" / "Webb_Specz_Feb2026_plus_COSMOS_field.fits"

FOOTPRINT_HULL_N_POINTS = 40000  # subsample for a fast, still-representative Delaunay hull


def main():
    parser = argparse.ArgumentParser(description="Combine Webb_Specz_Feb2026 (priority) with the full-COSMOS-field spec-z+photo-z catalog (true-margin fallback only)")
    parser.add_argument("--match-radius-arcsec", type=float, default=1.0,
                        help="Positional radius to consider a full-field galaxy already covered by Webb")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for footprint-hull subsampling")
    args = parser.parse_args()

    print(f"Loading Webb (priority): {WEBB_FITS}")
    webb = Table.read(WEBB_FITS).to_pandas()
    webb = webb.rename(columns={c: c for c in webb.columns})
    webb["source_catalog"] = "Webb"
    print(f"  {len(webb)} galaxies")

    print(f"Loading full-COSMOS-field compilation (fallback): {FULL_FIELD_FITS}")
    full_field = Table.read(FULL_FIELD_FITS).to_pandas()
    for col in full_field.select_dtypes(include=["object"]).columns:
        if full_field[col].map(lambda v: isinstance(v, bytes)).any():
            full_field[col] = full_field[col].str.decode("utf-8")
    full_field["source_catalog"] = "COSMOS_full_field"
    print(f"  {len(full_field)} galaxies")

    webb_coords = SkyCoord(ra=webb["RA"].values * u.deg, dec=webb["DEC"].values * u.deg)
    field_coords = SkyCoord(ra=full_field["RA"].values * u.deg, dec=full_field["DEC"].values * u.deg)
    _, sep2d, _ = field_coords.match_to_catalog_sky(webb_coords)
    already_in_webb = sep2d < args.match_radius_arcsec * u.arcsec

    print(f"\nFull-field galaxies already covered by Webb (within {args.match_radius_arcsec} arcsec): "
          f"{already_in_webb.sum()}/{len(full_field)}")

    unmatched = full_field[~already_in_webb].copy()

    print(f"Building Webb footprint hull (Delaunay over up to {FOOTPRINT_HULL_N_POINTS} Webb positions)...")
    webb_pts = webb[["RA", "DEC"]].dropna().values
    rng = np.random.default_rng(args.seed)
    if len(webb_pts) > FOOTPRINT_HULL_N_POINTS:
        webb_pts = webb_pts[rng.choice(len(webb_pts), FOOTPRINT_HULL_N_POINTS, replace=False)]
    hull = Delaunay(webb_pts)

    inside_footprint = hull.find_simplex(unmatched[["RA", "DEC"]].values) >= 0
    n_internal_gap = int(inside_footprint.sum())
    supplement = unmatched[~inside_footprint].copy()

    print(f"Unmatched full-field galaxies: {len(unmatched)}")
    print(f"  ...inside Webb footprint (internal gap, EXCLUDED -- not true margin): {n_internal_gap}")
    print(f"  ...outside Webb footprint (true margin, kept as supplement): {len(supplement)}")

    # Keep a common column schema: RA, DEC, zfin, dz, redshift_type, source_catalog
    webb_common = webb[["RA", "DEC", "zfin", "dz", "source_catalog"]].copy()
    webb_common["redshift_type"] = np.where(webb["dz"] <= 0.003, "spec-z", "photo-z")

    supplement_common = supplement[["RA", "DEC", "zfin", "dz", "redshift_type", "source_catalog"]].copy()

    combined = pd.concat([webb_common, supplement_common], ignore_index=True)
    print(f"\nCombined catalog: {len(combined)} galaxies "
          f"(Webb: {(combined['source_catalog']=='Webb').sum()}, "
          f"COSMOS_full_field supplement: {(combined['source_catalog']=='COSMOS_full_field').sum()})")
    print(combined.groupby(["source_catalog", "redshift_type"]).size())

    Table.from_pandas(combined).write(OUTPUT_FITS, format="fits", overwrite=True)
    combined.to_csv(OUTPUT_FITS.with_suffix(".csv"), index=False)
    print(f"\nSaved: {OUTPUT_FITS}")
    print(f"Saved: {OUTPUT_FITS.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
