#!/usr/bin/env python
"""
Build a full-COSMOS-field spec-z + photo-z catalog, using the same dz-tiering
convention as Webb_Specz_Feb2026.fits, but covering the *entire* COSMOS field
(COSMOS-Web is a sub-area within it) so it can supply real margin/border
coverage beyond the Webb footprint edge.

Sources:
  - data/speczcompilation/specz_compilation/specz_compilation_COSMOS_DR1.1_unique.fits
    (spectroscopic compilation incl. recent DESI spec-z; VVDS-style `flag` and
    `Confidence_level` per source; cross-matched IDs to COSMOS2020 Classic)
  - data/galaxy_catalog_cosmos2020/COSMOS2020_CLASSIC_R1_v2.2_p3.fits
    (photo-z + per-object error via lp_zPDF 68% interval, joined via
    Id_COS20_Classic -> COSMOS2020 ID)

Quality-tier convention (mirrors Webb_Specz_Feb2026.fits' own description):
  dz=0.001: highest spectroscopic quality (flag 3 or 4)
  dz=0.002: spectroscopic quality flag 2
  dz=0.003: flag 1 (lowest secure spec-z), only kept if consistent with photo-z
            (|specz - photoz| / (1+specz) < --specz-photz-consistency)
  otherwise: photometric redshift; dz = per-object photo-z error from COSMOS2020
             lp_zPDF (half of the 68% credible interval), falling back to
             the specz_compilation `photoz` value if not COSMOS2020-matched.

Usage:
    python build_full_cosmos_field_specz_photz.py
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
SPECZ_COMPILATION = (BASE_DIR / "data" / "speczcompilation" / "specz_compilation"
                     / "specz_compilation_COSMOS_DR1.1_unique.fits")
COSMOS2020_FITS = BASE_DIR / "data" / "galaxy_catalog_cosmos2020" / "COSMOS2020_CLASSIC_R1_v2.2_p3.fits"
OUTPUT_FITS = BASE_DIR / "data" / "specz" / "Full_COSMOS_Field_Specz_Photz_DR1.1.fits"

SPECZ_PHOTZ_CONSISTENCY = 0.05  # |specz-photoz|/(1+specz) threshold for flag==1 rescue
DEFAULT_PHOTZ_ERR_FRAC = 0.03   # fallback dz = this * (1+z) when no COSMOS2020 match


def main():
    print(f"Loading spec-z compilation: {SPECZ_COMPILATION}")
    comp = Table.read(SPECZ_COMPILATION).to_pandas()
    print(f"  {len(comp)} sources")

    print(f"Loading COSMOS2020 for photo-z errors: {COSMOS2020_FITS}")
    c20 = Table.read(COSMOS2020_FITS)[["ID", "lp_zPDF_l68", "lp_zPDF_u68"]].to_pandas()
    c20["photz_err"] = (c20["lp_zPDF_u68"] - c20["lp_zPDF_l68"]) / 2.0
    c20 = c20[np.isfinite(c20["photz_err"]) & (c20["photz_err"] > 0)]
    print(f"  {len(c20)} valid COSMOS2020 photo-z error entries")

    df = comp.merge(
        c20[["ID", "photz_err"]], left_on="Id_COS20_Classic", right_on="ID", how="left"
    )

    ra = df["ra_original"].values
    dec = df["dec_original"].values
    specz = df["specz"].values
    photoz = df["photoz"].values
    flag = df["flag"].values
    photz_err = df["photz_err"].values

    valid_specz = specz > -90
    valid_photoz = photoz > -90

    zfin = np.full(len(df), np.nan)
    dz = np.full(len(df), np.nan)
    source_flag = np.full(len(df), "photo-z", dtype=object)

    # Tier 1: flag 3 or 4 -> dz=0.001
    m = valid_specz & np.isin(flag, [3, 4])
    zfin[m] = specz[m]
    dz[m] = 0.001
    source_flag[m] = "spec-z_q3-4"

    # Tier 2: flag 2 -> dz=0.002
    m = valid_specz & (flag == 2) & np.isnan(zfin)
    zfin[m] = specz[m]
    dz[m] = 0.002
    source_flag[m] = "spec-z_q2"

    # Tier 3: flag 1, consistent with photo-z -> dz=0.003
    dz_norm_consistency = np.abs(specz - photoz) / (1 + specz)
    m = (valid_specz & valid_photoz & (flag == 1) & np.isnan(zfin)
         & (dz_norm_consistency < SPECZ_PHOTZ_CONSISTENCY))
    zfin[m] = specz[m]
    dz[m] = 0.003
    source_flag[m] = "spec-z_q1_photz_consistent"

    # Otherwise: photometric redshift, error from COSMOS2020 lp_zPDF if matched,
    # else a fixed fractional default
    m = np.isnan(zfin) & valid_photoz
    zfin[m] = photoz[m]
    has_err = m & np.isfinite(photz_err)
    dz[has_err] = photz_err[has_err]
    fallback = m & ~np.isfinite(photz_err)
    dz[fallback] = DEFAULT_PHOTZ_ERR_FRAC * (1 + photoz[fallback])
    source_flag[m] = "photo-z"

    out = pd.DataFrame({
        "RA": ra,
        "DEC": dec,
        "zfin": zfin,
        "dz": dz,
        "redshift_type": np.where(np.char.startswith(source_flag.astype(str), "spec-z"), "spec-z", "photo-z"),
        "source_detail": source_flag,
        "flag": flag,
        "Confidence_level": df["Confidence_level"].values,
        "survey": df["survey"].values,
        "GroupID": df["GroupID"].values,
    })
    out = out[np.isfinite(out["zfin"]) & np.isfinite(out["dz"]) & (out["zfin"] > 0)].reset_index(drop=True)

    print(f"\nFinal catalog: {len(out)} sources with valid zfin/dz")
    print(out["source_detail"].value_counts())
    print(f"\nspec-z: {(out['redshift_type']=='spec-z').sum()}, photo-z: {(out['redshift_type']=='photo-z').sum()}")

    OUTPUT_FITS.parent.mkdir(parents=True, exist_ok=True)
    Table.from_pandas(out).write(OUTPUT_FITS, format="fits", overwrite=True)
    out.to_csv(OUTPUT_FITS.with_suffix(".csv"), index=False)
    print(f"\nSaved: {OUTPUT_FITS}")
    print(f"Saved: {OUTPUT_FITS.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
