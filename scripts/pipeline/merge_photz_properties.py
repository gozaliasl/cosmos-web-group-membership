#!/usr/bin/env python
"""
Enrich membership member catalogs with galaxy properties from the full
COSMOS-Web photometric catalog (data/galaxy_catalog_photz/galaxy_catalog.fits):
stellar mass, multi-band magnitudes, and LePhare rest-frame absolute magnitudes
(including LP_MR_phys, the interpolated rest-frame R magnitude).

Matches by RA/DEC position (Webb_Specz_Feb2026 galaxies against
RA_MODEL/DEC_MODEL in the photo-z catalog) since the two catalogs don't share
a common ID. Unmatched members (positional offset > --match-radius-arcsec)
keep NaN for the new columns and are flagged via `photz_matched`.

Usage:
    python merge_photz_properties.py --members-file ../outputs/results/membership_dztier/cw_all_members_*.csv
    python merge_photz_properties.py --members-file <path> --match-radius-arcsec 0.5
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
PHOTZ_CATALOG = BASE_DIR / "data" / "galaxy_catalog_photz" / "galaxy_catalog.fits"

# All MAG_MODEL_* / MAG_ERR_MODEL_* bands available in the catalog (JWST + HST + ground-based)
_MAG_BANDS = [
    "F115W", "F150W", "F277W", "F444W", "HST-F814W", "F770W", "CFHT-u",
    "HSC-g", "HSC-r", "HSC-i", "HSC-z", "HSC-y", "HSC-NB0816", "HSC-NB0921", "HSC-NB1010",
    "UVISTA-Y", "UVISTA-J", "UVISTA-H", "UVISTA-Ks", "UVISTA-NB118",
]
_MAG_COLUMNS = [f"MAG_MODEL_{b}" for b in _MAG_BANDS] + [f"MAG_ERR_MODEL_{b}" for b in _MAG_BANDS]

# Columns to pull from the photo-z catalog into the membership catalog: full
# LePhare SED-fitting output (redshift/type/quality, stellar pop. params with
# uncertainties, rest-frame luminosities/absolute mags), Chandra X-ray flag,
# and every available magnitude band + its error.
PROPERTY_COLUMNS = [
    "LP_zfinal", "LP_warn_fl", "LP_type",
    "LP_zPDF", "LP_zPDF_l68", "LP_zPDF_u68", "LP_zChi2", "LP_chi2_best", "LP_NbFilt",
    "LP_zp_AGN", "LP_chi2_agn", "LP_mod_agn", "LP_mod_star", "LP_chi_star",
    "LP_mod_minchi2_phys", "LP_ebv_minchi2", "LP_law_minchi2",
    "LP_age_minchi2", "LP_age_l68_PDF", "LP_age_med_PDF", "LP_age_u68_PDF",
    "LP_mass_minchi2", "LP_mass_med_PDF", "LP_mass_l68_PDF", "LP_mass_u68_PDF",
    "LP_sfr_minchi2", "LP_sfr_med_PDF", "LP_sfr_l68_PDF", "LP_sfr_u68_PDF",
    "LP_ssfr_minchi2", "LP_ssfr_med_PDF", "LP_ssfr_l68_PDF", "LP_ssfr_u68_PDF",
    "LP_Lnuv", "LP_Lr", "LP_Lk",
    "LP_MNUV_phys", "LP_MR_phys", "LP_MJ_phys", "LP_MK_phys",
    "FLAG_CHANDRA",
] + _MAG_COLUMNS


def load_photz_catalog() -> pd.DataFrame:
    print(f"Loading photo-z property catalog: {PHOTZ_CATALOG}")
    cols = ["RA_MODEL", "DEC_MODEL"] + PROPERTY_COLUMNS
    t = Table.read(PHOTZ_CATALOG)
    available = [c for c in cols if c in t.colnames]
    missing = set(cols) - set(available)
    if missing:
        print(f"  Warning: columns not found, skipping: {missing}")
    df = t[available].to_pandas()
    print(f"  Loaded {len(df)} galaxies, {len(available)-2} property columns")
    return df


def merge_properties(members: pd.DataFrame, photz: pd.DataFrame, match_radius_arcsec: float) -> pd.DataFrame:
    valid = photz.dropna(subset=["RA_MODEL", "DEC_MODEL"])
    photz_coords = SkyCoord(ra=valid["RA_MODEL"].values * u.deg, dec=valid["DEC_MODEL"].values * u.deg)
    member_coords = SkyCoord(ra=members["RA"].values * u.deg, dec=members["DEC"].values * u.deg)

    idx, sep2d, _ = member_coords.match_to_catalog_sky(photz_coords)
    matched = sep2d < match_radius_arcsec * u.arcsec

    out = members.copy()
    out["photz_matched"] = matched
    out["photz_sep_arcsec"] = sep2d.arcsec

    prop_cols = [c for c in PROPERTY_COLUMNS if c in valid.columns]
    matched_props = valid.iloc[idx].reset_index(drop=True)[prop_cols]
    for col in prop_cols:
        vals = matched_props[col].values.astype(float)
        vals[~matched.value if hasattr(matched, "value") else ~matched] = np.nan
        out[col] = vals

    print(f"  Matched {matched.sum()}/{len(members)} members within {match_radius_arcsec} arcsec")
    return out


def main():
    parser = argparse.ArgumentParser(description="Merge COSMOS-Web photo-z catalog properties into a membership members file")
    parser.add_argument("--members-file", type=str, required=True,
                        help="Path to a membership *_members_*.csv or .fits file")
    parser.add_argument("--match-radius-arcsec", type=float, default=0.5,
                        help="Positional match radius in arcsec")
    parser.add_argument("--output", type=str, default=None,
                        help="Output path (default: <input>_with_photz.<ext> next to input)")
    args = parser.parse_args()

    members_path = Path(args.members_file)
    members = (Table.read(members_path).to_pandas() if members_path.suffix.lower() == ".fits"
               else pd.read_csv(members_path))
    print(f"Loaded {len(members)} members from {members_path.name}")

    photz = load_photz_catalog()
    enriched = merge_properties(members, photz, args.match_radius_arcsec)

    if args.output:
        out_stub = Path(args.output)
    else:
        out_stub = members_path.parent / f"{members_path.stem}_with_photz"

    enriched.to_csv(out_stub.with_suffix(".csv"), index=False)
    Table.from_pandas(enriched).write(out_stub.with_suffix(".fits"), format="fits", overwrite=True)
    print(f"Saved: {out_stub.with_suffix('.csv').name}, {out_stub.with_suffix('.fits').name}")


if __name__ == "__main__":
    main()
