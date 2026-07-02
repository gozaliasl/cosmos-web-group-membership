#!/usr/bin/env python
"""
Add AGN flags to a membership members file:

  - `is_sed_agn`: LePhare SED-fitting classification (LP_type == 2) from
    data/galaxy_catalog_photz/galaxy_catalog.fits, matched by RA/DEC. This is
    an SED-shape based AGN flag (not X-ray).
  - `is_xray_agn` / `xray_agn_sep_arcsec`: positional match (within
    --xray-match-radius-arcsec, default 1.5", the typical Chandra positional
    accuracy) to a nuclear X-ray point source in
    data/chandra_point_source_catalog/Chandra_COSMOS_Legacy_20151120_4d.fits.
    This identifies AGN via compact nuclear X-ray emission, distinct from the
    extended/diffuse group-scale X-ray emission measured by the group X-ray
    pipeline (cosmos-web-xray-igm).

Members flagged by either are also given a combined `agn_flag` summary string.

Usage:
    python add_agn_flags.py --members-file outputs/results/membership_dztier/fullfield_cw_all_members_*_with_photz.csv
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
CHANDRA_CATALOG = BASE_DIR / "data" / "chandra_point_source_catalog" / "Chandra_COSMOS_Legacy_20151120_4d.fits"
VLA_CATALOG = BASE_DIR / "data" / "radio-catalog" / "VLA_3GHz_counterpart_array_20170210_paper_smolcic_et_al.fits"


def _decode_bool_col(series: pd.Series) -> pd.Series:
    """VLA catalog stores true/false as padded byte-strings, not real booleans."""
    def _to_bool(v):
        if pd.isna(v):
            return False
        if isinstance(v, bytes):
            v = v.decode("utf-8")
        return str(v).strip().lower() == "true"
    return series.apply(_to_bool)


def add_sed_agn_flag(members: pd.DataFrame, match_radius_arcsec: float) -> pd.DataFrame:
    members = members.copy()
    if "LP_type" in members.columns:
        members["is_sed_agn"] = members["LP_type"] == 2
        return members

    print(f"  LP_type not present in members file; matching against {PHOTZ_CATALOG.name}")
    t = Table.read(PHOTZ_CATALOG)
    photz = t["RA_MODEL", "DEC_MODEL", "LP_type"].to_pandas().dropna(subset=["RA_MODEL", "DEC_MODEL"])
    photz_coords = SkyCoord(ra=photz["RA_MODEL"].values * u.deg, dec=photz["DEC_MODEL"].values * u.deg)
    member_coords = SkyCoord(ra=members["RA"].values * u.deg, dec=members["DEC"].values * u.deg)
    idx, sep2d, _ = member_coords.match_to_catalog_sky(photz_coords)
    matched = sep2d < match_radius_arcsec * u.arcsec

    lp_type = photz.iloc[idx]["LP_type"].values.astype(float)
    lp_type[~matched] = np.nan
    members["LP_type"] = lp_type
    members["is_sed_agn"] = lp_type == 2
    return members


def add_xray_agn_flag(members: pd.DataFrame, match_radius_arcsec: float) -> pd.DataFrame:
    members = members.copy()
    print(f"  Matching against nuclear X-ray point sources: {CHANDRA_CATALOG.name}")
    t = Table.read(CHANDRA_CATALOG)
    chandra = t["RA_x", "DEC_x", "flux_F", "flux_H", "flux_S"].to_pandas().dropna(subset=["RA_x", "DEC_x"])
    chandra_coords = SkyCoord(ra=chandra["RA_x"].values * u.deg, dec=chandra["DEC_x"].values * u.deg)
    member_coords = SkyCoord(ra=members["RA"].values * u.deg, dec=members["DEC"].values * u.deg)
    idx, sep2d, _ = member_coords.match_to_catalog_sky(chandra_coords)
    matched = sep2d < match_radius_arcsec * u.arcsec

    members["is_xray_agn"] = matched
    members["xray_agn_sep_arcsec"] = np.where(matched, sep2d.arcsec, np.nan)
    members["xray_agn_flux_full_erg_cm2_s"] = np.where(matched, chandra.iloc[idx]["flux_F"].values, np.nan)
    print(f"  X-ray AGN matches: {matched.sum()}/{len(members)}")
    return members


# VLA catalog's own true/false classification columns (Smolcic et al. 2017)
_VLA_BOOL_COLUMNS = [
    "Radio_excess", "HLAGN", "MLAGN", "Xray_AGN", "MIR_AGN", "SED_AGN",
    "Quiescent_MLAGN", "SFG", "Clean_SFG",
]


def add_vla_radio_agn_flag(members: pd.DataFrame, match_radius_arcsec: float) -> pd.DataFrame:
    members = members.copy()
    print(f"  Matching against VLA 3GHz counterparts: {VLA_CATALOG.name}")
    t = Table.read(VLA_CATALOG)
    cols = ["RA_CPT_J2000", "DEC_CPT_J2000", "FLUX_INT_3GHz", "Z_BEST"] + _VLA_BOOL_COLUMNS
    vla = t[cols].to_pandas().dropna(subset=["RA_CPT_J2000", "DEC_CPT_J2000"])
    for col in _VLA_BOOL_COLUMNS:
        vla[col] = _decode_bool_col(vla[col])

    vla_coords = SkyCoord(ra=vla["RA_CPT_J2000"].values * u.deg, dec=vla["DEC_CPT_J2000"].values * u.deg)
    member_coords = SkyCoord(ra=members["RA"].values * u.deg, dec=members["DEC"].values * u.deg)
    idx, sep2d, _ = member_coords.match_to_catalog_sky(vla_coords)
    matched = np.asarray(sep2d < match_radius_arcsec * u.arcsec)
    matched_vla = vla.iloc[idx].reset_index(drop=True)

    members["radio_is_vla_match"] = matched
    # All VLA-derived columns share the `radio_` prefix for easy grouping/grepping.
    # VLA's own X-ray/MIR/SED AGN classification (Smolcic et al.) is kept separate
    # from our independently-derived is_xray_agn / is_sed_agn columns.
    rename = {
        "Radio_excess": "radio_is_excess_agn",
        "HLAGN": "radio_is_hlagn",
        "MLAGN": "radio_is_mlagn",
        "Xray_AGN": "radio_is_xray_agn",
        "MIR_AGN": "radio_is_mir_agn",
        "SED_AGN": "radio_is_sed_agn",
        "Quiescent_MLAGN": "radio_is_quiescent_mlagn",
        "SFG": "radio_is_sfg",
        "Clean_SFG": "radio_is_clean_sfg",
    }
    for src, dst in rename.items():
        members[dst] = matched & matched_vla[src].values
    members["radio_flux_int_3ghz"] = np.where(matched, matched_vla["FLUX_INT_3GHz"].values, np.nan)
    members["radio_z_best"] = np.where(matched, matched_vla["Z_BEST"].values, np.nan)

    # Single categorical radio classification (VLA's own scheme, priority order
    # follows Smolcic et al. 2017: HLAGN > MLAGN > quiescent MLAGN > SFG)
    conditions = [
        matched & matched_vla["HLAGN"].values,
        matched & matched_vla["MLAGN"].values,
        matched & matched_vla["Quiescent_MLAGN"].values,
        matched & (matched_vla["SFG"].values | matched_vla["Clean_SFG"].values),
    ]
    choices = ["HLAGN", "MLAGN", "Quiescent_MLAGN", "SFG"]
    members["radio_type"] = np.select(conditions, choices, default=np.where(matched, "unclassified", "no_match"))

    n_radio_agn = (members["radio_is_excess_agn"] | members["radio_is_hlagn"]).sum()
    print(f"  VLA counterpart matches: {matched.sum()}/{len(members)}; radio-AGN (excess or HLAGN): {n_radio_agn}")
    print(f"  radio_type breakdown:\n{members['radio_type'].value_counts()}")
    return members


def main():
    parser = argparse.ArgumentParser(description="Add SED-AGN, X-ray-AGN, and VLA radio-AGN flags to a membership members file")
    parser.add_argument("--members-file", type=str, required=True)
    parser.add_argument("--sed-match-radius-arcsec", type=float, default=0.5)
    parser.add_argument("--xray-match-radius-arcsec", type=float, default=1.5,
                        help="Typical Chandra positional accuracy for nuclear point-source AGN matching")
    parser.add_argument("--vla-match-radius-arcsec", type=float, default=1.0,
                        help="Match radius against VLA 3GHz optical/IR counterpart position")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    members_path = Path(args.members_file)
    members = (Table.read(members_path).to_pandas() if members_path.suffix.lower() == ".fits"
               else pd.read_csv(members_path))
    print(f"Loaded {len(members)} members from {members_path.name}")

    members = add_sed_agn_flag(members, args.sed_match_radius_arcsec)
    members = add_xray_agn_flag(members, args.xray_match_radius_arcsec)
    members = add_vla_radio_agn_flag(members, args.vla_match_radius_arcsec)

    is_radio_agn = members["radio_is_excess_agn"] | members["radio_is_hlagn"]
    members["agn_flag"] = np.select(
        [members["is_sed_agn"] & members["is_xray_agn"] & is_radio_agn,
         members["is_xray_agn"] & is_radio_agn,
         members["is_sed_agn"] & is_radio_agn,
         members["is_sed_agn"] & members["is_xray_agn"],
         members["is_xray_agn"],
         members["is_sed_agn"],
         is_radio_agn],
        ["SED+X-ray+Radio_AGN", "X-ray+Radio_AGN", "SED+Radio_AGN",
         "SED+X-ray_AGN", "X-ray_AGN", "SED_AGN", "Radio_AGN"],
        default="none",
    )

    print(f"\nAGN flag breakdown:\n{members['agn_flag'].value_counts()}")

    out_stub = Path(args.output) if args.output else members_path.parent / f"{members_path.stem}_agn"
    members.to_csv(out_stub.with_suffix(".csv"), index=False)
    Table.from_pandas(members).write(out_stub.with_suffix(".fits"), format="fits", overwrite=True)
    print(f"\nSaved: {out_stub.with_suffix('.csv').name}, {out_stub.with_suffix('.fits').name}")


if __name__ == "__main__":
    main()
