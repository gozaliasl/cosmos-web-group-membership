#!/usr/bin/env python
"""
Fill stellar-mass / rest-frame-R-luminosity for members that came from the
COSMOS_full_field supplement (COSMOS2020, no JWST coverage) and therefore
can't match the Webb-native galaxy_catalog_photz.fits (JWST-only photometry).

Produces catalog-agnostic `mass_final` / `MR_final` / `Lr_final` columns:
  - Webb-sourced members: mass_final = LP_mass_med_PDF, MR_final = LP_MR_phys,
    Lr_final = LP_Lr (unchanged, straight passthrough).
  - COSMOS_full_field members without a Webb photz match: filled from
    COSMOS2020 CLASSIC's own LePhare fit (lp_mass_med, lp_MR), matched by
    RA/DEC. Lr_final is derived from lp_MR via L = 10^(-0.4*(M_R - M_sun,R))
    with M_sun,R = 4.42 (AB), so it's on the same log10(L/Lsun) scale as the
    Webb catalog's own LP_Lr.

Without this, BGG selection (select_bgg.py) and any stellar-mass analysis
silently drops border/margin members just because they lack JWST coverage,
even though COSMOS2020 has a perfectly usable independent mass/magnitude
estimate for them.

Usage:
    python fill_cosmos2020_properties.py --members-file <path with LP_mass_med_PDF/LP_MR_phys/LP_Lr columns>
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
COSMOS2020_FITS = BASE_DIR / "data" / "galaxy_catalog_cosmos2020" / "COSMOS2020_CLASSIC_R1_v2.2_p3.fits"
M_SUN_R_AB = 4.42  # AB rest-frame R-band solar absolute magnitude


def main():
    parser = argparse.ArgumentParser(description="Fill mass/luminosity for COSMOS_full_field members from COSMOS2020")
    parser.add_argument("--members-file", type=str, required=True)
    parser.add_argument("--match-radius-arcsec", type=float, default=1.0)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    members_path = Path(args.members_file)
    members = (Table.read(members_path).to_pandas() if members_path.suffix.lower() == ".fits"
               else pd.read_csv(members_path))
    print(f"Loaded {len(members)} members from {members_path.name}")

    for col in ["LP_mass_med_PDF", "LP_MR_phys", "LP_Lr", "source_catalog"]:
        if col not in members.columns:
            raise SystemExit(f"Required column '{col}' not found (run merge_photz_properties.py first)")

    members["mass_final"] = members["LP_mass_med_PDF"]
    members["MR_final"] = members["LP_MR_phys"]
    members["Lr_final"] = members["LP_Lr"]
    members["mass_source"] = np.where(members["LP_mass_med_PDF"].notna(), "Webb_LePhare", "none")

    needs_fill = members["mass_final"].isna() & (members["source_catalog"] == "COSMOS_full_field")
    n_needs_fill = int(needs_fill.sum())
    print(f"Members needing COSMOS2020 mass/magnitude fill: {n_needs_fill}")

    if n_needs_fill > 0:
        print(f"Loading COSMOS2020: {COSMOS2020_FITS}")
        t = Table.read(COSMOS2020_FITS)
        c20 = t["ALPHA_J2000", "DELTA_J2000", "lp_mass_med", "lp_MR"].to_pandas()
        c20 = c20.dropna(subset=["ALPHA_J2000", "DELTA_J2000"])

        c20_coords = SkyCoord(ra=c20["ALPHA_J2000"].values * u.deg, dec=c20["DELTA_J2000"].values * u.deg)
        fill_members = members[needs_fill]
        member_coords = SkyCoord(ra=fill_members["RA"].values * u.deg, dec=fill_members["DEC"].values * u.deg)
        idx, sep2d, _ = member_coords.match_to_catalog_sky(c20_coords)
        matched = np.asarray(sep2d < args.match_radius_arcsec * u.arcsec)

        matched_c20 = c20.iloc[idx].reset_index(drop=True)
        mass_fill = np.where(matched, matched_c20["lp_mass_med"].values, np.nan)
        mr_fill = np.where(matched, matched_c20["lp_MR"].values, np.nan)
        lr_fill = np.where(np.isfinite(mr_fill), 10 ** (-0.4 * (mr_fill - M_SUN_R_AB)), np.nan)
        lr_fill = np.log10(lr_fill, where=np.isfinite(lr_fill) & (lr_fill > 0),
                            out=np.full_like(lr_fill, np.nan))

        members.loc[needs_fill, "mass_final"] = mass_fill
        members.loc[needs_fill, "MR_final"] = mr_fill
        members.loc[needs_fill, "Lr_final"] = lr_fill
        members.loc[needs_fill & pd.Series(matched, index=fill_members.index), "mass_source"] = "COSMOS2020_LePhare"
        print(f"  Filled from COSMOS2020: {int(matched.sum())}/{n_needs_fill}")

    print(f"\nmass_source breakdown:\n{members['mass_source'].value_counts()}")
    print(f"mass_final still missing: {members['mass_final'].isna().sum()}/{len(members)}")

    out_stub = Path(args.output) if args.output else members_path.parent / f"{members_path.stem}_massfilled"
    members.to_csv(out_stub.with_suffix(".csv"), index=False)
    Table.from_pandas(members).write(out_stub.with_suffix(".fits"), format="fits", overwrite=True)
    print(f"\nSaved: {out_stub.with_suffix('.csv').name}, {out_stub.with_suffix('.fits').name}")


if __name__ == "__main__":
    main()
