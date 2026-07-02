#!/usr/bin/env python
"""
Build the final science-ready Master Group Catalog and Master Galaxy
Catalog for CW-All and CW-HCG, merging:
  - Project A production (X-ray properties, xray_release_v1.1)
  - Project B production (membership_release_v1.0: membership, BGG,
    dynamical mass, informational quality flags)
  - COSMOSWEB2025 (the JWST-based COSMOS-Web photometric/SED/morphology
    catalog, data/galaxy_catalog_photz/galaxy_catalog.fits) -- the PRIMARY
    source for every galaxy property below (photo-z/SED properties via
    merge_photz_properties.py; morphology -- Sersic index, effective
    radius, axis ratio -- via a direct positional cross-match performed in
    this script, since merge_photz_properties.py's column allowlist does
    not currently include the structural-fit columns)
  - AGN/radio catalogs (SED-AGN from COSMOSWEB2025, X-ray-AGN, VLA 3GHz
    radio-AGN/SFG, via add_agn_flags.py -- already run on the release
    members table)
  - COSMOS2020 -- used ONLY as the fallback for galaxies outside
    COSMOSWEB2025/JWST coverage (footprint-margin members; mass_source
    already reflects this upstream, per fill_cosmos2020_properties.py).
    COSMOS2020 morphology (ACS_MU_CLASS, a coarse compactness-based
    star/extended-source classifier, not a Sersic fit) is used only where
    COSMOSWEB2025 morphology is unmatched. NO independent COSMOS2020 AGN
    classification exists in the data holdings available to this pipeline
    -- see QUALITY_FLAGS.md, this is documented as unavailable, not
    fabricated.

Redshift/spectroscopic provenance: `redshift` (zfin), `redshift_uncertainty_dz`
(dz), and `redshift_type` (spec-z vs. photo-z dz-tier) already come, upstream
of this script, from `data/specz/Webb_Specz_Feb2026_plus_COSMOS_field.fits`
(built by `build_combined_specz_photz_catalog.py`) -- the default input of
`iterative_membership_bgg.py`, which produced the membership_release_v1.0
member tables this script reads. This is unchanged here; noted for
provenance clarity only.

Produces ONE Master Group Catalog and ONE Master Galaxy Catalog PER
INPUT CATALOG (CW-All, CW-HCG) -- these remain two distinct, independently
group-found catalogs (different aperture scales, different populations by
design; see membership_aperture_methodology.md) and are not merged into a
single cross-catalog table, since a shared "Group_ID" namespace does not
exist between them and forcing one would misrepresent the data. This
mirrors how Project A's xray_release_v1.1 and Project B's
membership_release_v1.0 are both already organized (CW-All / CW-HCG pairs).

No methodology modification: all group/membership/BGG/dynamical-mass
values are taken as-is from membership_release_v1.0 and
xray_release_v1.1. This script only merges, aggregates, and documents.

Usage:
    python build_master_catalog.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent.parent
XRAY_REPO = BASE_DIR.parent / "cosmos-web-xray-igm"
RELEASE_DIR = BASE_DIR / "outputs" / "release" / "membership_release_v1.0"
XRAY_RELEASE_DIR = XRAY_REPO / "outputs" / "release_v1.1"
MASTER_DIR = BASE_DIR / "outputs" / "release" / "master_catalog_v1.0"
MASTER_DIR.mkdir(parents=True, exist_ok=True)
COSMOSWEB2025_PATH = BASE_DIR / "data" / "galaxy_catalog_photz" / "galaxy_catalog.fits"
COSMOS2020_PATH = BASE_DIR.parent / "data" / "galaxy_catalog_cosmos2020" / "COSMOS2020_CLASSIC_R1_v2.2_p3.fits"
MATCH_RADIUS_ARCSEC = 0.5

sys.path.insert(0, str(Path(__file__).parent))

COSMOSWEB2025_MORPH_COLS = ["SERSIC", "SERSIC_err", "ANGLE_SERSIC", "ANGLE_SERSIC_err",
                            "RADIUS", "RADIUS_err", "AXRATIO", "AXRATIO_err"]


def load_cosmosweb2025_morphology() -> pd.DataFrame:
    """PRIMARY morphology source: direct positional cross-match against
    COSMOSWEB2025 (data/galaxy_catalog_photz/galaxy_catalog.fits, the
    JWST-based 2025 COSMOS-Web catalog) for its real structural-fit
    columns (Sersic index, effective radius, axis ratio -- SE++-based
    single-Sersic model fit, catalog-native RA_MODEL/DEC_MODEL positions,
    matching the same convention merge_photz_properties.py already uses
    for LP_* SED properties)."""
    t = Table.read(COSMOSWEB2025_PATH).to_pandas()
    cols = ["RA_MODEL", "DEC_MODEL"] + [c for c in COSMOSWEB2025_MORPH_COLS if c in t.columns]
    return t[cols].rename(columns={"RA_MODEL": "RA", "DEC_MODEL": "DEC"})


def load_cosmos2020_morphology() -> pd.DataFrame:
    """FALLBACK ONLY (used where COSMOSWEB2025 morphology is unmatched --
    footprint-margin galaxies outside JWST coverage): positional
    cross-match against COSMOS2020 for ACS_MU_CLASS (a compactness-based
    star/extended-source classifier -- NOT a Sersic-index morphology
    measurement). Also confirms, by direct column inspection, that
    COSMOS2020 has no independent AGN classification column, so
    COSMOS2020_AGN is reported as unavailable rather than fabricated (see
    QUALITY_FLAGS.md)."""
    t = Table.read(COSMOS2020_PATH).to_pandas()
    agn_cols = [c for c in t.columns if "AGN" in c.upper()]
    if agn_cols:
        print(f"  NOTE: COSMOS2020 unexpectedly has AGN-named column(s) {agn_cols} -- "
              f"not used here since this was assumed absent; revisit COSMOS2020_AGN handling.")
    return t[["ALPHA_J2000", "DELTA_J2000", "ACS_MU_CLASS"]].rename(
        columns={"ALPHA_J2000": "RA", "DELTA_J2000": "DEC"})


def load_group_release(slug: str) -> pd.DataFrame:
    stub = "CW_All" if slug == "cw_all" else "CW_HCG"
    g = pd.read_csv(RELEASE_DIR / f"{stub}_groups_release_v1.0.csv")
    g["Group_ID"] = g["Group_ID"].astype(str)
    return g


def load_xray_release(slug: str) -> pd.DataFrame:
    stub = "CW-All" if slug == "cw_all" else "CW-HCG"
    x = Table.read(XRAY_RELEASE_DIR / f"{stub}_xray_catalog_v1.1.fits").to_pandas()
    x["Group_ID"] = x["Group_ID"].astype(str)
    keep = ["Group_ID", "Is_Detected", "Flux_erg_cm2_s", "Flux_Error", "Luminosity_erg_s",
            "Luminosity_Error", "Temperature_keV", "Temperature_Error", "Significance_Sigma",
            "Background_Quality", "Is_Projected_Contaminated", "Contamination_Severity",
            "Is_Suspected_False_Positive", "LAMBDA_STAR"]
    return x[[c for c in keep if c in x.columns]]


def load_galaxy_enriched(tmp_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tmp_path)
    df["Group_ID"] = df["Group_ID"].astype(str)
    return df


def match_morphology(df: pd.DataFrame, ref: pd.DataFrame, value_cols: list, radius_arcsec: float) -> tuple:
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    valid_ref = ref.dropna(subset=["RA", "DEC"])
    gal_coord = SkyCoord(ra=df["RA"].values * u.deg, dec=df["DEC"].values * u.deg)
    ref_coord = SkyCoord(ra=valid_ref["RA"].values * u.deg, dec=valid_ref["DEC"].values * u.deg)
    idx, sep2d, _ = gal_coord.match_to_catalog_sky(ref_coord)
    matched = sep2d.arcsec <= radius_arcsec
    out = {}
    for col in value_cols:
        s = pd.Series(np.nan, index=df.index)
        s[matched] = valid_ref[col].values[idx[matched]]
        out[col] = s
    return out, matched


def build_master_galaxy_catalog(galaxies: pd.DataFrame, catalog_name: str,
                                 cw2025_morph: pd.DataFrame, c2020_morph: pd.DataFrame) -> pd.DataFrame:
    df = galaxies.copy()

    # --- AGN / radio flags (from add_agn_flags.py output; already computed, not re-derived) ---
    df["X_ray_AGN"] = df["is_xray_agn"].astype(bool)
    df["COSMOSWEB2025_AGN"] = df["is_sed_agn"].astype(bool)  # SED-AGN, fit against COSMOSWEB2025 photometry
    df["VLA_AGN"] = (df["radio_is_hlagn"].astype(bool) | df["radio_is_mlagn"].astype(bool)
                      | df["radio_is_excess_agn"].astype(bool))
    df["radio_galaxy_flag"] = df["radio_is_vla_match"].astype(bool)  # has a VLA 3GHz counterpart at all
    # No independent COSMOS2020 AGN classification exists in the data holdings
    # available to this pipeline (checked: COSMOS2020_CLASSIC_R1_v2.2_p3.fits
    # has no AGN column). Reported as unavailable, not fabricated -- see
    # QUALITY_FLAGS.md.
    df["COSMOS2020_AGN"] = np.nan
    df["any_AGN_flag"] = df["agn_flag"] != "none"

    # --- morphology: PRIMARY = COSMOSWEB2025 (real Sersic-fit structural
    # parameters), FALLBACK = COSMOS2020 ACS_MU_CLASS ONLY where COSMOSWEB2025
    # is unmatched (footprint-margin galaxies outside JWST coverage) ---
    cw2025_vals, cw2025_matched = match_morphology(
        df, cw2025_morph, [c for c in COSMOSWEB2025_MORPH_COLS if c in cw2025_morph.columns], MATCH_RADIUS_ARCSEC)
    df["morphology_sersic_index"] = cw2025_vals.get("SERSIC")
    df["morphology_sersic_index_err"] = cw2025_vals.get("SERSIC_err")
    df["morphology_radius_arcsec"] = cw2025_vals.get("RADIUS")
    df["morphology_radius_arcsec_err"] = cw2025_vals.get("RADIUS_err")
    df["morphology_axis_ratio"] = cw2025_vals.get("AXRATIO")
    df["morphology_axis_ratio_err"] = cw2025_vals.get("AXRATIO_err")
    df["morphology_source"] = np.where(cw2025_matched, "COSMOSWEB2025", "none")

    c2020_vals, c2020_matched = match_morphology(df, c2020_morph, ["ACS_MU_CLASS"], MATCH_RADIUS_ARCSEC)
    fallback = (~cw2025_matched) & c2020_matched
    df["morphology_class_ACS_MU_CLASS"] = c2020_vals["ACS_MU_CLASS"]
    df.loc[~fallback, "morphology_class_ACS_MU_CLASS"] = np.nan  # only keep as a genuine fallback value
    df.loc[fallback, "morphology_source"] = "COSMOS2020_fallback"

    print(f"  Morphology: COSMOSWEB2025 matched {int(cw2025_matched.sum())}/{len(df)}; "
          f"COSMOS2020 fallback used for {int(fallback.sum())} additional (footprint-margin) galaxies; "
          f"{int((~cw2025_matched & ~fallback).sum())} unmatched by either")

    # --- stellar mass / redshift / spectroscopic info ---
    df["stellar_mass_final"] = df["mass_final"]
    df["stellar_mass_LePhare_med_PDF"] = df.get("LP_mass_med_PDF", np.nan)  # LePhare fit against COSMOSWEB2025 photometry
    df["redshift"] = df["zfin"]
    df["redshift_uncertainty_dz"] = df["dz"]
    df["redshift_type"] = df["redshift_type"]  # spec-z or photo-z (dz-tier)

    # --- X-ray group flag: does this galaxy's group have an X-ray detection ---
    xray_detected_methods = {"XRAY", "XRAY+SPECZ"}
    df["X_ray_group_flag"] = df["Group_ID"].map(
        galaxies.drop_duplicates("Group_ID").set_index("Group_ID").get(
            "r200_method_placeholder", pd.Series(dtype=object))
    )  # placeholder overwritten by caller after group merge

    df["Catalog"] = catalog_name

    # --- resolve one-row-per-galaxy from possible cross-group multiplicity ---
    # (documented, not hidden -- see Membership_Catalog_Release_Report.md)
    dedup_key = ["RA", "DEC"]
    df = df.sort_values("membership_prob", ascending=False)
    n_memberships = df.groupby(dedup_key)["Group_ID"].transform("nunique")
    df["n_group_memberships"] = n_memberships
    other_groups = (df.groupby(dedup_key)["Group_ID"]
                     .apply(lambda s: ",".join(sorted(set(s.astype(str)))))
                     .rename("all_group_ids_this_catalog"))
    df = df.merge(other_groups, on=dedup_key, how="left")
    bgg_any = df.groupby(dedup_key)["is_bgg"].any().rename("is_bgg_any_group")
    df = df.merge(bgg_any, on=dedup_key, how="left")

    primary = df.drop_duplicates(subset=dedup_key, keep="first").copy()

    out_cols = [
        "Catalog", "Group_ID", "RA", "DEC", "redshift", "redshift_uncertainty_dz", "redshift_type",
        "is_bgg", "is_bgg_any_group", "n_group_memberships", "all_group_ids_this_catalog",
        "membership_prob", "sep_kpc", "dv",
        "stellar_mass_final", "stellar_mass_LePhare_med_PDF", "mass_source", "MR_final",
        "morphology_source", "morphology_sersic_index", "morphology_sersic_index_err",
        "morphology_radius_arcsec", "morphology_radius_arcsec_err",
        "morphology_axis_ratio", "morphology_axis_ratio_err", "morphology_class_ACS_MU_CLASS",
        "X_ray_AGN", "COSMOSWEB2025_AGN", "VLA_AGN", "radio_galaxy_flag", "COSMOS2020_AGN",
        "any_AGN_flag", "agn_flag", "radio_type",
    ]
    out_cols = [c for c in out_cols if c in primary.columns]
    return primary[out_cols].reset_index(drop=True)


def build_master_group_catalog(group_release: pd.DataFrame, xray_release: pd.DataFrame,
                                galaxy_master: pd.DataFrame, catalog_name: str) -> pd.DataFrame:
    g = group_release.merge(xray_release, on="Group_ID", how="left", suffixes=("", "_v11"))

    # BGG summary
    bgg = galaxy_master[galaxy_master["is_bgg"] == True]  # noqa: E712
    bgg_cols = bgg[["Group_ID", "RA", "DEC", "stellar_mass_final", "MR_final", "redshift", "redshift_type"]].rename(
        columns={"RA": "BGG_RA", "DEC": "BGG_DEC", "stellar_mass_final": "BGG_stellar_mass",
                 "MR_final": "BGG_MR", "redshift": "BGG_redshift", "redshift_type": "BGG_redshift_type"})
    g = g.merge(bgg_cols, on="Group_ID", how="left")

    # AGN / radio statistics per group
    agn_stats = galaxy_master.groupby("Group_ID").agg(
        n_members_total_master=("RA", "size"),
        n_X_ray_AGN=("X_ray_AGN", "sum"),
        n_COSMOSWEB2025_AGN=("COSMOSWEB2025_AGN", "sum"),
        n_VLA_AGN=("VLA_AGN", "sum"),
        n_radio_galaxies=("radio_galaxy_flag", "sum"),
        n_any_AGN=("any_AGN_flag", "sum"),
    ).reset_index()
    agn_stats["AGN_fraction"] = agn_stats["n_any_AGN"] / agn_stats["n_members_total_master"]
    agn_stats["radio_galaxy_fraction"] = agn_stats["n_radio_galaxies"] / agn_stats["n_members_total_master"]
    g = g.merge(agn_stats, on="Group_ID", how="left")

    g["Catalog_master"] = catalog_name
    return g


def main():
    print("Loading COSMOSWEB2025 (primary) and COSMOS2020 (fallback-only) for morphology cross-match...")
    cw2025_morph = load_cosmosweb2025_morphology()
    c2020_morph = load_cosmos2020_morphology()
    print(f"  COSMOSWEB2025: {len(cw2025_morph)} galaxies, "
          f"COSMOS2020 (fallback): {len(c2020_morph)} galaxies")

    for slug, catalog_name, tmp_prefix in [("cw_all", "CW-All", "_tmp_cw_all_photz_agn.csv"),
                                            ("cw_hcg", "CW-HCG", "_tmp_cw_hcg_photz_agn.csv")]:
        print(f"\n{'='*70}\n{catalog_name}\n{'='*70}")
        group_release = load_group_release(slug)
        xray_release = load_xray_release(slug)
        galaxies_enriched = load_galaxy_enriched(RELEASE_DIR / tmp_prefix)

        # attach r200_method to galaxy rows for the X-ray group flag before dedup
        method_map = group_release.set_index("Group_ID")["r200_method"]
        galaxies_enriched["r200_method_placeholder"] = galaxies_enriched["Group_ID"].map(method_map)

        galaxy_master = build_master_galaxy_catalog(galaxies_enriched, catalog_name, cw2025_morph, c2020_morph)
        xray_detected_ids = set(group_release.loc[group_release["r200_method"].isin(["XRAY", "XRAY+SPECZ"]), "Group_ID"])
        galaxy_master["X_ray_group_flag"] = galaxy_master["Group_ID"].isin(xray_detected_ids)

        group_master = build_master_group_catalog(group_release, xray_release, galaxy_master, catalog_name)

        gal_stub = f"Master_Galaxy_Catalog_{slug.upper()}"
        grp_stub = f"Master_Group_Catalog_{slug.upper()}"
        galaxy_master.to_csv(MASTER_DIR / f"{gal_stub}.csv", index=False)
        Table.from_pandas(galaxy_master).write(MASTER_DIR / f"{gal_stub}.fits", format="fits", overwrite=True)
        group_master.to_csv(MASTER_DIR / f"{grp_stub}.csv", index=False)
        Table.from_pandas(group_master).write(MASTER_DIR / f"{grp_stub}.fits", format="fits", overwrite=True)

        print(f"  Master Galaxy Catalog: {len(galaxy_master)} galaxies -> {gal_stub}.[csv,fits]")
        print(f"  Master Group Catalog: {len(group_master)} groups -> {grp_stub}.[csv,fits]")
        print(f"  AGN counts: X-ray={int(galaxy_master.X_ray_AGN.sum())}, "
              f"COSMOSWEB2025(SED)={int(galaxy_master.COSMOSWEB2025_AGN.sum())}, "
              f"VLA={int(galaxy_master.VLA_AGN.sum())}, "
              f"radio galaxies={int(galaxy_master.radio_galaxy_flag.sum())}, "
              f"any={int(galaxy_master.any_AGN_flag.sum())}")
        print(f"  n_group_memberships > 1 (cross-group multiplicity, unresolved by design): "
              f"{int((galaxy_master.n_group_memberships > 1).sum())}")

    print(f"\nAll master catalogs written to {MASTER_DIR}")


if __name__ == "__main__":
    main()
