#!/usr/bin/env python
"""
Spectroscopic Membership Determination using VRF Method

Determines spec-z membership for COSMOS-Web compact groups using the VRF
(Velocity-Radius Fitting) interloper rejection algorithm.

Features:
- VRF iterative interloper rejection (Mamon et al. 2010)
- NFW velocity dispersion profile
- Dynamical properties (σ_v, M_200, R_200)
- Quality flags for reliability assessment
- Support for both automatic and parameter-specified runs

Usage:
    # Process all groups with default parameters
    python determine_specz_membership.py --catalog both --method vrf
    
    # Test run with 10 groups
    python determine_specz_membership.py --catalog all --method vrf --test
    
    # Custom parameters
    python determine_specz_membership.py --catalog all --method vrf \
        --max-dz-norm 0.008 --max-velocity 1500 --radius 500

    # Default spec-z is data/specz/Webb_Specz_Feb2026.fits when present; override if needed:
    python determine_specz_membership.py --catalog both --method vrf --specz-catalog data/specz/OTHER.fits
"""

import numpy as np
import pandas as pd
import argparse
import warnings
from pathlib import Path
from typing import Optional
from tqdm import tqdm
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from membership_functions import find_specz_members
from astropy.cosmology import Planck18 as cosmo
from radius_optimizer import get_per_group_radius

warnings.filterwarnings('ignore')

# Paths (project root = .../membership_determination/scripts -> parent.parent.parent)
BASE_DIR = Path(__file__).resolve().parent.parent.parent
GROUP_CATALOG_DIR = BASE_DIR / 'data' / 'group-catalog'
# Final COSMOS Web galaxy spec-z catalog (use .csv or .fits)
SPECZ_DIR = BASE_DIR / 'data' / 'specz'
SPECZ_CATALOG_CSV = SPECZ_DIR / 'Webb_Specz.csv'
SPECZ_CATALOG_FITS = SPECZ_DIR / 'Webb_Specz.fits'
# Default full COSMOS-Web spec-z release for this pipeline
SPECZ_WEBB_FEB2026_FITS = SPECZ_DIR / 'Webb_Specz_Feb2026.fits'
# COMaGN/COSMOS-Web French spec-z release (small subset; pass explicitly if needed)
SPECZ_COMAGN_FITS = SPECZ_DIR / 'COMAGN_FR_Webb_Specz_Feb2026.fits'
# Final HCG group catalog (use .csv or .fits)
HCG_CATALOG_CSV = GROUP_CATALOG_DIR / 'Py18_Groups.csv'
HCG_CATALOG_FITS = GROUP_CATALOG_DIR / 'Py18_Groups.fits'
OUTPUT_DIR = BASE_DIR / 'membership_determination' / 'results' / 'specz'

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _normalize_specz_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Match column names expected by membership_functions (RA, DEC, zfin, ez)."""
    out = df.copy()
    if "RA" not in out.columns and "RA_web" in out.columns:
        out["RA"] = pd.to_numeric(out["RA_web"], errors="coerce")
    if "DEC" not in out.columns and "DEC_web" in out.columns:
        out["DEC"] = pd.to_numeric(out["DEC_web"], errors="coerce")
    if "zfin" not in out.columns and "z" in out.columns:
        out["zfin"] = pd.to_numeric(out["z"], errors="coerce")
    if "ez" not in out.columns and "dz" in out.columns:
        out["ez"] = pd.to_numeric(out["dz"], errors="coerce")
    elif "ez" not in out.columns:
        out["ez"] = np.nan
    return out


def _read_specz_file(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return _normalize_specz_dataframe(pd.read_csv(path))
    from astropy.table import Table
    return _normalize_specz_dataframe(Table.read(path).to_pandas())


def load_specz_catalog(specz_path: Optional[Path] = None):
    """
    Load COSMOS-Web spec-z catalog.

    Tries explicit path, then ``data/specz/Webb_Specz_Feb2026.fits`` (pinned default),
    then Webb_Specz.csv / Webb_Specz.fits, then newest Webb_Specz*.fits (excl. *with_photz*).
    """
    if specz_path is not None:
        p = Path(specz_path)
        if not p.exists():
            raise FileNotFoundError(f"Spec-z catalog not found: {p}")
        print(f"  Spec-z catalog file: {p.resolve()}")
        return _read_specz_file(p)

    if SPECZ_WEBB_FEB2026_FITS.is_file():
        print(f"  Spec-z catalog file: {SPECZ_WEBB_FEB2026_FITS.resolve()} (default)")
        return _read_specz_file(SPECZ_WEBB_FEB2026_FITS)

    if SPECZ_CATALOG_CSV.exists():
        print(f"  Spec-z catalog file: {SPECZ_CATALOG_CSV.resolve()}")
        return _read_specz_file(SPECZ_CATALOG_CSV)
    if SPECZ_CATALOG_FITS.exists():
        print(f"  Spec-z catalog file: {SPECZ_CATALOG_FITS.resolve()}")
        return _read_specz_file(SPECZ_CATALOG_FITS)

    candidates: list[Path] = []
    if SPECZ_DIR.is_dir():
        for p in SPECZ_DIR.glob("Webb_Specz*.fits"):
            if "with_photz" in p.name.lower():
                continue
            candidates.append(p)
        candidates.sort(key=lambda x: x.stat().st_mtime, reverse=True)

    if candidates:
        chosen = candidates[0]
        print(f"  Spec-z catalog file: {chosen.resolve()} (newest Webb_Specz*.fits, excl. *with_photz*)")
        return _read_specz_file(chosen)

    # Last resort: small COMaGN high-confidence subset (~10^2–10^3 sources), not the full Webb spec-z sample
    if SPECZ_COMAGN_FITS.is_file():
        print(
            f"  Spec-z catalog file: {SPECZ_COMAGN_FITS.resolve()} "
            "(fallback; prefer Webb_Specz_Feb2026.fits in data/specz/)"
        )
        return _read_specz_file(SPECZ_COMAGN_FITS)

    raise FileNotFoundError(
        "Spec-z catalog not found. Expected one of:\n"
        f"  {SPECZ_WEBB_FEB2026_FITS}\n"
        f"  {SPECZ_CATALOG_CSV}\n"
        f"  {SPECZ_CATALOG_FITS}\n"
        f"  {SPECZ_DIR / 'Webb_Specz*.fits'} (excluding *with_photz*)\n"
        "Or pass --specz-catalog /path/to/catalog.fits"
    )


def load_hcg_catalog(hcg_path: Optional[Path] = None):
    """Load final HCG group catalog (Py18_Groups.csv or Py18_Groups.fits)."""
    if hcg_path is not None:
        p = Path(hcg_path)
        if not p.exists():
            raise FileNotFoundError(f"HCG catalog not found: {p}")
        if p.suffix.lower() == ".csv":
            return pd.read_csv(p)
        from astropy.table import Table
        return Table.read(p).to_pandas()
    if HCG_CATALOG_CSV.exists():
        return pd.read_csv(HCG_CATALOG_CSV)
    if HCG_CATALOG_FITS.exists():
        from astropy.table import Table
        t = Table.read(HCG_CATALOG_FITS)
        return t.to_pandas()
    raise FileNotFoundError(
        f"HCG catalog not found. Expected {HCG_CATALOG_CSV} or {HCG_CATALOG_FITS}, "
        f"or pass --hcg-catalog /path/to/Py18_Groups.csv"
    )


def load_cw_all_catalog(cw_all_path: Optional[Path] = None):
    """
    Load CW-All group catalog. Tries CSV, then FITS, then refined-z FITS under
    data/group-catalog/ unless cw_all_path is given.
    """
    if cw_all_path is not None:
        p = Path(cw_all_path)
        if not p.exists():
            raise FileNotFoundError(f"CW-All group catalog not found: {p}")
        if p.suffix.lower() == ".csv":
            return pd.read_csv(p)
        from astropy.table import Table
        return Table.read(p).to_pandas()

    candidates = [
        GROUP_CATALOG_DIR / "cosmos_web_groups_catalog.csv",
        GROUP_CATALOG_DIR / "cosmos_web_groups_catalog.fits",
        GROUP_CATALOG_DIR / "cosmos_web_groups_catalog_refined_z.fits",
    ]
    for p in candidates:
        if not p.exists():
            continue
        if p.suffix.lower() == ".csv":
            return pd.read_csv(p)
        from astropy.table import Table
        return Table.read(p).to_pandas()

    raise FileNotFoundError(
        "CW-All group catalog not found. Expected one of:\n  "
        + "\n  ".join(str(c) for c in candidates)
        + "\nOr pass --cw-all-catalog /path/to/cosmos_web_groups_catalog.csv"
    )


def load_catalogs(
    cw_all_path: Optional[Path] = None,
    hcg_path: Optional[Path] = None,
    specz_path: Optional[Path] = None,
):
    """Load group and spec-z catalogs."""
    print("Loading catalogs...")
    
    # Group catalogs
    cw_all = load_cw_all_catalog(cw_all_path)
    cw_hcg = load_hcg_catalog(hcg_path)
    
    specz = load_specz_catalog(specz_path)
    
    print(f"  CW-All groups: {len(cw_all)}")
    print(f"  CW-HCG groups: {len(cw_hcg)}")
    print(f"  Spec-z galaxies: {len(specz)}")
    
    return cw_all, cw_hcg, specz


def get_xray_detected_group_ids(xray_path: Path, catalog_name: str) -> set:
    """Return set of Group_IDs that are X-ray detected in the given catalog."""
    if not xray_path or not xray_path.exists():
        return set()
    try:
        if xray_path.suffix.lower() == '.csv':
            df = pd.read_csv(xray_path)
        else:
            from astropy.table import Table
            df = Table.read(xray_path).to_pandas()
    except Exception:
        return set()
    if 'Group_ID' not in df.columns or 'Is_Detected' not in df.columns:
        return set()
    if 'Catalog_Name' in df.columns:
        mask = (
            (df['Catalog_Name'].astype(str).str.strip().str.upper() == catalog_name.upper()) &
            (df['Is_Detected'] == True)
        )
    else:
        mask = (df['Is_Detected'] == True)
    return set(df.loc[mask, 'Group_ID'].astype(str).unique())


def add_vrf_dynamical_properties(summary, specz_catalog):
    """
    Calculate dynamical properties from VRF membership.
    
    Adds:
    - sigma_v_kms: Velocity dispersion
    - sigma_v_err_kms: Error on velocity dispersion
    - n_vrf_members: Number of VRF members
    - vrf_quality: Quality flag (1-4)
    - z_vrf_refined: VRF-refined group redshift
    """
    c_kms = 299792.458  # Speed of light in km/s
    
    # Determine member column name
    member_col = 'vrf_member' if 'vrf_member' in specz_catalog.columns else 'gapper_member'
    # Group ID column in members: CW-All uses Group_Group_ID, CW-HCG uses Group_Grp
    group_id_col = 'Group_Group_ID' if 'Group_Group_ID' in specz_catalog.columns else 'Group_Grp'
    
    # Calculate properties for each group
    refined_data = []
    
    for _, row in tqdm(summary.iterrows(), total=len(summary), desc="Computing dynamical properties"):
        group_id = row['Group_ID']
        catalog_z = row['Group_z']
        
        # Get VRF members for this group
        group_members = specz_catalog[
            (specz_catalog[group_id_col] == group_id) &
            (specz_catalog[member_col] == True)
        ]
        
        n_members = len(group_members)
        
        if n_members >= 3:
            member_z = group_members['zfin'].values
            
            # Refined redshift (median of members)
            z_refined = np.median(member_z)
            z_std = np.std(member_z)
            
            # Velocity dispersion
            member_velocities = c_kms * (member_z - z_refined) / (1 + z_refined)
            sigma_v = np.std(member_velocities)
            
            # Error estimate
            if n_members >= 5:
                sigma_v_err = sigma_v / np.sqrt(2 * (n_members - 1))
            else:
                sigma_v_err = sigma_v / np.sqrt(n_members)
            
            # Quality flag
            if n_members >= 15:
                quality = 1  # Excellent
            elif n_members >= 8:
                quality = 2  # Good
            elif n_members >= 5:
                quality = 3  # Poor
            else:
                quality = 4  # Insufficient
            
            refined_data.append({
                'Group_ID': group_id,
                'z_vrf_refined': z_refined,
                'z_vrf_std': z_std,
                'sigma_v_kms': sigma_v,
                'sigma_v_err_kms': sigma_v_err,
                'n_vrf_members': n_members,
                'vrf_quality': quality,
                'n_specz_total': n_members  # Add total count
            })
        else:
            # Not enough members
            refined_data.append({
                'Group_ID': group_id,
                'z_vrf_refined': catalog_z,
                'z_vrf_std': np.nan,
                'sigma_v_kms': np.nan,
                'sigma_v_err_kms': np.nan,
                'n_vrf_members': n_members,
                'vrf_quality': 4,
                'n_specz_total': n_members
            })
    
    # Merge with summary
    refined_df = pd.DataFrame(refined_data)
    summary = summary.merge(refined_df, on='Group_ID', how='left')
    
    return summary


def process_catalog(groups, specz, catalog_name, radius_kpc, max_dz_norm,
                   max_velocity, method, sample_size=None,
                   radius_mode='fixed', xray_path=None, xray_scale=1.0,
                   radius_min=300.0, radius_max=1500.0, use_r500_fallback=True):
    """
    Process a group catalog to find spec-z members.
    
    Parameters:
    -----------
    groups : DataFrame
        Group catalog
    specz : DataFrame
        Spec-z catalog
    catalog_name : str
        'CW-All' or 'CW-HCG'
    radius_kpc : float
        Default search radius in kpc (fixed mode or fallback)
    max_dz_norm : float
        Maximum normalized redshift offset
    max_velocity : float
        Maximum velocity in km/s
    method : str
        'vrf' or 'gapper'
    sample_size : int, optional
        Number of groups to process (for testing)
    radius_mode : str
        'fixed' | 'xray' | 'redshift' | 'hybrid' for radius optimisation
    xray_path : Path, optional
        Path to X-ray catalog (for xray/hybrid mode)
    xray_scale : float
        Scale factor for X-ray R200/R500 (e.g. 1.2)
    radius_min, radius_max : float
        Clamp per-group radius to [radius_min, radius_max] kpc
    
    Returns:
    --------
    specz_members : DataFrame
        All spec-z members
    summary : DataFrame
        Summary statistics per group
    """
    
    if sample_size is not None:
        groups = groups.head(sample_size)
    
    # Determine group ID and position columns
    # CW-All: Group_ID; CW-HCG: Py18 has Group_ID, older catalogs use Grp
    if catalog_name == 'CW-All':
        id_col, ra_col, dec_col, z_col = 'Group_ID', 'Ra', 'Dec', 'z'
    else:
        if 'Group_ID' in groups.columns:
            id_col, ra_col, dec_col, z_col = 'Group_ID', 'Ra', 'Dec', 'z'
        else:
            id_col, ra_col, dec_col, z_col = 'Grp', 'Ra', 'Dec', 'z'
    
    # Per-group radius when optimised
    radius_arr = None
    if radius_mode != 'fixed':
        radius_arr, _ = get_per_group_radius(
            groups, catalog_name, id_col, z_col,
            radius_mode=radius_mode,
            fixed_radius_kpc=radius_kpc,
            xray_path=xray_path,
            xray_scale=xray_scale,
            xray_r_min=radius_min,
            xray_r_max=radius_max,
            use_r500_fallback=use_r500_fallback,
        )
    
    print(f"\nProcessing {catalog_name} catalog (N={len(groups)} groups)")
    if radius_arr is not None:
        print(f"  Radius: per-group ({radius_mode}), mean={np.nanmean(radius_arr):.0f} kpc, range=[{np.nanmin(radius_arr):.0f}, {np.nanmax(radius_arr):.0f}]")
    else:
        print(f"  Radius: {radius_kpc} kpc (fixed)")
    print(f"  Method: {method.upper()}")
    print(f"  max_dz_norm: {max_dz_norm}")
    print(f"  max_velocity: {max_velocity} km/s")
    
    if sample_size is not None:
        print(f"  TEST MODE: Processing only {sample_size} groups")
    
    specz_members_list = []
    summary_list = []
    
    for i, (idx, row) in enumerate(tqdm(groups.iterrows(), total=len(groups), desc=f"Processing {catalog_name}")):
        group_id = row[id_col]
        group_ra = row[ra_col]
        group_dec = row[dec_col]
        group_z = row[z_col]
        
        # Skip very high redshift groups
        if group_z > 4.0:
            continue
        
        # Per-group radius when optimised
        r_kpc = float(radius_arr[i]) if radius_arr is not None and i < len(radius_arr) else radius_kpc
        
        # Find spec-z members
        specz_mem = find_specz_members(
            group_ra, group_dec, group_z, specz,
            radius_kpc=r_kpc,
            max_dz_norm=max_dz_norm,
            max_velocity=max_velocity,
            use_gapper=True,
            remove_duplicates=True,
            method=method
        )
        
        # Add group info
        if len(specz_mem) > 0:
            specz_mem_copy = specz_mem.copy()
            
            # Handle duplicates
            if 'ez' in specz_mem_copy.columns:
                specz_mem_copy = specz_mem_copy.sort_values('ez')
                specz_mem_copy = specz_mem_copy.drop_duplicates(
                    subset=['RA', 'DEC'], keep='first'
                )
            
            # Add group properties
            for col in groups.columns:
                specz_mem_copy[f'Group_{col}'] = row[col]
            specz_mem_copy['catalog'] = catalog_name
            specz_mem_copy['redshift_type'] = 'spec-z'
            
            # Add membership quality flag
            member_col = 'vrf_member' if 'vrf_member' in specz_mem_copy.columns else 'gapper_member'
            
            specz_mem_copy['membership_flag'] = 3  # Default: low confidence
            
            # High confidence: Member + |Δv| < 500 km/s + |Δz/(1+z)| < 0.01
            high_conf_mask = (
                specz_mem_copy[member_col] & 
                (np.abs(specz_mem_copy['dv']) < 500) & 
                (np.abs(specz_mem_copy['dz_norm']) < 0.01)
            )
            specz_mem_copy.loc[high_conf_mask, 'membership_flag'] = 1
            
            # Medium confidence
            medium_conf_mask = (
                (specz_mem_copy[member_col]) | 
                ((np.abs(specz_mem_copy['dv']) < 1000) & (np.abs(specz_mem_copy['dz_norm']) < 0.015))
            ) & (~high_conf_mask)
            specz_mem_copy.loc[medium_conf_mask, 'membership_flag'] = 2
            
            specz_members_list.append(specz_mem_copy)
        
        # Summary statistics
        if len(specz_mem) > 0:
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            n_members = len(specz_mem[specz_mem[member_col]])
        else:
            n_members = 0
        
        summary_list.append({
            'Group_ID': group_id,
            'Group_Ra': group_ra,
            'Group_Dec': group_dec,
            'Group_z': group_z,
            'n_specz_members': n_members,
            'catalog': catalog_name
        })
    
    # Compile results
    if len(specz_members_list) > 0:
        specz_members = pd.concat(specz_members_list, ignore_index=True)
    else:
        specz_members = pd.DataFrame()
    
    if len(summary_list) > 0:
        summary = pd.DataFrame(summary_list)
    else:
        summary = pd.DataFrame()
    
    print(f"  Found spec-z members for {len(summary)} groups")
    print(f"  Total spec-z members: {len(specz_members)}")
    
    return specz_members, summary


def main():
    parser = argparse.ArgumentParser(
        description='Determine spectroscopic membership using VRF method',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog', choices=['all', 'hcg', 'both'], default='both',
                       help='Which catalog to process')
    parser.add_argument('--method', choices=['vrf', 'gapper'], default='vrf',
                       help='Membership determination method')
    parser.add_argument('--radius', type=float, default=500,
                       help='Physical search radius in kpc (fixed mode or fallback)')
    parser.add_argument('--radius-mode', choices=['fixed', 'xray', 'redshift', 'hybrid'], default='fixed',
                       help='Radius strategy: fixed | xray (R200/R500 from X-ray) | redshift (z-dependent) | hybrid (xray else redshift)')
    parser.add_argument('--xray-catalog', type=str, default=None,
                       help='Path to X-ray catalog (FITS/CSV) for radius optimisation; default: outputs/results/<catalog>/xray_catalog')
    parser.add_argument('--xray-scale', type=float, default=1.0,
                       help='Scale X-ray R200/R500 by this factor (e.g. 1.2 to include outskirts)')
    parser.add_argument('--radius-min', type=float, default=300,
                       help='Minimum radius in kpc when using xray/redshift/hybrid')
    parser.add_argument('--radius-max', type=float, default=1500,
                       help='Maximum radius in kpc when using xray/redshift/hybrid')
    parser.add_argument('--xray-detected-only', action='store_true',
                       help='Restrict to X-ray detected groups only (from xray_catalog Is_Detected)')
    parser.add_argument('--radius-use-r200-only', action='store_true',
                       help='When using xray radius: use only R200 (no R500 fallback)')
    parser.add_argument('--max-dz-norm', type=float, default=0.01,
                       help='Maximum normalized redshift offset |Δz/(1+z)|')
    parser.add_argument('--max-velocity', type=float, default=2000,
                       help='Maximum velocity offset in km/s')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory (default: results/specz)')
    parser.add_argument('--cw-all-catalog', type=str, default=None,
                       help='Path to CW-All group catalog (CSV or FITS); default: data/group-catalog/cosmos_web_groups_catalog.*')
    parser.add_argument('--hcg-catalog', type=str, default=None,
                       help='Path to HCG group catalog (CSV or FITS); default: data/group-catalog/Py18_Groups.*')
    parser.add_argument(
        '--specz-catalog',
        type=str,
        default=None,
        help=(
            'Path to spec-z catalog (CSV or FITS). Default: data/specz/Webb_Specz_Feb2026.fits if present, '
            'else Webb_Specz.csv / Webb_Specz.fits / newest Webb_Specz*.fits (excl. *with_photz*). '
            'COMaGN FR table is ~674 rows — pass it only when you want that subset.'
        ),
    )
    parser.add_argument('--test', action='store_true',
                       help='Test mode: process only 10 groups')
    
    args = parser.parse_args()
    
    # Set output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load catalogs
    cw_all_p = Path(args.cw_all_catalog) if args.cw_all_catalog else None
    hcg_p = Path(args.hcg_catalog) if args.hcg_catalog else None
    specz_p = Path(args.specz_catalog) if args.specz_catalog is not None else None
    cw_all, cw_hcg, specz = load_catalogs(
        cw_all_path=cw_all_p, hcg_path=hcg_p, specz_path=specz_p,
    )
    
    # Test mode
    sample_size = 10 if args.test else None
    
    # X-ray catalog path per catalog (for radius optimisation and xray-detected-only)
    xray_all = Path(args.xray_catalog) if args.xray_catalog else (BASE_DIR / 'outputs' / 'results' / 'cw_all' / 'xray_catalog.fits')
    if not xray_all.exists():
        xray_all = BASE_DIR / 'outputs' / 'results' / 'cw_all' / 'xray_catalog.csv'
    xray_hcg = Path(args.xray_catalog) if args.xray_catalog else (BASE_DIR / 'outputs' / 'results' / 'cw_hcg' / 'xray_catalog.fits')
    if not xray_hcg.exists():
        xray_hcg = BASE_DIR / 'outputs' / 'results' / 'cw_hcg' / 'xray_catalog.csv'
    
    # Restrict to X-ray detected groups only if requested
    use_r500_fallback = not getattr(args, 'radius_use_r200_only', False)
    if args.xray_detected_only:
        id_col_all = 'Group_ID'
        id_col_hcg = 'Group_ID' if 'Group_ID' in cw_hcg.columns else 'Grp'
        detected_all = get_xray_detected_group_ids(xray_all, 'CW-All')
        detected_hcg = get_xray_detected_group_ids(xray_hcg, 'CW-HCG')
        if args.catalog in ['all', 'both']:
            cw_all = cw_all[cw_all[id_col_all].astype(str).isin(detected_all)].copy()
            print(f"  X-ray detected only (CW-All): {len(cw_all)} groups")
        if args.catalog in ['hcg', 'both']:
            cw_hcg = cw_hcg[cw_hcg[id_col_hcg].astype(str).isin(detected_hcg)].copy()
            print(f"  X-ray detected only (CW-HCG): {len(cw_hcg)} groups")
    
    # Process catalogs
    if args.catalog in ['all', 'both']:
        members_all, summary_all = process_catalog(
            cw_all, specz, 'CW-All',
            args.radius, args.max_dz_norm, args.max_velocity, args.method,
            sample_size=sample_size,
            radius_mode=args.radius_mode, xray_path=xray_all,
            xray_scale=args.xray_scale, radius_min=args.radius_min, radius_max=args.radius_max,
            use_r500_fallback=use_r500_fallback,
        )
        
        # Add dynamical properties
        if len(summary_all) > 0 and args.method == 'vrf':
            summary_all = add_vrf_dynamical_properties(summary_all, members_all)
        
        # Save
        method_tag = f'_{args.method}' if args.method != 'gapper' else ''
        radius_tag = f'_r{args.radius_mode}' if args.radius_mode != 'fixed' else f'_r{int(args.radius)}kpc'
        if args.xray_detected_only:
            radius_tag += '_xrayonly'
        if args.radius_use_r200_only:
            radius_tag += '_r200only'
        if len(members_all) > 0:
            members_file = output_dir / f'cw_all_members_specz{method_tag}{radius_tag}.csv'
            members_all.to_csv(members_file, index=False)
            print(f"\nSaved: {members_file}")
        
        if len(summary_all) > 0:
            summary_file = output_dir / f'cw_all_summary{method_tag}{radius_tag}.csv'
            summary_all.to_csv(summary_file, index=False)
            print(f"Saved: {summary_file}")
    
    if args.catalog in ['hcg', 'both']:
        members_hcg, summary_hcg = process_catalog(
            cw_hcg, specz, 'CW-HCG',
            args.radius, args.max_dz_norm, args.max_velocity, args.method,
            sample_size=sample_size,
            radius_mode=args.radius_mode, xray_path=xray_hcg,
            xray_scale=args.xray_scale, radius_min=args.radius_min, radius_max=args.radius_max,
            use_r500_fallback=use_r500_fallback,
        )
        
        # Add dynamical properties
        if len(summary_hcg) > 0 and args.method == 'vrf':
            summary_hcg = add_vrf_dynamical_properties(summary_hcg, members_hcg)
        
        # Save
        method_tag = f'_{args.method}' if args.method != 'gapper' else ''
        radius_tag = f'_r{args.radius_mode}' if args.radius_mode != 'fixed' else f'_r{int(args.radius)}kpc'
        if args.xray_detected_only:
            radius_tag += '_xrayonly'
        if args.radius_use_r200_only:
            radius_tag += '_r200only'
        if len(members_hcg) > 0:
            members_file = output_dir / f'cw_hcg_members_specz{method_tag}{radius_tag}.csv'
            members_hcg.to_csv(members_file, index=False)
            print(f"\nSaved: {members_file}")
        
        if len(summary_hcg) > 0:
            summary_file = output_dir / f'cw_hcg_summary{method_tag}{radius_tag}.csv'
            summary_hcg.to_csv(summary_file, index=False)
            print(f"Saved: {summary_file}")
    
    print("\n" + "="*60)
    print("SPEC-Z MEMBERSHIP DETERMINATION COMPLETE")
    print("="*60)


if __name__ == '__main__':
    main()
