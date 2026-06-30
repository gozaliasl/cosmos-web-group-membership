#!/usr/bin/env python
"""
Photo-z Membership Determination Script

Determines galaxy group membership using photometric redshifts with probabilistic
approach and background contamination correction.

Features:
- Probabilistic photo-z membership (P_member calculation)
- Background density correction
- Configurable 2σ NMAD cylinder (compact groups)
- Handles both CW-All and CW-HCG catalogs
- Test mode for validation
- Saves results to CSV and FITS

Usage:
    python determine_photoz_membership.py --catalog all
    python determine_photoz_membership.py --catalog hcg --max-dist-error 2.5
    python determine_photoz_membership.py --catalog both --test

Requirements:
    conda activate astro-clean
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u
import argparse
import warnings
from pathlib import Path
from datetime import datetime

warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path('/Users/gozalig1/Projects/compact-groups-xray-analysis')
GROUP_CATALOG_DIR = BASE_DIR / 'data' / 'group-catalog'
OUTPUT_DIR = BASE_DIR / 'membership_determination' / 'results' / 'photoz'

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def calculate_background_density(ra_center, dec_center, z_center, catalog_df,
                                 inner_radius_kpc=500, outer_radius_kpc=1000,
                                 dz_norm=0.02):
    """
    Calculate background galaxy density in annulus around group.
    
    Parameters:
    -----------
    ra_center, dec_center : float
        Group center coordinates (degrees)
    z_center : float
        Group redshift
    catalog_df : DataFrame
        Full galaxy catalog
    inner_radius_kpc : float
        Inner radius of annulus (kpc)
    outer_radius_kpc : float
        Outer radius of annulus (kpc)
    dz_norm : float
        Normalized redshift window: |Δz/(1+z)| < dz_norm
    
    Returns:
    --------
    float : Background density (galaxies per kpc²)
    """
    # Convert radii to arcsec
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z_center).value / 60.0
    inner_arcsec = inner_radius_kpc / kpc_per_arcsec
    outer_arcsec = outer_radius_kpc / kpc_per_arcsec
    
    # Find galaxies in annulus
    coords_center = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg)
    coords_gals = SkyCoord(ra=catalog_df['Ra'].values*u.deg, 
                           dec=catalog_df['Dec'].values*u.deg)
    sep_arcsec = coords_center.separation(coords_gals).arcsec
    
    # Redshift cut
    dz_norm_vals = np.abs((catalog_df['z_phot'].values - z_center) / (1 + z_center))
    
    # Select annulus galaxies
    in_annulus = (sep_arcsec >= inner_arcsec) & (sep_arcsec < outer_arcsec)
    in_z_window = dz_norm_vals < dz_norm
    annulus_gals = np.sum(in_annulus & in_z_window)
    
    # Calculate area
    area_kpc2 = np.pi * (outer_radius_kpc**2 - inner_radius_kpc**2)
    
    # Density
    density = annulus_gals / area_kpc2 if area_kpc2 > 0 else 0
    
    return density


def find_photoz_members(ra_center, dec_center, z_center, catalog_df,
                        radius_kpc=500, max_dist_phot_z_error=2.0,
                        prob_threshold=0.5, use_background_correction=True):
    """
    Find photo-z group members using probabilistic approach.
    
    Parameters:
    -----------
    ra_center, dec_center : float
        Group center coordinates (degrees)
    z_center : float
        Group redshift
    catalog_df : DataFrame
        Galaxy catalog with photo-z information
    radius_kpc : float
        Search radius (kpc)
    max_dist_phot_z_error : float
        Maximum distance in units of photo-z NMAD (e.g., 2.0 = 2σ)
    prob_threshold : float
        Minimum P_member to be considered a member
    use_background_correction : bool
        Apply background contamination correction
    
    Returns:
    --------
    DataFrame : Photo-z members with membership probabilities
    """
    # Convert radius to arcsec
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z_center).value / 60.0
    radius_arcsec = radius_kpc / kpc_per_arcsec
    
    # Find galaxies within projected radius
    coords_center = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg)
    coords_gals = SkyCoord(ra=catalog_df['Ra'].values*u.deg, 
                           dec=catalog_df['Dec'].values*u.deg)
    sep_arcsec = coords_center.separation(coords_gals).arcsec
    sep_kpc = sep_arcsec * kpc_per_arcsec
    
    # Spatial cut
    spatial_cut = sep_arcsec < radius_arcsec
    
    if np.sum(spatial_cut) == 0:
        return pd.DataFrame()
    
    # Get photo-z errors (NMAD or similar)
    # Assuming photo-z error column exists; adjust as needed
    if 'z_phot_err' in catalog_df.columns:
        z_phot_err = catalog_df['z_phot_err'].values
    elif 'NMAD' in catalog_df.columns:
        z_phot_err = catalog_df['NMAD'].values
    else:
        # Default: assume 3% error if no error column
        z_phot_err = 0.03 * (1 + catalog_df['z_phot'].values)
    
    # Redshift distance in units of photo-z error
    z_phot = catalog_df['z_phot'].values
    dz = np.abs(z_phot - z_center)
    dz_sigma = dz / z_phot_err
    
    # Photo-z selection
    photoz_cut = dz_sigma < max_dist_phot_z_error
    
    # Combined cut
    candidates = spatial_cut & photoz_cut
    
    if np.sum(candidates) == 0:
        return pd.DataFrame()
    
    # Extract candidates
    candidate_gals = catalog_df[candidates].copy()
    candidate_gals['sep_kpc'] = sep_kpc[candidates]
    candidate_gals['dz_sigma'] = dz_sigma[candidates]
    
    # Calculate membership probability
    # Simple Gaussian model: P(z|member) vs P(z|background)
    # P_member = P(z|member) / [P(z|member) + P(z|background)]
    
    # Member probability (Gaussian centered at z_center)
    sigma_member = np.median(z_phot_err[candidates])
    prob_member = np.exp(-0.5 * (dz[candidates] / sigma_member)**2)
    prob_member /= (sigma_member * np.sqrt(2 * np.pi))
    
    # Background probability (uniform over wide z range)
    if use_background_correction:
        # Calculate background density
        bg_density = calculate_background_density(
            ra_center, dec_center, z_center, catalog_df,
            inner_radius_kpc=radius_kpc, 
            outer_radius_kpc=radius_kpc*2,
            dz_norm=max_dist_phot_z_error * sigma_member / (1 + z_center)
        )
        
        # Background probability (constant)
        prob_background = bg_density * np.ones(len(candidate_gals))
    else:
        prob_background = 0.1 * np.ones(len(candidate_gals))
    
    # Membership probability
    P_member = prob_member / (prob_member + prob_background)
    P_member = np.clip(P_member, 0, 1)  # Ensure [0, 1]
    
    candidate_gals['P_member'] = P_member
    candidate_gals['is_member'] = P_member >= prob_threshold
    
    return candidate_gals


def process_catalog(catalog_path, catalog_name, radius_kpc=500,
                   max_dist_phot_z_error=2.0, prob_threshold=0.5,
                   test_mode=False, groups_df=None):
    """
    Process entire catalog for photo-z membership.

    Parameters:
    -----------
    catalog_path : Path
        Path to group catalog (used for display; used to load if groups_df is None)
    catalog_name : str
        'cw-all' or 'cw-hcg'
    radius_kpc : float
        Search radius
    max_dist_phot_z_error : float
        Photo-z selection threshold (σ units)
    prob_threshold : float
        Minimum P_member threshold
    test_mode : bool
        Process only first 10 groups
    groups_df : pandas.DataFrame, optional
        Pre-loaded group catalog (e.g. from Py18_Groups.fits). If provided, catalog_path is only used for display.

    Returns:
    --------
    DataFrame : All groups with photo-z membership
    """
    print(f"\n{'='*60}")
    print(f"Processing {catalog_name.upper()} Catalog - Photo-z Membership")
    print(f"{'='*60}")
    print(f"  Catalog: {catalog_path.name}")
    print(f"  Search radius: {radius_kpc} kpc")
    print(f"  Photo-z threshold: {max_dist_phot_z_error}σ NMAD")
    print(f"  P_member threshold: {prob_threshold}")
    print(f"  Test mode: {test_mode}")
    print(f"{'='*60}\n")

    # Load catalog
    if groups_df is not None:
        groups = groups_df.copy()
    elif str(catalog_path).lower().endswith('.fits'):
        from astropy.table import Table
        groups = Table.read(catalog_path).to_pandas()
    else:
        groups = pd.read_csv(catalog_path)
    
    if test_mode:
        groups = groups.head(10)
        print(f"TEST MODE: Processing first {len(groups)} groups\n")
    
    # Determine columns
    if catalog_name == 'cw-all':
        id_col, ra_col, dec_col, z_col = 'Group_ID', 'Ra', 'Dec', 'z'
    else:
        id_col, ra_col, dec_col, z_col = 'Grp', 'Ra', 'Dec', 'z'
    
    # For photo-z, we use the same catalog as the source
    # (In real scenario, you'd use a separate galaxy catalog with photo-z)
    # Here assuming the group catalog has photo-z info or we use it as proxy
    
    # Add photo-z column if not present (placeholder)
    if 'z_phot' not in groups.columns:
        # Use spectroscopic z as proxy for demo
        groups['z_phot'] = groups[z_col]
        groups['z_phot_err'] = 0.03 * (1 + groups[z_col])
    
    results = []
    
    for idx, row in groups.iterrows():
        group_id = row[id_col]
        ra, dec, z = row[ra_col], row[dec_col], row[z_col]
        
        if idx % 100 == 0:
            print(f"  Processing group {idx+1}/{len(groups)} (ID={group_id})...")
        
        # Find photo-z members
        members = find_photoz_members(
            ra, dec, z, groups,
            radius_kpc=radius_kpc,
            max_dist_phot_z_error=max_dist_phot_z_error,
            prob_threshold=prob_threshold,
            use_background_correction=True
        )
        
        # Group summary
        n_candidates = len(members)
        n_members = np.sum(members['is_member']) if n_candidates > 0 else 0
        mean_P_member = members['P_member'].mean() if n_candidates > 0 else 0
        
        result = {
            'group_id': group_id,
            'catalog': catalog_name,
            'ra': ra,
            'dec': dec,
            'z': z,
            'n_candidates_photoz': n_candidates,
            'n_members_photoz': n_members,
            'mean_P_member': mean_P_member,
            'radius_kpc': radius_kpc,
            'max_dist_phot_z_error': max_dist_phot_z_error,
            'prob_threshold': prob_threshold,
            'timestamp': datetime.now().isoformat()
        }
        
        results.append(result)
    
    results_df = pd.DataFrame(results)
    
    # Summary statistics
    print(f"\n{'='*60}")
    print(f"PHOTO-Z MEMBERSHIP SUMMARY - {catalog_name.upper()}")
    print(f"{'='*60}")
    print(f"  Total groups: {len(results_df)}")
    print(f"  Groups with candidates: {np.sum(results_df['n_candidates_photoz'] > 0)}")
    print(f"  Groups with members (P>{prob_threshold}): {np.sum(results_df['n_members_photoz'] > 0)}")
    print(f"  Mean candidates per group: {results_df['n_candidates_photoz'].mean():.1f}")
    print(f"  Mean members per group: {results_df['n_members_photoz'].mean():.1f}")
    print(f"  Mean P_member: {results_df['mean_P_member'].mean():.3f}")
    print(f"{'='*60}\n")
    
    return results_df


def save_results(results_df, catalog_name, output_dir):
    """Save results to CSV and FITS."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # CSV
    csv_path = output_dir / f'photoz_membership_{catalog_name}_{timestamp}.csv'
    results_df.to_csv(csv_path, index=False)
    print(f"  Saved CSV: {csv_path}")
    
    # FITS
    fits_path = output_dir / f'photoz_membership_{catalog_name}_{timestamp}.fits'
    table = Table.from_pandas(results_df)
    table.write(fits_path, format='fits', overwrite=True)
    print(f"  Saved FITS: {fits_path}")
    
    return csv_path, fits_path


def main():
    parser = argparse.ArgumentParser(
        description='Determine photo-z group membership',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog', choices=['all', 'hcg', 'both'], required=True,
                       help='Which catalog to process')
    parser.add_argument('--radius', type=float, default=500,
                       help='Search radius in kpc')
    parser.add_argument('--max-dist-error', type=float, default=2.0,
                       help='Maximum distance in photo-z NMAD units (e.g., 2.0 = 2σ)')
    parser.add_argument('--prob-threshold', type=float, default=0.5,
                       help='Minimum P_member to be considered member')
    parser.add_argument('--test', action='store_true',
                       help='Test mode: process only first 10 groups')
    
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print("PHOTO-Z MEMBERSHIP DETERMINATION")
    print("="*60)
    print(f"\nConfiguration:")
    print(f"  Catalog: {args.catalog}")
    print(f"  Search radius: {args.radius} kpc")
    print(f"  Photo-z threshold: {args.max_dist_error}σ NMAD (2σ for compact groups)")
    print(f"  P_member threshold: {args.prob_threshold}")
    print(f"  Test mode: {args.test}")
    print(f"  Output directory: {OUTPUT_DIR}")
    print("="*60)
    
    # Process catalogs
    if args.catalog in ['all', 'both']:
        catalog_path = GROUP_CATALOG_DIR / 'cosmos_web_groups_catalog.csv'
        results_all = process_catalog(
            catalog_path, 'cw-all',
            radius_kpc=args.radius,
            max_dist_phot_z_error=args.max_dist_error,
            prob_threshold=args.prob_threshold,
            test_mode=args.test
        )
        save_results(results_all, 'cw-all', OUTPUT_DIR)
    
    if args.catalog in ['hcg', 'both']:
        # Final HCG catalog: Py18_Groups.csv or .fits
        hcg_csv = GROUP_CATALOG_DIR / 'Py18_Groups.csv'
        hcg_fits = GROUP_CATALOG_DIR / 'Py18_Groups.fits'
        if hcg_csv.exists():
            catalog_path = hcg_csv
            groups_hcg = pd.read_csv(hcg_csv)
        elif hcg_fits.exists():
            catalog_path = hcg_fits
            from astropy.table import Table
            groups_hcg = Table.read(hcg_fits).to_pandas()
        else:
            raise FileNotFoundError(f"HCG catalog not found. Expected {hcg_csv} or {hcg_fits}")
        results_hcg = process_catalog(
            catalog_path, 'cw-hcg',
            radius_kpc=args.radius,
            max_dist_phot_z_error=args.max_dist_error,
            prob_threshold=args.prob_threshold,
            test_mode=args.test,
            groups_df=groups_hcg
        )
        save_results(results_hcg, 'cw-hcg', OUTPUT_DIR)
    
    print("\n" + "="*60)
    print("PHOTO-Z MEMBERSHIP DETERMINATION COMPLETE")
    print("="*60 + "\n")


if __name__ == '__main__':
    main()
