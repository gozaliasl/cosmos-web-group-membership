"""
Main script to determine group membership for CW-All and CW-HCG catalogs.

Usage:
    python determine_membership.py --catalog all     # Process CW-All catalog
    python determine_membership.py --catalog hcg     # Process CW-HCG catalog
    python determine_membership.py --catalog both    # Process both catalogs
    python determine_membership.py --test            # Test on sample groups
"""

import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

from membership_functions import (
    find_specz_members,
    combine_specz_photoz_members,
    create_nmad_function
)

# Import improved photo-z membership function
from membership_funcs_improved import (
    find_photoz_members_improved,
    create_nmad_function_improved
)


# Paths
BASE_DIR = Path('/Users/gozalig1/Projects/compact-groups-xray-analysis')
GROUP_CATALOG_DIR = BASE_DIR / 'data' / 'group-catalog'
# Use the pre-matched catalog that has both spec-z and photo-z
SPECZ_CATALOG = BASE_DIR / 'data' / 'specz' / 'Webb_Specz_with_photz.csv'
PHOTOZ_CATALOG = BASE_DIR / 'data' / 'galaxy_catalog_photz' / 'galaxy_catalog_photz.csv'
OUTPUT_DIR = BASE_DIR / 'membership_determination' / 'results'
MATCHED_CATALOG = BASE_DIR / 'data' / 'specz' / 'Webb_Specz_with_photz.csv'

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def add_vrf_refined_redshifts(summary, specz_catalog):
    """
    Add VRF-refined group redshifts and dynamical properties to summary table.
    
    VRF iteratively determines the true group center, providing a more accurate
    redshift estimate than the initial catalog value. This also extracts the
    full dynamical properties calculated by VRF.
    
    Parameters:
    -----------
    summary : DataFrame
        Summary table with Group_ID and Group_z columns
    specz_catalog : DataFrame
        Spectroscopic member catalog with VRF membership flags
        
    Returns:
    --------
    summary : DataFrame
        Summary table with added columns:
        - z_vrf_refined: VRF-refined redshift (median of members)
        - z_vrf_std: Standard deviation of member redshifts
        - dz_vrf: Redshift offset (z_refined - z_catalog)
        - dv_vrf_kms: Velocity offset in km/s
        - sigma_v_kms: Velocity dispersion (km/s)
        - sigma_v_err_kms: Error on velocity dispersion (km/s)
        - n_vrf_members: Number of VRF members used
        - vrf_quality: Quality flag (1=excellent, 2=good, 3=poor, 4=insufficient)
    """
    c_kms = 299792.458  # Speed of light in km/s
    
    # Determine member column name
    member_col = 'vrf_member' if 'vrf_member' in specz_catalog.columns else 'gapper_member'
    
    # Calculate refined redshifts and dynamical properties for each group
    refined_data = []
    
    for _, row in summary.iterrows():
        group_id = row['Group_ID']
        catalog_z = row['Group_z']
        
        # Get VRF members for this group
        group_members = specz_catalog[
            (specz_catalog['Group_Group_ID'] == group_id) & 
            (specz_catalog[member_col] == True)
        ]
        
        n_members = len(group_members)
        
        if n_members >= 3:  # Need at least 3 members for basic estimates
            member_z = group_members['zfin'].values
            
            # Use median as robust estimator (resistant to outliers)
            z_refined = np.median(member_z)
            z_std = np.std(member_z)
            
            # Calculate offset
            dz = z_refined - catalog_z
            dv_kms = c_kms * dz / (1 + catalog_z)
            
            # Calculate velocity dispersion from member redshifts
            # Convert redshifts to velocities relative to group center
            member_velocities = c_kms * (member_z - z_refined) / (1 + z_refined)
            
            # Robust velocity dispersion estimate
            sigma_v = np.std(member_velocities)
            
            # Error estimate (bootstrapped approximation)
            if n_members >= 5:
                sigma_v_err = sigma_v / np.sqrt(2 * (n_members - 1))
            else:
                sigma_v_err = sigma_v / np.sqrt(n_members)
            
            # Quality flag based on number of members
            if n_members >= 15:
                quality = 1  # Excellent: well-constrained dynamics
            elif n_members >= 8:
                quality = 2  # Good: reasonable constraints
            elif n_members >= 5:
                quality = 3  # Poor: marginal constraints
            else:
                quality = 4  # Insufficient: unreliable
            
            refined_data.append({
                'Group_ID': group_id,
                'z_vrf_refined': z_refined,
                'z_vrf_std': z_std,
                'dz_vrf': dz,
                'dv_vrf_kms': dv_kms,
                'sigma_v_kms': sigma_v,
                'sigma_v_err_kms': sigma_v_err,
                'n_vrf_members': n_members,
                'vrf_quality': quality
            })
        else:
            # Not enough members for refinement
            refined_data.append({
                'Group_ID': group_id,
                'z_vrf_refined': catalog_z,
                'z_vrf_std': np.nan,
                'dz_vrf': 0.0,
                'dv_vrf_kms': 0.0,
                'sigma_v_kms': np.nan,
                'sigma_v_err_kms': np.nan,
                'n_vrf_members': n_members,
                'vrf_quality': 4  # Insufficient data
            })
    
    # Merge with summary
    refined_df = pd.DataFrame(refined_data)
    summary = summary.merge(refined_df, on='Group_ID', how='left')
    
    return summary


def reprocess_all_with_refined_redshifts(summary, specz, photoz, nmad_func, 
                                        radius_kpc=500, method='vrf'):
    """
    Reprocess ALL groups using VRF-refined redshifts to create complete refined catalog.
    
    This creates a parallel catalog where every group uses z_vrf_refined instead
    of the original catalog redshift, regardless of the offset magnitude.
    
    Parameters:
    -----------
    summary : DataFrame
        Summary table with z_vrf_refined and dv_vrf_kms columns
    specz, photoz : DataFrame
        Galaxy catalogs
    nmad_func : function
        NMAD function for photo-z
    radius_kpc : float
        Search radius in kpc
    method : str
        Method for spec-z membership
        
    Returns:
    --------
    vrf_members : DataFrame
        All members using refined redshift
    vrf_specz : DataFrame
        Spec-z members using refined redshift
    vrf_photoz : DataFrame
        Photo-z members using refined redshift
    vrf_summary : DataFrame
        Summary with refined membership counts
    """
    
    print(f"\nReprocessing ALL groups with VRF-refined redshifts...")
    print(f"  Creating complete refined catalog for {len(summary)} groups")
    
    all_members_list = []
    specz_members_list = []
    photoz_members_list = []
    summary_list = []
    
    for idx, row in tqdm(summary.iterrows(), 
                        total=len(summary), 
                        desc="Reprocessing all groups"):
        
        group_id = row['Group_ID']
        group_ra = row['Group_Ra']
        group_dec = row['Group_Dec']
        
        # Use refined redshift
        group_z = row['z_vrf_refined']
        original_z = row['Group_z']
        dv_offset = row['dv_vrf_kms']
        
        # Find spec-z members with refined redshift
        specz_mem = find_specz_members(
            group_ra, group_dec, group_z, specz,
            radius_kpc=radius_kpc,
            max_dz_norm=0.01,  # Tighter cut for compact groups: ~1800 km/s
            max_velocity=2000,  # Max 2000 km/s initial cut
            use_gapper=True,
            remove_duplicates=True,
            method=method
        )
        
        # Find photo-z members with refined redshift
        # Use VRF member count if available, otherwise use total candidates
        if method == 'vrf' and len(specz_mem) > 0:
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            n_specz_members = len(specz_mem[specz_mem[member_col]]) if member_col in specz_mem.columns else len(specz_mem)
        else:
            n_specz_members = len(specz_mem)
        
        if n_specz_members > 100:
            prob_thresh = 0.2
        elif n_specz_members > 50:
            prob_thresh = 0.15
        elif n_specz_members > 15:
            prob_thresh = 0.1
        else:
            prob_thresh = 0.05
        
        photoz_mem = find_photoz_members_improved(
            group_ra, group_dec, group_z, photoz,
            radius_kpc=radius_kpc,
            nmad_function=nmad_func,
            max_dist_phot_z_error=2.0,
            prob_threshold=prob_thresh,
            remove_duplicates=True,
            use_field_correction=True,
            use_color_weighting=False
        )
        
        # Process members
        if len(specz_mem) > 0:
            specz_mem_copy = specz_mem.copy()
            if 'ez' in specz_mem_copy.columns:
                specz_mem_copy = specz_mem_copy.sort_values('ez')
                specz_mem_copy = specz_mem_copy.drop_duplicates(
                    subset=['RA', 'DEC'], keep='first'
                )
            specz_mem_copy['Group_Group_ID'] = group_id
            specz_mem_copy['Group_Ra'] = group_ra
            specz_mem_copy['Group_Dec'] = group_dec
            specz_mem_copy['Group_z_catalog'] = original_z  # Original for reference
            specz_mem_copy['Group_z'] = group_z  # Refined z used
            specz_mem_copy['dv_vrf_kms'] = dv_offset
            specz_mem_copy['redshift_type'] = 'spec-z'
            specz_members_list.append(specz_mem_copy)
        
        if len(photoz_mem) > 0:
            photoz_mem_copy = photoz_mem.copy()
            photoz_mem_copy['Group_Group_ID'] = group_id
            photoz_mem_copy['Group_Ra'] = group_ra
            photoz_mem_copy['Group_Dec'] = group_dec
            photoz_mem_copy['Group_z_catalog'] = original_z
            photoz_mem_copy['Group_z'] = group_z
            photoz_mem_copy['dv_vrf_kms'] = dv_offset
            photoz_mem_copy['redshift_type'] = 'photo-z'
            photoz_members_list.append(photoz_mem_copy)
        
        # Combine members
        combined_mem = combine_specz_photoz_members(specz_mem, photoz_mem)
        
        if len(combined_mem) > 0:
            combined_mem_copy = combined_mem.copy()
            combined_mem_copy['Group_Group_ID'] = group_id
            combined_mem_copy['Group_Ra'] = group_ra
            combined_mem_copy['Group_Dec'] = group_dec
            combined_mem_copy['Group_z_catalog'] = original_z
            combined_mem_copy['Group_z'] = group_z
            combined_mem_copy['dv_vrf_kms'] = dv_offset
            all_members_list.append(combined_mem_copy)
        
        # Summary statistics
        n_specz = len(specz_mem)
        n_photoz = len(photoz_mem)
        n_total = len(combined_mem)
        mean_prob = photoz_mem['membership_prob'].mean() if len(photoz_mem) > 0 else 0.0
        
        # Calculate velocity dispersion for VRF members
        if n_specz > 0 and method == 'vrf':
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            vrf_members = specz_mem[specz_mem[member_col]] if member_col in specz_mem.columns else specz_mem
            
            if len(vrf_members) >= 3:
                c_kms = 299792.458
                member_z = vrf_members['zfin'].values
                z_center = np.median(member_z)
                member_velocities = c_kms * (member_z - z_center) / (1 + z_center)
                sigma_v = np.std(member_velocities)
                sigma_v_err = sigma_v / np.sqrt(2 * (len(vrf_members) - 1)) if len(vrf_members) >= 5 else sigma_v / np.sqrt(len(vrf_members))
                n_vrf = len(vrf_members)
                
                # Quality flag
                if n_vrf >= 15:
                    quality = 1
                elif n_vrf >= 8:
                    quality = 2
                elif n_vrf >= 5:
                    quality = 3
                else:
                    quality = 4
            else:
                sigma_v = np.nan
                sigma_v_err = np.nan
                n_vrf = len(vrf_members)
                quality = 4
        else:
            sigma_v = np.nan
            sigma_v_err = np.nan
            n_vrf = n_specz
            quality = 4
        
        summary_list.append({
            'Group_ID': group_id,
            'Group_Ra': group_ra,
            'Group_Dec': group_dec,
            'Group_z_catalog': original_z,
            'Group_z': group_z,
            'dv_vrf_kms': dv_offset,
            'z_vrf_std': row['z_vrf_std'],
            'sigma_v_kms': sigma_v,
            'sigma_v_err_kms': sigma_v_err,
            'n_vrf_members': n_vrf,
            'vrf_quality': quality,
            'n_specz_members': n_specz,
            'n_photoz_members': n_photoz,
            'n_total_members': n_total,
            'mean_photoz_prob': mean_prob
        })
    
    # Compile results
    if len(all_members_list) > 0:
        vrf_members = pd.concat(all_members_list, ignore_index=True)
    else:
        vrf_members = pd.DataFrame()
    
    if len(specz_members_list) > 0:
        vrf_specz = pd.concat(specz_members_list, ignore_index=True)
    else:
        vrf_specz = pd.DataFrame()
    
    if len(photoz_members_list) > 0:
        vrf_photoz = pd.concat(photoz_members_list, ignore_index=True)
    else:
        vrf_photoz = pd.DataFrame()
    
    if len(summary_list) > 0:
        vrf_summary = pd.DataFrame(summary_list)
    else:
        vrf_summary = pd.DataFrame()
    
    print(f"  VRF-refined catalog complete:")
    print(f"    Total members: {len(vrf_members)}")
    print(f"    Spec-z: {len(vrf_specz)}, Photo-z: {len(vrf_photoz)}")
    
    return vrf_members, vrf_specz, vrf_photoz, vrf_summary


def reprocess_with_refined_redshifts(summary, specz, photoz, nmad_func, 
                                     radius_kpc=500, method='vrf',
                                     dv_threshold=500):
    """
    Reprocess membership determination using VRF-refined redshifts.
    
    For groups where VRF provides a significantly different redshift
    (|Δv| > threshold), re-run membership determination using the refined z.
    This can improve both spec-z and photo-z membership.
    
    Parameters:
    -----------
    summary : DataFrame
        Summary table with z_vrf_refined and dv_vrf_kms columns
    specz, photoz : DataFrame
        Galaxy catalogs
    nmad_func : function
        NMAD function for photo-z
    radius_kpc : float
        Search radius in kpc
    method : str
        Method for spec-z membership
    dv_threshold : float
        Velocity offset threshold (km/s) for reprocessing
        Groups with |Δv| > threshold will be reprocessed
        
    Returns:
    --------
    vrf_members : DataFrame
        All members with refined redshift
    vrf_specz : DataFrame
        Spec-z members with refined redshift
    vrf_photoz : DataFrame
        Photo-z members with refined redshift
    vrf_summary : DataFrame
        Updated summary with refined membership counts
    """
    
    print(f"\nReprocessing groups with VRF-refined redshifts (|Δv| > {dv_threshold} km/s)...")
    
    # Identify groups needing reprocessing
    needs_reprocess = np.abs(summary['dv_vrf_kms']) > dv_threshold
    n_reprocess = needs_reprocess.sum()
    
    print(f"  Groups to reprocess: {n_reprocess} / {len(summary)}")
    if n_reprocess == 0:
        print("  No groups need reprocessing")
        return None, None, None, None
    
    all_members_list = []
    specz_members_list = []
    photoz_members_list = []
    summary_list = []
    
    for idx, row in tqdm(summary[needs_reprocess].iterrows(), 
                        total=n_reprocess, 
                        desc="Reprocessing"):
        
        group_id = row['Group_ID']
        group_ra = row['Group_Ra']
        group_dec = row['Group_Dec']
        
        # Use refined redshift
        group_z = row['z_vrf_refined']
        original_z = row['Group_z']
        
        # Find spec-z members with refined redshift
        specz_mem = find_specz_members(
            group_ra, group_dec, group_z, specz,
            radius_kpc=radius_kpc,
            max_dz_norm=0.01,  # Tighter cut for compact groups: ~1800 km/s
            max_velocity=2000,  # Max 2000 km/s initial cut
            use_gapper=True,
            remove_duplicates=True,
            method=method
        )
        
        # Find photo-z members with refined redshift
        # Use VRF member count if available, otherwise use total candidates
        if method == 'vrf' and len(specz_mem) > 0:
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            n_specz_members = len(specz_mem[specz_mem[member_col]]) if member_col in specz_mem.columns else len(specz_mem)
        else:
            n_specz_members = len(specz_mem)
        
        if n_specz_members > 100:
            prob_thresh = 0.2
        elif n_specz_members > 50:
            prob_thresh = 0.15
        elif n_specz_members > 15:
            prob_thresh = 0.1
        else:
            prob_thresh = 0.05
        
        photoz_mem = find_photoz_members_improved(
            group_ra, group_dec, group_z, photoz,
            radius_kpc=radius_kpc,
            nmad_function=nmad_func,
            max_dist_phot_z_error=2.0,
            prob_threshold=prob_thresh,
            remove_duplicates=True,
            use_field_correction=True,
            use_color_weighting=False
        )
        
        # Process members
        if len(specz_mem) > 0:
            specz_mem_copy = specz_mem.copy()
            if 'ez' in specz_mem_copy.columns:
                specz_mem_copy = specz_mem_copy.sort_values('ez')
                specz_mem_copy = specz_mem_copy.drop_duplicates(
                    subset=['RA', 'DEC'], keep='first'
                )
            specz_mem_copy['Group_Group_ID'] = group_id
            specz_mem_copy['Group_Ra'] = group_ra
            specz_mem_copy['Group_Dec'] = group_dec
            specz_mem_copy['Group_z'] = original_z  # Keep original for reference
            specz_mem_copy['Group_z_refined'] = group_z  # Add refined z
            specz_mem_copy['redshift_type'] = 'spec-z'
            specz_members_list.append(specz_mem_copy)
        
        if len(photoz_mem) > 0:
            photoz_mem_copy = photoz_mem.copy()
            photoz_mem_copy['Group_Group_ID'] = group_id
            photoz_mem_copy['Group_Ra'] = group_ra
            photoz_mem_copy['Group_Dec'] = group_dec
            photoz_mem_copy['Group_z'] = original_z
            photoz_mem_copy['Group_z_refined'] = group_z
            photoz_mem_copy['redshift_type'] = 'photo-z'
            photoz_members_list.append(photoz_mem_copy)
        
        # Combine members
        combined_mem = combine_specz_photoz_members(specz_mem, photoz_mem)
        
        if len(combined_mem) > 0:
            combined_mem_copy = combined_mem.copy()
            combined_mem_copy['Group_Group_ID'] = group_id
            combined_mem_copy['Group_Ra'] = group_ra
            combined_mem_copy['Group_Dec'] = group_dec
            combined_mem_copy['Group_z'] = original_z
            combined_mem_copy['Group_z_refined'] = group_z
            all_members_list.append(combined_mem_copy)
        
        # Summary statistics
        n_specz = len(specz_mem)
        n_photoz = len(photoz_mem)
        n_total = len(combined_mem)
        mean_prob = photoz_mem['membership_prob'].mean() if len(photoz_mem) > 0 else 0.0
        
        summary_list.append({
            'Group_ID': group_id,
            'Group_Ra': group_ra,
            'Group_Dec': group_dec,
            'Group_z': original_z,
            'Group_z_refined': group_z,
            'dv_vrf_kms': row['dv_vrf_kms'],
            'n_specz_members': n_specz,
            'n_photoz_members': n_photoz,
            'n_total_members': n_total,
            'mean_photoz_prob': mean_prob,
            'reprocessed': True
        })
    
    # Compile results
    if len(all_members_list) > 0:
        vrf_members = pd.concat(all_members_list, ignore_index=True)
    else:
        vrf_members = pd.DataFrame()
    
    if len(specz_members_list) > 0:
        vrf_specz = pd.concat(specz_members_list, ignore_index=True)
    else:
        vrf_specz = pd.DataFrame()
    
    if len(photoz_members_list) > 0:
        vrf_photoz = pd.concat(photoz_members_list, ignore_index=True)
    else:
        vrf_photoz = pd.DataFrame()
    
    if len(summary_list) > 0:
        vrf_summary = pd.DataFrame(summary_list)
    else:
        vrf_summary = pd.DataFrame()
    
    print(f"  Reprocessed members: {len(vrf_members)}")
    print(f"  Spec-z: {len(vrf_specz)}, Photo-z: {len(vrf_photoz)}")
    
    return vrf_members, vrf_specz, vrf_photoz, vrf_summary


def merge_vrf_refined_catalog(original_members, original_summary,
                              vrf_members, vrf_summary):
    """
    Create final VRF-refined catalog by merging original and reprocessed results.
    
    For groups that were reprocessed, use the VRF-refined membership.
    For groups that weren't reprocessed, use the original membership.
    
    Parameters:
    -----------
    original_members, original_summary : DataFrame
        Original membership results
    vrf_members, vrf_summary : DataFrame
        Reprocessed membership with refined redshifts
        
    Returns:
    --------
    final_members : DataFrame
        Combined catalog
    final_summary : DataFrame
        Combined summary
    """
    
    if vrf_summary is None or len(vrf_summary) == 0:
        print("  No reprocessing needed, using original catalog")
        return original_members, original_summary
    
    print("\nMerging VRF-refined catalog...")
    
    # Get list of reprocessed group IDs
    reprocessed_ids = set(vrf_summary['Group_ID'].values)
    
    # Determine column name for Group_ID in members catalog
    group_id_col = 'Group_Group_ID' if 'Group_Group_ID' in original_members.columns else 'Group_ID'
    
    # Keep original results for groups NOT reprocessed
    keep_original = ~original_members[group_id_col].isin(reprocessed_ids)
    original_kept = original_members[keep_original].copy()
    
    # Add Group_z_refined column to original (same as Group_z for these)
    if 'Group_z_refined' not in original_kept.columns:
        original_kept['Group_z_refined'] = original_kept['Group_z']
    
    # Combine
    final_members = pd.concat([original_kept, vrf_members], ignore_index=True)
    
    # Merge summaries
    original_summary_kept = original_summary[
        ~original_summary['Group_ID'].isin(reprocessed_ids)
    ].copy()
    
    # Add columns to match vrf_summary
    if 'Group_z_refined' not in original_summary_kept.columns:
        original_summary_kept['Group_z_refined'] = original_summary_kept['Group_z']
    if 'reprocessed' not in original_summary_kept.columns:
        original_summary_kept['reprocessed'] = False
    
    final_summary = pd.concat([original_summary_kept, vrf_summary], ignore_index=True)
    final_summary = final_summary.sort_values('Group_ID').reset_index(drop=True)
    
    n_reprocessed = (final_summary['reprocessed'] == True).sum()
    n_total = len(final_summary)
    
    print(f"  Final catalog: {len(final_members)} members")
    print(f"  Groups: {n_total} total, {n_reprocessed} reprocessed, " +
          f"{n_total - n_reprocessed} original")
    
    return final_members, final_summary


def _load_hcg_catalog():
    """Load final HCG group catalog (Py18_Groups.csv or .fits)."""
    hcg_csv = GROUP_CATALOG_DIR / 'Py18_Groups.csv'
    hcg_fits = GROUP_CATALOG_DIR / 'Py18_Groups.fits'
    if hcg_csv.exists():
        return pd.read_csv(hcg_csv)
    if hcg_fits.exists():
        from astropy.table import Table
        return Table.read(hcg_fits).to_pandas()
    raise FileNotFoundError(f"HCG catalog not found. Expected {hcg_csv} or {hcg_fits}")


def load_catalogs():
    """Load all necessary catalogs."""
    print("Loading catalogs...")
    
    # Group catalogs
    cw_all = pd.read_csv(GROUP_CATALOG_DIR / 'cosmos_web_groups_catalog.csv')
    cw_hcg = _load_hcg_catalog()
    
    # Galaxy catalogs
    specz = pd.read_csv(SPECZ_CATALOG)
    photoz = pd.read_csv(PHOTOZ_CATALOG)
    
    # Matched catalog for NMAD function
    if MATCHED_CATALOG.exists():
        matched = pd.read_csv(MATCHED_CATALOG)
    else:
        matched = None
    
    print(f"  CW-All groups: {len(cw_all)}")
    print(f"  CW-HCG groups: {len(cw_hcg)}")
    print(f"  Spec-z galaxies: {len(specz)}")
    print(f"  Photo-z galaxies: {len(photoz)}")
    
    return cw_all, cw_hcg, specz, photoz, matched


def process_cw_all_catalog(cw_all, specz, photoz, nmad_func, 
                           radius_kpc=500, sample_size=None, method='gapper'):
    """
    Process CW-All catalog to find members for all groups.
    
    Parameters:
    -----------
    cw_all : DataFrame
        CW-All group catalog
    specz : DataFrame
        Spec-z catalog
    photoz : DataFrame
        Photo-z catalog
    nmad_func : function
        NMAD function for photo-z
    radius_kpc : float
        Search radius in kpc
    sample_size : int or None
        Process only this many groups (for testing)
    method : str
        Method for spec-z membership: 'gapper' or 'vrf'
        
    Returns:
    --------
    all_members : DataFrame
        All members (combined)
    specz_only : DataFrame
        Spec-z members only
    photoz_only : DataFrame
        Photo-z members only
    summary : DataFrame
        Summary statistics per group
    """
    
    print(f"\nProcessing CW-All catalog (radius={radius_kpc} kpc)...")
    
    if sample_size is not None:
        cw_all = cw_all.head(sample_size)
        print(f"  Processing {sample_size} groups (test mode)")
    
    all_members_list = []
    specz_members_list = []
    photoz_members_list = []
    summary_list = []
    
    for idx, row in tqdm(cw_all.iterrows(), total=len(cw_all), desc="Processing groups"):
        group_id = row['Group_ID']
        group_ra = row['Ra']
        group_dec = row['Dec']
        group_z = row['z']
        
        # Skip very high redshift groups where photo-z might be unreliable
        if group_z > 4.0:
            continue
        
        # Find spec-z members
        specz_mem = find_specz_members(
            group_ra, group_dec, group_z, specz,
            radius_kpc=radius_kpc,
            max_dz_norm=0.01,  # Tighter cut for compact groups: |Δz/(1+z)| < 0.01 (~1800 km/s)
            max_velocity=2000,  # Max 2000 km/s for initial membership cut
            use_gapper=True,
            remove_duplicates=True,
            method=method  # Use specified method (gapper or vrf)
        )
        
        # Find photo-z members using improved method
        # Use 3σ redshift cylinder for cleaner selection
        # Probability threshold adapts to spec-z richness
        
        # Determine probability threshold based on spec-z richness
        # Use VRF member count if available, otherwise use total candidates
        if method == 'vrf' and len(specz_mem) > 0:
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            n_specz_members = len(specz_mem[specz_mem[member_col]]) if member_col in specz_mem.columns else len(specz_mem)
        else:
            n_specz_members = len(specz_mem)
        
        if n_specz_members > 100:  # Rich groups - be selective
            prob_thresh = 0.2
        elif n_specz_members > 50:  # Moderate richness
            prob_thresh = 0.15
        elif n_specz_members > 15:  # Small groups - adjusted for compact groups
            prob_thresh = 0.1
        else:  # Poor/compact groups (HCG typically 4-10 members) - be inclusive
            prob_thresh = 0.05
        
        photoz_mem = find_photoz_members_improved(
            group_ra, group_dec, group_z, photoz,
            radius_kpc=radius_kpc,  # Keep 500 kpc for wide net, filter by sep_kpc later
            nmad_function=nmad_func,
            max_dist_phot_z_error=2.0,  # 2σ cylinder for compact groups (was 3σ)
            prob_threshold=prob_thresh,  # Adaptive threshold based on richness
            remove_duplicates=True,
            use_field_correction=True,  # Enable background correction
            use_color_weighting=False  # Disable color (no color columns)
        )
        
        # Add group info to separate catalogs
        if len(specz_mem) > 0:
            specz_mem_copy = specz_mem.copy()
            
            # Handle spec-z duplicates: keep the one with lowest error
            if 'ez' in specz_mem_copy.columns:
                # Sort by redshift error and drop duplicates based on coordinates
                specz_mem_copy = specz_mem_copy.sort_values('ez')
                coords_str = specz_mem_copy['RA'].round(5).astype(str) + '_' + \
                            specz_mem_copy['DEC'].round(5).astype(str)
                specz_mem_copy = specz_mem_copy[~coords_str.duplicated(keep='first')]
            
            # Add all group properties
            for col in cw_all.columns:
                specz_mem_copy[f'Group_{col}'] = row[col]
            specz_mem_copy['catalog'] = 'CW-All'
            specz_mem_copy['redshift_type'] = 'spec-z'  # Flag for redshift type
            
            # Add membership quality flag for spec-z members
            # Handle both gapper and VRF methods
            member_col = 'vrf_member' if 'vrf_member' in specz_mem_copy.columns else 'gapper_member'
            
            specz_mem_copy['membership_flag'] = 3  # Default: low confidence
            
            # High confidence (flag=1): Member + |Δv| < 500 km/s + |Δz/(1+z)| < 0.01
            high_conf_mask = (
                specz_mem_copy[member_col] & 
                (np.abs(specz_mem_copy['dv']) < 500) & 
                (np.abs(specz_mem_copy['dz_norm']) < 0.01)
            )
            specz_mem_copy.loc[high_conf_mask, 'membership_flag'] = 1
            
            # Medium confidence (flag=2): Member OR reasonable offset
            medium_conf_mask = (
                (specz_mem_copy[member_col]) | 
                ((np.abs(specz_mem_copy['dv']) < 1000) & (np.abs(specz_mem_copy['dz_norm']) < 0.015))
            ) & (~high_conf_mask)
            specz_mem_copy.loc[medium_conf_mask, 'membership_flag'] = 2
            
            specz_members_list.append(specz_mem_copy)
        
        if len(photoz_mem) > 0:
            photoz_mem_copy = photoz_mem.copy()
            
            # Photo-z catalog is already clean (no duplicates expected)
            # But add duplicate check just in case
            if 'RA_MODEL' in photoz_mem_copy.columns:
                coords_str = photoz_mem_copy['RA_MODEL'].round(5).astype(str) + '_' + \
                            photoz_mem_copy['DEC_MODEL'].round(5).astype(str)
                photoz_mem_copy = photoz_mem_copy[~coords_str.duplicated(keep='first')]
            
            # Add all group properties
            for col in cw_all.columns:
                photoz_mem_copy[f'Group_{col}'] = row[col]
            photoz_mem_copy['catalog'] = 'CW-All'
            photoz_mem_copy['redshift_type'] = 'photo-z'  # Flag for redshift type
            
            # Add membership quality flag for photo-z members
            # Flag based on membership probability
            photoz_mem_copy['membership_flag'] = 3  # Default: low confidence
            
            # High confidence (flag=1): P > 0.5 (very likely member)
            photoz_mem_copy.loc[photoz_mem_copy['membership_prob'] > 0.5, 'membership_flag'] = 1
            
            # Medium confidence (flag=2): 0.2 < P <= 0.5 (probable member)
            medium_mask = (photoz_mem_copy['membership_prob'] > 0.2) & (photoz_mem_copy['membership_prob'] <= 0.5)
            photoz_mem_copy.loc[medium_mask, 'membership_flag'] = 2
            
            # Low confidence (flag=3): 0.05 < P <= 0.2 (possible member)
            # Already set as default
            
            photoz_members_list.append(photoz_mem_copy)
        
        # Combine for summary and create combined catalog with redshift_type flag
        if len(specz_mem) > 0 or len(photoz_mem) > 0:
            combined = combine_specz_photoz_members(specz_mem, photoz_mem)
            combined['Group_ID'] = group_id
            combined['Group_Ra'] = group_ra
            combined['Group_Dec'] = group_dec
            combined['Group_z'] = group_z
            combined['catalog'] = 'CW-All'
            
            # Add redshift_type flag based on source
            # 'source' column is added by combine_specz_photoz_members: 'specz' or 'photoz'
            if 'source' in combined.columns:
                combined['redshift_type'] = combined['source'].map({
                    'specz': 'spec-z',
                    'photoz': 'photo-z'
                })
            
            all_members_list.append(combined)
            
            # Summary - handle both VRF and gapper
            if len(specz_mem) > 0:
                member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
                n_specz = len(specz_mem[specz_mem[member_col]])
            else:
                n_specz = 0
            n_photoz = len(photoz_mem) if len(photoz_mem) > 0 else 0
            n_total = len(combined)
            mean_prob = photoz_mem['membership_prob'].mean() if len(photoz_mem) > 0 else np.nan
            
            summary_list.append({
                'Group_ID': group_id,
                'Group_Ra': group_ra,
                'Group_Dec': group_dec,
                'Group_z': group_z,
                'n_specz_members': n_specz,
                'n_photoz_members': n_photoz,
                'n_total_members': n_total,
                'mean_photoz_prob': mean_prob,
                'catalog': 'CW-All'
            })
    
    if len(all_members_list) > 0:
        all_members = pd.concat(all_members_list, ignore_index=True)
        summary = pd.DataFrame(summary_list)
    else:
        all_members = pd.DataFrame()
        summary = pd.DataFrame()
    
    if len(specz_members_list) > 0:
        specz_only = pd.concat(specz_members_list, ignore_index=True)
        
        # Add VRF-refined redshifts to summary (if we have members)
        if len(summary) > 0:
            summary = add_vrf_refined_redshifts(summary, specz_only)
    else:
        specz_only = pd.DataFrame()
    
    if len(photoz_members_list) > 0:
        photoz_only = pd.concat(photoz_members_list, ignore_index=True)
    else:
        photoz_only = pd.DataFrame()
    
    print(f"  Found members for {len(summary)} groups")
    print(f"  Total members: {len(all_members)}")
    print(f"  Spec-z members: {len(specz_only)}")
    print(f"  Photo-z members: {len(photoz_only)}")
    
    return all_members, specz_only, photoz_only, summary


def process_cw_hcg_catalog(cw_hcg, specz, photoz, nmad_func,
                           radius_kpc=500, sample_size=None, method='gapper'):
    """
    Process CW-HCG catalog to find members for all groups.
    
    Parameters:
    -----------
    cw_hcg : DataFrame
        CW-HCG group catalog
    specz : DataFrame
        Spec-z catalog
    photoz : DataFrame
        Photo-z catalog
    nmad_func : function
        NMAD function for photo-z
    radius_kpc : float
        Search radius in kpc
    sample_size : int or None
        Process only this many groups (for testing)
    method : str
        Method for spec-z membership: 'gapper' or 'vrf'
        
    Returns:
    --------
    all_members : DataFrame
        All members (combined)
    specz_only : DataFrame
        Spec-z members only
    photoz_only : DataFrame
        Photo-z members only
    summary : DataFrame
        Summary statistics per group
    """
    
    print(f"\nProcessing CW-HCG catalog (radius={radius_kpc} kpc)...")
    
    if sample_size is not None:
        cw_hcg = cw_hcg.head(sample_size)
        print(f"  Processing {sample_size} groups (test mode)")
    
    all_members_list = []
    specz_members_list = []
    photoz_members_list = []
    summary_list = []
    
    for idx, row in tqdm(cw_hcg.iterrows(), total=len(cw_hcg), desc="Processing groups"):
        group_id = row['Grp']
        group_ra = row['Ra']
        group_dec = row['Dec']
        group_z = row['z']
        
        # Skip very high redshift groups
        if group_z > 4.0:
            continue
        
        # Find spec-z members
        specz_mem = find_specz_members(
            group_ra, group_dec, group_z, specz,
            radius_kpc=radius_kpc,
            max_dz_norm=0.01,  # Tighter cut for compact groups: |Δz/(1+z)| < 0.01 (~1800 km/s)
            max_velocity=2000,  # Max 2000 km/s for initial membership cut
            use_gapper=True,
            remove_duplicates=True,
            method=method  # Use specified method (gapper or vrf)
        )
        
        # Find photo-z members using improved method
        # Use 3σ redshift cylinder for cleaner selection
        # Probability threshold adapts to spec-z richness
        
        # Determine probability threshold based on spec-z richness
        # Use VRF member count if available, otherwise use total candidates
        if method == 'vrf' and len(specz_mem) > 0:
            member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
            n_specz_members = len(specz_mem[specz_mem[member_col]]) if member_col in specz_mem.columns else len(specz_mem)
        else:
            n_specz_members = len(specz_mem)
        
        if n_specz_members > 100:  # Rich groups - be selective
            prob_thresh = 0.2
        elif n_specz_members > 50:  # Moderate richness
            prob_thresh = 0.15
        elif n_specz_members > 15:  # Small groups - adjusted for compact groups
            prob_thresh = 0.1
        else:  # Poor/compact groups (HCG typically 4-10 members) - be inclusive
            prob_thresh = 0.05
        
        photoz_mem = find_photoz_members_improved(
            group_ra, group_dec, group_z, photoz,
            radius_kpc=radius_kpc,  # Keep 500 kpc for wide net, filter by sep_kpc later
            nmad_function=nmad_func,
            max_dist_phot_z_error=2.0,  # 2σ cylinder for compact groups (was 3σ)
            prob_threshold=prob_thresh,  # Adaptive threshold based on richness
            remove_duplicates=True,
            use_field_correction=True,  # Enable background correction
            use_color_weighting=False  # Disable color (no color columns)
        )
        
        # Add group info to separate catalogs
        if len(specz_mem) > 0:
            specz_mem_copy = specz_mem.copy()
            
            # Handle spec-z duplicates: keep the one with lowest error
            if 'ez' in specz_mem_copy.columns:
                # Sort by redshift error and drop duplicates based on coordinates
                specz_mem_copy = specz_mem_copy.sort_values('ez')
                coords_str = specz_mem_copy['RA'].round(5).astype(str) + '_' + \
                            specz_mem_copy['DEC'].round(5).astype(str)
                specz_mem_copy = specz_mem_copy[~coords_str.duplicated(keep='first')]
            
            # Add all group properties
            for col in cw_hcg.columns:
                specz_mem_copy[f'Group_{col}'] = row[col]
            specz_mem_copy['catalog'] = 'CW-HCG'
            specz_mem_copy['redshift_type'] = 'spec-z'  # Flag for redshift type
            
            # Add membership quality flag for spec-z members
            # Handle both gapper and VRF methods
            member_col = 'vrf_member' if 'vrf_member' in specz_mem_copy.columns else 'gapper_member'
            
            specz_mem_copy['membership_flag'] = 3  # Default: low confidence
            
            # High confidence (flag=1): Member + |Δv| < 500 km/s + |Δz/(1+z)| < 0.01
            high_conf_mask = (
                specz_mem_copy[member_col] & 
                (np.abs(specz_mem_copy['dv']) < 500) & 
                (np.abs(specz_mem_copy['dz_norm']) < 0.01)
            )
            specz_mem_copy.loc[high_conf_mask, 'membership_flag'] = 1
            
            # Medium confidence (flag=2): Member OR reasonable offset
            medium_conf_mask = (
                (specz_mem_copy[member_col]) | 
                ((np.abs(specz_mem_copy['dv']) < 1000) & (np.abs(specz_mem_copy['dz_norm']) < 0.015))
            ) & (~high_conf_mask)
            specz_mem_copy.loc[medium_conf_mask, 'membership_flag'] = 2
            
            specz_members_list.append(specz_mem_copy)
        
        if len(photoz_mem) > 0:
            photoz_mem_copy = photoz_mem.copy()
            
            # Photo-z catalog duplicate check
            if 'RA_MODEL' in photoz_mem_copy.columns:
                coords_str = photoz_mem_copy['RA_MODEL'].round(5).astype(str) + '_' + \
                            photoz_mem_copy['DEC_MODEL'].round(5).astype(str)
                photoz_mem_copy = photoz_mem_copy[~coords_str.duplicated(keep='first')]
            
            # Add all group properties
            for col in cw_hcg.columns:
                photoz_mem_copy[f'Group_{col}'] = row[col]
            photoz_mem_copy['catalog'] = 'CW-HCG'
            photoz_mem_copy['redshift_type'] = 'photo-z'  # Flag for redshift type
            
            # Add membership quality flag for photo-z members
            photoz_mem_copy['membership_flag'] = 3  # Default: low confidence
            
            # High confidence (flag=1): P > 0.5
            photoz_mem_copy.loc[photoz_mem_copy['membership_prob'] > 0.5, 'membership_flag'] = 1
            
            # Medium confidence (flag=2): 0.2 < P <= 0.5
            medium_mask = (photoz_mem_copy['membership_prob'] > 0.2) & (photoz_mem_copy['membership_prob'] <= 0.5)
            photoz_mem_copy.loc[medium_mask, 'membership_flag'] = 2
            
            photoz_members_list.append(photoz_mem_copy)
        
        # Combine for summary and create combined catalog with redshift_type flag
        if len(specz_mem) > 0 or len(photoz_mem) > 0:
            combined = combine_specz_photoz_members(specz_mem, photoz_mem)
            combined['Group_ID'] = group_id
            combined['Group_Ra'] = group_ra
            combined['Group_Dec'] = group_dec
            combined['Group_z'] = group_z
            combined['catalog'] = 'CW-HCG'
            
            # Add redshift_type flag based on source
            if 'source' in combined.columns:
                combined['redshift_type'] = combined['source'].map({
                    'specz': 'spec-z',
                    'photoz': 'photo-z'
                })
            
            all_members_list.append(combined)
            
            # Summary - handle both VRF and gapper
            if len(specz_mem) > 0:
                member_col = 'vrf_member' if 'vrf_member' in specz_mem.columns else 'gapper_member'
                n_specz = len(specz_mem[specz_mem[member_col]])
            else:
                n_specz = 0
            n_photoz = len(photoz_mem) if len(photoz_mem) > 0 else 0
            n_total = len(combined)
            mean_prob = photoz_mem['membership_prob'].mean() if len(photoz_mem) > 0 else np.nan
            
            summary_list.append({
                'Group_ID': group_id,
                'Group_Ra': group_ra,
                'Group_Dec': group_dec,
                'Group_z': group_z,
                'n_specz_members': n_specz,
                'n_photoz_members': n_photoz,
                'n_total_members': n_total,
                'mean_photoz_prob': mean_prob,
                'catalog': 'CW-HCG'
            })
    
    if len(all_members_list) > 0:
        all_members = pd.concat(all_members_list, ignore_index=True)
        summary = pd.DataFrame(summary_list)
    else:
        all_members = pd.DataFrame()
        summary = pd.DataFrame()
    
    if len(specz_members_list) > 0:
        specz_only = pd.concat(specz_members_list, ignore_index=True)
        
        # Add VRF-refined redshifts to summary (if we have members)
        if len(summary) > 0:
            summary = add_vrf_refined_redshifts(summary, specz_only)
    else:
        specz_only = pd.DataFrame()
    
    if len(photoz_members_list) > 0:
        photoz_only = pd.concat(photoz_members_list, ignore_index=True)
    else:
        photoz_only = pd.DataFrame()
    
    print(f"  Found members for {len(summary)} groups")
    print(f"  Total members: {len(all_members)}")
    print(f"  Spec-z members: {len(specz_only)}")
    print(f"  Photo-z members: {len(photoz_only)}")
    
    return all_members, specz_only, photoz_only, summary


def plot_membership_statistics(summary_all, summary_hcg, output_dir):
    """Create summary plots of membership determination."""
    
    print("\nCreating summary plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Number of members vs redshift (CW-All)
    ax = axes[0, 0]
    if len(summary_all) > 0:
        ax.scatter(summary_all['Group_z'], summary_all['n_total_members'],
                  alpha=0.5, s=30, label='Total', c='blue')
        ax.scatter(summary_all['Group_z'], summary_all['n_specz_members'],
                  alpha=0.5, s=20, label='Spec-z', marker='s', c='red')
        if 'n_photoz_members' in summary_all.columns:
            photoz_data = summary_all[summary_all['n_photoz_members'] > 0]
            if len(photoz_data) > 0:
                ax.scatter(photoz_data['Group_z'], photoz_data['n_photoz_members'],
                          alpha=0.5, s=15, label='Photo-z', marker='^', c='green')
        ax.set_xlabel('Group Redshift', fontsize=11)
        ax.set_ylabel('Number of Members', fontsize=11)
        ax.set_title('CW-All: Members vs Redshift', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No CW-All data', ha='center', va='center', transform=ax.transAxes)
        ax.grid(alpha=0.3)
    
    # Plot 2: Number of members vs redshift (CW-HCG)
    ax = axes[0, 1]
    if len(summary_hcg) > 0:
        ax.scatter(summary_hcg['Group_z'], summary_hcg['n_total_members'],
                  alpha=0.5, s=30, label='Total', c='blue')
        ax.scatter(summary_hcg['Group_z'], summary_hcg['n_specz_members'],
                  alpha=0.5, s=20, label='Spec-z', marker='s', c='red')
        if 'n_photoz_members' in summary_hcg.columns:
            photoz_data = summary_hcg[summary_hcg['n_photoz_members'] > 0]
            if len(photoz_data) > 0:
                ax.scatter(photoz_data['Group_z'], photoz_data['n_photoz_members'],
                          alpha=0.5, s=15, label='Photo-z', marker='^', c='green')
        ax.set_xlabel('Group Redshift', fontsize=11)
        ax.set_ylabel('Number of Members', fontsize=11)
        ax.set_title('CW-HCG: Members vs Redshift', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No CW-HCG data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('CW-HCG: Members vs Redshift', fontsize=12, fontweight='bold')
    
    # Plot 3: Histogram of total members
    ax = axes[0, 2]
    has_data = False
    bins = np.arange(0, 150, 5)
    if len(summary_all) > 0:
        ax.hist(summary_all['n_total_members'], bins=bins, alpha=0.6,
               label=f'CW-All (N={len(summary_all)})', edgecolor='black')
        has_data = True
    if len(summary_hcg) > 0:
        ax.hist(summary_hcg['n_total_members'], bins=bins, alpha=0.6,
               label=f'CW-HCG (N={len(summary_hcg)})', edgecolor='black')
        has_data = True
    
    if has_data:
        ax.set_xlabel('Number of Members', fontsize=11)
        ax.set_ylabel('Number of Groups', fontsize=11)
        ax.set_title('Distribution of Member Counts', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
    
    # Plot 4: Spec-z vs Photo-z members (CW-All)
    ax = axes[1, 0]
    if len(summary_all) > 0 and 'n_photoz_members' in summary_all.columns:
        photoz_data = summary_all[summary_all['n_photoz_members'] > 0]
        if len(photoz_data) > 0:
            ax.scatter(summary_all['n_specz_members'], summary_all['n_photoz_members'],
                      alpha=0.5, s=30, c='blue')
            max_val = max(summary_all['n_specz_members'].max(),
                         summary_all['n_photoz_members'].max())
            ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='1:1 line')
            ax.set_xlabel('Spec-z Members', fontsize=11)
            ax.set_ylabel('Photo-z Members', fontsize=11)
            ax.set_title('CW-All: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No photo-z members', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('CW-All: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('CW-All: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
    
    # Plot 5: Spec-z vs Photo-z members (CW-HCG)
    ax = axes[1, 1]
    if len(summary_hcg) > 0 and 'n_photoz_members' in summary_hcg.columns:
        photoz_data = summary_hcg[summary_hcg['n_photoz_members'] > 0]
        if len(photoz_data) > 0:
            ax.scatter(summary_hcg['n_specz_members'], summary_hcg['n_photoz_members'],
                      alpha=0.5, s=30, c='blue')
            max_val = max(summary_hcg['n_specz_members'].max(),
                         summary_hcg['n_photoz_members'].max())
            ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='1:1 line')
            ax.set_xlabel('Spec-z Members', fontsize=11)
            ax.set_ylabel('Photo-z Members', fontsize=11)
            ax.set_title('CW-HCG: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
            ax.legend()
            ax.grid(alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No photo-z members', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('CW-HCG: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No CW-HCG data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('CW-HCG: Spec-z vs Photo-z', fontsize=12, fontweight='bold')
    
    # Plot 6: Mean photo-z probability distribution
    ax = axes[1, 2]
    has_data = False
    bins = np.linspace(0, 1, 20)
    
    if len(summary_all) > 0 and 'mean_photoz_prob' in summary_all.columns:
        all_probs = summary_all['mean_photoz_prob'].dropna()
        if len(all_probs) > 0:
            ax.hist(all_probs, bins=bins, alpha=0.6, label='CW-All', edgecolor='black')
            has_data = True
    
    if len(summary_hcg) > 0 and 'mean_photoz_prob' in summary_hcg.columns:
        hcg_probs = summary_hcg['mean_photoz_prob'].dropna()
        if len(hcg_probs) > 0:
            ax.hist(hcg_probs, bins=bins, alpha=0.6, label='CW-HCG', edgecolor='black')
            has_data = True
    
    if has_data:
        ax.set_xlabel('Mean Photo-z Probability', fontsize=11)
        ax.set_ylabel('Number of Groups', fontsize=11)
        ax.set_title('Distribution of Mean Probabilities', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No probability data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Distribution of Mean Probabilities', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'membership_summary.png', dpi=200, bbox_inches='tight')
    print(f"Saved: {output_dir / 'membership_summary.png'}")
    plt.close()


def clean_old_results(catalog_choice, method, radius_kpc):
    """
    Remove old result files before starting a new run.
    
    This prevents accidentally using stale results from previous runs
    with different parameters or buggy code.
    
    Parameters:
    -----------
    catalog_choice : str
        'all', 'hcg', or 'both'
    method : str
        'gapper' or 'vrf'
    radius_kpc : float
        Radius in kpc
    """
    import glob
    
    print("\n" + "="*80)
    print("CLEANING OLD RESULT FILES")
    print("="*80)
    
    method_tag = '_vrf' if method == 'vrf' else ''
    radius_tag = f'_r{radius_kpc:.0f}kpc'
    
    patterns_to_remove = []
    
    if catalog_choice in ['all', 'both']:
        patterns_to_remove.extend([
            f'cw_all_members_*{method_tag}{radius_tag}.csv',
            f'cw_all_summary*{method_tag}{radius_tag}.csv',
        ])
    
    if catalog_choice in ['hcg', 'both']:
        patterns_to_remove.extend([
            f'cw_hcg_members_*{method_tag}{radius_tag}.csv',
            f'cw_hcg_summary*{method_tag}{radius_tag}.csv',
        ])
    
    removed_count = 0
    for pattern in patterns_to_remove:
        full_pattern = str(OUTPUT_DIR / pattern)
        matching_files = glob.glob(full_pattern)
        
        for filepath in matching_files:
            try:
                Path(filepath).unlink()
                print(f"  Removed: {Path(filepath).name}")
                removed_count += 1
            except Exception as e:
                print(f"  Warning: Could not remove {Path(filepath).name}: {e}")
    
    if removed_count == 0:
        print("  No old result files found")
    else:
        print(f"  Total files removed: {removed_count}")
    
    print("="*80 + "\n")


def main():
    """Main execution function."""
    
    parser = argparse.ArgumentParser(description='Determine group membership')
    parser.add_argument('--catalog', choices=['all', 'hcg', 'both'], default='both',
                       help='Which catalog to process')
    parser.add_argument('--test', action='store_true',
                       help='Test mode: process only 10 groups')
    parser.add_argument('--radius', type=float, default=500,
                       help='Physical radius in kpc (default: 500)')
    parser.add_argument('--method', choices=['gapper', 'vrf'], default='gapper',
                       help='Method for spec-z membership: gapper (default) or vrf (VRF interloper rejection)')
    parser.add_argument('--vrf-refine', action='store_true',
                       help='Reprocess groups with refined VRF redshifts (requires --method vrf)')
    parser.add_argument('--dv-threshold', type=float, default=500,
                       help='Velocity threshold (km/s) for VRF refinement (default: 500)')
    parser.add_argument('--keep-old-results', action='store_true',
                       help='Keep old result files instead of removing them')
    
    args = parser.parse_args()
    
    # Clean old results unless --keep-old-results is specified
    if not args.keep_old_results:
        clean_old_results(args.catalog, args.method, args.radius)
    else:
        print("\n" + "="*80)
        print("KEEPING OLD RESULT FILES (--keep-old-results specified)")
        print("="*80 + "\n")
    
    # Display method info
    if args.method == 'vrf':
        print("\n" + "="*80)
        print("USING VRF INTERLOPER REJECTION FOR SPEC-Z MEMBERSHIP")
        print("="*80)
        print("Method: MBM10 algorithm (Mamon, Biviano & Murante 2010)")
        print("Features:")
        print("  - NFW-based velocity dispersion profile")
        print("  - Iterative interloper rejection")
        print("  - Radius-dependent velocity cuts")
        print("  - Self-consistent mass estimation")
        print("="*80 + "\n")
    else:
        print("\n" + "="*80)
        print("USING GAPPER METHOD FOR SPEC-Z MEMBERSHIP (DEFAULT)")
        print("="*80)
        print("Method: Gapper technique with velocity cuts")
        print("="*80 + "\n")
    
    # Load catalogs
    cw_all, cw_hcg, specz, photoz, matched = load_catalogs()
    
    # Create NMAD function using improved version
    print("\nCreating NMAD function...")
    nmad_func = create_nmad_function_improved(matched, use_magnitude=True)
    
    # Test mode
    sample_size = 10 if args.test else None
    
    # Process catalogs
    summary_all = pd.DataFrame()
    summary_hcg = pd.DataFrame()
    
    if args.catalog in ['all', 'both']:
        members_all, specz_all, photoz_all, summary_all = process_cw_all_catalog(
            cw_all, specz, photoz, nmad_func,
            radius_kpc=args.radius,
            sample_size=sample_size,
            method=args.method  # Pass method selection
        )
        
        # Create complete VRF-refined catalog (always for VRF method)
        if args.method == 'vrf' and len(summary_all) > 0:
            print("\n" + "="*80)
            print("CREATING VRF-REFINED CATALOG (ALL GROUPS)")
            print("="*80)
            print("This creates a parallel catalog using VRF-refined redshifts for ALL groups")
            
            # Process all groups with refined redshifts
            vrf_members_complete, vrf_specz_complete, vrf_photoz_complete, vrf_summary_complete = \
                reprocess_all_with_refined_redshifts(
                    summary_all, specz, photoz, nmad_func,
                    radius_kpc=args.radius,
                    method='vrf'
                )
            
            # Save complete VRF-refined catalog
            if len(vrf_members_complete) > 0:
                print("\nSaving VRF-refined catalogs (using z_vrf_refined for ALL groups)...")
                
                # Combined catalog
                refined_file = OUTPUT_DIR / f'cw_all_members_combined_vrf_refined_r{args.radius:.0f}kpc.csv'
                vrf_members_complete.to_csv(refined_file, index=False)
                print(f"  Saved: {refined_file}")
                
                # Spec-z only
                if len(vrf_specz_complete) > 0:
                    specz_file = OUTPUT_DIR / f'cw_all_members_specz_vrf_refined_r{args.radius:.0f}kpc.csv'
                    vrf_specz_complete.to_csv(specz_file, index=False)
                    print(f"  Saved: {specz_file}")
                
                # Photo-z only
                if len(vrf_photoz_complete) > 0:
                    photoz_file = OUTPUT_DIR / f'cw_all_members_photoz_vrf_refined_r{args.radius:.0f}kpc.csv'
                    vrf_photoz_complete.to_csv(photoz_file, index=False)
                    print(f"  Saved: {photoz_file}")
                
                # Summary
                summary_file = OUTPUT_DIR / f'cw_all_summary_vrf_refined_r{args.radius:.0f}kpc.csv'
                vrf_summary_complete.to_csv(summary_file, index=False)
                print(f"  Saved: {summary_file}")
                
                # Print comparison
                print("\n" + "-"*80)
                print("COMPARISON: Catalog Redshift vs VRF-Refined Redshift")
                print("-"*80)
                orig_total = summary_all['n_total_members'].sum()
                vrf_total = vrf_summary_complete['n_total_members'].sum()
                diff = vrf_total - orig_total
                pct_diff = 100 * diff / orig_total if orig_total > 0 else 0
                print(f"  Catalog redshift:     {orig_total:,} total members")
                print(f"  VRF-refined redshift: {vrf_total:,} total members")
                print(f"  Difference:           {diff:+,} ({pct_diff:+.1f}%)")
                print("-"*80)
        
        # Save original catalog (using catalog redshifts)
        if len(members_all) > 0:
            print("\nSaving original catalogs (using catalog redshifts)...")
            # Include method in filename
            method_tag = '_vrf' if args.method == 'vrf' else ''
            # Combined catalog
            members_file = OUTPUT_DIR / f'cw_all_members_combined{method_tag}_r{args.radius:.0f}kpc.csv'
            members_all.to_csv(members_file, index=False)
            print(f"  Saved: {members_file}")
            
            # Spec-z only catalog
            if len(specz_all) > 0:
                specz_file = OUTPUT_DIR / f'cw_all_members_specz{method_tag}_r{args.radius:.0f}kpc.csv'
                specz_all.to_csv(specz_file, index=False)
                print(f"Saved: {specz_file}")
            
            # Photo-z only catalog
            if len(photoz_all) > 0:
                photoz_file = OUTPUT_DIR / f'cw_all_members_photoz{method_tag}_r{args.radius:.0f}kpc.csv'
                photoz_all.to_csv(photoz_file, index=False)
                print(f"Saved: {photoz_file}")
            
            # Summary
            summary_file = OUTPUT_DIR / f'cw_all_summary{method_tag}_r{args.radius:.0f}kpc.csv'
            summary_all.to_csv(summary_file, index=False)
            print(f"Saved: {summary_file}")
    
    if args.catalog in ['hcg', 'both']:
        members_hcg, specz_hcg, photoz_hcg, summary_hcg = process_cw_hcg_catalog(
            cw_hcg, specz, photoz, nmad_func,
            radius_kpc=args.radius,
            sample_size=sample_size,
            method=args.method  # Pass method selection
        )
        
        # Create complete VRF-refined catalog (always for VRF method)
        if args.method == 'vrf' and len(summary_hcg) > 0:
            print("\n" + "="*80)
            print("CREATING VRF-REFINED CATALOG (HCG - ALL GROUPS)")
            print("="*80)
            print("This creates a parallel catalog using VRF-refined redshifts for ALL groups")
            
            # Process all groups with refined redshifts
            vrf_members_complete, vrf_specz_complete, vrf_photoz_complete, vrf_summary_complete = \
                reprocess_all_with_refined_redshifts(
                    summary_hcg, specz, photoz, nmad_func,
                    radius_kpc=args.radius,
                    method='vrf'
                )
            
            # Save complete VRF-refined catalog
            if len(vrf_members_complete) > 0:
                print("\nSaving VRF-refined catalogs (using z_vrf_refined for ALL groups)...")
                
                # Combined catalog
                refined_file = OUTPUT_DIR / f'cw_hcg_members_combined_vrf_refined_r{args.radius:.0f}kpc.csv'
                vrf_members_complete.to_csv(refined_file, index=False)
                print(f"  Saved: {refined_file}")
                
                # Spec-z only
                if len(vrf_specz_complete) > 0:
                    specz_file = OUTPUT_DIR / f'cw_hcg_members_specz_vrf_refined_r{args.radius:.0f}kpc.csv'
                    vrf_specz_complete.to_csv(specz_file, index=False)
                    print(f"  Saved: {specz_file}")
                
                # Photo-z only
                if len(vrf_photoz_complete) > 0:
                    photoz_file = OUTPUT_DIR / f'cw_hcg_members_photoz_vrf_refined_r{args.radius:.0f}kpc.csv'
                    vrf_photoz_complete.to_csv(photoz_file, index=False)
                    print(f"  Saved: {photoz_file}")
                
                # Summary
                summary_file = OUTPUT_DIR / f'cw_hcg_summary_vrf_refined_r{args.radius:.0f}kpc.csv'
                vrf_summary_complete.to_csv(summary_file, index=False)
                print(f"  Saved: {summary_file}")
                
                # Print comparison
                print("\n" + "-"*80)
                print("COMPARISON: Catalog Redshift vs VRF-Refined Redshift")
                print("-"*80)
                orig_total = summary_hcg['n_total_members'].sum()
                vrf_total = vrf_summary_complete['n_total_members'].sum()
                diff = vrf_total - orig_total
                pct_diff = 100 * diff / orig_total if orig_total > 0 else 0
                print(f"  Catalog redshift:     {orig_total:,} total members")
                print(f"  VRF-refined redshift: {vrf_total:,} total members")
                print(f"  Difference:           {diff:+,} ({pct_diff:+.1f}%)")
                print("-"*80)
        
        # Save original catalog (using catalog redshifts)
        if len(members_hcg) > 0:
            print("\nSaving original catalogs (using catalog redshifts)...")
            # Include method in filename
            method_tag = '_vrf' if args.method == 'vrf' else ''
            # Combined catalog
            members_file = OUTPUT_DIR / f'cw_hcg_members_combined{method_tag}_r{args.radius:.0f}kpc.csv'
            members_hcg.to_csv(members_file, index=False)
            print(f"Saved: {members_file}")
            
            # Spec-z only catalog
            if len(specz_hcg) > 0:
                specz_file = OUTPUT_DIR / f'cw_hcg_members_specz{method_tag}_r{args.radius:.0f}kpc.csv'
                specz_hcg.to_csv(specz_file, index=False)
                print(f"Saved: {specz_file}")
            
            # Photo-z only catalog
            if len(photoz_hcg) > 0:
                photoz_file = OUTPUT_DIR / f'cw_hcg_members_photoz{method_tag}_r{args.radius:.0f}kpc.csv'
                photoz_hcg.to_csv(photoz_file, index=False)
                print(f"Saved: {photoz_file}")
            
            # Summary
            summary_file = OUTPUT_DIR / f'cw_hcg_summary_r{args.radius:.0f}kpc.csv'
            summary_hcg.to_csv(summary_file, index=False)
            print(f"Saved: {summary_file}")
    
    # Create summary plots
    if len(summary_all) > 0 or len(summary_hcg) > 0:
        plot_membership_statistics(summary_all, summary_hcg, OUTPUT_DIR)
    
    # Print final summary
    print("\n" + "="*60)
    print("MEMBERSHIP DETERMINATION COMPLETE")
    print("="*60)
    
    if len(summary_all) > 0:
        print(f"\nCW-All:")
        print(f"  Groups with members: {len(summary_all)}")
        print(f"  Total members found: {summary_all['n_total_members'].sum():.0f}")
        print(f"  Mean members per group: {summary_all['n_total_members'].mean():.1f}")
        print(f"  Median members per group: {summary_all['n_total_members'].median():.0f}")
    
    if len(summary_hcg) > 0:
        print(f"\nCW-HCG:")
        print(f"  Groups with members: {len(summary_hcg)}")
        print(f"  Total members found: {summary_hcg['n_total_members'].sum():.0f}")
        print(f"  Mean members per group: {summary_hcg['n_total_members'].mean():.1f}")
        print(f"  Median members per group: {summary_hcg['n_total_members'].median():.0f}")


if __name__ == '__main__':
    main()
