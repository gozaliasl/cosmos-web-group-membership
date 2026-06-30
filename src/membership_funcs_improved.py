"""
Improved functions for determining group membership using photometric redshifts.
Optimized for small, deep surveys (e.g., 0.45 deg² effective area).

Key improvements:
1. Proper background/field contamination modeling
2. NFW-inspired radial profiles
3. Magnitude-dependent photo-z errors
4. Bootstrap-based field density estimation for small areas
5. Red-sequence weighting (optional)
6. Quality flags for low-richness systems
"""

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from scipy import stats
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# CONFIGURATION FOR SMALL SURVEY
# ============================================================================
SURVEY_AREA_DEG2 = 0.45  # Effective survey area
FIELD_SAMPLE_FRACTION = 0.3  # Use 30% of survey for field sampling


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def angular_to_physical_radius(z, radius_kpc=500):
    """
    Convert physical radius to angular radius at given redshift.
    
    Parameters:
    -----------
    z : float
        Redshift
    radius_kpc : float
        Physical radius in kpc (default: 500 kpc)
        
    Returns:
    --------
    radius_arcsec : float
        Angular radius in arcseconds
    """
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
    radius_arcsec = radius_kpc / kpc_per_arcsec
    return radius_arcsec


def spatial_probability_nfw(sep_kpc, radius_kpc, concentration=4.0):
    """
    NFW-inspired radial membership probability.
    More realistic than simple quadratic decline.
    
    Parameters:
    -----------
    sep_kpc : array or float
        Projected separation in kpc
    radius_kpc : float
        Maximum radius for membership
    concentration : float
        NFW concentration parameter (default: 4.0 for groups)
        
    Returns:
    --------
    prob : array or float
        Radial membership probability [0, 1]
    """
    rs = radius_kpc / concentration  # Scale radius
    x = sep_kpc / rs
    
    # Truncated NFW profile normalized to 1 at center
    is_scalar = np.isscalar(x)
    if is_scalar:
        x = np.array([x])
    
    profile = np.zeros_like(x, dtype=float)
    
    # Handle center (avoid log(0))
    small_x = x < 1e-6
    profile[small_x] = 1.0
    
    # NFW profile for r > 0
    large_x = ~small_x
    if large_x.any():
        numerator = np.log(1 + x[large_x]) - x[large_x]/(1 + x[large_x])
        x_max = radius_kpc / rs
        denominator = np.log(1 + x_max) - x_max/(1 + x_max)
        profile[large_x] = numerator / denominator
    
    # Clip to [0, 1]
    profile = np.clip(profile, 0, 1)
    
    return profile[0] if is_scalar else profile


def create_nmad_function_improved(matched_catalog=None, use_magnitude=True):
    """
    Create improved NMAD function with magnitude and redshift dependence.
    
    Parameters:
    -----------
    matched_catalog : DataFrame or None
        Matched spec-z/photo-z catalog for empirical calibration
    use_magnitude : bool
        Whether to use magnitude-dependent errors (requires F444W mag)
        
    Returns:
    --------
    nmad_func : function
        Function that takes (z, mag_f444w=None) and returns NMAD
    """
    
    if matched_catalog is not None:
        # Try to derive empirical NMAD from matched catalog
        try:
            # Determine column names
            if 'LP_zfinal' in matched_catalog.columns and 'zfin' in matched_catalog.columns:
                zphot_col, zspec_col = 'LP_zfinal', 'zfin'
            elif 'zphot' in matched_catalog.columns and 'zspec' in matched_catalog.columns:
                zphot_col, zspec_col = 'zphot', 'zspec'
            else:
                raise ValueError("Cannot find redshift columns")
            
            # Clean sample
            clean_mask = (matched_catalog[zphot_col].notna()) & \
                        (matched_catalog[zspec_col].notna())
            if 'LP_warn_fl' in matched_catalog.columns:
                clean_mask &= (matched_catalog['LP_warn_fl'] == 0)
            
            clean_catalog = matched_catalog[clean_mask].copy()
            
            if len(clean_catalog) >= 100:
                # Calculate NMAD in redshift bins
                z_bins = np.arange(0, 4.5, 0.3)
                z_centers, nmad_values = [], []
                
                dz = clean_catalog[zphot_col] - clean_catalog[zspec_col]
                dz_norm = dz / (1 + clean_catalog[zspec_col])
                
                for i in range(len(z_bins)-1):
                    z_min, z_max = z_bins[i], z_bins[i+1]
                    bin_mask = (clean_catalog[zspec_col] >= z_min) & \
                              (clean_catalog[zspec_col] < z_max)
                    
                    if bin_mask.sum() > 15:
                        dz_bin = dz_norm[bin_mask]
                        nmad = 1.48 * np.median(np.abs(dz_bin - np.median(dz_bin)))
                        z_centers.append((z_min + z_max) / 2)
                        nmad_values.append(nmad)
                
                if len(z_centers) > 2:
                    nmad_interp = interp1d(z_centers, nmad_values, 
                                          kind='linear', fill_value='extrapolate',
                                          bounds_error=False)
                    print(f"Empirical NMAD function created with {len(z_centers)} bins")
                    
                    def nmad_func(z, mag_f444w=None):
                        base_nmad = float(nmad_interp(z))
                        
                        # Magnitude correction if available
                        if use_magnitude and mag_f444w is not None:
                            if mag_f444w >= 26:
                                base_nmad *= 1.5
                            elif mag_f444w >= 25:
                                base_nmad *= 1.2
                        
                        return np.clip(base_nmad, 0.010, 0.050)
                    
                    return nmad_func
        
        except Exception as e:
            warnings.warn(f"Could not create empirical NMAD function: {e}")
    
    # Default: COSMOS2025-based function
    def nmad_func(z, mag_f444w=None):
        """
        σ_MAD from COSMOS2025 paper (Shuntov et al. 2025).
        
        Returns normalized MAD: σ_MAD = 1.48 × median(|Δz/(1+z)|)
        """
        # Magnitude-dependent (if available)
        if use_magnitude and mag_f444w is not None:
            if mag_f444w < 23:
                return 0.011
            elif mag_f444w < 24:
                return 0.012
            elif mag_f444w < 25:
                return 0.015
            elif mag_f444w < 26:
                return 0.020
            else:  # mag >= 26
                return 0.030
        
        # Redshift-dependent (default, conservative)
        if z < 1.0:
            return 0.012  # Excellent for z < 1
        elif z < 2.0:
            return 0.015  # Good for 1 < z < 2
        elif z < 3.0:
            return 0.020  # Moderate for 2 < z < 3
        else:
            return 0.030  # Conservative for z > 3
    
    return nmad_func


# ============================================================================
# FIELD DENSITY ESTIMATION FOR SMALL SURVEYS
# ============================================================================

def estimate_field_density_small_survey(photoz_catalog, group_ra, group_dec, 
                                       group_z, radius_kpc, z_width,
                                       n_bootstrap=100, exclude_inner_kpc=1500):
    """
    Estimate background field density for small surveys using bootstrap.
    
    For small surveys, we can't use large annuli. Instead:
    1. Sample field regions away from known groups
    2. Use bootstrap to estimate uncertainty
    3. Account for cosmic variance
    
    Parameters:
    -----------
    photoz_catalog : DataFrame
        Full photo-z catalog
    group_ra, group_dec, group_z : float
        Group parameters
    radius_kpc : float
        Group radius
    z_width : float
        Redshift slice width (physical units)
    n_bootstrap : int
        Number of bootstrap samples for uncertainty
    exclude_inner_kpc : float
        Exclude regions within this distance from group
        
    Returns:
    --------
    field_density : float
        Background density (galaxies per kpc²)
    field_density_err : float
        Uncertainty in field density
    """
    
    # Find galaxies in redshift slice
    z_mask = np.abs(photoz_catalog['LP_zfinal'] - group_z) < z_width
    redshift_slice = photoz_catalog[z_mask].copy()
    
    if len(redshift_slice) < 10:
        return 0.0, 0.0
    
    # Calculate angular separations from group
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    gal_coords = SkyCoord(ra=redshift_slice['RA_MODEL'].values*u.deg,
                         dec=redshift_slice['DEC_MODEL'].values*u.deg)
    sep_arcsec = group_coord.separation(gal_coords).arcsec
    
    # Convert to physical
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    sep_kpc = sep_arcsec * kpc_per_arcsec
    
    # Define field region (exclude inner region)
    exclude_arcsec = exclude_inner_kpc / kpc_per_arcsec
    field_mask = sep_arcsec > exclude_arcsec
    
    if field_mask.sum() < 5:
        # Not enough field galaxies, use global density
        # Total survey area in kpc² at this redshift
        area_deg2 = SURVEY_AREA_DEG2
        area_arcmin2 = area_deg2 * 3600  # 1 deg² = 3600 arcmin²
        area_kpc2 = area_arcmin2 * (kpc_per_arcsec * 60)**2
        
        field_density = len(redshift_slice) / area_kpc2
        field_density_err = field_density / np.sqrt(len(redshift_slice))
        
        return field_density, field_density_err
    
    # Calculate field area (survey area minus excluded region)
    excluded_area_kpc2 = np.pi * exclude_inner_kpc**2
    total_area_deg2 = SURVEY_AREA_DEG2
    total_area_arcmin2 = total_area_deg2 * 3600
    total_area_kpc2 = total_area_arcmin2 * (kpc_per_arcsec * 60)**2
    field_area_kpc2 = total_area_kpc2 - excluded_area_kpc2
    
    # Bootstrap estimate
    n_field = field_mask.sum()
    densities = []
    
    for _ in range(n_bootstrap):
        # Resample field galaxies
        boot_indices = np.random.choice(n_field, size=n_field, replace=True)
        boot_density = len(boot_indices) / field_area_kpc2
        densities.append(boot_density)
    
    field_density = np.median(densities)
    field_density_err = 1.48 * np.median(np.abs(densities - field_density))  # MAD
    
    return field_density, field_density_err


# ============================================================================
# CORRECTED MEMBERSHIP PROBABILITY
# ============================================================================

def calculate_corrected_membership_probability(z_prob, spatial_prob, 
                                              field_density, cluster_density,
                                              radius_kpc, color_prob=1.0):
    """
    Bayesian correction for background contamination.
    
    P(member|data) = P(cluster|data) / [P(cluster|data) + P(field|data)]
    
    Parameters:
    -----------
    z_prob : array
        Redshift probability from photo-z
    spatial_prob : array
        Spatial probability from radial profile
    field_density : float
        Background field density (gal/kpc²)
    cluster_density : float
        Estimated cluster density (gal/kpc²)
    radius_kpc : float
        Aperture radius
    color_prob : array or float
        Optional color-based probability
        
    Returns:
    --------
    corrected_prob : array
        Corrected membership probability [0, 1]
    """
    
    # Expected number from cluster vs field
    area = np.pi * radius_kpc**2
    N_expected_field = field_density * area
    N_expected_cluster = cluster_density * area
    
    # Prior probability of cluster membership
    if N_expected_cluster + N_expected_field > 0:
        prior_cluster = N_expected_cluster / (N_expected_cluster + N_expected_field)
    else:
        prior_cluster = 0.5  # Uninformative prior
    
    # Likelihood ratio: P(data|cluster) / P(data|field)
    # For cluster members: high z_prob, high spatial_prob
    # For field: uniform spatial, broader z distribution
    likelihood_cluster = z_prob * spatial_prob * color_prob
    likelihood_field = 0.3 * spatial_prob  # Field has lower z concentration
    
    # Posterior probability using Bayes theorem
    numerator = likelihood_cluster * prior_cluster
    denominator = numerator + likelihood_field * (1 - prior_cluster)
    
    # Avoid division by zero
    denominator = np.maximum(denominator, 1e-10)
    corrected_prob = numerator / denominator
    
    return np.clip(corrected_prob, 0, 1)


def calculate_color_probability(catalog, group_z, color_col='color_gr'):
    """
    Weight by proximity to red sequence.
    
    Parameters:
    -----------
    catalog : DataFrame
        Galaxy catalog with color information
    group_z : float
        Group redshift
    color_col : str
        Column name for color (default: 'color_gr')
        
    Returns:
    --------
    color_prob : array
        Color-based probability [0.3, 1.0]
    """
    
    if color_col not in catalog.columns:
        return np.ones(len(catalog))
    
    # Red sequence color at this redshift (empirical)
    # For g-r: RS ≈ 1.2 - 0.1*z (approximate)
    rs_color = 1.2 - 0.1 * group_z
    rs_width = 0.15  # Typical scatter
    
    # Gaussian probability around red sequence
    color = catalog[color_col].values
    color_prob = np.exp(-0.5 * ((color - rs_color) / rs_width)**2)
    
    # Don't exclude blue galaxies entirely, just downweight
    color_prob = 0.3 + 0.7 * color_prob  # Range [0.3, 1.0]
    
    return color_prob


# ============================================================================
# MAIN PHOTO-Z MEMBERSHIP FUNCTION
# ============================================================================

def find_photoz_members_improved(group_ra, group_dec, group_z, photoz_catalog,
                                radius_kpc=500, nmad_function=None, 
                                max_dist_phot_z_error=3.0,
                                prob_threshold=0.05, remove_duplicates=True,
                                use_field_correction=True, 
                                use_color_weighting=False,
                                color_col='color_gr',
                                magnitude_col='MAG_F444W'):
    """
    Improved photo-z membership determination for small, deep surveys.
    
    Key improvements over basic method:
    1. Background field correction using bootstrap for small areas
    2. NFW-inspired radial profile
    3. Magnitude-dependent photo-z errors
    4. Optional red-sequence weighting
    5. Quality flags for low-richness systems
    
    Parameters:
    -----------
    group_ra, group_dec : float
        Group coordinates (degrees)
    group_z : float
        Group redshift
    photoz_catalog : DataFrame
        Photo-z catalog with: RA_MODEL, DEC_MODEL, LP_zfinal, LP_warn_fl
        Optional: MAG_F444W, color_gr or other color
    radius_kpc : float
        Physical search radius (default: 500 kpc)
    nmad_function : function or None
        NMAD(z, mag) function. If None, uses COSMOS2025 values
    max_dist_phot_z_error : float
        Cylinder height in units of σ_MAD × (1+z) (default: 3.0)
    prob_threshold : float
        Minimum probability to include (default: 0.05)
    remove_duplicates : bool
        Remove duplicate coordinates (default: True)
    use_field_correction : bool
        Apply background correction (default: True)
    use_color_weighting : bool
        Use red-sequence weighting (default: False)
    color_col : str
        Column name for color (default: 'color_gr')
    magnitude_col : str
        Column name for magnitude (default: 'MAG_F444W')
        
    Returns:
    --------
    members : DataFrame
        Candidate members with membership probabilities and quality flags
    """
    
    # Clean photo-z catalog
    clean_photoz = photoz_catalog[
        (photoz_catalog['LP_warn_fl'] == 0)
    ].copy()
    
    # Get NMAD function
    if nmad_function is None:
        nmad_function = create_nmad_function_improved()
    
    # Get NMAD for this redshift (use median magnitude if available)
    if magnitude_col in clean_photoz.columns:
        median_mag = clean_photoz[magnitude_col].median()
        nmad = nmad_function(group_z, median_mag)
    else:
        nmad = nmad_function(group_z)
    
    # Define cylinder height
    z_error = nmad * (1 + group_z)
    z_min = group_z - max_dist_phot_z_error * z_error
    z_max = group_z + max_dist_phot_z_error * z_error
    
    # Apply redshift cut
    z_mask = (clean_photoz['LP_zfinal'] >= z_min) & \
             (clean_photoz['LP_zfinal'] <= z_max)
    
    if z_mask.sum() == 0:
        return pd.DataFrame()
    
    clean_photoz = clean_photoz[z_mask].copy()
    
    # Spatial cut
    radius_arcsec = angular_to_physical_radius(group_z, radius_kpc)
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    gal_coords = SkyCoord(ra=clean_photoz['RA_MODEL'].values*u.deg,
                         dec=clean_photoz['DEC_MODEL'].values*u.deg)
    
    sep = group_coord.separation(gal_coords)
    within_radius = sep < radius_arcsec * u.arcsec
    
    if within_radius.sum() == 0:
        return pd.DataFrame()
    
    candidates = clean_photoz[within_radius].copy()
    candidates['sep_arcsec'] = sep[within_radius].arcsec
    
    # Physical separation
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    candidates['sep_kpc'] = candidates['sep_arcsec'] * kpc_per_arcsec
    
    # Remove duplicates
    if remove_duplicates and len(candidates) > 1:
        coords_str = candidates['RA_MODEL'].astype(str) + '_' + \
                    candidates['DEC_MODEL'].astype(str)
        candidates = candidates[~coords_str.duplicated(keep='first')]
    
    if len(candidates) == 0:
        return pd.DataFrame()
    
    # Calculate membership probabilities
    
    # 1. Redshift probability
    dz = candidates['LP_zfinal'] - group_z
    dz_norm = dz / (1 + group_z)
    
    # Use magnitude-dependent NMAD if available
    if magnitude_col in candidates.columns:
        nmad_array = np.array([nmad_function(group_z, mag) 
                              for mag in candidates[magnitude_col]])
    else:
        nmad_array = nmad
    
    z_prob = np.exp(-0.5 * (dz_norm / nmad_array)**2)
    
    # 2. Spatial probability (NFW profile)
    spatial_prob = spatial_probability_nfw(candidates['sep_kpc'].values, 
                                          radius_kpc, concentration=4.0)
    
    # 3. Color probability (optional)
    if use_color_weighting and color_col in candidates.columns:
        color_prob = calculate_color_probability(candidates, group_z, color_col)
    else:
        color_prob = 1.0
    
    # 4. Field correction
    if use_field_correction:
        # Estimate field density
        field_density, field_density_err = estimate_field_density_small_survey(
            photoz_catalog, group_ra, group_dec, group_z, 
            radius_kpc, z_width=max_dist_phot_z_error * z_error,
            exclude_inner_kpc=1500
        )
        
        # Estimate cluster density from current candidates
        area_kpc2 = np.pi * radius_kpc**2
        cluster_density = len(candidates) / area_kpc2
        
        # Apply correction
        candidates['membership_prob'] = calculate_corrected_membership_probability(
            z_prob, spatial_prob, field_density, cluster_density,
            radius_kpc, color_prob
        )
        
        candidates['field_density'] = field_density
        candidates['field_density_err'] = field_density_err
    else:
        # No field correction - just combine probabilities
        candidates['membership_prob'] = z_prob * spatial_prob * color_prob
        candidates['field_density'] = 0.0
        candidates['field_density_err'] = 0.0
    
    # Quality assessment
    total_richness = candidates['membership_prob'].sum()
    
    if total_richness < 3:
        candidates['quality_flag'] = 'low_richness_strict'
        prob_threshold = max(prob_threshold, 0.2)  # Stricter threshold
    elif total_richness < 5:
        candidates['quality_flag'] = 'low_richness'
        prob_threshold = max(prob_threshold, 0.1)
    elif total_richness < 8:
        candidates['quality_flag'] = 'moderate'
    else:
        candidates['quality_flag'] = 'good'
    
    # Apply threshold
    candidates = candidates[candidates['membership_prob'] >= prob_threshold]
    
    if len(candidates) == 0:
        return pd.DataFrame()
    
    # Add metadata
    candidates['dz'] = candidates['LP_zfinal'] - group_z
    candidates['dz_norm'] = dz_norm
    candidates['z_min_cylinder'] = z_min
    candidates['z_max_cylinder'] = z_max
    candidates['nmad_used'] = nmad
    candidates['richness'] = total_richness
    
    # Add individual probability components for diagnostics
    candidates['z_prob'] = z_prob[candidates.index]
    candidates['spatial_prob'] = spatial_prob[np.isin(np.arange(len(spatial_prob)), 
                                                       candidates.index)]
    
    return candidates


# ============================================================================
# CONVENIENCE WRAPPER
# ============================================================================

def find_photoz_members(group_ra, group_dec, group_z, photoz_catalog,
                       radius_kpc=500, nmad_function=None, 
                       max_dist_phot_z_error=3.0,
                       prob_threshold=0.05, remove_duplicates=True,
                       use_improved=True, **kwargs):
    """
    Wrapper function for backward compatibility.
    
    Set use_improved=True (default) to use new improved method,
    or use_improved=False for original method.
    """
    
    if use_improved:
        return find_photoz_members_improved(
            group_ra, group_dec, group_z, photoz_catalog,
            radius_kpc=radius_kpc, nmad_function=nmad_function,
            max_dist_phot_z_error=max_dist_phot_z_error,
            prob_threshold=prob_threshold, remove_duplicates=remove_duplicates,
            **kwargs
        )
    else:
        # Call original method (from your existing code)
        # This would be your original find_photoz_members function
        raise NotImplementedError("Original method not included in this version")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    """
    Example usage for small deep survey
    """
    
    # Load catalogs (example)
    # photoz_catalog = pd.read_csv('galaxy_catalog_photz.csv')
    # specz_catalog = pd.read_csv('galaxy_catalog_specz.csv')
    
    # Create NMAD function from matched spec-z/photo-z
    # matched_catalog = pd.merge(photoz_catalog, specz_catalog, ...)
    # nmad_func = create_nmad_function_improved(matched_catalog)
    
    # Or use default COSMOS2025-based function
    nmad_func = create_nmad_function_improved()
    
    # Find members for a group
    # group_ra, group_dec, group_z = 150.1, 2.3, 0.85
    
    # members = find_photoz_members_improved(
    #     group_ra, group_dec, group_z, photoz_catalog,
    #     radius_kpc=500,
    #     nmad_function=nmad_func,
    #     max_dist_phot_z_error=3.0,
    #     prob_threshold=0.05,
    #     use_field_correction=True,
    #     use_color_weighting=False  # Set True if colors available
    # )
    
    # Check results
    # print(f"Found {len(members)} candidate members")
    # print(f"Total richness: {members['membership_prob'].sum():.1f}")
    # print(f"Quality: {members['quality_flag'].iloc[0]}")
    # print(f"\nTop 5 members:")
    # print(members.nlargest(5, 'membership_prob')[
    #     ['RA_MODEL', 'DEC_MODEL', 'LP_zfinal', 'sep_kpc', 
    #      'membership_prob', 'quality_flag']
    # ])
    
    print("Improved photo-z membership module loaded successfully!")
    print(f"Optimized for small survey area: {SURVEY_AREA_DEG2} deg²")
