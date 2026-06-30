"""
Functions for determining group membership using spectroscopic and photometric redshifts.
"""

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Import VRFEstimator for advanced spec-z membership
try:
    from VRFEstimator import VRFEstimator
    VRF_AVAILABLE = True
except ImportError:
    VRF_AVAILABLE = False
    warnings.warn("VRFEstimator not available. Will use gapper method only.")


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


def gapper_technique(z_values, v_escape_factor=2.5):
    """
    Apply gapper technique to identify velocity outliers.
    
    The gapper technique sorts redshifts and looks at gaps between consecutive values.
    Large gaps indicate potential interlopers.
    
    Parameters:
    -----------
    z_values : array
        Array of redshift values
    v_escape_factor : float
        Factor for velocity cut (default: 2.5 sigma)
        
    Returns:
    --------
    members : boolean array
        True for members, False for interlopers
    v_dispersion : float
        Velocity dispersion in km/s
    """
    
    if len(z_values) < 3:
        # Not enough galaxies for gapper technique, accept all
        return np.ones(len(z_values), dtype=bool), 0.0
    
    # Convert redshifts to velocities (cz in km/s)
    c = 299792.458  # speed of light in km/s
    v = z_values * c
    
    # Sort velocities
    v_sorted = np.sort(v)
    
    # Calculate gaps
    gaps = np.diff(v_sorted)
    
    # Gapper estimator for velocity dispersion
    # w_i are weights based on gaps
    n = len(v_sorted)
    w = np.arange(1, n) * (n - np.arange(1, n))
    # Beers et al. 1990, AJ 100:32, eq. 6 — note n*(n-1) in denominator
    v_dispersion = np.sqrt(np.pi) / (n * (n - 1)) * np.sum(w * gaps)
    
    # Identify members using velocity cut
    median_v = np.median(v)
    v_cut = v_escape_factor * v_dispersion
    
    members = np.abs(v - median_v) <= v_cut
    
    return members, v_dispersion


def find_specz_members_vrf(group_ra, group_dec, group_z, specz_catalog, 
                            radius_kpc=500, max_dz_norm=0.02, 
                            vrf_params=None, remove_duplicates=True):
    """
    Find spectroscopic redshift members using VRF interloper rejection algorithm.
    
    This uses the MBM10 (Mamon, Biviano & Murante 2010) iterative interloper 
    rejection with NFW profile and velocity dispersion estimation.
    
    Parameters:
    -----------
    group_ra : float
        Group RA in degrees
    group_dec : float
        Group Dec in degrees
    group_z : float
        Group redshift
    specz_catalog : DataFrame
        Spectroscopic redshift catalog with columns: RA, DEC, zfin, ez
    radius_kpc : float
        Physical radius for search in kpc (default: 500)
    max_dz_norm : float
        Maximum normalized redshift offset for initial cut (default: 0.02)
    vrf_params : dict, optional
        VRF algorithm parameters. If None, uses defaults:
        {
            'anis': 'ML',           # Anisotropy model (ML = Mamon-Lokas)
            'delta': 200,           # Overdensity definition
            'itmax': 20,            # Max iterations
            'evel': 150.0,          # Velocity error (km/s)
            'nscut': 2.7,           # Envelope width (sigma)
            'rvguess': 1000.0,      # Initial rvir guess (kpc)
            'H0': 70.0,             # Hubble constant
            'Om0': 0.3,             # Omega matter
            'Ode0': 0.7,            # Omega Lambda
            'widegap': 4.0,         # Weighted gap threshold
            'physgap': None,        # Physical gap (km/s), None to disable
            'nsigma': 2.7           # Sigma cut
        }
    remove_duplicates : bool
        Whether to remove duplicate coordinates (default: True)
        
    Returns:
    --------
    members : DataFrame
        Member galaxies with membership information, including:
        - All original galaxy properties
        - sep_arcsec, sep_kpc: Separations from group center
        - dv, dz_norm: Velocity offsets
        - vrf_member: VRF membership flag (1=member, 0=interloper)
        - v_dispersion: Group velocity dispersion
        - n_galaxies: Number of galaxies used in VRF
    """
    
    if not VRF_AVAILABLE:
        raise ImportError("VRFEstimator not available. Use find_specz_members() with use_gapper=True instead.")
    
    # Default VRF parameters
    default_vrf_params = {
        'anis': 'ML',
        'delta': 200,
        'itmax': 20,
        'evel': 150.0,
        'nscut': 2.7,
        'rvguess': 1000.0,
        'H0': 70.0,
        'Om0': 0.3,
        'Ode0': 0.7
    }
    
    if vrf_params is not None:
        default_vrf_params.update(vrf_params)
    vrf_params = default_vrf_params
    
    # STEP 1: Initial redshift cut
    dz_norm_all = np.abs(specz_catalog['zfin'] - group_z) / (1 + group_z)
    z_cut = dz_norm_all < max_dz_norm
    
    if z_cut.sum() == 0:
        return pd.DataFrame()
    
    specz_zcut = specz_catalog[z_cut].copy()
    
    # STEP 2: Spatial cut - find galaxies within radius
    radius_arcsec = angular_to_physical_radius(group_z, radius_kpc)
    
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    gal_coords = SkyCoord(ra=specz_zcut['RA'].values*u.deg, 
                          dec=specz_zcut['DEC'].values*u.deg)
    
    sep = group_coord.separation(gal_coords)
    within_radius = sep < radius_arcsec * u.arcsec
    
    if within_radius.sum() == 0:
        return pd.DataFrame()
    
    candidates = specz_zcut[within_radius].copy()
    candidates['sep_arcsec'] = sep[within_radius].arcsec
    
    # Calculate physical separation
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    candidates['sep_kpc'] = candidates['sep_arcsec'] * kpc_per_arcsec
    
    # Remove duplicates if requested
    if remove_duplicates and len(candidates) > 1:
        coords_str = candidates['RA'].astype(str) + '_' + candidates['DEC'].astype(str)
        duplicates = coords_str.duplicated(keep=False)
        
        if duplicates.any():
            unique_candidates = []
            for coord in coords_str[duplicates].unique():
                dup_mask = coords_str == coord
                best_idx = candidates.loc[dup_mask, 'sep_arcsec'].idxmin()
                unique_candidates.append(candidates.loc[best_idx])
            
            non_dup = candidates[~duplicates]
            candidates = pd.concat([non_dup] + unique_candidates, ignore_index=True)
    
    # STEP 3: Apply VRF interloper rejection
    if len(candidates) >= 5:  # Need at least 5 galaxies for VRF
        try:
            # Initialize VRF estimator
            vrf = VRFEstimator(vrf_params, None, None)
            
            # Prepare inputs for VRF
            rorig = candidates['sep_kpc'].values  # Projected radii in kpc
            zgal = candidates['zfin'].values      # Galaxy redshifts
            zcl = group_z                         # Group redshift
            
            # Run interloper rejection
            kwargs = {
                'widegap': vrf_params.get('widegap', 4.0),
                'nsigma': vrf_params.get('nsigma', 2.7),
                'noiter': False
            }
            
            # Add physgap if specified
            if 'physgap' in vrf_params and vrf_params['physgap'] is not None:
                kwargs['physgap'] = vrf_params['physgap']
            
            # Call VRF algorithm
            imembers = vrf.interlopmbmiter(rorig, zgal, zcl, **kwargs)
            
            # Mark members
            candidates['vrf_member'] = imembers.astype(bool)
            
            # Calculate velocity dispersion for members
            c = 299792.458  # km/s
            member_zs = candidates[candidates['vrf_member']]['zfin']
            if len(member_zs) > 1:
                v_disp = np.std((member_zs - group_z) * c)
            else:
                v_disp = 0.0
            
            candidates['v_dispersion'] = v_disp
            candidates['n_galaxies'] = len(member_zs)
            candidates['method'] = 'VRF'
            
        except Exception as e:
            # If VRF fails, fall back to all candidates and flag the issue
            warnings.warn(f"VRF failed for group at z={group_z:.3f}: {str(e)}. Marking all as members.")
            candidates['vrf_member'] = True
            candidates['v_dispersion'] = 0.0
            candidates['n_galaxies'] = len(candidates)
            candidates['method'] = 'VRF_failed'
    else:
        # Too few galaxies for VRF - accept all within radius
        candidates['vrf_member'] = True
        candidates['v_dispersion'] = 0.0
        candidates['n_galaxies'] = len(candidates)
        candidates['method'] = 'too_few'
    
    # Calculate velocity offsets
    c = 299792.458  # km/s
    candidates['dv'] = (candidates['zfin'] - group_z) * c
    candidates['dz_norm'] = (candidates['zfin'] - group_z) / (1 + group_z)
    
    return candidates


def find_specz_members(group_ra, group_dec, group_z, specz_catalog, 
                       radius_kpc=500, max_dz_norm=0.02, max_velocity=3000,
                       use_gapper=True, remove_duplicates=True, method='gapper',
                       vrf_params=None):
    """
    Find spectroscopic redshift members for a group.
    
    Parameters:
    -----------
    group_ra : float
        Group RA in degrees
    group_dec : float
        Group Dec in degrees
    group_z : float
        Group redshift
    specz_catalog : DataFrame
        Spectroscopic redshift catalog with columns: RA, DEC, zfin, ez
    radius_kpc : float
        Physical radius for search in kpc (default: 500)
    max_dz_norm : float
        Maximum normalized redshift offset |Δz/(1+z)| for initial cut (default: 0.02)
    max_velocity : float
        Maximum velocity offset in km/s for membership (default: 3000, used only for gapper method)
    use_gapper : bool
        Whether to use gapper technique (default: True, used only for gapper method)
    remove_duplicates : bool
        Whether to remove duplicate coordinates (default: True)
    method : str
        Method for membership determination: 'gapper' or 'vrf' (default: 'gapper')
        - 'gapper': Simple gapper technique with velocity cuts
        - 'vrf': VRF interloper rejection (MBM10 algorithm)
    vrf_params : dict, optional
        Parameters for VRF algorithm (used only if method='vrf')
        
    Returns:
    --------
    members : DataFrame
        Member galaxies with membership information
    """
    
    # Route to appropriate method
    if method == 'vrf':
        if not VRF_AVAILABLE:
            warnings.warn("VRF method requested but VRFEstimator not available. Falling back to gapper.")
            method = 'gapper'
        else:
            return find_specz_members_vrf(group_ra, group_dec, group_z, specz_catalog,
                                         radius_kpc=radius_kpc, max_dz_norm=max_dz_norm,
                                         vrf_params=vrf_params, remove_duplicates=remove_duplicates)
    
    # Continue with gapper method (original implementation)
    # STEP 1: Apply initial redshift cut to avoid including z=0 or z>>group_z galaxies
    # Use Δz/(1+z) < max_dz_norm as initial filter
    dz_norm_all = np.abs(specz_catalog['zfin'] - group_z) / (1 + group_z)
    z_cut = dz_norm_all < max_dz_norm
    
    if z_cut.sum() == 0:
        # No galaxies in redshift range
        return pd.DataFrame()
    
    # Apply redshift cut first
    specz_zcut = specz_catalog[z_cut].copy()
    
    # STEP 2: Spatial cut - find galaxies within radius
    # Convert physical radius to angular radius at group redshift
    radius_arcsec = angular_to_physical_radius(group_z, radius_kpc)
    
    # Create SkyCoord for group
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    
    # Create SkyCoord for galaxies
    gal_coords = SkyCoord(ra=specz_zcut['RA'].values*u.deg, 
                          dec=specz_zcut['DEC'].values*u.deg)
    
    # Find galaxies within radius
    sep = group_coord.separation(gal_coords)
    within_radius = sep < radius_arcsec * u.arcsec
    
    if within_radius.sum() == 0:
        # No galaxies found
        return pd.DataFrame()
    
    # Extract candidates
    candidates = specz_zcut[within_radius].copy()
    candidates['sep_arcsec'] = sep[within_radius].arcsec
    
    # Calculate physical separation in kpc
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    candidates['sep_kpc'] = candidates['sep_arcsec'] * kpc_per_arcsec
    
    # Remove duplicates if requested
    if remove_duplicates and len(candidates) > 1:
        # Keep the one closest to group center for each duplicate coordinate
        coords_str = candidates['RA'].astype(str) + '_' + candidates['DEC'].astype(str)
        duplicates = coords_str.duplicated(keep=False)
        
        if duplicates.any():
            # For duplicates, keep the one with smallest separation
            unique_candidates = []
            for coord in coords_str[duplicates].unique():
                dup_mask = coords_str == coord
                best_idx = candidates.loc[dup_mask, 'sep_arcsec'].idxmin()
                unique_candidates.append(candidates.loc[best_idx])
            
            # Add non-duplicates
            non_dup = candidates[~duplicates]
            candidates = pd.concat([non_dup] + unique_candidates, ignore_index=True)
    
    # STEP 3: Apply velocity/gapper technique if requested and enough galaxies
    if use_gapper and len(candidates) >= 5:
        # Use gapper technique but with maximum velocity limit
        member_mask, v_disp = gapper_technique(candidates['zfin'].values, v_escape_factor=2.5)
        
        # Additional check: velocity dispersion should be reasonable (<= 1500 km/s for groups)
        # If too high, it means we're including field galaxies - use stricter cut
        if v_disp > 1500:
            # Fall back to simple velocity cut with tighter limit
            c = 299792.458  # km/s
            dv_abs = np.abs(candidates['zfin'] - group_z) * c
            member_mask = dv_abs < min(max_velocity, 1500)  # Use 1500 km/s for high-z groups
            v_disp = np.std((candidates.loc[member_mask, 'zfin'] - group_z) * c) if member_mask.sum() > 1 else 0.0
            
        # Additional sanity check: if still too high, use median-based method
        if v_disp > 1200:
            c = 299792.458  # km/s
            z_median = np.median(candidates['zfin'])
            dv_from_median = np.abs(candidates['zfin'] - z_median) * c
            # Use robust scatter estimate
            mad_v = 1.48 * np.median(dv_from_median)  # MAD-based velocity dispersion
            if mad_v < 1200:
                # Accept galaxies within 3*MAD of median
                member_mask = dv_from_median < 3 * mad_v
                v_disp = mad_v
            else:
                # Very scattered - use tight cut around group redshift
                dv_abs = np.abs(candidates['zfin'] - group_z) * c
                member_mask = dv_abs < 1000  # Conservative 1000 km/s
                v_disp = np.std((candidates.loc[member_mask, 'zfin'] - group_z) * c) if member_mask.sum() > 1 else 0.0
        
        candidates['gapper_member'] = member_mask
        candidates['v_dispersion'] = v_disp
    else:
        # Simple velocity cut for few galaxies
        c = 299792.458  # km/s
        dv_abs = np.abs(candidates['zfin'] - group_z) * c
        candidates['gapper_member'] = dv_abs < max_velocity
        candidates['v_dispersion'] = np.std((candidates['zfin'] - group_z) * c) if len(candidates) > 1 else 0.0
    
    # Calculate velocity offset from group
    c = 299792.458  # km/s
    candidates['dv'] = (candidates['zfin'] - group_z) * c
    candidates['dz_norm'] = (candidates['zfin'] - group_z) / (1 + group_z)
    
    return candidates


def calculate_photoz_membership_probability(galaxy_ra, galaxy_dec, galaxy_z, galaxy_z_err,
                                            group_ra, group_dec, group_z,
                                            radius_kpc=500, nmad_redshift_dependent=None,
                                            stellar_mass=None, use_mass_weighting=False):
    """
    Calculate membership probability for a galaxy with photometric redshift.
    
    Uses a Bayesian approach considering:
    1. Redshift probability (based on photo-z error and accuracy)
    2. Spatial probability (based on offset from group center)
    3. Optionally: stellar mass weighting
    
    Parameters:
    -----------
    galaxy_ra, galaxy_dec : float or array
        Galaxy coordinates in degrees
    galaxy_z : float or array
        Galaxy photo-z
    galaxy_z_err : float or array
        Galaxy photo-z error (or NMAD if error not available)
    group_ra, group_dec : float
        Group coordinates in degrees
    group_z : float
        Group redshift
    radius_kpc : float
        Physical radius for membership consideration (default: 500 kpc)
    nmad_redshift_dependent : function or None
        Function that returns NMAD as function of redshift
    stellar_mass : float or array or None
        Stellar mass in solar masses (for optional weighting)
    use_mass_weighting : bool
        Whether to use stellar mass weighting (default: False)
        
    Returns:
    --------
    probability : float or array
        Membership probability (0 to 1)
    """
    
    # Convert to arrays if scalars
    scalar_input = np.isscalar(galaxy_ra)
    if scalar_input:
        galaxy_ra = np.array([galaxy_ra])
        galaxy_dec = np.array([galaxy_dec])
        galaxy_z = np.array([galaxy_z])
        galaxy_z_err = np.array([galaxy_z_err])
        if stellar_mass is not None:
            stellar_mass = np.array([stellar_mass])
    
    # 1. Redshift probability
    # P(z_gal | z_group) using Gaussian based on photo-z accuracy
    dz = galaxy_z - group_z
    dz_norm = dz / (1 + group_z)
    
    # Use σ_MAD from COSMOS2025 paper if available
    if nmad_redshift_dependent is not None:
        sigma_z = nmad_redshift_dependent(group_z)
    else:
        sigma_z = galaxy_z_err / (1 + group_z)
    
    # Gaussian probability: P(z) ~ exp(-0.5 * (Δz/σ)^2)
    # This gives P=1 at Δz=0, P=0.6 at Δz=σ, P=0.14 at Δz=2σ, P=0.01 at Δz=3σ
    z_prob = np.exp(-0.5 * (dz_norm / sigma_z)**2)
    
    # 2. Spatial probability
    # Use a flatter radial profile - most cluster studies use uniform weighting
    # or very gentle decline within the aperture
    # Reference: Rykoff et al. (2014) redMaPPer, Koester et al. (2007)
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    gal_coords = SkyCoord(ra=galaxy_ra*u.deg, dec=galaxy_dec*u.deg)
    sep_arcsec = group_coord.separation(gal_coords).arcsec
    
    # Convert to physical separation
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    sep_kpc = sep_arcsec * kpc_per_arcsec
    
    # Gentle radial weighting: P(r) = 1 - (r/R_max)^2
    # This keeps P > 0.5 for r < 0.7*R_max, and P > 0.25 for r < 0.87*R_max
    # Much less harsh than exponential decay
    spatial_prob = 1.0 - (sep_kpc / radius_kpc)**2
    spatial_prob = np.maximum(spatial_prob, 0.0)  # Clip to [0, 1]
    
    # Alternative: even simpler - constant weight within radius
    # spatial_prob = np.ones_like(sep_kpc)
    # spatial_prob[sep_kpc > radius_kpc] = 0.0
    
    # 3. Combined probability
    # The redshift probability dominates, spatial is gentle modifier
    combined_prob = z_prob * spatial_prob
    
    # 4. Optional: stellar mass weighting
    # More massive galaxies are more likely to be members
    if use_mass_weighting and stellar_mass is not None:
        # Normalize mass weight: typical mass ~ 10^10 M_sun
        mass_weight = np.log10(stellar_mass / 1e10 + 1) + 1
        mass_weight = np.clip(mass_weight, 0.5, 2.0)  # Keep reasonable range
        combined_prob *= mass_weight
    
    # Normalize to [0, 1]
    combined_prob = np.clip(combined_prob, 0, 1)
    
    if scalar_input:
        return combined_prob[0]
    return combined_prob


def find_photoz_members(group_ra, group_dec, group_z, photoz_catalog,
                       radius_kpc=500, nmad_function=None, max_dist_phot_z_error=3.0,
                       prob_threshold=0.05, remove_duplicates=True):
    """
    Find photometric redshift members for a group with membership probabilities.
    
    Uses cylinder selection with height based on photo-z accuracy from COSMOS2025 paper:
    z_range = z_group ± max_dist_phot_z_error × σ_MAD × (1 + z_group)
    
    where σ_MAD is the normalized median absolute deviation from the paper:
    σ_MAD = 1.48 × median(|Δz/(1+z)|)
    
    For typical groups (z ~ 0.5-1.5): σ_MAD ≈ 0.012-0.015
    For 3σ cut: Δz ≈ 0.04-0.06, capturing 99.7% of true members
    
    Parameters:
    -----------
    group_ra : float
        Group RA in degrees
    group_dec : float
        Group Dec in degrees
    group_z : float
        Group redshift
    photoz_catalog : DataFrame
        Photometric redshift catalog with columns: RA_MODEL, DEC_MODEL, LP_zfinal, LP_warn_fl
    radius_kpc : float
        Physical radius for search in kpc (default: 500)
    nmad_function : function or None
        Function that returns σ_MAD as function of redshift (from COSMOS2025 paper)
    max_dist_phot_z_error : float
        Cylinder height in units of σ_MAD × (1 + z) (default: 3.0 for 3σ cut)
        - 3.0: captures 99.7% of true members (recommended)
        - 2.5: more conservative, 98.8% completeness
        - 5.0: very loose, includes more contaminants
    prob_threshold : float
        Minimum probability to include (default: 0.05)
    remove_duplicates : bool
        Whether to remove duplicate coordinates (default: True)
        
    Returns:
    --------
    members : DataFrame
        Member galaxies with membership probabilities
    """
    
    # Clean photo-z catalog - use high-quality photo-z only (LP_warn_fl==0)
    # This gives 694,074 clean galaxies from galaxy_catalog_photz.csv
    clean_photoz = photoz_catalog[
        (photoz_catalog['LP_warn_fl'] == 0)
    ].copy()
    
    if len(photoz_catalog) > 0:
        print(f"Photo-z catalog filtering: {len(photoz_catalog)} → {len(clean_photoz)} galaxies (LP_warn_fl==0)")
    
    # Get NMAD for this redshift
    if nmad_function is not None:
        nmad = nmad_function(group_z)
    else:
        # Default: σ_MAD from COSMOS2025 paper (Shuntov et al. 2025)
        # Paper reports σ_MAD = 1.48 × median(|Δz/(1+z)|)
        if group_z < 1.0:
            nmad = 0.012  # Excellent performance for z < 1
        elif group_z < 2.0:
            nmad = 0.015  # Good performance for 1 < z < 2
        elif group_z < 3.0:
            nmad = 0.020  # Moderate performance for 2 < z < 3
        else:
            nmad = 0.030  # Conservative for z > 3
    
    # Define cylinder height based on photo-z accuracy
    z_error = nmad * (1 + group_z)
    z_min = group_z - max_dist_phot_z_error * z_error
    z_max = group_z + max_dist_phot_z_error * z_error
    
    # Apply redshift cut
    z_mask = (clean_photoz['LP_zfinal'] >= z_min) & (clean_photoz['LP_zfinal'] <= z_max)
    
    if z_mask.sum() == 0:
        return pd.DataFrame()
    
    clean_photoz = clean_photoz[z_mask].copy()
    
    # Convert physical radius to angular radius at group redshift
    radius_arcsec = angular_to_physical_radius(group_z, radius_kpc)
    
    # Create SkyCoord for group
    group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
    
    # Create SkyCoord for galaxies
    gal_coords = SkyCoord(ra=clean_photoz['RA_MODEL'].values*u.deg, 
                          dec=clean_photoz['DEC_MODEL'].values*u.deg)
    
    # Find galaxies within radius
    sep = group_coord.separation(gal_coords)
    within_radius = sep < radius_arcsec * u.arcsec
    
    if within_radius.sum() == 0:
        # No galaxies found
        return pd.DataFrame()
    
    # Extract candidates
    candidates = clean_photoz[within_radius].copy()
    candidates['sep_arcsec'] = sep[within_radius].arcsec
    
    # Calculate physical separation in kpc
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(group_z).value / 60.0
    candidates['sep_kpc'] = candidates['sep_arcsec'] * kpc_per_arcsec
    
    # Remove duplicates if requested
    if remove_duplicates and len(candidates) > 1:
        coords_str = candidates['RA_MODEL'].astype(str) + '_' + candidates['DEC_MODEL'].astype(str)
        candidates = candidates[~coords_str.duplicated(keep='first')]
    
    # Calculate membership probabilities
    candidates['membership_prob'] = calculate_photoz_membership_probability(
        candidates['RA_MODEL'].values,
        candidates['DEC_MODEL'].values,
        candidates['LP_zfinal'].values,
        nmad,  # Use NMAD as error estimate
        group_ra, group_dec, group_z,
        radius_kpc=radius_kpc,
        nmad_redshift_dependent=nmad_function
    )
    
    # Filter by probability threshold
    candidates = candidates[candidates['membership_prob'] >= prob_threshold]
    
    # Calculate redshift offset
    candidates['dz'] = candidates['LP_zfinal'] - group_z
    candidates['dz_norm'] = candidates['dz'] / (1 + group_z)
    
    # Add cylinder limits for reference
    candidates['z_min_cylinder'] = z_min
    candidates['z_max_cylinder'] = z_max
    candidates['nmad_used'] = nmad
    
    return candidates


def create_nmad_function(matched_catalog=None):
    """
    Create a function that returns NMAD as function of redshift.
    
    Parameters:
    -----------
    matched_catalog : DataFrame or None
        Matched spec-z/photo-z catalog with columns: zspec, zphot
        If None, uses default values from exploration
        
    Returns:
    --------
    nmad_func : function
        Function that takes redshift and returns NMAD
    """
    
    if matched_catalog is not None:
        # Calculate NMAD in bins using clean sample (LP_warn_fl==0)
        z_bins = np.arange(0, 4.5, 0.2)
        z_centers = []
        nmad_values = []
        
        # Use the actual column names from the matched catalog
        # Check which column names are present
        if 'zphot' in matched_catalog.columns and 'zspec' in matched_catalog.columns:
            zphot_col, zspec_col = 'zphot', 'zspec'
        elif 'LP_zfinal' in matched_catalog.columns and 'zfin' in matched_catalog.columns:
            zphot_col, zspec_col = 'LP_zfinal', 'zfin'
        else:
            print("Warning: Could not find redshift columns, using default NMAD")
            matched_catalog = None
            
        if matched_catalog is not None:
            # Filter for clean sample: LP_warn_fl==0 for reliable NMAD calculation
            clean_catalog = matched_catalog.copy()
            if 'LP_warn_fl' in clean_catalog.columns:
                clean_mask = (clean_catalog['LP_warn_fl'] == 0) & \
                            (clean_catalog[zphot_col].notna()) & \
                            (clean_catalog[zspec_col].notna())
                clean_catalog = clean_catalog[clean_mask]
                print(f"NMAD calculation: Using {len(clean_catalog)}/{len(matched_catalog)} clean galaxies (LP_warn_fl==0)")
            
            if len(clean_catalog) < 100:
                print(f"Warning: Only {len(clean_catalog)} clean galaxies for NMAD, using default")
                matched_catalog = None
            else:
                # Calculate normalized redshift difference
                dz = clean_catalog[zphot_col] - clean_catalog[zspec_col]
                dz_norm = dz / (1 + clean_catalog[zspec_col])
            
                for i in range(len(z_bins)-1):
                    z_min, z_max = z_bins[i], z_bins[i+1]
                    bin_mask = (clean_catalog[zspec_col] >= z_min) & (clean_catalog[zspec_col] < z_max)
                
                    if bin_mask.sum() > 10:
                        dz_bin = dz_norm[bin_mask]
                        nmad = 1.48 * np.median(np.abs(dz_bin - np.median(dz_bin)))
                        z_centers.append((z_min + z_max) / 2)
                        nmad_values.append(nmad)
            
                # Interpolate
                if len(z_centers) > 1:
                    from scipy.interpolate import interp1d
                    nmad_func = interp1d(z_centers, nmad_values, 
                                        kind='linear', fill_value='extrapolate')
                    print(f"NMAD function created with {len(z_centers)} redshift bins")
                else:
                    print("Warning: Not enough redshift bins for NMAD interpolation, using default")
                    matched_catalog = None
    
    if matched_catalog is None:
        # Default function from COSMOS2025 paper (Shuntov et al. 2025, arXiv:2506.03243)
        def nmad_func(z, mag_f444w=None):
            """
            σ_MAD from COSMOS2025 paper.
            
            Paper reports σ_MAD = 1.48 × median(|Δz/(1+z)|):
            - Overall: σ_MAD = 0.012 at mF444W < 28
            - Magnitude-dependent: 0.011 (bright) to 0.030 (faint)
            - Color-dependent: 0.010 (blue) to 0.035 (red/dusty)
            
            This is the NORMALIZED MAD, not the raw photo-z error.
            """
            # Magnitude-dependent (if F444W magnitude available)
            if mag_f444w is not None:
                if mag_f444w < 23:
                    return 0.011
                elif mag_f444w < 24:
                    return 0.014
                elif mag_f444w < 25:
                    return 0.015
                elif mag_f444w < 26:
                    return 0.020
                else:  # mag >= 26
                    return 0.030
            
            # Redshift-dependent (default, conservative)
            if z < 1.0:
                return 0.012  # Excellent performance for z < 1
            elif z < 2.0:
                return 0.015  # Good performance for 1 < z < 2
            elif z < 3.0:
                return 0.020  # Moderate performance for 2 < z < 3
            else:
                return 0.030  # Conservative for z > 3
    
    return nmad_func


def combine_specz_photoz_members(specz_members, photoz_members, match_radius_arcsec=1.0):
    """
    Combine spec-z and photo-z member catalogs, avoiding duplicates.
    
    Gives priority to spec-z when galaxies are matched.
    
    Parameters:
    -----------
    specz_members : DataFrame
        Spec-z members
    photoz_members : DataFrame
        Photo-z members
    match_radius_arcsec : float
        Matching radius to identify same galaxy (default: 1 arcsec)
        
    Returns:
    --------
    combined : DataFrame
        Combined catalog with source flag
    """
    
    if len(specz_members) == 0:
        photoz_members['source'] = 'photoz'
        return photoz_members
    
    if len(photoz_members) == 0:
        specz_members['source'] = 'specz'
        return specz_members
    
    # Remove any NaN coordinates before matching
    specz_valid = specz_members.dropna(subset=['RA', 'DEC'])
    photoz_valid = photoz_members.dropna(subset=['RA_MODEL', 'DEC_MODEL'])
    
    if len(specz_valid) == 0:
        photoz_members['source'] = 'photoz'
        return photoz_members
    
    if len(photoz_valid) == 0:
        specz_members['source'] = 'specz'
        return specz_members
    
    # Match catalogs
    specz_coords = SkyCoord(ra=specz_valid['RA'].values*u.deg,
                           dec=specz_valid['DEC'].values*u.deg)
    photoz_coords = SkyCoord(ra=photoz_valid['RA_MODEL'].values*u.deg,
                            dec=photoz_valid['DEC_MODEL'].values*u.deg)
    
    idx, sep2d, _ = photoz_coords.match_to_catalog_sky(specz_coords)
    matched = sep2d < match_radius_arcsec * u.arcsec
    
    # Keep unmatched photo-z members (only from the valid ones we matched)
    unmatched_photoz = photoz_valid[~matched].copy()
    
    # Also add back any photo-z members that were dropped due to NaN coordinates
    # (they can't match spec-z anyway, so they should be included)
    photoz_invalid_coords = photoz_members[
        photoz_members['RA_MODEL'].isna() | photoz_members['DEC_MODEL'].isna()
    ].copy()
    
    # Combine all photo-z members (unmatched valid + invalid coords)
    all_photoz = pd.concat([unmatched_photoz, photoz_invalid_coords], ignore_index=True)
    all_photoz['source'] = 'photoz'
    
    # Add spec-z members (use original, not filtered)
    specz_members_copy = specz_members.copy()
    specz_members_copy['source'] = 'specz'
    
    # Combine
    combined = pd.concat([specz_members_copy, all_photoz], ignore_index=True)
    
    return combined
