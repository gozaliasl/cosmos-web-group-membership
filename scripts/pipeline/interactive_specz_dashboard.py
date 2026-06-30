#!/usr/bin/env python
"""
Interactive Spectroscopic Membership Dashboard

Real-time visualization and parameter adjustment for spec-z membership determination
using the VRF method. Allows users to inspect membership, adjust parameters, and
navigate through groups interactively.

Features:
- 4-panel visualization (spatial, velocity-radius, cumulative, VRF history)
- Interactive parameter adjustment (sliders)
- Real-time VRF recomputation
- Member/non-member highlighting
- Physical property display (σ_v, M_200, R_200)
- Save and navigate controls

Usage:
    python interactive_specz_dashboard.py --catalog cw-all --group-id 1
    python interactive_specz_dashboard.py --catalog cw-hcg --group-id Py12_65 --radius 500

Requirements:
    conda activate astro-clean
    matplotlib backend with interactive support (Qt5Agg, TkAgg)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.patches import Circle, Ellipse
from scipy.stats import norm
import argparse
from pathlib import Path
import sys
import warnings

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from membership_functions import find_specz_members
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord
from astropy import units as u

warnings.filterwarnings('ignore')

# Use interactive backend (try TkAgg first, fallback to default)
try:
    import matplotlib
    matplotlib.use('TkAgg')
except:
    pass  # Use default backend

# Paths
BASE_DIR = Path('/Users/gozalig1/Projects/compact-groups-xray-analysis')
GROUP_CATALOG_DIR = BASE_DIR / 'data' / 'group-catalog'
# Final HCG group catalog (Py18_Groups.csv or .fits)
HCG_CATALOG_CSV = GROUP_CATALOG_DIR / 'Py18_Groups.csv'
HCG_CATALOG_FITS = GROUP_CATALOG_DIR / 'Py18_Groups.fits'
SPECZ_CATALOG = BASE_DIR / 'data' / 'specz' / 'Webb_Specz_with_photz.csv'
FULL_GALAXY_CATALOG = BASE_DIR / 'data' / 'galaxy_catalog_photz' / 'galaxy_catalog_photz.csv'
OUTPUT_DIR = BASE_DIR / 'membership_determination' / 'results' / 'interactive_review'

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_hcg_catalog():
    """Load final HCG group catalog (Py18_Groups.csv or Py18_Groups.fits)."""
    if HCG_CATALOG_CSV.exists():
        return pd.read_csv(HCG_CATALOG_CSV)
    if HCG_CATALOG_FITS.exists():
        from astropy.table import Table
        t = Table.read(HCG_CATALOG_FITS)
        return t.to_pandas()
    raise FileNotFoundError(f"HCG catalog not found. Expected {HCG_CATALOG_CSV} or {HCG_CATALOG_FITS}")


def gapper_estimator(velocities):
    """
    Calculate velocity dispersion using the Gapper estimator.
    
    The Gapper estimator is robust and doesn't assume a Gaussian distribution.
    It uses gaps between sorted velocities weighted by position.
    
    Parameters:
    -----------
    velocities : array-like
        Velocity values (in km/s)
    
    Returns:
    --------
    sigma_gapper : float
        Velocity dispersion estimate
    
    Reference:
    ----------
    Wainer & Thissen (1976), Beers et al. (1990)
    """
    v = np.sort(velocities)
    n = len(v)
    
    if n < 2:
        return np.nan
    
    # Calculate gaps (differences between consecutive sorted velocities)
    gaps = np.diff(v)
    
    # Weight by position: w_i = i * (n - i)
    weights = np.arange(1, n) * np.arange(n-1, 0, -1)
    
    # Gapper estimator: sqrt(π) * Σ(w_i * g_i) / (n * (n-1))
    sigma_gapper = np.sqrt(np.pi) * np.sum(weights * gaps) / (n * (n - 1))
    
    return sigma_gapper


class InteractiveSpeczDashboard:
    """Interactive dashboard for spec-z membership review and adjustment."""
    
    def __init__(self, catalog_name, start_group_id, radius_kpc=500):
        """
        Initialize dashboard.
        
        Parameters:
        -----------
        catalog_name : str
            'cw-all' or 'cw-hcg'
        start_group_id : int or str
            Starting group ID
        radius_kpc : float
            Search radius in kpc
        """
        self.catalog_name = catalog_name
        self.current_group_idx = 0
        self.radius_kpc = radius_kpc
        
        # Load catalogs
        print("Loading catalogs...")
        if catalog_name == 'cw-all':
            self.groups = pd.read_csv(GROUP_CATALOG_DIR / 'cosmos_web_groups_catalog.csv')
            self.id_col = 'Group_ID'
            self.ra_col, self.dec_col, self.z_col = 'Ra', 'Dec', 'z'
            
            # Find starting index
            start_idx = self.groups[self.groups[self.id_col] == start_group_id].index
            if len(start_idx) > 0:
                self.current_group_idx = start_idx[0]
        else:
            self.groups = load_hcg_catalog()
            self.id_col = 'Grp'
            self.ra_col, self.dec_col, self.z_col = 'Ra', 'Dec', 'z'
            
            # Find starting index
            start_idx = self.groups[self.groups[self.id_col] == start_group_id].index
            if len(start_idx) > 0:
                self.current_group_idx = start_idx[0]
        
        self.specz = pd.read_csv(SPECZ_CATALOG)
        
        # Load full galaxy catalog (spec-z + photo-z) for BGG identification
        print(f"  Loading full galaxy catalog...")
        self.full_catalog = pd.read_csv(FULL_GALAXY_CATALOG)
        
        print(f"  Loaded {len(self.groups)} groups")
        print(f"  Loaded {len(self.specz)} spec-z galaxies")
        print(f"  Loaded {len(self.full_catalog)} total galaxies (for BGG search)")
        print(f"  Starting at group index {self.current_group_idx}")
        
        # Default parameters
        self.max_dz_norm = 0.01
        self.max_velocity = 2000
        self.group_z_original = None  # Will be set on first plot
        self.group_z_adjusted = None  # Adjusted by user
        
        # VRF results storage
        self.vrf_M200 = None
        self.vrf_Rvir = None
        self.vrf_sigma_v = None
        self.vrf_z_median = None
        
        # Results storage
        self.reviewed_groups = []
        
        # Create figure and axes
        self.setup_figure()
        
        # Initial plot
        self.update_plot()
        
    def setup_figure(self):
        """Create figure layout with 4 panels and controls."""
        self.fig = plt.figure(figsize=(18, 12))
        
        # Create grid: 2x2 for plots, bottom for controls
        gs = self.fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.3], hspace=0.3, wspace=0.3)
        
        # Plot axes
        self.ax_spatial = self.fig.add_subplot(gs[0, 0])
        self.ax_velocity = self.fig.add_subplot(gs[0, 1])
        self.ax_cumulative = self.fig.add_subplot(gs[1, 0])
        self.ax_nfw = self.fig.add_subplot(gs[1, 1])
        
        # Control axes
        self.ax_controls = self.fig.add_subplot(gs[2, :])
        self.ax_controls.axis('off')
        
        # Sliders with input boxes
        slider_x_start = 0.12
        slider_width = 0.18
        slider_height = 0.025
        textbox_width = 0.04
        textbox_offset = 0.005
        
        # Group redshift adjustment slider + textbox
        ax_z_slider = self.fig.add_axes([slider_x_start, 0.20, slider_width, slider_height])
        self.slider_z = Slider(
            ax_z_slider, 'Group z',
            0.0, 2.0, valinit=0.71,
            valstep=0.0001, valfmt='%.4f'  # Finer step: 0.0001
        )
        self.slider_z.on_changed(self.on_slider_change)
        
        ax_z_textbox = self.fig.add_axes([slider_x_start + slider_width + textbox_offset, 0.20, textbox_width, slider_height])
        self.textbox_z = TextBox(ax_z_textbox, '', initial=f'{0.71:.4f}')
        self.textbox_z.on_submit(self.on_z_textbox_submit)
        
        # Radius slider (kpc) + textbox
        ax_radius_slider = self.fig.add_axes([slider_x_start, 0.16, slider_width, slider_height])
        self.slider_radius = Slider(
            ax_radius_slider, 'Radius (kpc)',
            200, 1000, valinit=self.radius_kpc,
            valstep=10  # Finer step: 10 kpc
        )
        self.slider_radius.on_changed(self.on_slider_change)
        
        ax_radius_textbox = self.fig.add_axes([slider_x_start + slider_width + textbox_offset, 0.16, textbox_width, slider_height])
        self.textbox_radius = TextBox(ax_radius_textbox, '', initial=f'{self.radius_kpc:.0f}')
        self.textbox_radius.on_submit(self.on_radius_textbox_submit)
        
        # max_dz_norm slider + textbox
        ax_dz_slider = self.fig.add_axes([slider_x_start, 0.12, slider_width, slider_height])
        self.slider_dz = Slider(
            ax_dz_slider, 'max_dz_norm',
            0.005, 0.03, valinit=self.max_dz_norm,
            valstep=0.0001  # Finer step: 0.0001
        )
        self.slider_dz.on_changed(self.on_slider_change)
        
        ax_dz_textbox = self.fig.add_axes([slider_x_start + slider_width + textbox_offset, 0.12, textbox_width, slider_height])
        self.textbox_dz = TextBox(ax_dz_textbox, '', initial=f'{self.max_dz_norm:.4f}')
        self.textbox_dz.on_submit(self.on_dz_textbox_submit)
        
        # max_velocity slider + textbox
        ax_vel_slider = self.fig.add_axes([slider_x_start, 0.08, slider_width, slider_height])
        self.slider_vel = Slider(
            ax_vel_slider, 'max_velocity (km/s)',
            1000, 4000, valinit=self.max_velocity,
            valstep=50  # Finer step: 50 km/s
        )
        self.slider_vel.on_changed(self.on_slider_change)
        
        ax_vel_textbox = self.fig.add_axes([slider_x_start + slider_width + textbox_offset, 0.08, textbox_width, slider_height])
        self.textbox_vel = TextBox(ax_vel_textbox, '', initial=f'{self.max_velocity:.0f}')
        self.textbox_vel.on_submit(self.on_vel_textbox_submit)
        
        # Buttons
        button_width = 0.10
        button_height = 0.035
        button_y = 0.02
        button_x_start = 0.45
        
        # Rerun button
        ax_rerun = self.fig.add_axes([button_x_start, button_y, button_width, button_height])
        self.btn_rerun = Button(ax_rerun, 'Rerun VRF', color='lightblue', hovercolor='skyblue')
        self.btn_rerun.on_clicked(self.on_rerun_click)
        
        # Accept button
        ax_accept = self.fig.add_axes([button_x_start + 0.12, button_y, button_width, button_height])
        self.btn_accept = Button(ax_accept, 'Accept & Next', color='lightgreen', hovercolor='lightgreen')
        self.btn_accept.on_clicked(self.on_accept_click)
        
        # Previous button
        ax_prev = self.fig.add_axes([button_x_start + 0.24, button_y, button_width, button_height])
        self.btn_prev = Button(ax_prev, 'Previous', color='lightyellow', hovercolor='yellow')
        self.btn_prev.on_clicked(self.on_prev_click)
        
        # Save button
        ax_save = self.fig.add_axes([button_x_start + 0.36, button_y, button_width, button_height])
        self.btn_save = Button(ax_save, 'Save Results', color='lightcoral', hovercolor='coral')
        self.btn_save.on_clicked(self.on_save_click)
        
        # Auto-Optimize button (above other buttons)
        ax_optimize = self.fig.add_axes([button_x_start, button_y + 0.045, button_width * 1.2, button_height])
        self.btn_optimize = Button(ax_optimize, 'Auto-Optimize', color='orange', hovercolor='darkorange')
        self.btn_optimize.on_clicked(self.on_optimize_click)
        
        # Recenter on BGG button
        ax_recenter = self.fig.add_axes([button_x_start + 0.14, button_y + 0.045, button_width * 1.2, button_height])
        self.btn_recenter = Button(ax_recenter, 'Center on BGG', color='gold', hovercolor='yellow')
        self.btn_recenter.on_clicked(self.on_recenter_bgg_click)
        
        # Jump to group textbox
        ax_jump = self.fig.add_axes([0.12, button_y, 0.08, button_height])
        self.textbox_jump = TextBox(ax_jump, 'Jump to:', initial=str(self.current_group_idx))
        self.textbox_jump.on_submit(self.on_jump_submit)
        
        self.fig.suptitle('Spectroscopic Membership Interactive Dashboard', 
                         fontsize=16, fontweight='bold')
        
    def get_current_group(self):
        """Get current group information."""
        row = self.groups.iloc[self.current_group_idx]
        
        return {
            'id': row[self.id_col],
            'ra': row[self.ra_col],
            'dec': row[self.dec_col],
            'z': row[self.z_col],
            'idx': self.current_group_idx
        }
    
    def update_plot(self):
        """Update all plots for current group."""
        group = self.get_current_group()
        
        # Initialize/update group redshift
        if self.group_z_original is None or self.group_z_original != group['z']:
            self.group_z_original = group['z']
            self.group_z_adjusted = group['z']
            # Update slider range and value
            self.slider_z.valmin = max(0.0, group['z'] - 0.05)
            self.slider_z.valmax = group['z'] + 0.05
            self.slider_z.set_val(group['z'])
        
        # Use adjusted redshift
        working_z = self.group_z_adjusted
        
        print(f"\n{'='*60}")
        print(f"Group {group['id']} (index {group['idx']}/{len(self.groups)-1})")
        print(f"z_catalog = {group['z']:.4f}, z_working = {working_z:.4f}, RA = {group['ra']:.4f}, Dec = {group['dec']:.4f}")
        print(f"Parameters: radius={self.radius_kpc:.0f} kpc, max_dz_norm={self.max_dz_norm:.3f}, max_velocity={self.max_velocity:.0f} km/s")
        
        # Find spec-z members using adjusted redshift
        self.current_members = find_specz_members(
            group['ra'], group['dec'], working_z, self.specz,
            radius_kpc=self.radius_kpc,
            max_dz_norm=self.max_dz_norm,
            max_velocity=self.max_velocity,
            use_gapper=True,
            remove_duplicates=True,
            method='vrf'
        )
        
        if len(self.current_members) == 0:
            print("  No spec-z members found!")
            self.clear_all_axes()
            return
        
        # Separate members and non-members
        member_col = 'vrf_member' if 'vrf_member' in self.current_members.columns else 'gapper_member'
        members = self.current_members[self.current_members[member_col] == True]
        non_members = self.current_members[self.current_members[member_col] == False]
        
        n_total = len(self.current_members)
        n_members = len(members)
        n_rejected = len(non_members)
        rejection_rate = 100 * n_rejected / n_total if n_total > 0 else 0
        
        # Identify BGG (Most Massive Galaxy) - search ALL galaxies in radius (spec-z + photo-z)
        self.bgg_data = None
        self.bgg_idx = None
        self.bgg_has_specz = False
        
        # Search in FULL catalog (not just spec-z) within the search radius
        group_ra = self.groups.loc[self.current_group_idx, self.ra_col]
        group_dec = self.groups.loc[self.current_group_idx, self.dec_col]
        group_z = self.group_z_adjusted if self.group_z_adjusted is not None else self.groups.loc[self.current_group_idx, self.z_col]
        
        # Calculate projected distance to all galaxies in full catalog
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        
        group_coord = SkyCoord(ra=group_ra*u.deg, dec=group_dec*u.deg)
        full_coords = SkyCoord(ra=self.full_catalog['RA_DETEC'].values*u.deg, 
                               dec=self.full_catalog['DEC_DETEC'].values*u.deg)
        
        angular_sep = group_coord.separation(full_coords)
        D_A = cosmo.angular_diameter_distance(group_z)
        proj_distance = (angular_sep.radian * D_A).to(u.kpc).value
        
        # Find galaxies within search radius
        in_radius = proj_distance < self.radius_kpc
        
        if in_radius.sum() > 0:
            nearby_galaxies = self.full_catalog[in_radius].copy()
            nearby_galaxies['proj_distance'] = proj_distance[in_radius]
            
            # Import membership probability function
            import sys
            sys.path.insert(0, str(BASE_DIR / 'membership_determination'))
            from membership_functions import calculate_photoz_membership_probability
            
            # Filter by valid photo-z and mass
            has_photoz = ~nearby_galaxies['LP_zfinal'].isna()
            has_mass = ~nearby_galaxies['LP_mass_med_PDF'].isna()
            valid_for_bgg = has_photoz & has_mass
            
            if valid_for_bgg.sum() > 0:
                candidates = nearby_galaxies[valid_for_bgg].copy()
                
                # Calculate membership probability for each candidate
                # Use photo-z uncertainty (LP_zPDF_u68 - LP_zPDF_l68)/2 if available
                if 'LP_zPDF_u68' in candidates.columns and 'LP_zPDF_l68' in candidates.columns:
                    z_err = (candidates['LP_zPDF_u68'] - candidates['LP_zPDF_l68']) / 2.0
                    # Use 0.05 as default for missing values (typical COSMOS photo-z uncertainty)
                    z_err = z_err.fillna(0.05)
                else:
                    z_err = np.full(len(candidates), 0.05)  # Default NMAD for COSMOS
                
                # Calculate membership probability
                mem_prob = calculate_photoz_membership_probability(
                    galaxy_ra=candidates['RA_DETEC'].values,
                    galaxy_dec=candidates['DEC_DETEC'].values,
                    galaxy_z=candidates['LP_zfinal'].values,
                    galaxy_z_err=z_err.values,
                    group_ra=group_ra,
                    group_dec=group_dec,
                    group_z=group_z,
                    radius_kpc=self.radius_kpc,
                    stellar_mass=None,  # Don't use mass weighting (we're selecting by max mass)
                    use_mass_weighting=False
                )
                
                candidates['membership_prob'] = mem_prob
                
                # Select candidates with P > 0.5 (likely members)
                likely_members = candidates[candidates['membership_prob'] > 0.5]
                
                if len(likely_members) > 0:
                    # Find most massive among likely members
                    bgg_idx_in_likely = np.argmax(likely_members['LP_mass_med_PDF'].values)
                    self.bgg_data = likely_members.iloc[bgg_idx_in_likely]
                    self.bgg_has_specz = False  # From full catalog (may or may not have spec-z)
                    bgg_mem_prob = self.bgg_data['membership_prob']
                    
                    # Check if this BGG also has spec-z and is in VRF members
                    bgg_ra = self.bgg_data['RA_DETEC']
                    bgg_dec = self.bgg_data['DEC_DETEC']
                    
                    # Match to spec-z catalog (within 1 arcsec)
                    if len(self.current_members) > 0:
                        ra_match = np.abs(self.current_members['RA'] - bgg_ra) < (1/3600)
                        dec_match = np.abs(self.current_members['DEC'] - bgg_dec) < (1/3600)
                        match = ra_match & dec_match
                        
                        if match.sum() > 0:
                            self.bgg_has_specz = True
                            matched_galaxy = self.current_members[match].iloc[0]
                            is_member = matched_galaxy[member_col] if member_col in matched_galaxy else False
                            bgg_status = "VRF MEMBER" if is_member else "REJECTED"
                            bgg_z = matched_galaxy['zfin']
                        else:
                            bgg_status = "NO SPEC-Z"
                            bgg_z = self.bgg_data['LP_zfinal']
                    else:
                        bgg_status = "NO SPEC-Z"
                        bgg_z = self.bgg_data['LP_zfinal']
                    
                    # Calculate offset from catalog center
                    bgg_coord = SkyCoord(ra=bgg_ra*u.deg, dec=bgg_dec*u.deg)
                    offset_angle = group_coord.separation(bgg_coord)
                    offset_kpc = (offset_angle.radian * D_A).to(u.kpc).value
                    
                    # Display in log10 format (e.g., 10^11.5 instead of 1.15e+11)
                    mass_log = self.bgg_data['LP_mass_med_PDF']
                    bgg_photoz = self.bgg_data['LP_zfinal']
                    delta_z = abs(bgg_z - group_z)
                    
                    print(f"  BGG (Most Massive Galaxy):")
                    print(f"    Stellar Mass: 10^{mass_log:.2f} M☉ (log10 M* = {mass_log:.2f})")
                    print(f"    RA: {bgg_ra:.5f}, Dec: {bgg_dec:.5f}")
                    print(f"    z = {bgg_z:.4f} (photo-z = {bgg_photoz:.4f}, Δz = {delta_z:.4f})")
                    print(f"    Membership probability: {bgg_mem_prob:.3f}")
                    print(f"    Offset from catalog center: {offset_kpc:.1f} kpc")
                    print(f"    Status: {bgg_status}")
                    
                    # Note if BGG is significantly offset (>100 kpc) - suggest recentering
                    if offset_kpc > 100:
                        print(f"\n  ⚠️  BGG offset = {offset_kpc:.1f} kpc > 100 kpc threshold")
                        print(f"  💡 Recommendation: Click 'Center on BGG' button to recenter and rerun VRF")
                        print(f"     (This typically improves membership and mass estimates)")
                else:
                    print(f"  No BGG identified (no likely members with P > 0.5)!")
                    print(f"    Searched {len(candidates)} galaxies in radius with mass & photo-z")
                    print(f"    Max membership probability: {candidates['membership_prob'].max():.3f}")
            else:
                print(f"  No BGG identified (no galaxies with mass + photo-z in radius)!")
                print(f"    Found {has_photoz.sum()} with photo-z, {has_mass.sum()} with mass")
        else:
            print(f"  No BGG identified (no galaxies in radius)!")
        
        # Reset auto-recenter flag if it exists
        if hasattr(self, '_auto_recenter_done'):
            delattr(self, '_auto_recenter_done')
        
        # Calculate velocity dispersion and extract VRF properties
        c_kms = 299792.458
        if n_members >= 3:
            member_z = members['zfin'].values
            self.vrf_z_median = np.median(member_z)
            member_velocities = c_kms * (member_z - self.vrf_z_median) / (1 + self.vrf_z_median)
            
            # Use robust biweight scale estimator for velocity dispersion
            # Better than np.std() for groups with small N or outliers
            from astropy.stats import biweight_scale
            from scipy.stats import bootstrap
            
            if n_members >= 5:
                # Biweight scale (robust, efficient for N ≥ 5)
                self.vrf_sigma_v = biweight_scale(member_velocities, c=9.0)
                
                # Bootstrap uncertainty estimation (1000 resamples)
                if n_members >= 10:  # Only compute uncertainty for reasonable N
                    try:
                        rng = np.random.default_rng(seed=42)
                        res = bootstrap(
                            (member_velocities,),
                            lambda v: biweight_scale(v[0], c=9.0),
                            n_resamples=1000,
                            random_state=rng,
                            method='percentile'
                        )
                        self.vrf_sigma_v_err = res.standard_error
                    except:
                        self.vrf_sigma_v_err = None
                else:
                    self.vrf_sigma_v_err = None
            else:
                # Fall back to std for very small N (biweight needs ≥5)
                self.vrf_sigma_v = np.std(member_velocities)
                self.vrf_sigma_v_err = None
            
            # Also calculate alternative methods for comparison
            self.vrf_sigma_v_std = np.std(member_velocities)  # Standard deviation
            self.vrf_sigma_v_gapper = gapper_estimator(member_velocities)  # Gapper estimator
            
            if n_members >= 5:
                self.vrf_sigma_v_biweight = self.vrf_sigma_v
            else:
                self.vrf_sigma_v_biweight = None
            
            # Try to extract VRF properties if available
            if 'vrf_M200' in members.columns and len(members) > 0:
                self.vrf_M200 = members['vrf_M200'].iloc[0]
            if 'vrf_Rvir' in members.columns and len(members) > 0:
                self.vrf_Rvir = members['vrf_Rvir'].iloc[0]
        else:
            self.vrf_sigma_v = np.nan
            self.vrf_sigma_v_err = None
            self.vrf_sigma_v_std = np.nan
            self.vrf_sigma_v_gapper = np.nan
            self.vrf_sigma_v_biweight = None
            self.vrf_z_median = working_z
            self.vrf_M200 = None
            self.vrf_Rvir = None
        
        print(f"  Candidates: {n_total}")
        print(f"  VRF Members: {n_members} ({100*n_members/n_total:.1f}%)")
        print(f"  Rejected: {n_rejected} ({rejection_rate:.1f}%)")
        
        # Display velocity dispersion with uncertainty if available
        if not np.isnan(self.vrf_sigma_v):
            if self.vrf_sigma_v_err is not None:
                print(f"  σ_v = {self.vrf_sigma_v:.0f} ± {self.vrf_sigma_v_err:.0f} km/s (biweight)")
            else:
                print(f"  σ_v = {self.vrf_sigma_v:.0f} km/s")
            
            # Display multi-method comparison if available (show all 3 methods)
            if n_members >= 3:
                print(f"  σ_v Methods:")
                print(f"    Biweight: {self.vrf_sigma_v:.0f} km/s (robust, recommended)")
                if not np.isnan(self.vrf_sigma_v_gapper):
                    print(f"    Gapper:   {self.vrf_sigma_v_gapper:.0f} km/s (robust, no distribution assumption)")
                print(f"    Std Dev:  {self.vrf_sigma_v_std:.0f} km/s")
                
                # Show max difference between methods
                methods = [self.vrf_sigma_v, self.vrf_sigma_v_gapper, self.vrf_sigma_v_std]
                methods = [m for m in methods if not np.isnan(m)]
                if len(methods) > 1:
                    max_diff_pct = (max(methods) - min(methods)) / np.mean(methods) * 100
                    if max_diff_pct > 10:
                        print(f"    → Max spread: {max_diff_pct:.1f}% (outliers likely present)")
        else:
            print(f"  σ_v = N/A")
        
        if self.vrf_z_median is not None:
            print(f"  z_VRF = {self.vrf_z_median:.4f}")
        if self.vrf_M200 is not None:
            print(f"  M_halo (M_200) = {self.vrf_M200:.2e} M☉")
        if self.vrf_Rvir is not None:
            print(f"  R_vir = {self.vrf_Rvir:.1f} kpc")
        
        # Print halo-to-stellar mass ratio if both available
        if self.vrf_M200 is not None and self.bgg_data is not None:
            bgg_mass = self.bgg_data['LP_mass_med_PDF']
            mass_ratio = bgg_mass / self.vrf_M200  # Stellar-to-halo mass ratio
            print(f"  M*_BGG / M_halo = {mass_ratio:.4f} ({mass_ratio*100:.2f}%)")
        
        # Clear axes
        self.clear_all_axes()
        
        # Panel 1: Spatial distribution (RA-Dec)
        self.plot_spatial(group, members, non_members)
        
        # Panel 2: Velocity vs Radius
        self.plot_velocity_radius(group, members, non_members)
        
        # Panel 3: Velocity distribution
        self.plot_cumulative_velocity(self.current_members, members)
        
        # Panel 4: NFW profile
        self.plot_nfw_profile(members, working_z)
        
        # Update title with statistics - emphasize halo mass and BGG
        title_line1 = f'Group {group["id"]} | z_cat={group["z"]:.4f} z_VRF={self.vrf_z_median:.4f} | N={n_members}/{n_total} ({100*n_members/n_total:.0f}%)'
        
        title_line2 = ''
        if self.vrf_M200 is not None:
            title_line2 += f'M_halo={self.vrf_M200:.2e} M☉'
        if self.vrf_Rvir is not None:
            title_line2 += f' | R_vir={self.vrf_Rvir:.0f} kpc'
        
        # Display all three σ_v estimates side-by-side for comparison
        title_line3 = ''
        if not np.isnan(self.vrf_sigma_v):
            if self.vrf_sigma_v_err is not None:
                title_line3 += f'σ_v: Biweight={self.vrf_sigma_v:.0f}±{self.vrf_sigma_v_err:.0f}'
            else:
                title_line3 += f'σ_v: Biweight={self.vrf_sigma_v:.0f}'
            
            # Add Gapper and std deviation comparison
            if not np.isnan(self.vrf_sigma_v_gapper):
                title_line3 += f' | Gapper={self.vrf_sigma_v_gapper:.0f}'
            if not np.isnan(self.vrf_sigma_v_std):
                title_line3 += f' | Std={self.vrf_sigma_v_std:.0f} km/s'
            
            # Add difference indicator if methods disagree significantly
            if not np.isnan(self.vrf_sigma_v_std) and not np.isnan(self.vrf_sigma_v_gapper):
                methods = [self.vrf_sigma_v, self.vrf_sigma_v_gapper, self.vrf_sigma_v_std]
                max_diff_pct = (max(methods) - min(methods)) / np.mean(methods) * 100
                if max_diff_pct > 15:
                    title_line3 += f' (Δ={max_diff_pct:.0f}%)'
        
        if self.bgg_data is not None:
            bgg_mass = self.bgg_data['LP_mass_med_PDF']
            if title_line3:
                title_line3 += f' | BGG M*={bgg_mass:.2e} M☉'
            else:
                title_line3 = f'BGG M*={bgg_mass:.2e} M☉'
        
        # Combine title lines
        if title_line3:
            full_title = title_line1 + '\n' + title_line2 + '\n' + title_line3
        else:
            full_title = title_line1 + '\n' + title_line2
        
        self.fig.suptitle(full_title, fontsize=11, fontweight='bold')
        
        self.fig.canvas.draw_idle()
    
    def plot_spatial(self, group, members, non_members):
        """Plot spatial distribution in RA-Dec with BGG highlighted."""
        ax = self.ax_spatial
        
        # Plot members (with size scaled by mass if available)
        if len(members) > 0:
            if 'LP_mass_med_PDF' in members.columns:
                # Size galaxies by stellar mass
                masses = members['LP_mass_med_PDF'].values
                valid_mass = ~np.isnan(masses)
                if valid_mass.sum() > 0:
                    # Normalize sizes: 30-150 based on log(mass)
                    log_masses = np.log10(masses[valid_mass])
                    mass_min, mass_max = log_masses.min(), log_masses.max()
                    if mass_max > mass_min:
                        sizes = 30 + 120 * (log_masses - mass_min) / (mass_max - mass_min)
                    else:
                        sizes = 50
                    
                    # Plot members with valid mass
                    ax.scatter(members.loc[valid_mass, 'RA'], members.loc[valid_mass, 'DEC'], 
                              c='blue', s=sizes, alpha=0.6, 
                              label=f'Members (N={len(members)})', 
                              edgecolors='darkblue', linewidths=1, zorder=2)
                    
                    # Plot members without mass data
                    if (~valid_mass).sum() > 0:
                        ax.scatter(members.loc[~valid_mass, 'RA'], members.loc[~valid_mass, 'DEC'],
                                  c='lightblue', s=30, alpha=0.5, edgecolors='blue', linewidths=0.5, zorder=2)
                else:
                    ax.scatter(members['RA'], members['DEC'], c='blue', s=50, alpha=0.7, 
                              label=f'Members (N={len(members)})', edgecolors='darkblue', linewidths=1, zorder=2)
            else:
                ax.scatter(members['RA'], members['DEC'], c='blue', s=50, alpha=0.7, 
                          label=f'Members (N={len(members)})', edgecolors='darkblue', linewidths=1, zorder=2)
        
        # Plot BGG (most massive galaxy) with star marker
        if self.bgg_data is not None:
            bgg_mass = self.bgg_data['LP_mass_med_PDF']
            bgg_ra = self.bgg_data['RA_DETEC']
            bgg_dec = self.bgg_data['DEC_DETEC']
            
            ax.scatter(bgg_ra, bgg_dec, 
                      c='gold', marker='*', s=500, 
                      edgecolors='darkorange', linewidths=3,
                      label=f'BGG: M*=10^{bgg_mass:.2f} M☉', zorder=10)
            
            # Add annotation for BGG with mass in log10 format
            mass_text = f'BGG\nM*=10^{bgg_mass:.2f}'
            ax.annotate(mass_text, xy=(bgg_ra, bgg_dec),
                       xytext=(15, 15), textcoords='offset points',
                       fontsize=9, fontweight='bold', color='darkorange',
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8, edgecolor='darkorange', linewidth=2),
                       arrowprops=dict(arrowstyle='->', color='darkorange', lw=2.5))
        
        # Plot non-members
        if len(non_members) > 0:
            ax.scatter(non_members['RA'], non_members['DEC'],
                      c='red', marker='x', s=50, alpha=0.5,
                      label=f'Rejected (N={len(non_members)})',
                      linewidths=1.5, zorder=1)
        
        # Plot group center (catalog position)
        ax.scatter(group['ra'], group['dec'], 
                  c='black', marker='+', s=300, 
                  linewidths=3,
                  label='Catalog Center', zorder=9)
        
        # Plot search radius circle
        radius_arcsec = self.radius_kpc / (cosmo.kpc_proper_per_arcmin(group['z']).value / 60.0)
        radius_deg = radius_arcsec / 3600.0
        
        circle = Circle((group['ra'], group['dec']), radius_deg,
                       fill=False, edgecolor='black', linestyle='--',
                       linewidth=1.5, alpha=0.5, label=f'{self.radius_kpc} kpc', zorder=0)
        ax.add_patch(circle)
        
        ax.set_xlabel('RA (deg)', fontsize=11)
        ax.set_ylabel('Dec (deg)', fontsize=11)
        ax.set_title('Spatial Distribution (size ∝ stellar mass)', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()  # RA increases to the left
        ax.set_aspect('equal', adjustable='box')
    
    def plot_velocity_radius(self, group, members, non_members):
        """Plot velocity offset vs projected radius with BGG highlighted."""
        ax = self.ax_velocity
        
        # Plot members (with size scaled by mass)
        if len(members) > 0:
            if 'LP_mass_med_PDF' in members.columns:
                masses = members['LP_mass_med_PDF'].values
                valid_mass = ~np.isnan(masses)
                if valid_mass.sum() > 0:
                    log_masses = np.log10(masses[valid_mass])
                    mass_min, mass_max = log_masses.min(), log_masses.max()
                    if mass_max > mass_min:
                        sizes = 30 + 120 * (log_masses - mass_min) / (mass_max - mass_min)
                    else:
                        sizes = 50
                    
                    ax.scatter(members.loc[valid_mass, 'sep_kpc'], members.loc[valid_mass, 'dv'],
                              c='blue', s=sizes, alpha=0.6,
                              label=f'Members (N={len(members)})',
                              edgecolors='darkblue', linewidths=1, zorder=2)
                    
                    if (~valid_mass).sum() > 0:
                        ax.scatter(members.loc[~valid_mass, 'sep_kpc'], members.loc[~valid_mass, 'dv'],
                                  c='lightblue', s=30, alpha=0.5, edgecolors='blue', linewidths=0.5, zorder=2)
                else:
                    ax.scatter(members['sep_kpc'], members['dv'], c='blue', s=50, alpha=0.7,
                              label=f'Members (N={len(members)})', edgecolors='darkblue', linewidths=1, zorder=2)
            else:
                ax.scatter(members['sep_kpc'], members['dv'], c='blue', s=50, alpha=0.7,
                          label=f'Members (N={len(members)})', edgecolors='darkblue', linewidths=1, zorder=2)
        
        # Highlight BGG (only if it has spec-z and is in the candidate list)
        if self.bgg_data is not None and self.bgg_has_specz and 'dv' in self.bgg_data:
            bgg_dv = self.bgg_data['dv']
            bgg_sep = self.bgg_data['sep_kpc']
            ax.scatter(bgg_sep, bgg_dv,
                      c='gold', marker='*', s=500,
                      edgecolors='darkorange', linewidths=3,
                      label='BGG (has spec-z)', zorder=10)
            
            # Add vertical/horizontal lines to BGG
            ax.axvline(bgg_sep, color='gold', linestyle=':', linewidth=1.5, alpha=0.5, zorder=1)
            ax.axhline(bgg_dv, color='gold', linestyle=':', linewidth=1.5, alpha=0.5, zorder=1)
        
        # Plot non-members
        if len(non_members) > 0:
            ax.scatter(non_members['sep_kpc'], non_members['dv'],
                      c='red', marker='x', s=50, alpha=0.5,
                      label=f'Rejected (N={len(non_members)})',
                      linewidths=1.5, zorder=1)
        
        # Plot velocity limits
        ax.axhline(self.max_velocity, color='gray', linestyle='--', 
                  linewidth=1.5, alpha=0.7, label=f'±{self.max_velocity} km/s')
        ax.axhline(-self.max_velocity, color='gray', linestyle='--', 
                  linewidth=1.5, alpha=0.7)
        ax.axhline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        
        ax.set_xlabel('Projected Radius (kpc)', fontsize=11)
        ax.set_ylabel('Velocity Offset (km/s)', fontsize=11)
        ax.set_title('Velocity vs Radius (size ∝ stellar mass)', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, self.radius_kpc * 1.1)
    
    def plot_cumulative_velocity(self, all_candidates, members):
        """Plot velocity distribution with normal distribution overlay."""
        ax = self.ax_cumulative
        
        member_vel = members['dv'].values if len(members) > 0 else np.array([])
        
        if len(member_vel) > 0:
            # Calculate histogram
            n_bins = min(30, max(10, len(member_vel) // 3))
            counts, bins, patches = ax.hist(member_vel, bins=n_bins, 
                                           density=True, alpha=0.6, 
                                           color='blue', edgecolor='black',
                                           label='VRF Members')
            
            # Calculate mean and std for normal distribution
            mean_vel = np.mean(member_vel)
            std_vel = np.std(member_vel)
            
            # Plot normal distribution overlay
            vel_range = np.linspace(bins[0], bins[-1], 200)
            normal_dist = norm.pdf(vel_range, mean_vel, std_vel)
            ax.plot(vel_range, normal_dist, 'r-', linewidth=2.5, 
                   label=f'Normal(μ={mean_vel:.0f}, σ={std_vel:.0f})', alpha=0.8)
            
            # Add mean line
            ax.axvline(mean_vel, color='darkblue', linestyle='--', 
                      linewidth=2, alpha=0.7, label=f'Mean={mean_vel:.0f} km/s')
            
            # Add ±1σ and ±2σ lines
            ax.axvline(mean_vel + std_vel, color='orange', linestyle=':', 
                      linewidth=1.5, alpha=0.6, label=f'±1σ')
            ax.axvline(mean_vel - std_vel, color='orange', linestyle=':', 
                      linewidth=1.5, alpha=0.6)
            ax.axvline(mean_vel + 2*std_vel, color='red', linestyle=':', 
                      linewidth=1.5, alpha=0.5, label=f'±2σ')
            ax.axvline(mean_vel - 2*std_vel, color='red', linestyle=':', 
                      linewidth=1.5, alpha=0.5)
        
        # Mark velocity limits
        ax.axvline(self.max_velocity, color='gray', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'Cut: ±{self.max_velocity:.0f} km/s')
        ax.axvline(-self.max_velocity, color='gray', linestyle='--', 
                  linewidth=2, alpha=0.7)
        ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
        
        ax.set_xlabel('Velocity Offset (km/s)', fontsize=11)
        ax.set_ylabel('Probability Density', fontsize=11)
        ax.set_title('Velocity Distribution vs Normal', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
    
    def plot_nfw_profile(self, members, group_z):
        """Plot NFW velocity dispersion profile with data and theoretical model."""
        ax = self.ax_nfw
        
        if len(members) < 3:
            ax.text(0.5, 0.5, 'Too few members\nfor NFW fit', 
                   ha='center', va='center', fontsize=14, transform=ax.transAxes)
            ax.set_title('NFW Velocity Dispersion Profile', fontsize=12, fontweight='bold')
            return
        
        # Get member data
        radii = members['sep_kpc'].values
        c_kms = 299792.458
        velocities = np.abs((members['zfin'].values - self.vrf_z_median) * c_kms / (1 + self.vrf_z_median))
        
        # Sort by radius for better visualization
        sort_idx = np.argsort(radii)
        radii_sorted = radii[sort_idx]
        velocities_sorted = velocities[sort_idx]
        
        # === PLOT 1: Individual Members (gray points, size by mass) ===
        if 'LP_mass_med_PDF' in members.columns:
            masses = members['LP_mass_med_PDF'].values
            valid_mass = ~np.isnan(masses)
            if valid_mass.sum() > 0:
                log_masses = np.log10(masses[valid_mass])
                mass_min, mass_max = log_masses.min(), log_masses.max()
                if mass_max > mass_min:
                    sizes = 20 + 80 * (log_masses - mass_min) / (mass_max - mass_min)
                else:
                    sizes = 30
                
                ax.scatter(radii[valid_mass], velocities[valid_mass], 
                          s=sizes, color='lightgray', alpha=0.4,
                          edgecolors='gray', linewidths=0.5, label='Individual Members', zorder=1)
                
                if (~valid_mass).sum() > 0:
                    ax.scatter(radii[~valid_mass], velocities[~valid_mass],
                              s=20, color='lightgray', alpha=0.3, edgecolors='gray', linewidths=0.5, zorder=1)
            else:
                ax.scatter(radii, velocities, alpha=0.4, s=30, color='lightgray', 
                          edgecolors='gray', linewidths=0.5, label='Individual Members', zorder=1)
        else:
            ax.scatter(radii, velocities, alpha=0.4, s=30, color='lightgray', 
                      edgecolors='gray', linewidths=0.5, label='Individual Members', zorder=1)
        
        # Highlight BGG in NFW plot (only if it has spec-z)
        if self.bgg_data is not None and self.bgg_has_specz and 'zfin' in self.bgg_data:
            c_kms = 299792.458
            bgg_vel = np.abs((self.bgg_data['zfin'] - self.vrf_z_median) * c_kms / (1 + self.vrf_z_median))
            bgg_radius = self.bgg_data['sep_kpc']
            ax.scatter(bgg_radius, bgg_vel,
                      c='gold', marker='*', s=500,
                      edgecolors='darkorange', linewidths=3,
                      label='BGG (has spec-z)', zorder=10)
        
        # === PLOT 2: Binned Data (blue points with error bars) ===
        n_bins = min(8, max(3, len(members) // 6))
        bin_edges = np.percentile(radii, np.linspace(0, 100, n_bins + 1))
        bin_centers = []
        bin_sigmas = []
        bin_errors = []
        bin_counts = []
        
        for i in range(n_bins):
            mask = (radii >= bin_edges[i]) & (radii < bin_edges[i+1])
            n_in_bin = mask.sum()
            if n_in_bin >= 2:
                bin_centers.append(np.mean(radii[mask]))
                bin_sigmas.append(np.std(velocities[mask]))
                bin_errors.append(bin_sigmas[-1] / np.sqrt(n_in_bin))
                bin_counts.append(n_in_bin)
        
        if len(bin_centers) > 0:
            # Scale marker size by number of galaxies in bin
            marker_sizes = [50 + 20 * np.sqrt(n) for n in bin_counts]
            ax.errorbar(bin_centers, bin_sigmas, yerr=bin_errors,
                       fmt='o', markersize=10, color='dodgerblue', ecolor='navy',
                       markeredgecolor='navy', markeredgewidth=1.5,
                       capsize=5, capthick=2, alpha=0.9, label='Binned σ_v', zorder=3)
        
        # === PLOT 3: Theoretical NFW Profile (if Rvir available) ===
        if self.vrf_Rvir is not None and self.vrf_Rvir > 0 and not np.isnan(self.vrf_sigma_v):
            # Create smooth radius array for NFW profile
            r_model = np.linspace(0, self.radius_kpc * 1.1, 100)
            
            # NFW velocity dispersion profile (simplified)
            # σ(r) decreases with radius approximately as:
            # σ(r) ≈ σ_vir * sqrt(f(r/R_vir)) where f is NFW function
            r_norm = r_model / self.vrf_Rvir
            
            # Approximate NFW profile (decreasing with radius)
            # Using empirical fit: σ(r) ≈ σ_0 / sqrt(1 + (r/R_vir)^2)
            sigma_profile = self.vrf_sigma_v / np.sqrt(1 + 0.5 * r_norm**2)
            
            ax.plot(r_model, sigma_profile, color='green', linestyle='-', 
                   linewidth=2.5, alpha=0.7, label='NFW Profile (Expected)', zorder=2)
        
        # === PLOT 4: Average σ_v line (dashed red) ===
        if not np.isnan(self.vrf_sigma_v):
            ax.axhline(self.vrf_sigma_v, color='orangered', linestyle='--', 
                      linewidth=2.5, alpha=0.8, label=f'Average σ_v = {self.vrf_sigma_v:.0f} km/s', zorder=2)
        
        # === PLOT 5: R_vir line (vertical red - optimal result location) ===
        if self.vrf_Rvir is not None and self.vrf_Rvir > 0:
            ax.axvline(self.vrf_Rvir, color='red', linestyle='-', 
                      linewidth=3, alpha=0.85, label=f'R_vir = {self.vrf_Rvir:.0f} kpc (Optimal)', zorder=4)
            
            # Add shaded region around R_vir (±20%)
            ax.axvspan(self.vrf_Rvir * 0.8, self.vrf_Rvir * 1.2, 
                      color='red', alpha=0.1, zorder=0)
        
        # === Add text annotation for halo, BGG, and σ_v properties ===
        textstr = f'N_members: {len(members)}\n'
        
        # Add all three σ_v estimates (Biweight, Gapper, Std)
        if not np.isnan(self.vrf_sigma_v):
            textstr += '\nVelocity Dispersion:\n'
            if self.vrf_sigma_v_err is not None:
                textstr += f'  Biweight: {self.vrf_sigma_v:.0f} ± {self.vrf_sigma_v_err:.0f} km/s\n'
            else:
                textstr += f'  Biweight: {self.vrf_sigma_v:.0f} km/s\n'
            
            if not np.isnan(self.vrf_sigma_v_gapper):
                textstr += f'  Gapper:   {self.vrf_sigma_v_gapper:.0f} km/s\n'
            
            if not np.isnan(self.vrf_sigma_v_std):
                textstr += f'  Std Dev:  {self.vrf_sigma_v_std:.0f} km/s\n'
                
                # Show spread between methods
                methods = [self.vrf_sigma_v, self.vrf_sigma_v_gapper, self.vrf_sigma_v_std]
                methods = [m for m in methods if not np.isnan(m)]
                if len(methods) > 1:
                    max_diff_pct = (max(methods) - min(methods)) / np.mean(methods) * 100
                    if max_diff_pct > 10:
                        textstr += f'  (Spread: {max_diff_pct:.1f}%)\n'
        
        if self.vrf_M200 is not None:
            textstr += f'\nM_halo: {self.vrf_M200:.2e} M☉\n'
        if self.bgg_data is not None:
            bgg_mass = self.bgg_data['LP_mass_med_PDF']
            textstr += f'BGG M*: {bgg_mass:.2e} M☉\n'
            if self.vrf_M200 is not None:
                mass_ratio = bgg_mass / self.vrf_M200  # Stellar-to-halo
                textstr += f'M*/M_halo: {mass_ratio:.4f} ({mass_ratio*100:.2f}%)'
        
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, 
               fontsize=8, verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=1.5))
        
        ax.set_xlabel('Projected Radius (kpc)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Velocity Dispersion (km/s)', fontsize=11, fontweight='bold')
        ax.set_title('NFW Velocity Dispersion Profile', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle=':', linewidth=1)
        ax.set_xlim(0, self.radius_kpc * 1.1)
        
        # Set y-axis limits intelligently
        if len(velocities) > 0:
            y_max = max(np.max(velocities), self.vrf_sigma_v * 1.5 if not np.isnan(self.vrf_sigma_v) else 0)
            ax.set_ylim(0, y_max * 1.15)
    
    def clear_all_axes(self):
        """Clear all plot axes."""
        for ax in [self.ax_spatial, self.ax_velocity, self.ax_cumulative, self.ax_nfw]:
            ax.clear()
    
    def on_slider_change(self, val):
        """Handle slider changes."""
        self.max_dz_norm = self.slider_dz.val
        self.max_velocity = self.slider_vel.val
        self.radius_kpc = self.slider_radius.val
        self.group_z_adjusted = self.slider_z.val
        
        # Update textboxes to match sliders
        self.textbox_z.set_val(f'{self.group_z_adjusted:.4f}')
        self.textbox_radius.set_val(f'{self.radius_kpc:.0f}')
        self.textbox_dz.set_val(f'{self.max_dz_norm:.4f}')
        self.textbox_vel.set_val(f'{self.max_velocity:.0f}')
    
    def on_z_textbox_submit(self, text):
        """Handle group redshift textbox input."""
        try:
            val = float(text)
            if 0.0 <= val <= 2.0:
                self.group_z_adjusted = val
                self.slider_z.set_val(val)
                self.update_plot()
            else:
                print(f"Warning: Redshift {val:.4f} out of range [0.0, 2.0]")
        except ValueError:
            print(f"Warning: Invalid redshift input '{text}'")
    
    def on_radius_textbox_submit(self, text):
        """Handle radius textbox input."""
        try:
            val = float(text)
            if 200 <= val <= 1000:
                self.radius_kpc = val
                self.slider_radius.set_val(val)
                self.update_plot()
            else:
                print(f"Warning: Radius {val:.0f} kpc out of range [200, 1000]")
        except ValueError:
            print(f"Warning: Invalid radius input '{text}'")
    
    def on_dz_textbox_submit(self, text):
        """Handle max_dz_norm textbox input."""
        try:
            val = float(text)
            if 0.005 <= val <= 0.03:
                self.max_dz_norm = val
                self.slider_dz.set_val(val)
                self.update_plot()
            else:
                print(f"Warning: max_dz_norm {val:.4f} out of range [0.005, 0.03]")
        except ValueError:
            print(f"Warning: Invalid max_dz_norm input '{text}'")
    
    def on_vel_textbox_submit(self, text):
        """Handle max_velocity textbox input."""
        try:
            val = float(text)
            if 1000 <= val <= 4000:
                self.max_velocity = val
                self.slider_vel.set_val(val)
                self.update_plot()
            else:
                print(f"Warning: max_velocity {val:.0f} km/s out of range [1000, 4000]")
        except ValueError:
            print(f"Warning: Invalid max_velocity input '{text}'")
    
    def on_rerun_click(self, event):
        """Handle Rerun VRF button click."""
        print("\nRerunning VRF with updated parameters...")
        self.update_plot()
    
    def on_accept_click(self, event):
        """Handle Accept & Next button click."""
        group = self.get_current_group()
        
        # Save current group results
        result = {
            'group_id': group['id'],
            'group_z_catalog': group['z'],
            'group_z_adjusted': self.group_z_adjusted,
            'VRF_group_z': self.vrf_z_median if self.vrf_z_median is not None else self.group_z_adjusted,
            'VRF_M200': self.vrf_M200,
            'VRF_Rvir': self.vrf_Rvir,
            'VRF_sigma_v': self.vrf_sigma_v,
            'radius_kpc': self.radius_kpc,
            'max_dz_norm': self.max_dz_norm,
            'max_velocity': self.max_velocity,
            'n_candidates': len(self.current_members),
            'n_members': len(self.current_members[self.current_members['vrf_member']]) if 'vrf_member' in self.current_members.columns else 0,
            'status': 'accepted'
        }
        self.reviewed_groups.append(result)
        
        print(f"\n✓ Accepted group {group['id']}")
        print(f"  z_catalog={group['z']:.4f}, z_adjusted={self.group_z_adjusted:.4f}, z_VRF={result['VRF_group_z']:.4f}")
        print(f"  Total reviewed: {len(self.reviewed_groups)}")
        
        # Move to next group
        if self.current_group_idx < len(self.groups) - 1:
            self.current_group_idx += 1
            # Reset redshift for next group
            self.group_z_original = None
            self.group_z_adjusted = None
            self.update_plot()
        else:
            print("\n  Last group reached!")
            self.on_save_click(None)
    
    def on_prev_click(self, event):
        """Handle Previous button click."""
        if self.current_group_idx > 0:
            self.current_group_idx -= 1
            # Reset redshift for previous group
            self.group_z_original = None
            self.group_z_adjusted = None
            self.update_plot()
        else:
            print("\n  Already at first group!")
    
    def on_jump_submit(self, text):
        """Handle jump to group textbox submission."""
        try:
            idx = int(text)
            if 0 <= idx < len(self.groups):
                self.current_group_idx = idx
                # Reset redshift for jumped-to group
                self.group_z_original = None
                self.group_z_adjusted = None
                self.update_plot()
            else:
                print(f"\n  Invalid index: {idx} (must be 0-{len(self.groups)-1})")
        except ValueError:
            print(f"\n  Invalid input: {text}")
    
    def on_save_click(self, event):
        """Handle Save Results button click."""
        if len(self.reviewed_groups) == 0:
            print("\n  No groups reviewed yet!")
            return
        
        # Save reviewed groups
        df_reviewed = pd.DataFrame(self.reviewed_groups)
        
        timestamp = pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')
        output_file = OUTPUT_DIR / f'reviewed_groups_{self.catalog_name}_{timestamp}.csv'
        
        df_reviewed.to_csv(output_file, index=False)
        
        print(f"\n{'='*60}")
        print(f"SAVED RESULTS")
        print(f"{'='*60}")
        print(f"  File: {output_file}")
        print(f"  Groups reviewed: {len(self.reviewed_groups)}")
        print(f"  Catalog: {self.catalog_name}")
        print(f"{'='*60}")
    
    def on_recenter_bgg_click(self, event):
        """Handle Center on BGG button - recenter group on most massive galaxy."""
        if self.bgg_data is None:
            print("\n  No BGG identified (no mass data or no members)!")
            return
        
        group = self.get_current_group()
        bgg_ra = self.bgg_data['RA_DETEC']
        bgg_dec = self.bgg_data['DEC_DETEC']
        bgg_mass = self.bgg_data['LP_mass_med_PDF']
        
        print(f"\n{'='*60}")
        print("RECENTERING ON BGG (MOST MASSIVE GALAXY)")
        print(f"{'='*60}")
        print(f"Original Center: RA={group['ra']:.5f}, Dec={group['dec']:.5f}")
        print(f"BGG Position:    RA={bgg_ra:.5f}, Dec={bgg_dec:.5f}")
        print(f"BGG Mass: 10^{bgg_mass:.2f} M☉ (log10 M* = {bgg_mass:.2f})")
        
        # Update group center in catalog temporarily
        old_ra = group['ra']
        old_dec = group['dec']
        
        self.groups.loc[self.current_group_idx, self.ra_col] = bgg_ra
        self.groups.loc[self.current_group_idx, self.dec_col] = bgg_dec
        
        # Rerun membership with new center
        print(f"\nRerunning VRF with BGG-centered coordinates...")
        self.update_plot()
        
        print(f"{'='*60}")
        print("Note: Center changed temporarily. Use 'Accept' to save or 'Previous/Next' to revert.")
        print(f"{'='*60}")
    
    def on_optimize_click(self, event):
        """Handle Auto-Optimize button click - intelligently find best parameters using coarse-to-fine search."""
        print("\n" + "="*70)
        print("STARTING SMART AUTO-OPTIMIZATION (Coarse-to-Fine Strategy)")
        print("="*70)
        
        group = self.get_current_group()
        print(f"Group: {group['id']}")
        print(f"Original z_catalog: {group['z']:.4f}")
        
        # Store original parameters
        original_z = self.group_z_adjusted
        original_radius = self.radius_kpc
        original_dz = self.max_dz_norm
        original_vel = self.max_velocity
        
        best_config = None
        best_score = -np.inf
        
        # PHASE 1: Coarse grid search (fast exploration)
        print(f"\nPhase 1: Coarse grid search...")
        z_offsets_coarse = np.linspace(-0.02, 0.02, 5)  # 5 z points
        radii_coarse = [400, 500, 600, 700]  # 4 radii
        dz_norms_coarse = [0.010, 0.012, 0.015]  # 3 dz values
        velocities_coarse = [1500, 2000, 2500]  # 3 velocities
        
        total_phase1 = len(z_offsets_coarse) * len(radii_coarse) * len(dz_norms_coarse) * len(velocities_coarse)
        print(f"Testing {total_phase1} combinations (vs 900 in old method)...")
        
        iteration = 0
        for z_offset in z_offsets_coarse:
            test_z = group['z'] + z_offset
            if test_z < 0.0 or test_z > 2.0:
                continue
                
            for radius in radii_coarse:
                for dz_norm in dz_norms_coarse:
                    for velocity in velocities_coarse:
                        iteration += 1
                        if iteration % 20 == 0:
                            print(f"  {iteration}/{total_phase1}...", end='', flush=True)
                        
                        # Set parameters
                        self.group_z_adjusted = test_z
                        self.radius_kpc = radius
                        self.max_dz_norm = dz_norm
                        self.max_velocity = velocity
                        
                        # Run VRF
                        candidates, members, vrf_result = self.run_vrf_membership(group)
                        
                        if vrf_result is None:
                            continue
                        
                        # Calculate quality score
                        score = self.calculate_quality_score(
                            vrf_result, members, candidates, test_z, group['z']
                        )
                        
                        # Track best configuration
                        if score > best_score:
                            best_score = score
                            best_config = {
                                'group_z': test_z,
                                'radius_kpc': radius,
                                'max_dz_norm': dz_norm,
                                'max_velocity': velocity,
                                'n_members': len(members),
                                'n_candidates': len(candidates),
                                'sigma_v': vrf_result.get('sigma_v', np.nan),
                                'M200': vrf_result.get('M200', np.nan),
                                'Rvir': vrf_result.get('Rvir', np.nan),
                                'algorithm_flag': vrf_result.get('algorithm_flag', -1),
                                'z_offset': z_offset,
                                'score': score
                            }
                            
                        # Early stopping: if we find an excellent solution (score > 250), stop Phase 1
                        if best_score > 250:
                            print(f"\n  Early stopping! Found excellent solution (score={best_score:.1f})")
                            break
                    if best_score > 250:
                        break
                if best_score > 250:
                    break
            if best_score > 250:
                break
        
        print(f"\n  Phase 1 complete! Best score: {best_score:.1f}")
        
        # PHASE 2: Fine grid search around best result (skip if Phase 1 found excellent solution)
        if best_config is not None and best_score < 250:
            print(f"\nPhase 2: Fine grid refinement around best parameters...")
            
            # Define fine search around best config
            z_best = best_config['group_z']
            z_offsets_fine = np.linspace(
                max(z_best - 0.005, group['z'] - 0.02),
                min(z_best + 0.005, group['z'] + 0.02),
                5
            )
            
            r_best = best_config['radius_kpc']
            radii_fine = [max(300, r_best - 100), r_best, min(800, r_best + 100)]
            
            dz_best = best_config['max_dz_norm']
            dz_norms_fine = [max(0.008, dz_best - 0.002), dz_best, min(0.020, dz_best + 0.002)]
            
            vel_best = best_config['max_velocity']
            velocities_fine = [max(1000, vel_best - 500), vel_best, min(3000, vel_best + 500)]
            
            total_phase2 = len(z_offsets_fine) * len(radii_fine) * len(dz_norms_fine) * len(velocities_fine)
            print(f"Testing {total_phase2} combinations around best result...")
            
            phase2_count = 0
            for test_z in z_offsets_fine:
                if test_z < 0.0 or test_z > 2.0:
                    continue
                    
                for radius in radii_fine:
                    for dz_norm in dz_norms_fine:
                        for velocity in velocities_fine:
                            phase2_count += 1
                            if phase2_count % 20 == 0:
                                print(f"  {phase2_count}/{total_phase2}...", end='', flush=True)
                            
                            # Set parameters
                            self.group_z_adjusted = test_z
                            self.radius_kpc = radius
                            self.max_dz_norm = dz_norm
                            self.max_velocity = velocity
                            
                            # Run VRF
                            candidates, members, vrf_result = self.run_vrf_membership(group)
                            
                            if vrf_result is None:
                                continue
                            
                            # Calculate quality score
                            score = self.calculate_quality_score(
                                vrf_result, members, candidates, test_z, group['z']
                            )
                            
                            # Update best configuration if better
                            if score > best_score:
                                best_score = score
                                z_offset = test_z - group['z']
                                best_config = {
                                    'group_z': test_z,
                                    'radius_kpc': radius,
                                    'max_dz_norm': dz_norm,
                                    'max_velocity': velocity,
                                    'n_members': len(members),
                                    'n_candidates': len(candidates),
                                    'sigma_v': vrf_result.get('sigma_v', np.nan),
                                    'M200': vrf_result.get('M200', np.nan),
                                    'Rvir': vrf_result.get('Rvir', np.nan),
                                    'algorithm_flag': vrf_result.get('algorithm_flag', -1),
                                    'z_offset': z_offset,
                                    'score': score
                                }
            
            print(f"\n  Phase 2 complete! Final score: {best_score:.1f}")
            total_tests = total_phase1 + total_phase2
        else:
            total_tests = iteration
            if best_score > 250:
                print("  (Skipped Phase 2 due to excellent Phase 1 result)")
        
        print(f"\n  Total VRF runs: {total_tests} (vs 900 in old method = {100*total_tests/900:.0f}% of original)")
        print("  Done!")
        
        
        print(f"\n  Total VRF runs: {total_tests} (vs 900 in old method = {100*total_tests/900:.0f}% of original)")
        print("  Done!")
        
        # Apply best configuration
        if best_config is not None:
            print(f"\n{'='*70}")
            print("OPTIMIZATION RESULTS")
            print(f"{'='*70}")
            print(f"Best Score: {best_score:.2f}")
            print(f"\nOptimal Parameters:")
            print(f"  Group z: {best_config['group_z']:.4f} (offset: {best_config['z_offset']:+.4f})")
            print(f"  Radius: {best_config['radius_kpc']:.0f} kpc")
            print(f"  max_dz_norm: {best_config['max_dz_norm']:.4f}")
            print(f"  max_velocity: {best_config['max_velocity']:.0f} km/s")
            print(f"\nVRF Results:")
            print(f"  Members: {best_config['n_members']} / {best_config['n_candidates']} candidates")
            print(f"  σ_v: {best_config['sigma_v']:.0f} km/s")
            print(f"  M_200: {best_config['M200']:.2e} M☉")
            print(f"  R_vir: {best_config['Rvir']:.0f} kpc")
            print(f"  Algorithm Flag: {best_config['algorithm_flag']}")
            print(f"{'='*70}")
            
            # Update UI with best parameters
            self.group_z_adjusted = best_config['group_z']
            self.radius_kpc = best_config['radius_kpc']
            self.max_dz_norm = best_config['max_dz_norm']
            self.max_velocity = best_config['max_velocity']
            
            # Update sliders and textboxes
            self.slider_z.set_val(best_config['group_z'])
            self.slider_radius.set_val(best_config['radius_kpc'])
            self.slider_dz.set_val(best_config['max_dz_norm'])
            self.slider_vel.set_val(best_config['max_velocity'])
            
            self.textbox_z.set_val(f"{best_config['group_z']:.4f}")
            self.textbox_radius.set_val(f"{best_config['radius_kpc']:.0f}")
            self.textbox_dz.set_val(f"{best_config['max_dz_norm']:.4f}")
            self.textbox_vel.set_val(f"{best_config['max_velocity']:.0f}")
            
            # Update plot
            self.update_plot()
            
        else:
            print("\n  Optimization failed - no valid configuration found!")
            print("  Restoring original parameters...")
            
            # Restore original parameters
            self.group_z_adjusted = original_z
            self.radius_kpc = original_radius
            self.max_dz_norm = original_dz
            self.max_velocity = original_vel
            
            self.slider_z.set_val(original_z)
            self.slider_radius.set_val(original_radius)
            self.slider_dz.set_val(original_dz)
            self.slider_vel.set_val(original_vel)
    
    def run_vrf_membership(self, group):
        """
        Run VRF membership determination with current parameters.
        
        Returns:
            candidates: DataFrame of all candidates
            members: DataFrame of confirmed members
            vrf_result: Dictionary with VRF properties
        """
        # Find spec-z members
        candidates = find_specz_members(
            group['ra'], group['dec'], self.group_z_adjusted, self.specz,
            radius_kpc=self.radius_kpc,
            max_dz_norm=self.max_dz_norm,
            max_velocity=self.max_velocity,
            use_gapper=True,
            remove_duplicates=True,
            method='vrf'
        )
        
        if len(candidates) == 0:
            return None, None, None
        
        # Separate members and non-members
        member_col = 'vrf_member' if 'vrf_member' in candidates.columns else 'gapper_member'
        members = candidates[candidates[member_col] == True]
        
        if len(members) < 3:
            return candidates, members, None
        
        # Extract VRF properties
        vrf_result = {}
        
        c_kms = 299792.458
        member_z = members['zfin'].values
        z_median = np.median(member_z)
        member_velocities = c_kms * (member_z - z_median) / (1 + z_median)
        vrf_result['sigma_v'] = np.std(member_velocities)
        vrf_result['z_median'] = z_median
        
        # Extract VRF properties if available
        if 'vrf_M200' in members.columns:
            vrf_result['M200'] = members['vrf_M200'].iloc[0]
        else:
            vrf_result['M200'] = np.nan
            
        if 'vrf_Rvir' in members.columns:
            vrf_result['Rvir'] = members['vrf_Rvir'].iloc[0]
        else:
            vrf_result['Rvir'] = np.nan
            
        if 'vrf_algorithm_flag' in members.columns:
            vrf_result['algorithm_flag'] = members['vrf_algorithm_flag'].iloc[0]
        else:
            vrf_result['algorithm_flag'] = 0  # Assume success if not provided
        
        return candidates, members, vrf_result
    
    def calculate_quality_score(self, vrf_result, members, candidates, test_z, catalog_z):
        """
        Calculate quality score for VRF result.
        
        Higher score = better configuration
        
        Scoring criteria:
        - Algorithm success (flag=0)
        - Reasonable number of members (10-150)
        - Reasonable velocity dispersion (200-1500 km/s)
        - Reasonable M200 (1e13 - 1e15 M☉)
        - Member fraction not too low/high
        - Redshift close to catalog value
        """
        score = 0.0
        
        # Algorithm success
        algorithm_flag = vrf_result.get('algorithm_flag', -1)
        if algorithm_flag == 0:
            score += 100  # Major bonus for convergence
        elif algorithm_flag == 1:
            score += 20   # Some credit for reasonable result
        else:
            return -1000  # Heavily penalize failed runs
        
        # Number of members
        n_members = len(members)
        if 10 <= n_members <= 150:
            score += 50
            # Bonus for sweet spot (20-80 members)
            if 20 <= n_members <= 80:
                score += 20
        elif n_members < 5:
            score -= 100  # Too few members
        elif n_members > 200:
            score -= 50   # Too many members
        
        # Member fraction
        member_fraction = n_members / len(candidates) if len(candidates) > 0 else 0
        if 0.15 <= member_fraction <= 0.85:
            score += 30
            # Optimal range
            if 0.25 <= member_fraction <= 0.70:
                score += 20
        elif member_fraction > 0.95:
            score -= 50  # Likely not selective enough
        
        # Velocity dispersion
        sigma_v = vrf_result.get('sigma_v', np.nan)
        if not np.isnan(sigma_v):
            if 200 <= sigma_v <= 1500:
                score += 40
                # Optimal range for compact groups
                if 300 <= sigma_v <= 800:
                    score += 30
            elif sigma_v > 2000:
                score -= 50  # Too high
            elif sigma_v < 100:
                score -= 50  # Too low
        
        # Mass reasonability
        M200 = vrf_result.get('M200', np.nan)
        if not np.isnan(M200):
            if 1e13 <= M200 <= 1e15:
                score += 30
                # Optimal for compact groups
                if 1e13 <= M200 <= 5e14:
                    score += 20
            elif M200 > 1e16:
                score -= 40  # Unrealistically high
            elif M200 < 1e12:
                score -= 40  # Too low
        
        # Redshift proximity to catalog value
        z_diff = abs(test_z - catalog_z)
        if z_diff < 0.005:
            score += 20
        elif z_diff < 0.010:
            score += 10
        elif z_diff > 0.030:
            score -= 30
        
        # R_vir reasonability
        Rvir = vrf_result.get('Rvir', np.nan)
        if not np.isnan(Rvir):
            if 300 <= Rvir <= 2000:
                score += 20
            elif Rvir > 3000:
                score -= 30
            elif Rvir < 100:
                score -= 30
        
        # BGG presence and mass (bonus for having massive galaxies)
        if 'LP_mass_med_PDF' in members.columns:
            masses = members['LP_mass_med_PDF'].values
            valid_masses = masses[~np.isnan(masses)]
            
            if len(valid_masses) > 0:
                max_mass = np.max(valid_masses)
                median_mass = np.median(valid_masses)
                
                # Bonus for having a massive BGG (>10^11 M☉)
                if max_mass > 1e11:
                    score += 15
                    if max_mass > 5e11:  # Very massive BGG
                        score += 10
                
                # Bonus for good mass range (indicates complete sampling)
                mass_ratio = max_mass / median_mass if median_mass > 0 else 1
                if 3 <= mass_ratio <= 20:  # Reasonable mass hierarchy
                    score += 10
        
        return score
    
    def show(self):
        """Show the dashboard."""
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Interactive spectroscopic membership dashboard',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog', choices=['cw-all', 'cw-hcg'], required=True,
                       help='Which catalog to use')
    parser.add_argument('--group-id', required=True,
                       help='Starting group ID (integer for cw-all, string for cw-hcg)')
    parser.add_argument('--radius', type=float, default=500,
                       help='Search radius in kpc')
    
    args = parser.parse_args()
    
    # Convert group ID to appropriate type
    if args.catalog == 'cw-all':
        try:
            group_id = int(args.group_id)
        except ValueError:
            print(f"Error: group-id must be an integer for cw-all catalog")
            sys.exit(1)
    else:
        group_id = args.group_id
    
    print("\n" + "="*60)
    print("INTERACTIVE SPECTROSCOPIC MEMBERSHIP DASHBOARD")
    print("="*60)
    print(f"\nCatalog: {args.catalog}")
    print(f"Starting Group ID: {group_id}")
    print(f"Search Radius: {args.radius} kpc")
    print("\nControls:")
    print("  - Adjust sliders to change parameters")
    print("  - Click 'Rerun VRF' to recompute with new parameters")
    print("  - Click 'Accept & Next' when satisfied")
    print("  - Click 'Previous' to go back")
    print("  - Click 'Save Results' to save reviewed groups")
    print("  - Enter index in textbox to jump to specific group")
    print("="*60 + "\n")
    
    # Create and show dashboard
    dashboard = InteractiveSpeczDashboard(args.catalog, group_id, args.radius)
    dashboard.show()


if __name__ == '__main__':
    main()
