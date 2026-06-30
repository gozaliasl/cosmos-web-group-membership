#!/usr/bin/env python
"""
Publication-Quality Analysis Plots

Creates high-quality figures for publication from membership determination results:
- Member statistics (N vs z, histograms)
- Velocity dispersion distributions
- VRF performance metrics
- Mass-richness relations
- Comparison plots (VRF vs gapper, spec-z vs photo-z)

Features:
- 300 DPI PDF output
- Consistent styling
- Publication-ready formatting
- Multi-panel layouts

Usage:
    python create_publication_plots.py --catalog all
    python create_publication_plots.py --catalog both --output-dir figures/

Requirements:
    conda activate astro-clean
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import argparse
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# Publication settings
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.framealpha'] = 0.9

# Paths
BASE_DIR = Path('/Users/gozalig1/Projects/compact-groups-xray-analysis/membership_determination')
RESULTS_DIR = BASE_DIR / 'results'
OUTPUT_DIR = BASE_DIR / 'figures' / 'publication'

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_latest_results(catalog_name):
    """Load most recent results for a catalog."""
    # Combined results
    combined_dir = RESULTS_DIR / 'combined'
    combined_files = list(combined_dir.glob(f'combined_membership_{catalog_name}_*.csv'))
    
    if len(combined_files) == 0:
        print(f"  Warning: No combined results found for {catalog_name}")
        return None
    
    latest_file = max(combined_files, key=lambda p: p.stat().st_mtime)
    print(f"  Loading: {latest_file.name}")
    
    return pd.read_csv(latest_file)


def plot_membership_statistics(df, catalog_name, output_dir):
    """
    Create multi-panel figure showing membership statistics.
    
    Panels:
    1. N_members vs redshift (spec-z and photo-z)
    2. Histograms of N_members
    3. Rejection rate vs redshift
    4. N_specz vs N_photoz comparison
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Membership Statistics: {catalog_name.upper()}', 
                 fontsize=14, fontweight='bold')
    
    # Panel 1: N_members vs redshift
    ax = axes[0, 0]
    
    if 'n_members_specz' in df.columns and 'z_specz' in df.columns:
        ax.scatter(df['z_specz'], df['n_members_specz'], 
                  alpha=0.5, s=30, c='blue', label='Spec-z', edgecolors='darkblue')
    
    if 'n_members_photoz' in df.columns and 'z_photoz' in df.columns:
        ax.scatter(df['z_photoz'], df['n_members_photoz'], 
                  alpha=0.5, s=30, c='red', marker='s', label='Photo-z', edgecolors='darkred')
    
    ax.set_xlabel('Redshift', fontsize=12)
    ax.set_ylabel('Number of Members', fontsize=12)
    ax.set_title('Members vs Redshift', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Histograms
    ax = axes[0, 1]
    
    if 'n_members_specz' in df.columns:
        specz_members = df['n_members_specz'][df['n_members_specz'] > 0]
        ax.hist(specz_members, bins=20, alpha=0.6, color='blue', 
               label=f'Spec-z (N={len(specz_members)})', edgecolor='darkblue')
    
    if 'n_members_photoz' in df.columns:
        photoz_members = df['n_members_photoz'][df['n_members_photoz'] > 0]
        ax.hist(photoz_members, bins=20, alpha=0.6, color='red',
               label=f'Photo-z (N={len(photoz_members)})', edgecolor='darkred')
    
    ax.set_xlabel('Number of Members', fontsize=12)
    ax.set_ylabel('Number of Groups', fontsize=12)
    ax.set_title('Member Distribution', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 3: Rejection rate vs redshift
    ax = axes[1, 0]
    
    if 'n_candidates_specz' in df.columns and 'n_members_specz' in df.columns:
        mask = df['n_candidates_specz'] > 0
        rejection_rate = 100 * (1 - df['n_members_specz'] / df['n_candidates_specz'])
        
        ax.scatter(df['z_specz'][mask], rejection_rate[mask],
                  alpha=0.5, s=30, c='blue', edgecolors='darkblue')
        
        # Running median
        z_bins = np.linspace(df['z_specz'].min(), df['z_specz'].max(), 10)
        z_centers = (z_bins[:-1] + z_bins[1:]) / 2
        median_rejection = []
        
        for i in range(len(z_bins)-1):
            in_bin = (df['z_specz'] >= z_bins[i]) & (df['z_specz'] < z_bins[i+1]) & mask
            if np.sum(in_bin) > 0:
                median_rejection.append(np.median(rejection_rate[in_bin]))
            else:
                median_rejection.append(np.nan)
        
        ax.plot(z_centers, median_rejection, 'r-', linewidth=2, label='Running Median')
    
    ax.set_xlabel('Redshift', fontsize=12)
    ax.set_ylabel('VRF Rejection Rate (%)', fontsize=12)
    ax.set_title('Rejection Rate vs Redshift', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 100)
    
    # Panel 4: Spec-z vs Photo-z comparison
    ax = axes[1, 1]
    
    if 'n_members_specz' in df.columns and 'n_members_photoz' in df.columns:
        mask = (df['n_members_specz'] > 0) | (df['n_members_photoz'] > 0)
        
        ax.scatter(df['n_members_specz'][mask], df['n_members_photoz'][mask],
                  alpha=0.5, s=30, c='purple', edgecolors='darkviolet')
        
        # 1:1 line
        max_val = max(df['n_members_specz'].max(), df['n_members_photoz'].max())
        ax.plot([0, max_val], [0, max_val], 'k--', linewidth=2, alpha=0.5, label='1:1')
        
        # Calculate correlation
        corr = np.corrcoef(df['n_members_specz'][mask], df['n_members_photoz'][mask])[0, 1]
        ax.text(0.05, 0.95, f'r = {corr:.2f}', transform=ax.transAxes,
               fontsize=11, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax.set_xlabel('N$_{members}$ (Spec-z)', fontsize=12)
    ax.set_ylabel('N$_{members}$ (Photo-z)', fontsize=12)
    ax.set_title('Spec-z vs Photo-z Members', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / f'membership_statistics_{catalog_name}.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file.name}")
    
    plt.close()


def plot_velocity_dispersion(df, catalog_name, output_dir):
    """
    Create figure showing velocity dispersion properties.
    
    Panels:
    1. σ_v distribution (histogram)
    2. σ_v vs redshift
    3. σ_v vs N_members
    4. σ_v vs richness (if available)
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Velocity Dispersion Analysis: {catalog_name.upper()}',
                 fontsize=14, fontweight='bold')
    
    if 'sigma_v' not in df.columns:
        print(f"  Warning: No sigma_v column for {catalog_name}")
        plt.close()
        return
    
    # Filter valid values
    valid = (df['sigma_v'] > 0) & (df['sigma_v'] < 5000)
    sigma_v = df['sigma_v'][valid]
    
    if len(sigma_v) == 0:
        print(f"  Warning: No valid sigma_v values for {catalog_name}")
        plt.close()
        return
    
    # Panel 1: Histogram
    ax = axes[0, 0]
    
    ax.hist(sigma_v, bins=30, alpha=0.7, color='steelblue', edgecolor='darkblue')
    
    # Statistics
    median_sigma = np.median(sigma_v)
    mean_sigma = np.mean(sigma_v)
    
    ax.axvline(median_sigma, color='red', linestyle='--', linewidth=2, label=f'Median: {median_sigma:.0f} km/s')
    ax.axvline(mean_sigma, color='orange', linestyle='--', linewidth=2, label=f'Mean: {mean_sigma:.0f} km/s')
    
    ax.set_xlabel('σ$_v$ (km/s)', fontsize=12)
    ax.set_ylabel('Number of Groups', fontsize=12)
    ax.set_title('Velocity Dispersion Distribution', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 2: σ_v vs redshift
    ax = axes[0, 1]
    
    if 'z_specz' in df.columns:
        ax.scatter(df['z_specz'][valid], sigma_v, alpha=0.5, s=30, 
                  c='steelblue', edgecolors='darkblue')
        
        # Running median
        z_bins = np.linspace(df['z_specz'][valid].min(), df['z_specz'][valid].max(), 10)
        z_centers = (z_bins[:-1] + z_bins[1:]) / 2
        median_sigma = []
        
        for i in range(len(z_bins)-1):
            in_bin = (df['z_specz'][valid] >= z_bins[i]) & (df['z_specz'][valid] < z_bins[i+1])
            if np.sum(in_bin) > 0:
                median_sigma.append(np.median(sigma_v[in_bin.values]))
            else:
                median_sigma.append(np.nan)
        
        ax.plot(z_centers, median_sigma, 'r-', linewidth=2, label='Running Median')
    
    ax.set_xlabel('Redshift', fontsize=12)
    ax.set_ylabel('σ$_v$ (km/s)', fontsize=12)
    ax.set_title('σ$_v$ vs Redshift', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel 3: σ_v vs N_members
    ax = axes[1, 0]
    
    if 'n_members_specz' in df.columns:
        n_members = df['n_members_specz'][valid]
        
        ax.scatter(n_members, sigma_v, alpha=0.5, s=30,
                  c='steelblue', edgecolors='darkblue')
        
        # Binned median
        n_bins = np.arange(5, max(n_members)+5, 5)
        n_centers = (n_bins[:-1] + n_bins[1:]) / 2
        median_sigma = []
        
        for i in range(len(n_bins)-1):
            in_bin = (n_members >= n_bins[i]) & (n_members < n_bins[i+1])
            if np.sum(in_bin) > 0:
                median_sigma.append(np.median(sigma_v[in_bin.values]))
            else:
                median_sigma.append(np.nan)
        
        ax.plot(n_centers, median_sigma, 'r-', linewidth=2, label='Binned Median')
    
    ax.set_xlabel('N$_{members}$ (Spec-z)', fontsize=12)
    ax.set_ylabel('σ$_v$ (km/s)', fontsize=12)
    ax.set_title('σ$_v$ vs Richness', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel 4: Quality distribution
    ax = axes[1, 1]
    
    if 'vrf_quality' in df.columns:
        quality_counts = df['vrf_quality'][valid].value_counts().sort_index()
        
        colors = ['green', 'yellowgreen', 'orange', 'red']
        labels = ['Excellent (N≥15)', 'Good (8≤N<15)', 'Poor (5≤N<8)', 'Insufficient (N<5)']
        
        bars = ax.bar(quality_counts.index, quality_counts.values,
                     color=[colors[int(q)-1] if q <= 4 else 'gray' for q in quality_counts.index],
                     edgecolor='black', alpha=0.7)
        
        # Add counts on bars
        for bar, count in zip(bars, quality_counts.values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(count)}', ha='center', va='bottom', fontsize=10)
        
        ax.set_xlabel('VRF Quality Flag', fontsize=12)
        ax.set_ylabel('Number of Groups', fontsize=12)
        ax.set_title('VRF Quality Distribution', fontsize=12, fontweight='bold')
        ax.set_xticks(range(1, 5))
        ax.set_xticklabels(['1\nExcellent', '2\nGood', '3\nPoor', '4\nInsuff.'])
        ax.grid(True, alpha=0.3, axis='y')
    else:
        ax.text(0.5, 0.5, 'Quality data not available', transform=ax.transAxes,
               ha='center', va='center', fontsize=12)
        ax.axis('off')
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / f'velocity_dispersion_{catalog_name}.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file.name}")
    
    plt.close()


def plot_mass_richness(df, catalog_name, output_dir):
    """
    Create figure showing halo mass vs richness relations.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'Mass-Richness Relations: {catalog_name.upper()}',
                 fontsize=14, fontweight='bold')
    
    if 'M_200' not in df.columns:
        print(f"  Warning: No M_200 column for {catalog_name}")
        plt.close()
        return
    
    # Filter valid masses
    valid = (df['M_200'] > 0) & (df['M_200'] < 1e16)
    
    if np.sum(valid) == 0:
        print(f"  Warning: No valid M_200 values for {catalog_name}")
        plt.close()
        return
    
    log_mass = np.log10(df['M_200'][valid])
    
    # Panel 1: M_200 vs N_members
    ax = axes[0]
    
    if 'n_members_specz' in df.columns:
        n_members = df['n_members_specz'][valid]
        
        ax.scatter(n_members, log_mass, alpha=0.5, s=30,
                  c='steelblue', edgecolors='darkblue')
        
        # Fit power law
        coeffs = np.polyfit(np.log10(n_members[n_members > 0]), 
                           log_mass[n_members > 0], 1)
        
        n_fit = np.linspace(n_members.min(), n_members.max(), 100)
        log_mass_fit = coeffs[0] * np.log10(n_fit) + coeffs[1]
        
        ax.plot(n_fit, log_mass_fit, 'r--', linewidth=2,
               label=f'log(M) = {coeffs[0]:.2f}·log(N) + {coeffs[1]:.2f}')
    
    ax.set_xlabel('N$_{members}$ (Spec-z)', fontsize=12)
    ax.set_ylabel('log(M$_{200}$ / M$_\\odot$)', fontsize=12)
    ax.set_title('Halo Mass vs Richness', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Panel 2: M_200 vs σ_v
    ax = axes[1]
    
    if 'sigma_v' in df.columns:
        sigma_v_valid = df['sigma_v'][valid]
        mask = sigma_v_valid > 0
        
        ax.scatter(sigma_v_valid[mask], log_mass[mask], alpha=0.5, s=30,
                  c='steelblue', edgecolors='darkblue')
        
        # Theoretical M-σ relation: M ∝ σ^3
        sigma_range = np.linspace(sigma_v_valid[mask].min(), sigma_v_valid[mask].max(), 100)
        # Normalize to median
        median_sigma = np.median(sigma_v_valid[mask])
        median_mass = np.median(log_mass[mask])
        log_mass_theory = median_mass + 3 * np.log10(sigma_range / median_sigma)
        
        ax.plot(sigma_range, log_mass_theory, 'r--', linewidth=2,
               label='M $\\propto$ σ$^3$')
    
    ax.set_xlabel('σ$_v$ (km/s)', fontsize=12)
    ax.set_ylabel('log(M$_{200}$ / M$_\\odot$)', fontsize=12)
    ax.set_title('Halo Mass vs Velocity Dispersion', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / f'mass_richness_{catalog_name}.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file.name}")
    
    plt.close()


def create_all_plots(catalog_name, output_dir):
    """Create all publication plots for a catalog."""
    print(f"\n{'='*60}")
    print(f"Creating publication plots: {catalog_name.upper()}")
    print(f"{'='*60}")
    
    # Load data
    df = load_latest_results(catalog_name)
    
    if df is None:
        print(f"  Skipping {catalog_name} - no data")
        return
    
    print(f"  Groups: {len(df)}")
    
    # Create plots
    print("\n  Generating figures...")
    
    plot_membership_statistics(df, catalog_name, output_dir)
    plot_velocity_dispersion(df, catalog_name, output_dir)
    plot_mass_richness(df, catalog_name, output_dir)
    
    print(f"\n  ✓ All plots created for {catalog_name}")


def main():
    parser = argparse.ArgumentParser(
        description='Create publication-quality analysis plots',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog', choices=['all', 'hcg', 'both'], required=True,
                       help='Which catalog to plot')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory (default: figures/publication/)')
    
    args = parser.parse_args()
    
    # Determine output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = OUTPUT_DIR
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*60)
    print("PUBLICATION PLOT GENERATION")
    print("="*60)
    print(f"\nConfiguration:")
    print(f"  Catalog: {args.catalog}")
    print(f"  Output directory: {output_dir}")
    print(f"  Format: PDF (300 DPI)")
    print("="*60)
    
    # Create plots
    if args.catalog in ['all', 'both']:
        create_all_plots('cw-all', output_dir)
    
    if args.catalog in ['hcg', 'both']:
        create_all_plots('cw-hcg', output_dir)
    
    print("\n" + "="*60)
    print("PLOT GENERATION COMPLETE")
    print("="*60)
    print(f"  Output directory: {output_dir}")
    print("="*60 + "\n")


if __name__ == '__main__':
    main()
