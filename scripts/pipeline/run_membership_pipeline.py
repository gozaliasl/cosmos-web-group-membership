#!/usr/bin/env python
"""
Membership Pipeline Orchestrator

Manages the complete workflow for galaxy group membership determination:
1. Spectroscopic membership (VRF method)
2. Photometric membership (probabilistic)
3. Combined membership catalog
4. Optional interactive review

Features:
- Automated batch processing
- Optional interactive QA/QC after automated run
- Progress tracking and error handling
- Combines spec-z and photo-z results
- Generates summary reports

Usage:
    # Run both spec-z and photo-z
    python run_membership_pipeline.py --catalog both --method vrf
    
    # Run with interactive review
    python run_membership_pipeline.py --catalog all --interactive
    
    # Test mode
    python run_membership_pipeline.py --catalog both --test

Requirements:
    conda activate astro-clean
"""

import subprocess
import sys
import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BASE_DIR = PROJECT_ROOT / 'membership_determination'
RESULTS_DIR = BASE_DIR / 'results'


class MembershipPipeline:
    """Orchestrator for membership determination pipeline."""
    
    def __init__(self, catalog, method='vrf', radius=500, 
                 max_dz_norm=0.01, max_velocity=2000,
                 max_dist_phot_z_error=2.0, prob_threshold=0.5,
                 test_mode=False, interactive=False,
                 specz_catalog=None):
        """
        Initialize pipeline.
        
        Parameters:
        -----------
        catalog : str
            'all', 'hcg', or 'both'
        method : str
            'vrf' or 'gapper' for spec-z
        radius : float
            Search radius in kpc
        max_dz_norm : float
            Normalized redshift cut for spec-z
        max_velocity : float
            Velocity cut for spec-z (km/s)
        max_dist_phot_z_error : float
            Photo-z threshold (σ units)
        prob_threshold : float
            P_member threshold for photo-z
        test_mode : bool
            Process only 10 groups per catalog
        interactive : bool
            Launch interactive dashboard after automated run
        """
        self.catalog = catalog
        self.method = method
        self.radius = radius
        self.max_dz_norm = max_dz_norm
        self.max_velocity = max_velocity
        self.max_dist_phot_z_error = max_dist_phot_z_error
        self.prob_threshold = prob_threshold
        self.test_mode = test_mode
        self.interactive = interactive
        self.specz_catalog = specz_catalog
        
        self.results = {
            'specz': {'all': None, 'hcg': None},
            'photoz': {'all': None, 'hcg': None}
        }
        
        self.start_time = datetime.now()
    
    def run(self):
        """Execute complete pipeline."""
        print("\n" + "="*80)
        print("MEMBERSHIP DETERMINATION PIPELINE")
        print("="*80)
        print(f"\nConfiguration:")
        print(f"  Catalog: {self.catalog}")
        print(f"  Method: {self.method}")
        print(f"  Radius: {self.radius} kpc")
        print(f"  Spec-z cuts: max_dz_norm={self.max_dz_norm}, max_velocity={self.max_velocity} km/s")
        print(f"  Photo-z cuts: {self.max_dist_phot_z_error}σ NMAD, P_member>={self.prob_threshold}")
        print(f"  Test mode: {self.test_mode}")
        print(f"  Interactive: {self.interactive}")
        print(f"  Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*80 + "\n")
        
        # Step 1: Spectroscopic membership
        print("\n" + "="*80)
        print("STEP 1: SPECTROSCOPIC MEMBERSHIP (VRF)")
        print("="*80)
        self.run_specz_membership()
        
        # Step 2: Photometric membership
        print("\n" + "="*80)
        print("STEP 2: PHOTOMETRIC MEMBERSHIP")
        print("="*80)
        self.run_photoz_membership()
        
        # Step 3: Combine results
        print("\n" + "="*80)
        print("STEP 3: COMBINE RESULTS")
        print("="*80)
        self.combine_results()
        
        # Step 4: Generate summary
        print("\n" + "="*80)
        print("STEP 4: GENERATE SUMMARY")
        print("="*80)
        self.generate_summary()
        
        # Step 5: Interactive review (optional)
        if self.interactive:
            print("\n" + "="*80)
            print("STEP 5: INTERACTIVE REVIEW")
            print("="*80)
            self.launch_interactive_dashboard()
        
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds()
        
        print("\n" + "="*80)
        print("PIPELINE COMPLETE")
        print("="*80)
        print(f"  Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Duration: {duration/60:.1f} minutes")
        print(f"  Results directory: {RESULTS_DIR}")
        print("="*80 + "\n")
    
    def run_specz_membership(self):
        """Run spectroscopic membership determination."""
        script_path = SCRIPT_DIR / 'determine_specz_membership.py'
        
        if not script_path.exists():
            print(f"  ERROR: Script not found: {script_path}")
            return
        
        cmd = [
            'python', str(script_path),
            '--catalog', self.catalog,
            '--method', self.method,
            '--radius', str(self.radius),
            '--max-dz-norm', str(self.max_dz_norm),
            '--max-velocity', str(self.max_velocity)
        ]
        
        if self.test_mode:
            cmd.append('--test')
        
        if self.specz_catalog:
            cmd.extend(['--specz-catalog', self.specz_catalog])
        
        print(f"\nRunning: {' '.join(cmd)}\n")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            print("\n  ✓ Spec-z membership complete")
        except subprocess.CalledProcessError as e:
            print(f"\n  ✗ Spec-z membership failed: {e}")
            sys.exit(1)
    
    def run_photoz_membership(self):
        """Run photometric membership determination."""
        script_path = SCRIPT_DIR / 'determine_photoz_membership.py'
        
        if not script_path.exists():
            print(f"  ERROR: Script not found: {script_path}")
            return
        
        cmd = [
            'python', str(script_path),
            '--catalog', self.catalog,
            '--radius', str(self.radius),
            '--max-dist-error', str(self.max_dist_phot_z_error),
            '--prob-threshold', str(self.prob_threshold)
        ]
        
        if self.test_mode:
            cmd.append('--test')
        
        print(f"\nRunning: {' '.join(cmd)}\n")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            print("\n  ✓ Photo-z membership complete")
        except subprocess.CalledProcessError as e:
            print(f"\n  ✗ Photo-z membership failed: {e}")
            sys.exit(1)
    
    def combine_results(self):
        """Combine spec-z and photo-z results into unified catalog."""
        print("\nCombining spec-z and photo-z results...")
        
        combined_dir = RESULTS_DIR / 'combined'
        combined_dir.mkdir(parents=True, exist_ok=True)
        
        # Get latest result files
        catalogs = []
        if self.catalog in ['all', 'both']:
            catalogs.append('cw-all')
        if self.catalog in ['hcg', 'both']:
            catalogs.append('cw-hcg')
        
        for cat in catalogs:
            # Find latest spec-z results
            specz_dir = RESULTS_DIR / 'specz'
            specz_files = list(specz_dir.glob(f'specz_membership_{cat}_*.csv'))
            
            # Find latest photo-z results
            photoz_dir = RESULTS_DIR / 'photoz'
            photoz_files = list(photoz_dir.glob(f'photoz_membership_{cat}_*.csv'))
            
            if len(specz_files) == 0 or len(photoz_files) == 0:
                print(f"  WARNING: Missing results for {cat}")
                continue
            
            # Get most recent files
            specz_file = max(specz_files, key=lambda p: p.stat().st_mtime)
            photoz_file = max(photoz_files, key=lambda p: p.stat().st_mtime)
            
            print(f"\n  Combining {cat}:")
            print(f"    Spec-z: {specz_file.name}")
            print(f"    Photo-z: {photoz_file.name}")
            
            # Load results
            specz_df = pd.read_csv(specz_file)
            photoz_df = pd.read_csv(photoz_file)
            
            # Merge on group_id
            combined_df = specz_df.merge(
                photoz_df, 
                on='group_id', 
                how='outer',
                suffixes=('_specz', '_photoz')
            )
            
            # Save combined catalog
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_file = combined_dir / f'combined_membership_{cat}_{timestamp}.csv'
            combined_df.to_csv(output_file, index=False)
            
            print(f"    Saved: {output_file.name}")
            print(f"    Groups: {len(combined_df)}")
        
        print("\n  ✓ Results combined")
    
    def generate_summary(self):
        """Generate summary statistics and report."""
        print("\nGenerating summary report...")
        
        summary_lines = []
        summary_lines.append("="*80)
        summary_lines.append("MEMBERSHIP DETERMINATION SUMMARY REPORT")
        summary_lines.append("="*80)
        summary_lines.append(f"\nDate: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        summary_lines.append(f"Catalog: {self.catalog}")
        summary_lines.append(f"Test mode: {self.test_mode}")
        summary_lines.append("\nPARAMETERS:")
        summary_lines.append(f"  Spec-z method: {self.method}")
        summary_lines.append(f"  Search radius: {self.radius} kpc")
        summary_lines.append(f"  Spec-z cuts: max_dz_norm={self.max_dz_norm}, max_velocity={self.max_velocity} km/s")
        summary_lines.append(f"  Photo-z cuts: {self.max_dist_phot_z_error}σ NMAD, P_member>={self.prob_threshold}")
        
        # Analyze combined results
        combined_dir = RESULTS_DIR / 'combined'
        combined_files = list(combined_dir.glob('combined_membership_*.csv'))
        
        if len(combined_files) > 0:
            summary_lines.append("\nRESULTS:")
            
            for cat_name in ['cw-all', 'cw-hcg']:
                cat_files = [f for f in combined_files if cat_name in f.name]
                if len(cat_files) == 0:
                    continue
                
                # Get most recent
                latest_file = max(cat_files, key=lambda p: p.stat().st_mtime)
                df = pd.read_csv(latest_file)
                
                summary_lines.append(f"\n  {cat_name.upper()}:")
                summary_lines.append(f"    Total groups: {len(df)}")
                
                if 'n_candidates_specz' in df.columns:
                    summary_lines.append(f"    Groups with spec-z candidates: {np.sum(df['n_candidates_specz'] > 0)}")
                    summary_lines.append(f"    Mean spec-z members: {df['n_members_specz'].mean():.1f}")
                
                if 'n_candidates_photoz' in df.columns:
                    summary_lines.append(f"    Groups with photo-z candidates: {np.sum(df['n_candidates_photoz'] > 0)}")
                    summary_lines.append(f"    Mean photo-z members: {df['n_members_photoz'].mean():.1f}")
                
                if 'sigma_v' in df.columns:
                    valid_sigma = df['sigma_v'][df['sigma_v'] > 0]
                    if len(valid_sigma) > 0:
                        summary_lines.append(f"    Mean σ_v: {valid_sigma.mean():.0f} km/s (N={len(valid_sigma)})")
        
        summary_lines.append("\n" + "="*80)
        
        summary_text = "\n".join(summary_lines)
        print(summary_text)
        
        # Save summary
        summary_file = RESULTS_DIR / f'pipeline_summary_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt'
        with open(summary_file, 'w') as f:
            f.write(summary_text)
        
        print(f"\n  ✓ Summary saved: {summary_file}")
    
    def launch_interactive_dashboard(self):
        """Launch interactive dashboard for manual review."""
        script_path = SCRIPT_DIR / 'interactive_specz_dashboard.py'
        
        if not script_path.exists():
            print(f"  ERROR: Interactive dashboard script not found: {script_path}")
            return
        
        # Determine which catalog to review
        if self.catalog == 'both':
            review_catalog = 'cw-all'  # Default to all
            print(f"  Note: Multiple catalogs available. Launching dashboard for: {review_catalog}")
        else:
            review_catalog = 'cw-all' if self.catalog == 'all' else 'cw-hcg'
        
        cmd = [
            'python', str(script_path),
            '--catalog', review_catalog,
            '--group-id', '0' if review_catalog == 'cw-all' else 'Py12_1',
            '--radius', str(self.radius)
        ]
        
        print(f"\nLaunching interactive dashboard...")
        print(f"  Command: {' '.join(cmd)}")
        print(f"  Catalog: {review_catalog}")
        print(f"\n  Instructions:")
        print(f"    - Review membership visually")
        print(f"    - Adjust parameters with sliders")
        print(f"    - Click 'Accept & Next' to save and move on")
        print(f"    - Click 'Save Results' when finished\n")
        
        try:
            # Launch in separate process (non-blocking)
            subprocess.Popen(cmd)
            print("  ✓ Interactive dashboard launched")
        except Exception as e:
            print(f"  ✗ Failed to launch dashboard: {e}")


def main():
    parser = argparse.ArgumentParser(
        description='Run complete membership determination pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--catalog', choices=['all', 'hcg', 'both'], required=True,
                       help='Which catalog(s) to process')
    parser.add_argument('--method', choices=['vrf', 'gapper'], default='vrf',
                       help='Spec-z membership method')
    parser.add_argument('--radius', type=float, default=500,
                       help='Search radius in kpc')
    parser.add_argument('--max-dz-norm', type=float, default=0.01,
                       help='Normalized redshift cut for spec-z')
    parser.add_argument('--max-velocity', type=float, default=2000,
                       help='Velocity cut for spec-z (km/s)')
    parser.add_argument('--max-dist-error', type=float, default=2.0,
                       help='Photo-z threshold (σ NMAD units)')
    parser.add_argument('--prob-threshold', type=float, default=0.5,
                       help='P_member threshold for photo-z')
    parser.add_argument('--test', action='store_true',
                       help='Test mode: process only 10 groups')
    parser.add_argument('--interactive', action='store_true',
                       help='Launch interactive dashboard after automated run')
    parser.add_argument(
        '--specz-catalog',
        type=str,
        default=None,
        help=(
            'Forwarded to determine_specz_membership.py. Default: auto (Webb spec-z). '
            'Pass COMMAGN path only for the small subset.'
        ),
    )
    
    args = parser.parse_args()
    
    # Create and run pipeline
    pipeline = MembershipPipeline(
        catalog=args.catalog,
        method=args.method,
        radius=args.radius,
        max_dz_norm=args.max_dz_norm,
        max_velocity=args.max_velocity,
        max_dist_phot_z_error=args.max_dist_error,
        prob_threshold=args.prob_threshold,
        test_mode=args.test,
        interactive=args.interactive,
        specz_catalog=args.specz_catalog,
    )
    
    pipeline.run()


if __name__ == '__main__':
    main()
