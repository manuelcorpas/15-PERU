#!/usr/bin/env python3
"""
ADMIXTURE Analysis Pipeline - Modified to preserve all samples
Never drops individuals during preparation (maintains all 736 samples)
Expects input from: ANALYSIS/00-06-IBD/paper936k/merged_936k_final.*
"""

import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

def check_dependencies():
    """Check if required tools are available"""
    tools = {
        'plink': 'PLINK is required for LD pruning',
        'admixture': 'ADMIXTURE is required for ancestry analysis'
    }
    
    missing = []
    for tool, msg in tools.items():
        try:
            subprocess.run([tool, '--help'], capture_output=True, check=False)
        except FileNotFoundError:
            missing.append(f"  - {msg}")
    
    if missing:
        print("ERROR: Missing required tools:")
        for msg in missing:
            print(msg)
        print("\nPlease install missing tools and ensure they are in your PATH")
        sys.exit(1)

def count_samples(fam_file):
    """Count number of samples in a .fam file"""
    with open(fam_file, 'r') as f:
        return sum(1 for line in f)

def prepare_for_admixture(plink_path, input_prefix, output_dir, 
                         window_kb=200, window_step=25, r2_threshold=0.5,
                         maf_threshold=0.01, geno_threshold=0.05):
    """
    Prepare data for ADMIXTURE: LD pruning while preserving all samples
    
    Parameters:
    -----------
    window_kb : int, default 200
        Window size in kb for LD pruning
    window_step : int, default 25 
        Step size (number of variants) for LD pruning
    r2_threshold : float, default 0.5
        r^2 threshold for LD pruning
    maf_threshold : float, default 0.01
        Minor allele frequency threshold (can be adjusted)
    geno_threshold : float, default 0.05
        Maximum missing rate per SNP (can be adjusted)
    """
    output_prefix = f'{output_dir}/merged_936k_admixture'
    
    print("\n=== Preparing data for ADMIXTURE ===")
    print(f"Input: {input_prefix}")
    print(f"LD pruning parameters: window={window_kb}kb, step={window_step}, r²<{r2_threshold}")
    print(f"Variant filters: MAF>{maf_threshold}, missing<{geno_threshold}")
    
    # Check input files exist
    for ext in ['.bed', '.bim', '.fam']:
        if not os.path.exists(f'{input_prefix}{ext}'):
            raise FileNotFoundError(f"Required input file not found: {input_prefix}{ext}")
    
    # Count initial samples and SNPs
    n_input_samples = count_samples(f'{input_prefix}.fam')
    with open(f'{input_prefix}.bim', 'r') as f:
        initial_snps = sum(1 for line in f)
    print(f"\nInitial dataset: {n_input_samples} samples, {initial_snps:,} SNPs")
    
    # Check if pruned files already exist and have correct sample count
    pruned_bed = f'{output_prefix}_pruned.bed'
    pruned_fam = f'{output_prefix}_pruned.fam'
    prune_in = f'{output_prefix}.prune.in'
    
    if os.path.exists(pruned_bed) and os.path.exists(pruned_fam):
        existing_samples = count_samples(pruned_fam)
        if existing_samples == n_input_samples:
            print(f"\nExisting pruned files found with correct sample count ({existing_samples})")
            print("Reusing existing pruned dataset...")
            return f'{output_prefix}_pruned'
        else:
            print(f"\nExisting pruned files have wrong sample count ({existing_samples} vs {n_input_samples})")
            print("Rebuilding pruned dataset...")
    
    # Step 1: LD pruning to identify independent SNPs
    prune_cmd = [
        plink_path,
        '--bfile', input_prefix,
        '--indep-pairwise', str(window_kb), str(window_step), str(r2_threshold),
        '--out', output_prefix,
        '--threads', str(mp.cpu_count())
    ]
    
    print(f"\nStep 1: Identifying LD-independent SNPs...")
    result = subprocess.run(prune_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in LD pruning: {result.stderr}")
        sys.exit(1)
    
    # Count pruned SNPs
    with open(prune_in, 'r') as f:
        pruned_count = sum(1 for line in f)
    print(f"  SNPs passing LD filter: {pruned_count:,}")
    
    # Step 2: Extract pruned SNPs WITH variant filters
    print(f"\nStep 2: Extracting pruned SNPs with variant filters...")
    extract_cmd = [
        plink_path,
        '--bfile', input_prefix,
        '--extract', prune_in,
        '--make-bed',
        '--maf', str(maf_threshold),
        '--geno', str(geno_threshold),
        # NO --mind flag here!
        '--out', f'{output_prefix}_pruned',
        '--threads', str(mp.cpu_count())
    ]
    
    result = subprocess.run(extract_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in extraction: {result.stderr}")
        sys.exit(1)
    
    # Verify sample count is preserved
    n_output_samples = count_samples(f'{output_prefix}_pruned.fam')
    with open(f'{output_prefix}_pruned.bim', 'r') as f:
        final_snps = sum(1 for line in f)
    
    print(f"  After variant filters: {n_output_samples} samples, {final_snps:,} SNPs")
    
    # If samples were lost (shouldn't happen without --mind), rebuild without variant filters
    if n_output_samples != n_input_samples:
        print(f"\nWARNING: Sample count mismatch ({n_output_samples} != {n_input_samples})")
        print("Rebuilding without MAF/GENO filters to preserve all samples...")
        
        # Fallback: Extract only LD-pruned SNPs without any other filters
        fallback_cmd = [
            plink_path,
            '--bfile', input_prefix,
            '--extract', prune_in,
            '--make-bed',
            # No MAF, no GENO, no MIND
            '--out', f'{output_prefix}_pruned',
            '--threads', str(mp.cpu_count())
        ]
        
        result = subprocess.run(fallback_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error in fallback extraction: {result.stderr}")
            sys.exit(1)
        
        # Re-verify sample count
        n_output_samples = count_samples(f'{output_prefix}_pruned.fam')
        with open(f'{output_prefix}_pruned.bim', 'r') as f:
            final_snps = sum(1 for line in f)
        
        if n_output_samples != n_input_samples:
            print(f"\nERROR: Unable to preserve all samples!")
            print(f"  Input samples: {n_input_samples}")
            print(f"  Output samples: {n_output_samples}")
            print("This should not happen. Please check your input files.")
            sys.exit(1)
        
        print(f"  Fallback successful: {n_output_samples} samples, {final_snps:,} SNPs")
    
    # Final assertion
    assert n_output_samples == n_input_samples, \
        f"Sample count mismatch: {n_output_samples} != {n_input_samples}"
    
    print(f"\n✓ Sample count preserved: {n_output_samples} samples")
    print(f"✓ Final SNP count: {final_snps:,}")
    
    # Check we have enough SNPs
    if final_snps < 10000:
        print(f"WARNING: Only {final_snps} SNPs remain after pruning. Consider adjusting parameters.")
    
    return f'{output_prefix}_pruned'

def run_single_admixture(K, bed_basename, output_dir, use_cv=False, cv_folds=5):
    """
    Run ADMIXTURE for a single K value
    """
    print(f"  Starting K={K}...")
    
    # Save current directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        if use_cv:
            cmd = [
                'admixture',
                '--cv=' + str(cv_folds),
                '-j' + str(mp.cpu_count()),
                '-s', str(np.random.randint(1000)),
                f'{bed_basename}.bed',
                str(K)
            ]
        else:
            cmd = [
                'admixture',
                '-j' + str(mp.cpu_count()),
                '-s', str(np.random.randint(1000)),
                f'{bed_basename}.bed',
                str(K)
            ]
        
        log_file = f'admixture_K{K}.log'
        with open(log_file, 'w') as log:
            result = subprocess.run(cmd, capture_output=True, text=True)
            log.write(result.stdout)
            log.write(result.stderr)
            
            if result.returncode != 0:
                print(f"  ERROR in K={K}: Check {log_file}")
                return K, None
        
        # Parse CV error if applicable
        cv_error = None
        if use_cv:
            for line in result.stdout.split('\n'):
                if 'CV error' in line:
                    try:
                        cv_error = float(line.split('=')[-1].strip())
                    except:
                        pass
                    break
        
        print(f"  Completed K={K}" + (f", CV error={cv_error:.6f}" if cv_error else ""))
        
    finally:
        os.chdir(original_dir)
    
    return K, cv_error

def run_admixture_parallel(bed_file_prefix, output_dir, K_range=(2, 8), use_cv=False, cv_folds=5):
    """
    Run ADMIXTURE in parallel for multiple K values
    """
    bed_basename = os.path.basename(bed_file_prefix)
    
    print(f"\nRunning ADMIXTURE for K={K_range[0]} to K={K_range[1]}...")
    cv_errors = {}
    
    # Use ProcessPoolExecutor for parallel execution
    max_workers = min(mp.cpu_count() // 2, K_range[1] - K_range[0] + 1)
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for K in range(K_range[0], K_range[1] + 1):
            future = executor.submit(
                run_single_admixture, 
                K, bed_basename, output_dir, use_cv, cv_folds
            )
            futures.append(future)
        
        # Collect results as they complete
        for future in as_completed(futures):
            K, cv_error = future.result()
            if cv_error is not None:
                cv_errors[K] = cv_error
    
    # Save CV errors if applicable
    if cv_errors:
        cv_file = f'{output_dir}/cv_errors.txt'
        with open(cv_file, 'w') as f:
            f.write("K\tCV_error\n")
            for k, error in sorted(cv_errors.items()):
                f.write(f"{k}\t{error:.6f}\n")
        print(f"\nCV errors saved to {cv_file}")
    
    return cv_errors

def quick_run_mode(bed_file_prefix, output_dir, K_values=[2, 3, 4, 5]):
    """
    Quick run for specific K values without CV
    """
    bed_basename = os.path.basename(bed_file_prefix)
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        for K in K_values:
            print(f"  Quick run for K={K}...")
            cmd = [
                'admixture',
                '-j' + str(mp.cpu_count()),
                '-s', str(np.random.randint(1000)),
                f'{bed_basename}.bed',
                str(K)
            ]
            
            log_file = f'admixture_K{K}_quick.log'
            with open(log_file, 'w') as log:
                result = subprocess.run(cmd, capture_output=True, text=True)
                log.write(result.stdout)
                log.write(result.stderr)
                
                if result.returncode != 0:
                    print(f"  ERROR in K={K}: Check {log_file}")
                else:
                    print(f"  K={K} complete")
    finally:
        os.chdir(original_dir)

def plot_cv_errors(cv_errors, output_dir):
    """
    Plot cross-validation errors to determine optimal K
    """
    if not cv_errors:
        return None
    
    K_values = sorted(cv_errors.keys())
    errors = [cv_errors[k] for k in K_values]
    
    plt.figure(figsize=(10, 6))
    plt.plot(K_values, errors, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('K (number of ancestral populations)', fontsize=12)
    plt.ylabel('Cross-validation error', fontsize=12)
    plt.title('ADMIXTURE Cross-Validation Error', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Mark the minimum
    min_k = min(cv_errors, key=cv_errors.get)
    min_error = cv_errors[min_k]
    plt.axvline(x=min_k, color='r', linestyle='--', alpha=0.5)
    plt.annotate(f'Optimal K={min_k}\nCV={min_error:.6f}', 
                xy=(min_k, min_error), 
                xytext=(min_k + 0.5, min_error + (max(errors) - min(errors)) * 0.05),
                fontsize=11, color='red',
                arrowprops=dict(arrowstyle='->', color='red', alpha=0.5))
    
    plt.tight_layout()
    plot_file = f'{output_dir}/cv_errors.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"CV error plot saved to {plot_file}")
    return min_k

def load_population_map(pop_file):
    """
    Load population mapping file with columns FID IID POP
    """
    if not os.path.exists(pop_file):
        return None
    
    try:
        # Try to read with different delimiters
        for sep in ['\t', ' ', ',']:
            try:
                pop_df = pd.read_csv(pop_file, sep=sep)
                # Normalize column names to uppercase
                pop_df.columns = [col.upper() for col in pop_df.columns]
                
                if all(col in pop_df.columns for col in ['FID', 'IID', 'POP']):
                    population_map = dict(zip(pop_df['IID'], pop_df['POP']))
                    return population_map
            except:
                continue
        
        # If no header, try as headerless file
        pop_df = pd.read_csv(pop_file, sep=None, header=None, engine='python')
        if len(pop_df.columns) >= 3:
            pop_df.columns = ['FID', 'IID', 'POP'] + list(pop_df.columns[3:])
            population_map = dict(zip(pop_df['IID'], pop_df['POP']))
            return population_map
            
    except Exception as e:
        print(f"  Warning: Could not load population file {pop_file}: {e}")
    
    return None

def plot_admixture_results(Q_file, fam_file, K, output_dir, population_map=None, cv_error=None):
    """
    Create ADMIXTURE bar plot
    """
    # Read Q matrix
    Q = np.loadtxt(Q_file)
    n_samples = Q.shape[0]
    
    # Read sample IDs from .fam file
    fam = pd.read_csv(fam_file, sep=' ', header=None, usecols=[0, 1])
    fam.columns = ['FID', 'IID']
    
    # Create DataFrame
    Q_df = pd.DataFrame(Q, columns=[f'Ancestry_{i+1}' for i in range(K)])
    Q_df['IID'] = fam['IID'].values
    Q_df['FID'] = fam['FID'].values
    
    # Add population labels if available
    n_with_pop = 0
    if population_map:
        Q_df['Population'] = Q_df['IID'].map(population_map)
        n_with_pop = Q_df['Population'].notna().sum()
        
        if n_with_pop > 0:
            # Sort by population and then by dominant ancestry
            Q_df['max_ancestry'] = Q_df[[f'Ancestry_{i+1}' for i in range(K)]].max(axis=1)
            Q_df['Population'] = Q_df['Population'].fillna('Unknown')
            Q_df = Q_df.sort_values(['Population', 'max_ancestry'], ascending=[True, False])
        else:
            # No valid population assignments
            Q_df = Q_df.sort_values('IID')
    else:
        # Sort by IID if no population map
        Q_df = Q_df.sort_values('IID')
    
    # Create plot
    fig, ax = plt.subplots(figsize=(20, 6))
    
    # Use distinct colors
    if K <= 8:
        colors = plt.cm.Set1(np.linspace(0, 0.9, K))
    else:
        colors = plt.cm.tab20(np.linspace(0, 1, K))
    
    # Plot stacked bars
    ind = np.arange(len(Q_df))
    bottom = np.zeros(len(Q_df))
    
    for k in range(K):
        values = Q_df[f'Ancestry_{k+1}'].values
        ax.bar(ind, values, bottom=bottom, width=1.0, 
               color=colors[k], label=f'Ancestry {k+1}', linewidth=0)
        bottom += values
    
    # Add population boundaries if we have populations
    if population_map and n_with_pop > 0:
        pop_counts = Q_df.groupby('Population').size()
        pop_boundaries = pop_counts.cumsum()
        pop_centers = pop_boundaries - pop_counts/2
        
        for boundary in pop_boundaries[:-1]:
            ax.axvline(x=boundary, color='white', linewidth=1)
        
        ax.set_xticks(pop_centers)
        ax.set_xticklabels(pop_counts.index, rotation=45, ha='right')
        
        print(f"  Samples with population assignment: {n_with_pop}/{n_samples}")
    else:
        ax.set_xticks([])
        ax.set_xlabel(f'Samples (n={n_samples})', fontsize=12)
    
    ax.set_ylabel('Ancestry Proportion', fontsize=12)
    
    # Title with CV if available
    title = f'ADMIXTURE Analysis (K={K}, n={n_samples}'
    if cv_error is not None:
        title += f', CV={cv_error:.6f}'
    title += ')'
    ax.set_title(title, fontsize=14)
    
    ax.set_ylim([0, 1])
    ax.set_xlim([0, len(Q_df)])
    
    # Place legend outside plot area
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', ncol=1)
    
    plt.tight_layout()
    plot_file = f'{output_dir}/admixture_K{K}.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Plot saved: {plot_file}")

def main():
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(
        description='Run ADMIXTURE analysis preserving all samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick run without CV (K=2,3,4,5)
  python %(prog)s --quick
  
  # Full run with cross-validation (K=2-8)
  python %(prog)s --cv
  
  # Custom K range with CV and population file
  python %(prog)s --cv --k-range 2 10 --pop-file ANALYSIS/00-08-PCA/pop_736.updated.tsv
        """
    )
    parser.add_argument('--quick', action='store_true', 
                       help='Quick mode: K=2,3,4,5 without CV')
    parser.add_argument('--cv', action='store_true', 
                       help='Run with cross-validation to find optimal K')
    parser.add_argument('--k-range', nargs=2, type=int, default=[2, 8], 
                       help='K range (default: 2 8)')
    parser.add_argument('--cv-folds', type=int, default=5, 
                       help='Number of CV folds (default: 5)')
    parser.add_argument('--input', type=str, 
                       default='ANALYSIS/00-06-IBD/paper936k/merged_936k_final',
                       help='Input prefix from IBD pipeline')
    parser.add_argument('--output-dir', type=str,
                       default='ANALYSIS/00-07-ADMIXTURE',
                       help='Output directory')
    parser.add_argument('--pop-file', type=str,
                       default='ANALYSIS/00-08-PCA/pop_736.updated.tsv',
                       help='Population mapping file with columns FID IID POP')
    args = parser.parse_args()
    
    # Check dependencies
    print("=" * 60)
    print("ADMIXTURE ANALYSIS PIPELINE")
    print("=" * 60)
    print("\nChecking dependencies...")
    check_dependencies()
    
    # Setup
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    plink_path = 'plink'
    input_prefix = args.input
    
    # Check input files exist
    print(f"\nChecking input files...")
    missing_files = []
    for ext in ['.bed', '.bim', '.fam']:
        filepath = f'{input_prefix}{ext}'
        if not os.path.exists(filepath):
            missing_files.append(filepath)
    
    if missing_files:
        print(f"ERROR: Required input files not found:")
        for f in missing_files:
            print(f"  - {f}")
        print(f"\nExpected input from IBD pipeline: {input_prefix}.*")
        print("Please run the IBD de-duplication pipeline first")
        sys.exit(1)
    
    # Count samples
    with open(f'{input_prefix}.fam', 'r') as f:
        n_samples = sum(1 for line in f)
    print(f"✓ Input dataset: {n_samples} samples")
    
    if n_samples != 736:
        print(f"WARNING: Expected 736 samples but found {n_samples}")
        print("Continuing with sample preservation...")
    
    print(f"Using {mp.cpu_count()} CPU cores")
    print(f"Output directory: {output_dir}")
    
    # Step 1: LD pruning (preserving all samples)
    print("\n" + "=" * 40)
    print("STEP 1: LD PRUNING")
    print("=" * 40)
    pruned_prefix = prepare_for_admixture(plink_path, input_prefix, output_dir)
    
    # Verify sample count one more time
    final_samples = count_samples(f'{pruned_prefix}.fam')
    assert final_samples == n_samples, f"Sample loss detected: {final_samples} != {n_samples}"
    print(f"✓ Confirmed: All {final_samples} samples preserved")
    
    # Step 2: Run ADMIXTURE
    print("\n" + "=" * 40)
    print("STEP 2: ADMIXTURE ANALYSIS")
    print("=" * 40)
    
    if args.quick:
        print("Quick mode - Running K=2,3,4,5 without CV...")
        quick_run_mode(pruned_prefix, output_dir, K_values=[2, 3, 4, 5])
        optimal_k = 3  # Default for quick mode
        cv_errors = {}
    else:
        print(f"Running ADMIXTURE K={args.k_range[0]}-{args.k_range[1]}" + 
              (" with CV" if args.cv else " without CV"))
        cv_errors = run_admixture_parallel(
            pruned_prefix, output_dir, 
            K_range=tuple(args.k_range),
            use_cv=args.cv,
            cv_folds=args.cv_folds
        )
        
        if cv_errors:
            optimal_k = plot_cv_errors(cv_errors, output_dir)
            print(f"\n✓ Optimal K = {optimal_k} (lowest CV error)")
        else:
            optimal_k = 3  # Default if no CV
    
    # Step 3: Create plots
    print("\n" + "=" * 40)
    print("STEP 3: VISUALIZATION")
    print("=" * 40)
    
    # Load population labels if available
    population_map = None
    if args.pop_file:
        print(f"Loading population labels from: {args.pop_file}")
        population_map = load_population_map(args.pop_file)
        if population_map:
            print(f"✓ Loaded {len(population_map)} population assignments")
        else:
            print("  Population file not found or invalid format")
    
    # Plot results for each K
    fam_file = f'{pruned_prefix}.fam'
    plots_created = []
    
    for k in range(args.k_range[0], args.k_range[1] + 1):
        Q_file = f'{output_dir}/{os.path.basename(pruned_prefix)}.{k}.Q'
        if os.path.exists(Q_file):
            print(f"Creating plot for K={k}...")
            cv_error = cv_errors.get(k) if cv_errors else None
            plot_admixture_results(Q_file, fam_file, k, output_dir, 
                                 population_map, cv_error)
            plots_created.append(k)
        elif args.quick and k in [2, 3, 4, 5]:
            print(f"  WARNING: Q matrix not found for K={k}")
    
    # Summary
    print("\n" + "=" * 60)
    print("ADMIXTURE ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Results directory: {output_dir}/")
    print(f"Samples preserved: {final_samples}")
    print(f"Plots created for K: {plots_created}")
    if cv_errors:
        print(f"Optimal K: {optimal_k} (based on CV)")
        print(f"CV errors: {output_dir}/cv_errors.txt")
        print(f"CV plot: {output_dir}/cv_errors.png")
    print(f"\nKey outputs:")
    print(f"  - Pruned dataset: {pruned_prefix}.*")
    print(f"  - Q matrices: {output_dir}/*.Q")
    print(f"  - P matrices: {output_dir}/*.P")
    print(f"  - Plots: {output_dir}/admixture_K*.png")
    print(f"  - Logs: {output_dir}/admixture_K*.log")

if __name__ == "__main__":
    main()