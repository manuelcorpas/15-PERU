#!/usr/bin/env python3
"""
Complete Pairwise Fst Analysis Pipeline with Clustering
Calculates pairwise Fst between all 28 populations and creates clustered visualizations
All outputs to ANALYSIS/00-09-FST/
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from pathlib import Path
import multiprocessing as mp
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

# Configuration
BASE_DIR = 'ANALYSIS/00-09-FST'
INPUT_PREFIX = 'ANALYSIS/00-06-IBD/paper936k/merged_936k_final'
POP_FILE = 'ANALYSIS/00-08-PCA/pop_736.updated.tsv'
THREADS = mp.cpu_count()

# ============================================================================
# PART 1: FST CALCULATION
# ============================================================================

def check_dependencies():
    """Check required tools are available"""
    tools = {'plink': 'PLINK v1.9+ is required for Fst calculation'}
    missing = []
    
    for tool, msg in tools.items():
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=False)
            print(f"✓ {tool} found")
        except FileNotFoundError:
            missing.append(msg)
    
    if missing:
        print("\nERROR: Missing required tools:")
        for msg in missing:
            print(f"  - {msg}")
        sys.exit(1)

def prepare_directories():
    """Create output directory structure"""
    dirs = [
        BASE_DIR,
        f"{BASE_DIR}/pairwise_fst",
        f"{BASE_DIR}/population_subsets",
        f"{BASE_DIR}/logs",
        f"{BASE_DIR}/results"
    ]
    
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    
    print(f"✓ Created output directory: {BASE_DIR}")

def load_population_data():
    """Load and validate population assignments"""
    print("\nLoading population data...")
    
    # Check files exist
    for ext in ['.bed', '.bim', '.fam']:
        if not os.path.exists(f"{INPUT_PREFIX}{ext}"):
            print(f"ERROR: Required file not found: {INPUT_PREFIX}{ext}")
            sys.exit(1)
    
    if not os.path.exists(POP_FILE):
        print(f"ERROR: Population file not found: {POP_FILE}")
        sys.exit(1)
    
    # Load population assignments
    pop_df = pd.read_csv(POP_FILE, sep='\t')
    
    # Standardize column names
    pop_df.columns = [col.upper() for col in pop_df.columns]
    
    # Ensure required columns exist
    if 'IID' not in pop_df.columns or 'POP' not in pop_df.columns:
        print("ERROR: Population file must have IID and POP columns")
        sys.exit(1)
    
    # Add FID if missing
    if 'FID' not in pop_df.columns:
        pop_df['FID'] = pop_df['IID']
    
    # Clean data
    pop_df['IID'] = pop_df['IID'].astype(str)
    pop_df['FID'] = pop_df['FID'].astype(str)
    pop_df['POP'] = pop_df['POP'].astype(str)
    
    # Get population counts
    pop_counts = pop_df['POP'].value_counts().sort_index()
    
    print(f"✓ Loaded {len(pop_df)} samples from {len(pop_counts)} populations")
    print("\nPopulation counts:")
    for pop, count in pop_counts.items():
        print(f"  {pop}: {count} samples")
    
    return pop_df, pop_counts

def create_population_pairs(populations):
    """Generate all pairwise combinations of populations"""
    pairs = list(combinations(sorted(populations), 2))
    print(f"\n✓ Generated {len(pairs)} population pairs for Fst calculation")
    return pairs

def calculate_fst_pair(pop1, pop2, pop_df, pair_idx, total_pairs):
    """Calculate Fst for a single pair of populations"""
    print(f"  [{pair_idx}/{total_pairs}] Calculating Fst: {pop1} vs {pop2}")
    
    # Create keep files for each population
    keep1_file = f"{BASE_DIR}/population_subsets/{pop1}.keep"
    keep2_file = f"{BASE_DIR}/population_subsets/{pop2}.keep"
    
    # Write samples for population 1
    pop1_samples = pop_df[pop_df['POP'] == pop1][['FID', 'IID']]
    pop1_samples.to_csv(keep1_file, sep='\t', header=False, index=False)
    
    # Write samples for population 2
    pop2_samples = pop_df[pop_df['POP'] == pop2][['FID', 'IID']]
    pop2_samples.to_csv(keep2_file, sep='\t', header=False, index=False)
    
    # Combine both populations
    combined_file = f"{BASE_DIR}/population_subsets/{pop1}_{pop2}.keep"
    combined_samples = pd.concat([pop1_samples, pop2_samples])
    combined_samples.to_csv(combined_file, sep='\t', header=False, index=False)
    
    # Create cluster file for these two populations
    cluster_file = f"{BASE_DIR}/population_subsets/{pop1}_{pop2}.clust"
    with open(cluster_file, 'w') as f:
        for _, row in pop1_samples.iterrows():
            f.write(f"{row['FID']}\t{row['IID']}\t1\n")
        for _, row in pop2_samples.iterrows():
            f.write(f"{row['FID']}\t{row['IID']}\t2\n")
    
    # Output prefix
    output_prefix = f"{BASE_DIR}/pairwise_fst/{pop1}_{pop2}"
    
    # Run PLINK Fst calculation
    cmd = [
        'plink',
        '--bfile', INPUT_PREFIX,
        '--keep', combined_file,
        '--fst',
        '--within', cluster_file,
        '--out', output_prefix
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse Fst from log file
        log_file = f"{output_prefix}.log"
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log_content = f.read()
                
            # Look for Fst value in log
            for line in log_content.split('\n'):
                if 'Fst estimate:' in line:
                    fst_value = float(line.split(':')[1].strip())
                    return pop1, pop2, fst_value
                elif 'weighted Fst estimate:' in line:
                    fst_value = float(line.split(':')[1].strip())
                    return pop1, pop2, fst_value
        
        # If not found in log, try reading .fst file
        fst_file = f"{output_prefix}.fst"
        if os.path.exists(fst_file):
            fst_df = pd.read_csv(fst_file, sep=r'\s+')
            if 'FST' in fst_df.columns:
                fst_value = fst_df['FST'].mean()
                return pop1, pop2, fst_value
        
        print(f"    WARNING: Could not parse Fst value for {pop1} vs {pop2}")
        return pop1, pop2, np.nan
        
    except subprocess.CalledProcessError as e:
        print(f"    ERROR calculating Fst for {pop1} vs {pop2}: {e}")
        return pop1, pop2, np.nan

def calculate_all_fst(pop_df, populations):
    """Calculate all pairwise Fst values"""
    pairs = create_population_pairs(populations)
    total_pairs = len(pairs)
    
    print(f"\nCalculating pairwise Fst...")
    
    # Store results
    fst_results = []
    
    # Sequential processing (more stable for file I/O)
    for idx, (pop1, pop2) in enumerate(pairs, 1):
        result = calculate_fst_pair(pop1, pop2, pop_df, idx, total_pairs)
        fst_results.append(result)
    
    return fst_results

def create_fst_matrix(fst_results, populations):
    """Create symmetric matrix from pairwise Fst results"""
    n_pops = len(populations)
    pop_list = sorted(populations)
    
    # Initialize matrix
    fst_matrix = np.zeros((n_pops, n_pops))
    
    # Fill matrix
    for pop1, pop2, fst in fst_results:
        if not np.isnan(fst):
            idx1 = pop_list.index(pop1)
            idx2 = pop_list.index(pop2)
            fst_matrix[idx1, idx2] = fst
            fst_matrix[idx2, idx1] = fst  # Symmetric
    
    # Create DataFrame
    fst_df = pd.DataFrame(fst_matrix, index=pop_list, columns=pop_list)
    
    return fst_df

def save_fst_results(fst_df, fst_results, pop_counts):
    """Save Fst results in multiple formats"""
    print("\nSaving Fst results...")
    
    # Save matrix as CSV
    matrix_file = f"{BASE_DIR}/results/fst_matrix.csv"
    fst_df.to_csv(matrix_file)
    print(f"✓ Saved Fst matrix: {matrix_file}")
    
    # Save pairwise list
    pairwise_file = f"{BASE_DIR}/results/fst_pairwise.tsv"
    pairwise_df = pd.DataFrame(fst_results, columns=['Population1', 'Population2', 'Fst'])
    
    # Add sample sizes
    pairwise_df['N1'] = pairwise_df['Population1'].map(pop_counts)
    pairwise_df['N2'] = pairwise_df['Population2'].map(pop_counts)
    
    # Sort by Fst value
    pairwise_df = pairwise_df.sort_values('Fst', ascending=False)
    pairwise_df.to_csv(pairwise_file, sep='\t', index=False)
    print(f"✓ Saved pairwise Fst: {pairwise_file}")
    
    # Calculate summary statistics
    fst_values = pairwise_df['Fst'].dropna()
    
    summary = {
        'n_populations': len(fst_df),
        'n_pairs': len(pairwise_df),
        'mean_fst': fst_values.mean(),
        'median_fst': fst_values.median(),
        'min_fst': fst_values.min(),
        'max_fst': fst_values.max(),
        'std_fst': fst_values.std()
    }
    
    summary_file = f"{BASE_DIR}/results/fst_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("FST ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        for key, value in summary.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.6f}\n")
            else:
                f.write(f"{key}: {value}\n")
        
        f.write("\nTop 10 highest Fst pairs:\n")
        f.write("-" * 30 + "\n")
        for _, row in pairwise_df.head(10).iterrows():
            f.write(f"{row['Population1']} vs {row['Population2']}: {row['Fst']:.6f}\n")
        
        f.write("\nTop 10 lowest Fst pairs:\n")
        f.write("-" * 30 + "\n")
        for _, row in pairwise_df.tail(10).iterrows():
            f.write(f"{row['Population1']} vs {row['Population2']}: {row['Fst']:.6f}\n")
    
    print(f"✓ Saved summary: {summary_file}")
    
    return summary

# ============================================================================
# PART 2: CLUSTERING AND VISUALIZATION
# ============================================================================

def perform_clustering(fst_df, method='average'):
    """Perform hierarchical clustering on Fst matrix"""
    # Convert Fst to distance matrix
    # Handle negative Fst values by setting them to 0
    fst_dist = fst_df.values.copy()
    fst_dist[fst_dist < 0] = 0
    
    # Ensure it's a proper distance matrix
    np.fill_diagonal(fst_dist, 0)
    fst_dist = (fst_dist + fst_dist.T) / 2
    
    # Convert to condensed distance matrix for clustering
    condensed_dist = squareform(fst_dist)
    
    # Perform hierarchical clustering
    Z = linkage(condensed_dist, method=method)
    
    # Get the order of populations after clustering
    order = leaves_list(Z)
    
    return Z, order

def plot_alphabetical_heatmap(fst_df, output_prefix):
    """Create traditional alphabetical heatmap"""
    print("  Creating alphabetical heatmap...")
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Create mask for upper triangle
    mask = np.triu(np.ones_like(fst_df, dtype=bool))
    
    sns.heatmap(fst_df, 
                mask=mask,
                cmap='YlOrRd',
                vmin=0,
                vmax=fst_df.max().max(),
                annot=True,
                fmt='.3f',
                annot_kws={'size': 7},
                square=True,
                linewidths=0.5,
                cbar_kws={'label': 'Fst'},
                ax=ax)
    
    ax.set_title('Pairwise Fst between Peruvian Populations (Alphabetical)', 
                 fontsize=16, pad=20)
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_prefix}.{ext}"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"    ✓ Saved: {output_file}")
    
    plt.close()

def plot_clustered_heatmap(fst_df, order, output_prefix):
    """Create clustered heatmap"""
    print("  Creating clustered heatmap...")
    
    # Reorder the matrix
    clustered_df = fst_df.iloc[order, order]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Create mask for upper triangle
    mask = np.triu(np.ones_like(clustered_df, dtype=bool))
    
    sns.heatmap(clustered_df,
                mask=mask,
                cmap='RdYlBu_r',
                center=0,
                vmin=clustered_df.min().min(),
                vmax=clustered_df.max().max(),
                annot=False,  # No values for cleaner look
                square=True,
                linewidths=0.5,
                cbar_kws={'label': 'Fst', 'shrink': 0.8},
                ax=ax)
    
    ax.set_title('Pairwise Fst - Populations Clustered by Genetic Similarity', 
                 fontsize=16, pad=20, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)
    
    fig.text(0.02, 0.02, 
             'Populations ordered by hierarchical clustering (UPGMA) based on Fst distances',
             fontsize=8, style='italic', alpha=0.7)
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_prefix}.{ext}"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"    ✓ Saved: {output_file}")
    
    plt.close()

def plot_dendrogram_with_heatmap(fst_df, Z, order, output_prefix):
    """Create combined dendrogram and heatmap"""
    print("  Creating dendrogram with heatmap...")
    
    # Reorder the matrix
    clustered_df = fst_df.iloc[order, order]
    
    # Create figure with GridSpec
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.3, 1], height_ratios=[0.3, 1],
                          hspace=0.02, wspace=0.02)
    
    # Dendrogram on top
    ax_dendro_top = fig.add_subplot(gs[0, 1])
    dendro_top = dendrogram(Z, ax=ax_dendro_top, orientation='top',
                           labels=None, no_labels=True,
                           color_threshold=0, above_threshold_color='black')
    ax_dendro_top.set_xticks([])
    ax_dendro_top.set_yticks([])
    
    # Dendrogram on left
    ax_dendro_left = fig.add_subplot(gs[1, 0])
    dendro_left = dendrogram(Z, ax=ax_dendro_left, orientation='left',
                            labels=fst_df.index.tolist(),
                            color_threshold=0, above_threshold_color='black')
    ax_dendro_left.set_xticks([])
    ax_dendro_left.set_yticks([])
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1, 1])
    
    mask = np.triu(np.ones_like(clustered_df, dtype=bool))
    
    sns.heatmap(clustered_df,
                mask=mask,
                cmap='RdYlBu_r',
                center=0,
                vmin=clustered_df.min().min(),
                vmax=clustered_df.max().max(),
                annot=False,
                square=True,
                linewidths=0.3,
                cbar_kws={'label': 'Fst', 'shrink': 0.5, 'pad': 0.02},
                xticklabels=clustered_df.columns,
                yticklabels=False,
                ax=ax_heatmap)
    
    plt.setp(ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    
    fig.suptitle('Hierarchical Clustering of Peruvian Populations based on Fst',
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_prefix}.{ext}"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"    ✓ Saved: {output_file}")
    
    plt.close()

def plot_fst_distribution(fst_results, output_prefix):
    """Plot distribution of Fst values"""
    print("  Creating distribution plots...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Extract Fst values
    fst_values = [fst for _, _, fst in fst_results if not np.isnan(fst)]
    
    # Histogram
    ax = axes[0]
    ax.hist(fst_values, bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(np.mean(fst_values), color='red', linestyle='--', 
               label=f'Mean = {np.mean(fst_values):.4f}')
    ax.axvline(np.median(fst_values), color='green', linestyle='--',
               label=f'Median = {np.median(fst_values):.4f}')
    ax.set_xlabel('Fst', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of Pairwise Fst Values', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Box plot
    ax = axes[1]
    ax.boxplot([fst_values], labels=['All pairs'])
    ax.set_ylabel('Fst', fontsize=12)
    ax.set_title('Fst Value Range', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_prefix}.{ext}"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"    ✓ Saved: {output_file}")
    
    plt.close()

def identify_population_clusters(fst_df, Z, threshold=0.03):
    """Identify major population clusters based on Fst threshold"""
    # Get cluster assignments
    clusters = fcluster(Z, threshold, criterion='distance')
    
    # Create cluster DataFrame
    cluster_df = pd.DataFrame({
        'Population': fst_df.index,
        'Cluster': clusters
    })
    
    # Save cluster assignments
    cluster_file = f"{BASE_DIR}/results/population_clusters.tsv"
    cluster_df.to_csv(cluster_file, sep='\t', index=False)
    
    # Count clusters
    cluster_counts = cluster_df['Cluster'].value_counts().sort_index()
    
    print(f"\n=== POPULATION CLUSTERS (Fst threshold: {threshold}) ===")
    print(f"Number of clusters: {len(cluster_counts)}\n")
    
    for cluster_id in sorted(cluster_df['Cluster'].unique()):
        pops = cluster_df[cluster_df['Cluster'] == cluster_id]['Population'].tolist()
        print(f"Cluster {cluster_id} ({len(pops)} populations):")
        print(f"  {', '.join(pops)}\n")
    
    print(f"✓ Saved cluster assignments: {cluster_file}")
    
    return cluster_df

def save_clustered_order(fst_df, order):
    """Save the clustered population order"""
    reordered_pops = fst_df.index[order].tolist()
    order_file = f"{BASE_DIR}/results/population_order_clustered.txt"
    with open(order_file, 'w') as f:
        for pop in reordered_pops:
            f.write(f"{pop}\n")
    print(f"✓ Saved clustered order: {order_file}")
    
    # Also save clustered matrix
    clustered_matrix = fst_df.iloc[order, order]
    matrix_file = f"{BASE_DIR}/results/fst_matrix_clustered.csv"
    clustered_matrix.to_csv(matrix_file)
    print(f"✓ Saved clustered matrix: {matrix_file}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    print("=" * 70)
    print("COMPLETE PAIRWISE FST ANALYSIS WITH CLUSTERING")
    print("=" * 70)
    print(f"Using {THREADS} CPU cores")
    
    # Check dependencies
    check_dependencies()
    
    # Create output directories
    prepare_directories()
    
    # Load population data
    pop_df, pop_counts = load_population_data()
    
    # Get list of populations
    populations = set(pop_df['POP'].unique())
    
    # ========== PART 1: CALCULATE FST ==========
    print("\n" + "=" * 50)
    print("PART 1: FST CALCULATION")
    print("=" * 50)
    
    # Calculate all pairwise Fst
    fst_results = calculate_all_fst(pop_df, populations)
    
    # Create matrix
    fst_matrix = create_fst_matrix(fst_results, populations)
    
    # Save results
    summary = save_fst_results(fst_matrix, fst_results, pop_counts)
    
    # ========== PART 2: CLUSTERING ANALYSIS ==========
    print("\n" + "=" * 50)
    print("PART 2: CLUSTERING ANALYSIS")
    print("=" * 50)
    
    print("\nPerforming hierarchical clustering...")
    Z, order = perform_clustering(fst_matrix, method='average')
    print("✓ Clustering complete (UPGMA method)")
    
    # Save clustered order and matrix
    save_clustered_order(fst_matrix, order)
    
    # Identify clusters
    cluster_df = identify_population_clusters(fst_matrix, Z, threshold=0.03)
    
    # ========== PART 3: VISUALIZATIONS ==========
    print("\n" + "=" * 50)
    print("PART 3: CREATING VISUALIZATIONS")
    print("=" * 50)
    
    # 1. Traditional alphabetical heatmap
    plot_alphabetical_heatmap(fst_matrix, f"{BASE_DIR}/results/fst_heatmap_alphabetical")
    
    # 2. Clustered heatmap
    plot_clustered_heatmap(fst_matrix, order, f"{BASE_DIR}/results/fst_heatmap_clustered")
    
    # 3. Dendrogram with heatmap
    plot_dendrogram_with_heatmap(fst_matrix, Z, order, f"{BASE_DIR}/results/fst_dendrogram_heatmap")
    
    # 4. Distribution plots
    plot_fst_distribution(fst_results, f"{BASE_DIR}/results/fst_distribution")
    
    # ========== FINAL SUMMARY ==========
    print("\n" + "=" * 70)
    print("FST ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Analyzed {summary['n_populations']} populations")
    print(f"Calculated {summary['n_pairs']} pairwise Fst values")
    print(f"Mean Fst: {summary['mean_fst']:.6f}")
    print(f"Range: {summary['min_fst']:.6f} - {summary['max_fst']:.6f}")
    print(f"\nAll results saved to: {BASE_DIR}/")
    print("\nKey outputs:")
    print("\nData files:")
    print(f"  - Fst matrix: {BASE_DIR}/results/fst_matrix.csv")
    print(f"  - Fst matrix (clustered): {BASE_DIR}/results/fst_matrix_clustered.csv")
    print(f"  - Pairwise values: {BASE_DIR}/results/fst_pairwise.tsv")
    print(f"  - Population clusters: {BASE_DIR}/results/population_clusters.tsv")
    print(f"  - Clustered order: {BASE_DIR}/results/population_order_clustered.txt")
    print(f"  - Summary stats: {BASE_DIR}/results/fst_summary.txt")
    print("\nVisualizations:")
    print(f"  - Alphabetical heatmap: {BASE_DIR}/results/fst_heatmap_alphabetical.png")
    print(f"  - Clustered heatmap: {BASE_DIR}/results/fst_heatmap_clustered.png")
    print(f"  - Dendrogram+heatmap: {BASE_DIR}/results/fst_dendrogram_heatmap.png")
    print(f"  - Distribution plots: {BASE_DIR}/results/fst_distribution.png")
    print("=" * 70)

if __name__ == "__main__":
    main()