#!/usr/bin/env python3
"""
Order populations by similarity in ADMIXTURE and generate publication-ready figures.
Uses mean ADMIXTURE profiles to compute population similarity and applies consistent
ordering and colors across all K values.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os
import re
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
ORDERING_MODE = "global"  # "global" or "within-region" 
PANEL_K_MIN = 2
PANEL_K_MAX = 8

# Try scipy for optimal methods
try:
    from scipy.optimize import linear_sum_assignment
    from scipy.cluster.hierarchy import linkage, dendrogram, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import pdist, squareform, cosine
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not found. Using fallback methods for clustering and matching.")

def parse_population_file(pop_file):
    """Parse population labels with optional REGION column."""
    if not os.path.exists(pop_file):
        raise FileNotFoundError(f"Population file not found: {pop_file}")
    
    # Try different delimiters
    for sep in ['\t', ' ', ',']:
        try:
            # Try with header
            df = pd.read_csv(pop_file, sep=sep)
            df.columns = [col.upper() for col in df.columns]
            
            if all(col in df.columns for col in ['FID', 'IID', 'POP']):
                result = df[['FID', 'IID', 'POP']]
                if 'REGION' in df.columns:
                    result['REGION'] = df['REGION']
                print(f"  Successfully parsed population file with {len(df)} entries")
                return result
        except:
            pass
    
    # Try whitespace delimiter
    try:
        df = pd.read_csv(pop_file, sep=r'\s+', engine='python')
        df.columns = [col.upper() for col in df.columns]
        if all(col in df.columns for col in ['FID', 'IID', 'POP']):
            result = df[['FID', 'IID', 'POP']]
            if 'REGION' in df.columns:
                result['REGION'] = df['REGION']
            print(f"  Successfully parsed population file with {len(df)} entries")
            return result
    except:
        pass
    
    raise ValueError(f"Could not parse population file: {pop_file}")

def parse_cv_error(log_file):
    """Extract CV error from ADMIXTURE log."""
    if not os.path.exists(log_file):
        return None
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        patterns = [
            r'CV error \(K=\d+\):\s*([\d.]+)',
            r'CV error.*?=\s*([\d.]+)',
            r'CV:\s*([\d.]+)'
        ]
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                return float(match.group(1))
    except:
        pass
    return None

def get_okabe_ito_palette(n):
    """Okabe-Ito colorblind-safe palette."""
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', 
              '#0072B2', '#D55E00', '#CC79A7', '#999999',
              '#000000', '#B19CD9', '#FF6961', '#77DD77']
    
    if n > len(colors):
        # Generate additional colors
        import colorsys
        extra = []
        for i in range(n - len(colors)):
            hue = (i * 360 / (n - len(colors))) / 360
            rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.8)
            extra.append('#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)))
        colors.extend(extra)
    
    return colors[:n]

def align_ancestry_components(Q_matrices, K_values):
    """Align ancestry components across K for consistent colors."""
    if len(K_values) < 2:
        return {K_values[0]: list(range(Q_matrices[K_values[0]].shape[1]))}
    
    alignments = {}
    alignments[K_values[0]] = list(range(Q_matrices[K_values[0]].shape[1]))
    
    for i in range(1, len(K_values)):
        prev_K, curr_K = K_values[i-1], K_values[i]
        Q_prev, Q_curr = Q_matrices[prev_K], Q_matrices[curr_K]
        n_prev, n_curr = Q_prev.shape[1], Q_curr.shape[1]
        
        # Correlation matrix
        corr_matrix = np.zeros((n_curr, max(n_prev, n_curr)))
        for j in range(n_curr):
            for k in range(n_prev):
                if alignments[prev_K][k] < corr_matrix.shape[1]:
                    corr = np.corrcoef(Q_curr[:, j], Q_prev[:, k])[0, 1]
                    if not np.isnan(corr):
                        corr_matrix[j, alignments[prev_K][k]] = corr
        
        if HAS_SCIPY:
            # Hungarian algorithm for optimal matching
            row_ind, col_ind = linear_sum_assignment(-corr_matrix[:, :max(n_prev, n_curr)])
            alignment = col_ind[:n_curr].tolist()
        else:
            # Greedy matching fallback
            alignment = []
            used = set()
            for j in range(n_curr):
                correlations = [(k, corr_matrix[j, k]) for k in range(max(n_prev, n_curr)) if k not in used]
                if correlations:
                    best_k = max(correlations, key=lambda x: x[1])[0]
                    alignment.append(best_k)
                    used.add(best_k)
                else:
                    for k in range(max(n_prev, n_curr)):
                        if k not in used:
                            alignment.append(k)
                            used.add(k)
                            break
        
        alignments[curr_K] = alignment
    
    return alignments

def compute_cosine_distance_matrix(pop_means):
    """Compute pairwise cosine distance between population mean Q vectors."""
    n_pops = len(pop_means)
    distance_matrix = np.zeros((n_pops, n_pops))
    
    for i in range(n_pops):
        for j in range(i+1, n_pops):
            # Cosine distance = 1 - cosine similarity
            vec_i = pop_means.iloc[i].values
            vec_j = pop_means.iloc[j].values
            
            # Compute cosine similarity
            dot_product = np.dot(vec_i, vec_j)
            norm_i = np.linalg.norm(vec_i)
            norm_j = np.linalg.norm(vec_j)
            
            if norm_i > 0 and norm_j > 0:
                cos_sim = dot_product / (norm_i * norm_j)
                distance = 1 - cos_sim
            else:
                distance = 1.0
            
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    
    return distance_matrix

def seriation_fallback(distance_matrix, pop_names):
    """Simple seriation using nearest neighbor traveling salesman approach."""
    n = len(pop_names)
    if n <= 2:
        return list(range(n))
    
    # Start with the population with minimum average distance to others
    avg_distances = distance_matrix.mean(axis=1)
    current = np.argmin(avg_distances)
    
    visited = {current}
    order = [current]
    
    # Nearest neighbor approach
    while len(visited) < n:
        distances = distance_matrix[current].copy()
        # Set visited distances to infinity
        for v in visited:
            distances[v] = np.inf
        
        nearest = np.argmin(distances)
        visited.add(nearest)
        order.append(nearest)
        current = nearest
    
    return order

def compute_population_similarity_order(Q_matrix, pop_df, K_star, output_dir):
    """Compute population similarity order based on mean ADMIXTURE profiles."""
    print(f"\nComputing population similarity order using K={K_star}")
    
    # Create Q dataframe
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K_star)])
    Q_df['POP'] = pop_df['POP'].values
    
    # Check for REGION column
    has_region = 'REGION' in pop_df.columns
    if has_region:
        Q_df['REGION'] = pop_df['REGION'].values
        print(f"  Found REGION column, using ordering mode: {ORDERING_MODE}")
    
    # Calculate population means
    pop_means = Q_df.groupby('POP')[[f'Anc_{i+1}' for i in range(K_star)]].mean()
    pop_counts = Q_df.groupby('POP').size()
    
    # Warn about small populations
    small_pops = pop_counts[pop_counts < 5]
    if len(small_pops) > 0:
        print(f"  WARNING: {len(small_pops)} populations with <5 individuals: {list(small_pops.index)}")
    
    print(f"  Number of populations: {len(pop_means)}")
    print(f"  Distance metric: cosine distance")
    
    # Compute distance matrix
    distance_matrix = compute_cosine_distance_matrix(pop_means)
    
    if has_region and ORDERING_MODE == "within-region":
        # Get region for each population
        pop_regions = Q_df.groupby('POP')['REGION'].first()
        
        # Define region order
        region_order = ['Coastal', 'Highland', 'Andes', 'Amazon']
        
        # Order populations within each region
        final_order = []
        for region in region_order:
            if region not in pop_regions.values:
                continue
            
            region_pops = pop_regions[pop_regions == region].index.tolist()
            if len(region_pops) <= 1:
                final_order.extend(region_pops)
                continue
            
            # Get indices for this region
            region_indices = [list(pop_means.index).index(pop) for pop in region_pops]
            
            # Extract sub-matrix for this region
            region_dist = distance_matrix[np.ix_(region_indices, region_indices)]
            
            if HAS_SCIPY and len(region_pops) > 2:
                # Hierarchical clustering for this region
                condensed_dist = squareform(region_dist)
                Z = linkage(condensed_dist, method='average')
                Z_ordered = optimal_leaf_ordering(Z, condensed_dist)
                region_order_idx = leaves_list(Z_ordered)
            else:
                # Fallback seriation
                region_order_idx = seriation_fallback(region_dist, region_pops)
            
            # Add to final order
            ordered_region_pops = [region_pops[i] for i in region_order_idx]
            final_order.extend(ordered_region_pops)
        
        ordered_pops = final_order
        print(f"  Ordering: within-region (regions: {', '.join(region_order)})")
        
    else:
        # Global clustering
        if HAS_SCIPY and len(pop_means) > 2:
            # Hierarchical clustering with optimal leaf ordering
            condensed_dist = squareform(distance_matrix)
            Z = linkage(condensed_dist, method='average')
            Z_ordered = optimal_leaf_ordering(Z, condensed_dist)
            order_indices = leaves_list(Z_ordered)
            
            # Plot dendrogram
            fig, ax = plt.subplots(figsize=(12, 8))
            dendrogram(Z_ordered, labels=pop_means.index.tolist(), 
                      leaf_rotation=90, leaf_font_size=10, ax=ax)
            ax.set_title(f'Population Clustering Dendrogram (K={K_star})', fontsize=14, fontweight='bold')
            ax.set_xlabel('Population', fontsize=12)
            ax.set_ylabel('Cosine Distance', fontsize=12)
            plt.tight_layout()
            dendro_file = output_dir / f'admixture_population_dendrogram_K{K_star}.png'
            plt.savefig(dendro_file, dpi=400, bbox_inches='tight')
            plt.close()
            print(f"  Saved dendrogram: {dendro_file}")
            print(f"  Clustering method: scipy hierarchical with optimal leaf ordering")
            
        else:
            # Fallback seriation
            order_indices = seriation_fallback(distance_matrix, pop_means.index.tolist())
            print(f"  Clustering method: seriation fallback (nearest neighbor)")
        
        ordered_pops = [pop_means.index[i] for i in order_indices]
        print(f"  Ordering: global clustering")
    
    # Save population order
    order_file = output_dir / 'admixture_population_order.txt'
    with open(order_file, 'w') as f:
        for pop in ordered_pops:
            f.write(f'{pop}\n')
    print(f"  Saved population order: {order_file}")
    
    # Save population counts
    counts_file = output_dir / 'counts_by_population.tsv'
    counts_df = pd.DataFrame({
        'Population': ordered_pops,
        'N': [pop_counts[pop] for pop in ordered_pops]
    })
    counts_df.to_csv(counts_file, sep='\t', index=False)
    print(f"  Saved population counts: {counts_file}")
    
    return ordered_pops

def plot_admixture_grouped(Q_matrix, pop_df, K, colors, cv_error, output_file, pop_order):
    """Create population-grouped stacked bar plot with specified order."""
    n_samples = len(pop_df)
    
    # Create Q dataframe
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K)])
    Q_df = pd.concat([pop_df.reset_index(drop=True), Q_df], axis=1)
    
    # Sort by population order, then by dominant ancestry
    Q_df['POP'] = pd.Categorical(Q_df['POP'], categories=pop_order, ordered=True)
    Q_df['max_anc'] = Q_df[[f'Anc_{i+1}' for i in range(K)]].max(axis=1)
    Q_df = Q_df.sort_values(['POP', 'max_anc'], ascending=[True, False])
    
    # Dynamic figure width
    fig_width = max(24, n_samples / 25)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    
    # Plot stacked bars
    ind = np.arange(len(Q_df))
    bottom = np.zeros(len(Q_df))
    
    for k in range(K):
        values = Q_df[f'Anc_{k+1}'].values
        ax.bar(ind, values, bottom=bottom, width=1.0,
               color=colors[k], linewidth=0, edgecolor='none')
        bottom += values
    
    # Population boundaries and labels
    pop_groups = Q_df.groupby('POP', sort=False)
    pop_boundaries = []
    pop_centers = []
    pop_labels = []
    
    current_pos = 0
    for pop_name, group in pop_groups:
        group_size = len(group)
        pop_centers.append(current_pos + group_size / 2)
        pop_labels.append(f'{pop_name} (n={group_size})')
        if current_pos > 0:
            pop_boundaries.append(current_pos)
        current_pos += group_size
    
    # White separators
    for boundary in pop_boundaries:
        ax.axvline(x=boundary, color='white', linewidth=0.5, zorder=3)
    
    # X-axis
    ax.set_xticks(pop_centers)
    ax.set_xticklabels(pop_labels, rotation=45, ha='right', fontsize=10)
    
    # Title and labels
    title = f'ADMIXTURE (K={K}; n={n_samples}'
    if cv_error is not None:
        title += f'; CV={cv_error:.6f}'
    title += ')'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_ylabel('Ancestry Proportion', fontsize=12)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, len(Q_df)])
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    legend_elements = [mpatches.Patch(facecolor=colors[i], label=f'Ancestry {i+1}')
                      for i in range(K)]
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
             frameon=False, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=400, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def plot_population_heatmap(Q_matrix, pop_df, K, output_png, output_csv, pop_order):
    """Create population mean heatmap with specified order."""
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Ancestry_{i+1}' for i in range(K)])
    Q_df['POP'] = pop_df['POP'].values
    
    # Calculate population means
    pop_means = Q_df.groupby('POP').mean()
    pop_counts = Q_df.groupby('POP').size()
    
    # Reorder populations
    pop_means = pop_means.reindex(pop_order)
    pop_counts = pop_counts.reindex(pop_order)
    
    # Save to CSV
    pop_means_with_counts = pop_means.copy()
    pop_means_with_counts['Sample_Count'] = pop_counts
    pop_means_with_counts.to_csv(output_csv)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(max(8, K*0.5), len(pop_means) * 0.4 + 1))
    
    im = ax.imshow(pop_means.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
    
    ax.set_xticks(np.arange(K))
    ax.set_xticklabels([f'Ancestry {i+1}' for i in range(K)])
    ax.set_yticks(np.arange(len(pop_means)))
    ax.set_yticklabels([f'{pop} (n={pop_counts[pop]})' for pop in pop_means.index])
    
    # Annotate cells
    for i in range(len(pop_means)):
        for j in range(K):
            value = pop_means.iloc[i, j]
            color = 'white' if value > 0.5 else 'black'
            ax.text(j, i, f'{value:.2f}', ha='center', va='center',
                   color=color, fontsize=8, fontweight='bold')
    
    plt.colorbar(im, ax=ax, label='Mean Ancestry Proportion')
    ax.set_title(f'Population Mean Ancestry Proportions (K={K})', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_png, dpi=400, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_png}")
    print(f"  Saved: {output_csv}")

def plot_multi_panel(Q_matrices, pop_df_dict, K_values, colors_dict, cv_errors, output_file, pop_order):
    """Create multi-panel comparison with specified order."""
    n_panels = len(K_values)
    fig = plt.figure(figsize=(24, n_panels * 2.5))
    gs = GridSpec(n_panels, 1, figure=fig, hspace=0.05)
    
    for idx, K in enumerate(K_values):
        ax = fig.add_subplot(gs[idx, 0])
        
        Q_matrix = Q_matrices[K]
        pop_df = pop_df_dict[K]
        colors = colors_dict[K]
        cv_error = cv_errors.get(K)
        
        # Prepare data
        Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K)])
        Q_df = pd.concat([pop_df.reset_index(drop=True), Q_df], axis=1)
        Q_df['POP'] = pd.Categorical(Q_df['POP'], categories=pop_order, ordered=True)
        Q_df['max_anc'] = Q_df[[f'Anc_{i+1}' for i in range(K)]].max(axis=1)
        Q_df = Q_df.sort_values(['POP', 'max_anc'], ascending=[True, False])
        
        # Plot
        ind = np.arange(len(Q_df))
        bottom = np.zeros(len(Q_df))
        
        for k in range(K):
            values = Q_df[f'Anc_{k+1}'].values
            ax.bar(ind, values, bottom=bottom, width=1.0,
                   color=colors[k], linewidth=0, edgecolor='none')
            bottom += values
        
        # Population boundaries
        pop_groups = Q_df.groupby('POP', sort=False)
        pop_boundaries = []
        pop_centers = []
        pop_labels = []
        
        current_pos = 0
        for pop_name, group in pop_groups:
            group_size = len(group)
            pop_centers.append(current_pos + group_size / 2)
            pop_labels.append(f'{pop_name} (n={group_size})')
            if current_pos > 0:
                pop_boundaries.append(current_pos)
            current_pos += group_size
        
        for boundary in pop_boundaries:
            ax.axvline(x=boundary, color='white', linewidth=0.5, zorder=3)
        
        # Title
        title = f'K={K}'
        if cv_error is not None:
            title += f' (CV={cv_error:.6f})'
        ax.set_title(title, fontsize=10, loc='left', pad=2, fontweight='bold')
        
        ax.set_ylim([0, 1])
        ax.set_xlim([0, len(Q_df)])
        ax.set_ylabel('Ancestry', fontsize=9)
        
        # X-axis labels only on bottom
        if idx == n_panels - 1:
            ax.set_xticks(pop_centers)
            ax.set_xticklabels(pop_labels, rotation=45, ha='right', fontsize=8)
        else:
            ax.set_xticks([])
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0', '0.5', '1'], fontsize=8)
    
    fig.suptitle(f'ADMIXTURE Analysis: K={K_values[0]} to K={K_values[-1]}', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=400, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def plot_cv_curve(cv_errors, output_file):
    """Plot CV error curve."""
    if not cv_errors:
        return
    
    K_values = sorted(cv_errors.keys())
    cv_values = [cv_errors[k] for k in K_values]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(K_values, cv_values, 'o-', linewidth=2, markersize=8, color='#0072B2')
    
    # Highlight minimum
    min_k = min(cv_errors, key=cv_errors.get)
    min_cv = cv_errors[min_k]
    ax.plot(min_k, min_cv, 'ro', markersize=12)
    ax.annotate(f'Optimal K={min_k}\nCV={min_cv:.6f}',
                xy=(min_k, min_cv), xytext=(min_k + 0.5, min_cv),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red')
    
    ax.set_xlabel('K (Number of Ancestral Populations)', fontsize=12)
    ax.set_ylabel('Cross-Validation Error', fontsize=12)
    ax.set_title('ADMIXTURE Cross-Validation Error', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks(K_values)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=400, bbox_inches='tight')
    plt.close()
    print(f"  Saved CV curve: {output_file}")

def main():
    # Set matplotlib defaults
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    # Setup paths
    base_dir = Path('ANALYSIS/00-07-ADMIXTURE')
    pop_file = Path('ANALYSIS/00-08-PCA/pop_736.updated.tsv')
    output_dir = base_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("ADMIXTURE POPULATION SIMILARITY ORDERING AND VISUALIZATION")
    print("=" * 70)
    print(f"ADMIXTURE directory: {base_dir}")
    print(f"Population file: {pop_file}")
    print(f"Output directory: {output_dir}")
    
    # Find FAM file
    fam_candidates = [
        base_dir / 'merged_936k_admixture_pruned.fam',
        base_dir / 'merged_936k_pruned.fam',
        base_dir / 'admixture_pruned.fam'
    ]
    
    fam_file = None
    for candidate in fam_candidates:
        if candidate.exists():
            fam_file = candidate
            break
    
    if not fam_file:
        print("ERROR: Could not find pruned FAM file!")
        return
    
    print(f"FAM file: {fam_file}")
    
    # Read FAM file
    print("\nLoading data...")
    fam_df = pd.read_csv(fam_file, sep=' ', header=None, usecols=[0, 1])
    fam_df.columns = ['FID', 'IID']
    print(f"  FAM file: {len(fam_df)} samples")
    
    # Read population labels
    pop_df = parse_population_file(pop_file)
    has_region = 'REGION' in pop_df.columns
    
    # Merge FAM with population data
    merged_df = fam_df.merge(pop_df, on='IID', how='left')
    
    # Check samples with population assignments
    has_pop = merged_df['POP'].notna()
    n_with_pop = has_pop.sum()
    n_without_pop = (~has_pop).sum()
    n_total = len(merged_df)
    
    print(f"  Samples with population: {n_with_pop}/{n_total}")
    if n_without_pop > 0:
        print(f"  WARNING: {n_without_pop} samples without population (will be excluded)")
    
    if n_with_pop == 0:
        print("\nERROR: No samples with population assignments!")
        return
    
    # Filter to samples with population
    merged_df = merged_df[has_pop].copy()
    
    # Detect available K values
    available_K = []
    for K in range(2, 31):
        q_candidates = [
            base_dir / f'merged_936k_admixture_pruned.{K}.Q',
            base_dir / f'merged_936k_pruned.{K}.Q',
            base_dir / f'admixture_pruned.{K}.Q'
        ]
        for candidate in q_candidates:
            if candidate.exists():
                available_K.append(K)
                break
    
    if not available_K:
        print("\nERROR: No Q matrices found!")
        return
    
    print(f"  Available K values: {available_K}")
    
    # Load Q matrices and CV errors
    print("\nLoading ADMIXTURE results...")
    Q_matrices = {}
    Q_matrices_filtered = {}
    pop_df_dict = {}
    cv_errors = {}
    
    for K in available_K:
        q_candidates = [
            base_dir / f'merged_936k_admixture_pruned.{K}.Q',
            base_dir / f'merged_936k_pruned.{K}.Q',
            base_dir / f'admixture_pruned.{K}.Q'
        ]
        
        q_file = None
        for candidate in q_candidates:
            if candidate.exists():
                q_file = candidate
                break
        
        if q_file:
            Q_full = np.loadtxt(q_file)
            Q_matrices[K] = Q_full
            Q_matrices_filtered[K] = Q_full[has_pop]
            pop_df_dict[K] = merged_df.copy()
            print(f"  Loaded K={K} ({q_file.name})")
            
            # Try to find CV error
            log_file = base_dir / f'admixture_K{K}.log'
            cv_error = parse_cv_error(log_file)
            if cv_error is not None:
                cv_errors[K] = cv_error
                print(f"    CV error: {cv_error:.6f}")
    
    # Choose K* for similarity computation
    if cv_errors:
        K_star = min(cv_errors, key=cv_errors.get)
        print(f"\nChosen K* for similarity: {K_star} (lowest CV error)")
    else:
        K_star = 7 if 7 in available_K else available_K[len(available_K)//2]
        print(f"\nChosen K* for similarity: {K_star} (default)")
    
    # Align ancestry components
    print("\nAligning ancestry components across K values...")
    alignments = align_ancestry_components(Q_matrices_filtered, available_K)
    
    # Create color dictionaries
    max_K = max(available_K)
    base_colors = get_okabe_ito_palette(max_K)
    
    colors_dict = {}
    for K in available_K:
        alignment = alignments[K]
        colors_dict[K] = [base_colors[alignment[i] if alignment[i] < len(base_colors) else i]
                         for i in range(K)]
    
    # Compute population similarity order
    pop_order = compute_population_similarity_order(
        Q_matrices_filtered[K_star], pop_df_dict[K_star], K_star, output_dir
    )
    
    # Generate plots with new order
    print("\nGenerating figures with similarity-based ordering...")
    
    # Individual K plots
    for K in available_K:
        # Grouped bar plot
        output_file = output_dir / f'admixture_grouped_K{K}_ordered.png'
        plot_admixture_grouped(Q_matrices_filtered[K], pop_df_dict[K], K, 
                              colors_dict[K], cv_errors.get(K), output_file, pop_order)
        
        # Population heatmap
        output_png = output_dir / f'admixture_population_means_K{K}_ordered.png'
        output_csv = output_dir / f'admixture_population_means_K{K}_ordered.csv'
        plot_population_heatmap(Q_matrices_filtered[K], pop_df_dict[K], K, 
                               output_png, output_csv, pop_order)
    
    # Multi-panel plot
    panel_K = [K for K in available_K if PANEL_K_MIN <= K <= PANEL_K_MAX]
    if panel_K:
        print(f"\nCreating multi-panel plot for K={panel_K[0]}-{panel_K[-1]}...")
        output_file = output_dir / f'admixture_grouped_K{panel_K[0]}_K{panel_K[-1]}_ordered_panel.png'
        plot_multi_panel({K: Q_matrices_filtered[K] for K in panel_K},
                        {K: pop_df_dict[K] for K in panel_K},
                        panel_K, {K: colors_dict[K] for K in panel_K},
                        cv_errors, output_file, pop_order)
    
    # CV curve
    if cv_errors:
        print("\nCreating CV error curve...")
        output_file = output_dir / 'admixture_cv_curve.png'
        plot_cv_curve(cv_errors, output_file)
        
        # Save CV errors
        cv_file = output_dir / 'cv_errors.txt'
        with open(cv_file, 'w') as f:
            f.write("K\tCV_Error\n")
            for K in sorted(cv_errors.keys()):
                f.write(f"{K}\t{cv_errors[K]:.6f}\n")
        print(f"  Saved CV errors: {cv_file}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Chosen K* for similarity: {K_star}")
    print(f"Number of populations: {merged_df['POP'].nunique()}")
    print(f"Distance metric: cosine distance")
    if HAS_SCIPY:
        print(f"Clustering method: hierarchical with optimal leaf ordering")
    else:
        print(f"Clustering method: seriation fallback (nearest neighbor)")
    print(f"Samples analyzed: {n_with_pop}")
    print(f"K values processed: {available_K}")
    print(f"Figures saved to: {output_dir}")
    print("\nAll processing completed successfully!")

if __name__ == "__main__":
    main()