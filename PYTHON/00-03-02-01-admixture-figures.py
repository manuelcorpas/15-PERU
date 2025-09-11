#!/usr/bin/env python3
"""
Polish ADMIXTURE figures with similarity-ordered populations (publication-ready).
Generates high-quality visualizations with consistent colors, no label overlaps,
and both PNG and PDF outputs.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os
import re
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
POPULATION_LABELS = {
    "AFRODESCENDIENTES": "Afro-Peruvian",
    "CHOPCCAS": "Chopccas",
    "JAQARUS": "Jaqaru",
    "ASHANINKA": "Ashaninka",
    "MACHIGUENGA": "Machiguenga",
    "AWAJUN": "Awajun",
    "CANDOSHI": "Candoshi",
    "MATZES": "Matses",
    "SHIPIBO": "Shipibo",
    "NAHUA": "Nahua",
    "LAMAS": "Lamas",
    "IQUITOS": "Iquitos",
    "CHACHAPOYAS": "Chachapoyas",
    "LAMBAYEQUE": "Lambayeque",
    "MOCHES": "Moches",
    "TALLAN": "Tallan",
    "TUMBES": "Tumbes",
    "LIMA": "Lima",
    "AREQUIPA": "Arequipa",
    "MOQUEGUA": "Moquegua",
    "TACNA": "Tacna",
    "PUNO": "Puno",
    "UROS": "Uros",
    "QUEROS": "Queros",
    "CUSCO": "Cusco",
    "AYACUCHO": "Ayacucho",
    "HUARAZ": "Huaraz"
}

# Region order for display
REGION_ORDER = ["Amazon", "Highlands", "Coast"]

# Try scipy for optimal methods
try:
    from scipy.optimize import linear_sum_assignment
    from scipy.cluster.hierarchy import linkage, dendrogram, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import pdist, squareform
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("INFO: scipy not found. Using fallback methods.")

def parse_arguments():
    """Parse command line arguments."""
    recompute_order = '--recompute-order' in sys.argv
    return recompute_order

def parse_population_file(pop_file):
    """Parse population labels with optional REGION column."""
    if not os.path.exists(pop_file):
        raise FileNotFoundError(f"Population file not found: {pop_file}")
    
    # Try different delimiters
    for sep in ['\t', ' ', ',']:
        try:
            df = pd.read_csv(pop_file, sep=sep)
            df.columns = [col.upper() for col in df.columns]
            
            if all(col in df.columns for col in ['FID', 'IID', 'POP']):
                result = df[['FID', 'IID', 'POP']]
                if 'REGION' in df.columns:
                    result['REGION'] = df['REGION']
                print(f"  Loaded population file: {len(df)} entries")
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
            print(f"  Loaded population file: {len(df)} entries")
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
    """Okabe-Ito colorblind-safe palette, extended deterministically."""
    colors = [
        '#E69F00',  # Orange
        '#56B4E9',  # Sky blue
        '#009E73',  # Bluish green
        '#F0E442',  # Yellow
        '#0072B2',  # Blue
        '#D55E00',  # Vermillion
        '#CC79A7',  # Reddish purple
        '#999999'   # Gray
    ]
    
    if n > len(colors):
        # Extend deterministically
        import colorsys
        extra = []
        for i in range(n - len(colors)):
            hue = (i * 360 / (n - len(colors))) / 360
            rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.75)
            extra.append('#{:02x}{:02x}{:02x}'.format(
                int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)))
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
            # Hungarian algorithm
            row_ind, col_ind = linear_sum_assignment(-corr_matrix[:, :max(n_prev, n_curr)])
            alignment = col_ind[:n_curr].tolist()
        else:
            # Greedy matching
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
            vec_i = pop_means.iloc[i].values
            vec_j = pop_means.iloc[j].values
            
            # Cosine similarity
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

def seriation_fallback(distance_matrix):
    """Simple seriation using nearest neighbor."""
    n = len(distance_matrix)
    if n <= 2:
        return list(range(n))
    
    # Start with population with minimum average distance
    avg_distances = distance_matrix.mean(axis=1)
    current = np.argmin(avg_distances)
    
    visited = {current}
    order = [current]
    
    # Nearest neighbor
    while len(visited) < n:
        distances = distance_matrix[current].copy()
        for v in visited:
            distances[v] = np.inf
        
        nearest = np.argmin(distances)
        visited.add(nearest)
        order.append(nearest)
        current = nearest
    
    return order

def compute_population_order(Q_matrix, pop_df, K_star, output_dir):
    """Compute population similarity order."""
    print(f"\nComputing population similarity order using K={K_star}")
    
    # Create Q dataframe
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K_star)])
    Q_df['POP'] = pop_df['POP'].values
    
    # Calculate population means
    pop_means = Q_df.groupby('POP')[[f'Anc_{i+1}' for i in range(K_star)]].mean()
    pop_counts = Q_df.groupby('POP').size()
    
    # Warn about small populations
    small_pops = pop_counts[pop_counts < 5]
    if len(small_pops) > 0:
        print(f"  WARNING: {len(small_pops)} populations with <5 individuals: {list(small_pops.index)}")
    
    print(f"  Number of populations: {len(pop_means)}")
    print(f"  Distance metric: 1 - cosine similarity")
    
    # Compute distance matrix
    distance_matrix = compute_cosine_distance_matrix(pop_means)
    
    if HAS_SCIPY and len(pop_means) > 2:
        # Hierarchical clustering
        condensed_dist = squareform(distance_matrix)
        Z = linkage(condensed_dist, method='average')
        Z_ordered = optimal_leaf_ordering(Z, condensed_dist)
        order_indices = leaves_list(Z_ordered)
        
        # Plot dendrogram
        fig, ax = plt.subplots(figsize=(14, 8))
        dendrogram(Z_ordered, labels=pop_means.index.tolist(), 
                  leaf_rotation=45, leaf_font_size=10, ax=ax)
        ax.set_title(f'Population Clustering Dendrogram (K={K_star})', fontsize=14, fontweight='bold')
        ax.set_xlabel('Population', fontsize=12)
        ax.set_ylabel('Cosine Distance', fontsize=12)
        plt.tight_layout()
        dendro_file = output_dir / f'admixture_population_dendrogram_K{K_star}.png'
        plt.savefig(dendro_file, dpi=450, bbox_inches='tight')
        plt.close()
        print(f"  Saved dendrogram: {dendro_file}")
        print(f"  Clustering: hierarchical with optimal leaf ordering")
    else:
        # Fallback seriation
        order_indices = seriation_fallback(distance_matrix)
        print(f"  Clustering: seriation fallback (nearest neighbor)")
    
    ordered_pops = [pop_means.index[i] for i in order_indices]
    
    # Save population order
    order_file = output_dir / 'admixture_population_order.txt'
    with open(order_file, 'w') as f:
        for pop in ordered_pops:
            f.write(f'{pop}\n')
    print(f"  Saved population order: {order_file}")
    
    return ordered_pops, pop_counts

def get_adaptive_fontsize(width_pixels):
    """Get adaptive font size based on population block width."""
    if width_pixels < 30:
        return 8
    elif width_pixels < 60:
        return 9
    else:
        return 10

def plot_admixture_grouped_nice(Q_matrix, pop_df, K, colors, cv_error, output_base, pop_order, pop_counts, dpi=450):
    """Create polished population-grouped stacked bar plot."""
    n_samples = len(pop_df)
    
    # Create Q dataframe
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K)])
    Q_df = pd.concat([pop_df.reset_index(drop=True), Q_df], axis=1)
    
    # Check for REGION column
    has_region = 'REGION' in pop_df.columns
    
    # Sort by population order, then by dominant ancestry
    Q_df['POP'] = pd.Categorical(Q_df['POP'], categories=pop_order, ordered=True)
    Q_df['max_anc'] = Q_df[[f'Anc_{i+1}' for i in range(K)]].max(axis=1)
    Q_df = Q_df.sort_values(['POP', 'max_anc'], ascending=[True, False])
    
    # Figure setup
    fig_width = max(24, n_samples / 25)
    fig_height = 6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
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
    pop_widths = []
    
    current_pos = 0
    for pop_name, group in pop_groups:
        group_size = len(group)
        pop_centers.append(current_pos + group_size / 2)
        pop_widths.append(group_size)
        
        # Get short label
        short_name = POPULATION_LABELS.get(pop_name, pop_name)
        count = pop_counts[pop_name]
        
        # Add asterisk for small populations
        if count < 5:
            short_name += '*'
        
        label = f'{short_name}\n(n={count})'
        pop_labels.append(label)
        
        if current_pos > 0:
            pop_boundaries.append(current_pos)
        current_pos += group_size
    
    # Population separators
    for boundary in pop_boundaries:
        ax.axvline(x=boundary, color='black', alpha=0.4, linewidth=0.6, zorder=3)
    
    # Adaptive font sizes for labels
    for center, label, width in zip(pop_centers, pop_labels, pop_widths):
        # Convert width to approximate pixels (assuming ~1000 pixels for full width)
        width_pixels = width * (1000 / n_samples)
        fontsize = get_adaptive_fontsize(width_pixels)
        ax.text(center, -0.02, label, ha='center', va='top', 
                rotation=40, fontsize=fontsize, transform=ax.get_xaxis_transform())
    
    # Region super-labels if available
    if has_region:
        region_groups = []
        current_region = None
        start_idx = 0
        
        for i, (_, row) in enumerate(Q_df.iterrows()):
            if row.get('REGION') != current_region:
                if current_region is not None:
                    region_groups.append((current_region, start_idx, i))
                current_region = row.get('REGION')
                start_idx = i
        
        if current_region is not None:
            region_groups.append((current_region, start_idx, len(Q_df)))
        
        # Draw region labels
        for region, start, end in region_groups:
            if region in REGION_ORDER:
                center = (start + end) / 2
                ax.text(center, 1.02, region, ha='center', va='bottom',
                       fontsize=11, fontweight='bold', transform=ax.get_xaxis_transform())
    
    # Title and labels
    title = f'ADMIXTURE (K={K}; n={n_samples}'
    if cv_error is not None:
        title += f'; CV={cv_error:.6f}'
    title += ')'
    ax.set_title(title, fontsize=14, fontweight='bold', loc='left')
    ax.set_ylabel('Ancestry Proportion', fontsize=12)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, len(Q_df)])
    ax.set_xticks([])
    
    # Y-axis
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(['0', '0.5', '1'], fontsize=10)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    legend_elements = [mpatches.Patch(facecolor=colors[i], label=f'Ancestry {i+1}')
                      for i in range(K)]
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
             frameon=False, fontsize=10)
    
    plt.subplots_adjust(bottom=0.20, right=0.85)
    
    # Save PNG and PDF
    for ext in ['png', 'pdf']:
        output_file = f"{output_base}.{ext}"
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=ext)
        print(f"  Saved: {output_file}")
    
    plt.close()

def plot_multi_panel_nice(Q_matrices, pop_df_dict, K_values, colors_dict, cv_errors, output_base, pop_order, pop_counts, dpi=450):
    """Create polished multi-panel comparison."""
    n_panels = len(K_values)
    n_samples = len(pop_df_dict[K_values[0]])
    
    # Calculate figure dimensions
    fig_width = max(24, n_samples / 25)
    fig_height = 1.6 * n_panels + 1.5
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create subplots
    axes = []
    for idx, K in enumerate(K_values):
        ax = plt.subplot(n_panels, 1, idx + 1)
        axes.append(ax)
        
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
        
        # Plot stacked bars
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
        pop_widths = []
        
        current_pos = 0
        for pop_name, group in pop_groups:
            group_size = len(group)
            pop_centers.append(current_pos + group_size / 2)
            pop_widths.append(group_size)
            
            # Get short label
            short_name = POPULATION_LABELS.get(pop_name, pop_name)
            count = pop_counts[pop_name]
            if count < 5:
                short_name += '*'
            label = f'{short_name}\n(n={count})'
            pop_labels.append(label)
            
            if current_pos > 0:
                pop_boundaries.append(current_pos)
            current_pos += group_size
        
        # Population separators
        for boundary in pop_boundaries:
            ax.axvline(x=boundary, color='black', alpha=0.4, linewidth=0.6, zorder=3)
        
        # Title (left-aligned)
        title = f'K={K}'
        if cv_error is not None:
            title += f' (CV={cv_error:.6f})'
        ax.set_title(title, fontsize=11, loc='left', pad=2, fontweight='bold')
        
        ax.set_ylim([0, 1])
        ax.set_xlim([0, len(Q_df)])
        ax.set_ylabel('', fontsize=9)
        
        # Y-ticks
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(['0', '0.5', '1'], fontsize=9)
        
        # X-axis labels only on bottom panel
        if idx == n_panels - 1:
            for center, label, width in zip(pop_centers, pop_labels, pop_widths):
                width_pixels = width * (1000 / n_samples)
                fontsize = get_adaptive_fontsize(width_pixels) - 1  # Slightly smaller
                ax.text(center, -0.02, label, ha='center', va='top',
                       rotation=40, fontsize=fontsize, transform=ax.get_xaxis_transform())
            ax.set_xticks([])
        else:
            ax.set_xticks([])
        
        # Remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Global title
    fig.suptitle(f'ADMIXTURE (n={n_samples}; K={K_values[0]}..{K_values[-1]})',
                fontsize=14, fontweight='bold', y=0.98)
    
    # Global legend
    max_K = max(K_values)
    legend_elements = [mpatches.Patch(facecolor=colors_dict[max_K][i], label=f'Ancestry {i+1}')
                      for i in range(max_K)]
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.87, 0.5),
              frameon=False, fontsize=10, ncol=1)
    
    # Adjust layout
    plt.subplots_adjust(bottom=0.20, right=0.85, hspace=0.10, top=0.92)
    
    # Save PNG and PDF
    for ext in ['png', 'pdf']:
        output_file = f"{output_base}.{ext}"
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=ext)
        print(f"  Saved: {output_file}")
    
    plt.close()

def plot_population_heatmap(Q_matrix, pop_df, K, output_base, pop_order, pop_counts):
    """Create population mean heatmap."""
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Ancestry_{i+1}' for i in range(K)])
    Q_df['POP'] = pop_df['POP'].values
    
    # Calculate population means
    pop_means = Q_df.groupby('POP').mean()
    pop_means = pop_means.reindex(pop_order)
    
    # Save CSV
    pop_means_with_counts = pop_means.copy()
    pop_means_with_counts['Sample_Count'] = pop_counts.reindex(pop_order)
    pop_means_with_counts.to_csv(f"{output_base}.csv")
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(max(8, K*0.6), len(pop_means) * 0.4 + 2))
    
    im = ax.imshow(pop_means.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
    
    ax.set_xticks(np.arange(K))
    ax.set_xticklabels([f'Ancestry {i+1}' for i in range(K)], fontsize=10)
    ax.set_yticks(np.arange(len(pop_means)))
    
    # Use short labels for y-axis
    y_labels = []
    for pop in pop_means.index:
        short_name = POPULATION_LABELS.get(pop, pop)
        count = pop_counts[pop]
        if count < 5:
            short_name += '*'
        y_labels.append(f'{short_name} (n={count})')
    ax.set_yticklabels(y_labels, fontsize=9)
    
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
    
    # Save PNG and PDF
    for ext in ['png', 'pdf']:
        output_file = f"{output_base}.{ext}"
        plt.savefig(output_file, dpi=450, bbox_inches='tight', format=ext)
        print(f"  Saved: {output_file}")
    
    plt.close()

def plot_cv_curve(cv_errors, output_base):
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
    plt.savefig(f"{output_base}.png", dpi=450, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.png")

def main():
    # Set matplotlib defaults
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    # Parse arguments
    recompute_order = parse_arguments()
    
    # Setup paths
    base_dir = Path('ANALYSIS/00-07-ADMIXTURE')
    pop_file = Path('ANALYSIS/00-08-PCA/pop_736.updated.tsv')
    output_dir = base_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("POLISH ADMIXTURE FIGURES WITH SIMILARITY-ORDERED POPULATIONS")
    print("=" * 80)
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
    print("\n1) Loading and joining data...")
    fam_df = pd.read_csv(fam_file, sep=' ', header=None, usecols=[0, 1])
    fam_df.columns = ['FID', 'IID']
    print(f"  FAM: {len(fam_df)} samples")
    
    # Read population labels
    pop_df = parse_population_file(pop_file)
    has_region = 'REGION' in pop_df.columns
    if has_region:
        print(f"  Found REGION column in population file")
    
    # Merge FAM with population data
    merged_df = fam_df.merge(pop_df, on='IID', how='left')
    
    # Check samples
    has_pop = merged_df['POP'].notna()
    n_with_pop = has_pop.sum()
    n_without_pop = (~has_pop).sum()
    
    if n_without_pop > 0:
        print(f"  Dropped {n_without_pop} samples without POP assignment")
    
    if n_with_pop == 0:
        print("ERROR: No samples with population assignments!")
        return
    
    # Filter to samples with population
    merged_df = merged_df[has_pop].copy()
    print(f"  Final sample count: {n_with_pop}")
    
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
        print("ERROR: No Q matrices found!")
        return
    
    print(f"  Detected K values: {available_K}")
    
    # Load Q matrices and CV errors
    print("\n2) Loading ADMIXTURE results and aligning components...")
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
            
            # Parse CV error
            log_file = base_dir / f'admixture_K{K}.log'
            cv_error = parse_cv_error(log_file)
            if cv_error is not None:
                cv_errors[K] = cv_error
    
    print(f"  Loaded {len(Q_matrices_filtered)} Q matrices")
    if cv_errors:
        print(f"  Found CV errors for K={list(cv_errors.keys())}")
    
    # Align ancestry components
    alignments = align_ancestry_components(Q_matrices_filtered, available_K)
    
    # Create color palette
    max_K = max(available_K)
    base_colors = get_okabe_ito_palette(max_K)
    
    # Save palette
    palette_file = output_dir / 'palette_ancestry_colors.tsv'
    with open(palette_file, 'w') as f:
        f.write("Ancestry\tHex_Color\n")
        for i in range(max_K):
            f.write(f"Ancestry_{i+1}\t{base_colors[i]}\n")
    print(f"  Saved color palette: {palette_file}")
    
    # Create aligned color dictionaries
    colors_dict = {}
    for K in available_K:
        alignment = alignments[K]
        colors_dict[K] = [base_colors[alignment[i] if alignment[i] < len(base_colors) else i]
                         for i in range(K)]
    
    # Check for precomputed population order
    print("\n3) Computing population order by similarity...")
    precomputed_order_file = output_dir / 'admixture_population_order.txt'
    
    if precomputed_order_file.exists() and not recompute_order:
        print(f"  Using precomputed order from: {precomputed_order_file}")
        with open(precomputed_order_file, 'r') as f:
            pop_order = [line.strip() for line in f if line.strip()]
        
        # Get population counts
        Q_df_temp = pd.DataFrame(Q_matrices_filtered[available_K[0]])
        Q_df_temp['POP'] = pop_df_dict[available_K[0]]['POP'].values
        pop_counts = Q_df_temp.groupby('POP').size()
        
        # Validate order
        existing_pops = set(pop_counts.index)
        order_pops = set(pop_order)
        if existing_pops != order_pops:
            print(f"  WARNING: Population mismatch, recomputing order...")
            pop_order = None
    else:
        pop_order = None
    
    if pop_order is None:
        # Choose K* for similarity
        if cv_errors:
            K_star = min(cv_errors, key=cv_errors.get)
            print(f"  K* = {K_star} (lowest CV error)")
        else:
            K_star = 5 if 5 in available_K else available_K[len(available_K)//2]
            print(f"  K* = {K_star} (default)")
        
        # Compute order
        pop_order, pop_counts = compute_population_order(
            Q_matrices_filtered[K_star], pop_df_dict[K_star], K_star, output_dir
        )
    
    # Generate figures
    print("\n4-8) Generating polished figures...")
    
    # Individual K plots
    for K in available_K:
        output_base = output_dir / f'admixture_grouped_K{K}_ordered_nice'
        plot_admixture_grouped_nice(Q_matrices_filtered[K], pop_df_dict[K], K,
                                   colors_dict[K], cv_errors.get(K), 
                                   output_base, pop_order, pop_counts)
        
        # Population heatmap
        output_base = output_dir / f'admixture_population_means_K{K}_ordered'
        plot_population_heatmap(Q_matrices_filtered[K], pop_df_dict[K], K,
                               output_base, pop_order, pop_counts)
    
    # Multi-K panel (K=2..8 or available subset)
    panel_K = [K for K in available_K if 2 <= K <= 8]
    if panel_K:
        print(f"\n7) Creating multi-K panel (K={panel_K[0]}..{panel_K[-1]})...")
        output_base = output_dir / f'admixture_grouped_K{panel_K[0]}_K{panel_K[-1]}_ordered_panel_nice'
        plot_multi_panel_nice({K: Q_matrices_filtered[K] for K in panel_K},
                            {K: pop_df_dict[K] for K in panel_K},
                            panel_K, colors_dict, cv_errors, 
                            output_base, pop_order, pop_counts)
    
    # CV curve
    if cv_errors:
        print("\n9) Creating CV error curve...")
        output_base = output_dir / 'admixture_cv_curve'
        plot_cv_curve(cv_errors, output_base)
        
        # Save CV errors
        cv_file = output_dir / 'cv_errors.txt'
        with open(cv_file, 'w') as f:
            f.write("K\tCV_Error\n")
            for K in sorted(cv_errors.keys()):
                f.write(f"{K}\t{cv_errors[K]:.6f}\n")
        print(f"  Saved: {cv_file}")
    
    # Save population counts
    print("\n10) Saving diagnostics...")
    counts_file = output_dir / 'counts_by_population.tsv'
    counts_df = pd.DataFrame({
        'Population': pop_order,
        'N': [pop_counts[pop] for pop in pop_order]
    })
    counts_df.to_csv(counts_file, sep='\t', index=False)
    print(f"  Saved: {counts_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Number of populations: {len(pop_order)}")
    print(f"Distance metric: 1 - cosine similarity")
    if HAS_SCIPY:
        print(f"Clustering method: hierarchical with optimal leaf ordering")
    else:
        print(f"Clustering method: seriation fallback")
    print(f"Samples analyzed: {n_with_pop}")
    if n_without_pop > 0:
        print(f"Samples dropped (no POP): {n_without_pop}")
    print(f"K values processed: {available_K}")
    print(f"Output directory: {output_dir}")
    print("\nAll figures generated successfully!")

if __name__ == "__main__":
    main()