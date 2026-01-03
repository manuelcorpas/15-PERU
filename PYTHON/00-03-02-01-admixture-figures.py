#!/usr/bin/env python3
"""
ADMIXTURE plots with improved labels for Nature Health submission.

Modifications from original 00-03-02-01:
- Meaningful regional ancestry labels for K=5 instead of generic "Ancestry 1-5"
- Legend with frame and title for K=5 figures
- Excluded populations with very small sample sizes: Nahua (n=2), Machiguenga (n=3), Ashaninka (n=1)
- Small Amazonian populations (n<10) have bars displayed but no labels
  (consistent with manuscript stating 28 populations analysed)

Ancestry labels (K=5):
  Component 1: African (orange)
  Component 2: Northern Andean (light blue)
  Component 3: Amazonian Native (teal)
  Component 4: European/Coastal (yellow)
  Component 5: Altiplano Andean (dark blue)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.transforms import Bbox
import os
import re
import glob
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Short label mapping
SHORT_LABEL_MAP = {
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

# Abbreviated labels for very small populations (n<=3) to reduce crowding
ABBREVIATED_LABEL_MAP = {
    "NAHUA": "Nah.",
    "MACHIGUENGA": "Mach.",
    "ASHANINKA": "Ash.",
    "MATZES": "Mats.",
    "CANDOSHI": "Cand.",
    "AWAJUN": "Awaj.",
    "SHIPIBO": "Ship.",
    "LAMAS": "Lam."
}

# Fallback Amazonian populations (case-insensitive)
AMAZONIAN_FALLBACK = {
    'NAHUA', 'MACHIGUENGA', 'ASHANINKA', 'MATSES', 'MATZES', 
    'SHIPIBO', 'CANDOSHI', 'AWAJUN', 'LAMAS', 'IQUITOS'
}

# Populations to exclude from plots entirely (too small for analysis)
EXCLUDED_POPULATIONS = {'NAHUA', 'MACHIGUENGA', 'ASHANINKA'}

# Meaningful ancestry labels for K=5 based on population distribution
# Maps component index to regional/ancestral interpretation
ANCESTRY_LABELS_K5 = {
    0: 'African',               # Orange - Afrodescendientes
    1: 'Northern Andean',       # Light blue - highland populations
    2: 'Amazonian Native',      # Teal - dominant in Amazonian Indigenous groups
    3: 'European/Coastal',      # Yellow - admixed coastal/urban populations
    4: 'Altiplano Andean'       # Dark blue - southern highlands (Puno, Uros, etc.)
}

def get_ancestry_label(k_index, K):
    """Return meaningful ancestry label for K=5, generic for other K values."""
    if K == 5 and k_index in ANCESTRY_LABELS_K5:
        return ANCESTRY_LABELS_K5[k_index]
    return f'Ancestry {k_index + 1}'

# Try scipy for component alignment
try:
    from scipy.optimize import linear_sum_assignment
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("INFO: scipy not found. Using fallback methods.")

def parse_population_file(pop_file):
    """Parse population labels with optional REGION column."""
    if not os.path.exists(pop_file):
        raise FileNotFoundError(f"Population file not found: {pop_file}")
    
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

def parse_cv_from_logs(log_dir):
    """Parse CV errors from ADMIXTURE log files."""
    cv_errors = {}
    log_pattern = os.path.join(log_dir, "admixture_K*.log")
    log_files = glob.glob(log_pattern)
    
    for log_file in log_files:
        # Extract K from filename
        match = re.search(r'admixture_K(\d+)\.log', os.path.basename(log_file))
        if not match:
            continue
        k = int(match.group(1))
        
        # Parse CV error from log
        try:
            with open(log_file, 'r') as f:
                content = f.read()
                # Look for CV error line - standard ADMIXTURE format
                cv_match = re.search(r'CV error \(K=\d+\):\s*([\d.]+)', content)
                if cv_match:
                    cv_errors[k] = float(cv_match.group(1))
        except Exception as e:
            print(f"  Warning: Could not parse {log_file}: {e}")
    
    return cv_errors
def plot_and_save_cv(cv_errors, out_dir, fname="admixture_cv_curve.png"):
    """Plot and save CV error curve."""
    if not cv_errors:
        print("No CV errors found; skipping CV plot.")
        return
    
    # Sort by K
    k_values = sorted(cv_errors.keys())
    cv_values = [cv_errors[k] for k in k_values]
    
    # Find minimum
    min_k = min(cv_errors, key=cv_errors.get)
    min_cv = cv_errors[min_k]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(k_values, cv_values, 'o-', linewidth=2, markersize=8, color='#0072B2')
    
    # Highlight minimum
    ax.plot(min_k, min_cv, 'ro', markersize=12)
    ax.annotate(f'Lowest CV\nK={min_k}\nCV={min_cv:.6f}',
                xy=(min_k, min_cv), xytext=(min_k + 0.5, min_cv),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red')
    
    # Formatting
    ax.set_xlabel('K (Number of Ancestral Populations)', fontsize=12)
    ax.set_ylabel('Cross-Validation Error', fontsize=12)
    ax.set_title('ADMIXTURE Cross-Validation Error', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks(k_values)
    
    plt.tight_layout()
    
    # Save plot
    output_path = os.path.join(out_dir, fname)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved CV curve: {output_path}")
    
    # Save text table
    txt_path = os.path.join(out_dir, "cv_errors.txt")
    with open(txt_path, 'w') as f:
        f.write("K\tCV_error\n")
        for k in k_values:
            f.write(f"{k}\t{cv_errors[k]:.6f}\n")
    print(f"  Saved CV table: {txt_path}")

def get_okabe_ito_palette(n):
    """Okabe-Ito colorblind-safe palette."""
    colors = [
        '#E69F00', '#56B4E9', '#009E73', '#F0E442',
        '#0072B2', '#D55E00', '#CC79A7', '#999999',
        '#000000', '#B19CD9', '#FF6961', '#77DD77'
    ]
    
    if n > len(colors):
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
        
        corr_matrix = np.zeros((n_curr, max(n_prev, n_curr)))
        for j in range(n_curr):
            for k in range(n_prev):
                if alignments[prev_K][k] < corr_matrix.shape[1]:
                    corr = np.corrcoef(Q_curr[:, j], Q_prev[:, k])[0, 1]
                    if not np.isnan(corr):
                        corr_matrix[j, alignments[prev_K][k]] = corr
        
        if HAS_SCIPY:
            row_ind, col_ind = linear_sum_assignment(-corr_matrix[:, :max(n_prev, n_curr)])
            alignment = col_ind[:n_curr].tolist()
        else:
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

def get_block_width_pixels(ax, x_start, x_end, fig_width_inches, dpi=500):
    """Calculate block width in pixels."""
    xlim = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    block_width_data = x_end - x_start
    block_width_fraction = block_width_data / x_range
    width_pixels = block_width_fraction * fig_width_inches * dpi
    return width_pixels

def get_adaptive_fontsize(width_pixels):
    """Get adaptive font size based on block width in pixels."""
    if width_pixels < 22:
        return 7.5
    elif width_pixels < 40:
        return 8.5
    else:
        return 10

def check_overlap(bbox1, bbox2):
    """Check if two bounding boxes overlap."""
    return not (bbox1.x1 < bbox2.x0 or bbox2.x1 < bbox1.x0 or 
                bbox1.y1 < bbox2.y0 or bbox2.y1 < bbox1.y0)

def adjust_overlapping_labels(fig, annotations):
    """Adjust positions of overlapping micro-labels."""
    if not annotations:
        return []
    
    renderer = fig.canvas.get_renderer()
    bboxes = []
    for ann in annotations:
        bbox = ann.get_window_extent(renderer)
        bboxes.append(bbox)
    
    max_iterations = 15  # IMPROVED: More iterations
    y_increment = 0.03   # IMPROVED: Larger vertical adjustment
    max_y = 1.35         # IMPROVED: Allow labels to go higher
    
    for iteration in range(max_iterations):
        overlap_found = False
        
        for i in range(len(annotations)):
            for j in range(i + 1, len(annotations)):
                if check_overlap(bboxes[i], bboxes[j]):
                    overlap_found = True
                    
                    ann_i_y = annotations[i].get_position()[1]
                    ann_j_y = annotations[j].get_position()[1]
                    
                    if ann_i_y >= ann_j_y:
                        new_y = min(ann_i_y + y_increment, max_y)
                        annotations[i].set_position((annotations[i].get_position()[0], new_y))
                        annotations[i].xy = (annotations[i].xy[0], 1.0)
                    else:
                        new_y = min(ann_j_y + y_increment, max_y)
                        annotations[j].set_position((annotations[j].get_position()[0], new_y))
                        annotations[j].xy = (annotations[j].xy[0], 1.0)
                    
                    bboxes[i] = annotations[i].get_window_extent(renderer)
                    bboxes[j] = annotations[j].get_window_extent(renderer)
        
        if not overlap_found:
            break
    
    # IMPROVED: If still overlapping, reduce font size more aggressively
    if overlap_found:
        for i in range(len(annotations)):
            for j in range(i + 1, len(annotations)):
                if check_overlap(bboxes[i], bboxes[j]):
                    current_size_i = annotations[i].get_fontsize()
                    current_size_j = annotations[j].get_fontsize()
                    annotations[i].set_fontsize(max(5.5, current_size_i - 1.5))
                    annotations[j].set_fontsize(max(5.5, current_size_j - 1.5))
    
    return annotations

def plot_admixture_with_fixed_labels(Q_matrix, pop_df, K, colors, cv_error, output_base, 
                                     pop_order, pop_counts, fig_width, dpi=500):
    """Create ADMIXTURE plot with fixed labels for small Amazonian populations."""
    n_samples = len(pop_df)
    
    Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K)])
    Q_df = pd.concat([pop_df.reset_index(drop=True), Q_df], axis=1)
    
    has_region = 'REGION' in pop_df.columns
    
    # Filter out excluded populations (too small for analysis)
    Q_df = Q_df[~Q_df['POP'].str.upper().isin(EXCLUDED_POPULATIONS)].copy()
    
    Q_df['POP'] = pd.Categorical(Q_df['POP'], categories=[p for p in pop_order if p.upper() not in EXCLUDED_POPULATIONS], ordered=True)
    Q_df['max_anc'] = Q_df[[f'Anc_{i+1}' for i in range(K)]].max(axis=1)
    Q_df = Q_df.sort_values(['POP', 'max_anc'], ascending=[True, False])
    
    fig_height = 6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    ind = np.arange(len(Q_df))
    bottom = np.zeros(len(Q_df))
    
    for k in range(K):
        values = Q_df[f'Anc_{k+1}'].values
        ax.bar(ind, values, bottom=bottom, width=1.0,
               color=colors[k], linewidth=0, edgecolor='none')
        bottom += values
    
    amazonian_pops = set()
    if has_region:
        region_map = Q_df.groupby('POP')['REGION'].first()
        amazonian_pops = set(region_map[region_map == 'Amazon'].index)
    else:
        for pop in pop_counts.index:
            if pop.upper() in AMAZONIAN_FALLBACK:
                amazonian_pops.add(pop)
    
    pop_groups = Q_df.groupby('POP', sort=False)
    pop_boundaries = []
    label_info = []
    micro_annotations = []
    
    current_pos = 0
    # IMPROVED: Wider vertical spacing and more stagger levels to prevent overlap
    stagger_cycle = [1.05, 1.12, 1.19, 1.26]
    stagger_idx = 0
    # Track positions of small population labels for horizontal offset
    small_pop_positions = []
    
    for pop_name, group in pop_groups:
        group_size = len(group)
        center = current_pos + group_size / 2
        count = pop_counts[pop_name]
        
        # Use abbreviated name for very small populations (n<=3)
        if count <= 3 and pop_name in ABBREVIATED_LABEL_MAP:
            short_name = ABBREVIATED_LABEL_MAP[pop_name]
        else:
            short_name = SHORT_LABEL_MAP.get(pop_name, pop_name)
        label = f'{short_name}\n(n={count})'
        
        width_pixels = get_block_width_pixels(ax, current_pos, current_pos + group_size, 
                                             fig_width, dpi)
        
        is_small_amazonian = (pop_name in amazonian_pops and count < 10)
        
        if is_small_amazonian:
            # MODIFIED: Skip labels for small Amazonian populations (n<10)
            # These populations are not included in the 28 analysed populations
            # They still appear in the bars but without labels to avoid overlap
            small_pop_positions.append(center)
            stagger_idx += 1
            
            label_info.append({
                'POP': pop_name,
                'n': count,
                'x_center': center,
                'y_final': None,  # No label
                'fontsize': None,
                'stagger_level': 0
            })
        else:
            fontsize = 9 if width_pixels < 30 else (10 if width_pixels < 60 else 11)
            ax.text(center, -0.02, label, ha='center', va='top',
                   rotation=40, fontsize=fontsize, transform=ax.get_xaxis_transform())
            
            label_info.append({
                'POP': pop_name,
                'n': count,
                'x_center': center,
                'y_final': -0.02,
                'fontsize': fontsize,
                'stagger_level': 0
            })
        
        if current_pos > 0:
            pop_boundaries.append(current_pos)
        current_pos += group_size
    
    if micro_annotations:
        micro_annotations = adjust_overlapping_labels(fig, micro_annotations)
        
        for i, ann in enumerate(micro_annotations):
            for info in label_info:
                if info['stagger_level'] > 0:
                    if abs(info['x_center'] - ann.get_position()[0]) < 0.1:
                        info['y_final'] = ann.get_position()[1]
                        info['fontsize'] = ann.get_fontsize()
                        break
    
    for boundary in pop_boundaries:
        ax.axvline(x=boundary, color='black', alpha=0.4, linewidth=0.6, zorder=3)
    
    if has_region:
        region_starts = {}
        region_ends = {}
        current_region = None
        
        for i, (_, row) in enumerate(Q_df.iterrows()):
            if row.get('REGION') != current_region:
                if current_region is not None:
                    region_ends[current_region] = i
                current_region = row.get('REGION')
                region_starts[current_region] = i
        
        if current_region is not None:
            region_ends[current_region] = len(Q_df)
        
        for region in ['Amazon', 'Highlands', 'Coast']:
            if region in region_starts:
                start = region_starts[region]
                end = region_ends[region]
                center = (start + end) / 2
                ax.text(center, 1.01, region, ha='center', va='bottom',
                       fontsize=11, fontweight='bold', transform=ax.get_xaxis_transform())
    
    title = f'K={K}'
    if cv_error is not None:
        title += f' (CV={cv_error:.6f})'
    ax.set_title(title, fontsize=12, loc='left', fontweight='bold')
    ax.set_ylabel('Ancestry Proportion', fontsize=11)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, len(Q_df)])
    ax.set_xticks([])
    
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(['0', '0.5', '1'], fontsize=10)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Create legend with meaningful labels for K=5, generic for other K
    legend_elements = [mpatches.Patch(facecolor=colors[i], 
                                      edgecolor='black' if K == 5 else 'none',
                                      linewidth=0.5,
                                      label=get_ancestry_label(i, K))
                      for i in range(K)]
    legend = ax.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(1.12, 0.5),
                       frameon=True if K == 5 else False, fontsize=10,
                       title='Ancestry' if K == 5 else None,
                       edgecolor='0.8')
    if K == 5:
        legend.get_frame().set_linewidth(0.8)
    
    plt.subplots_adjust(bottom=0.20, right=0.85, top=0.92)
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_base}.{ext}"
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=ext)
        print(f"  Saved: {output_file}")
    
    plt.close()
    
    return label_info

def plot_multi_panel_fixed(Q_matrices, pop_df_dict, K_values, colors_dict, cv_errors, 
                          output_base, pop_order, pop_counts, fig_width, dpi=500):
    """Create multi-panel plot with fixed labels."""
    n_panels = len(K_values)
    n_samples = len(pop_df_dict[K_values[0]])
    
    fig_height = 1.6 * n_panels + 1.5
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    all_label_info = {}
    
    for idx, K in enumerate(K_values):
        ax = plt.subplot(n_panels, 1, idx + 1)
        
        Q_matrix = Q_matrices[K]
        pop_df = pop_df_dict[K]
        colors = colors_dict[K]
        cv_error = cv_errors.get(K)
        
        Q_df = pd.DataFrame(Q_matrix, columns=[f'Anc_{i+1}' for i in range(K)])
        Q_df = pd.concat([pop_df.reset_index(drop=True), Q_df], axis=1)
        
        has_region = 'REGION' in pop_df.columns
        
        # Filter out excluded populations (too small for analysis)
        Q_df = Q_df[~Q_df['POP'].str.upper().isin(EXCLUDED_POPULATIONS)].copy()
        
        Q_df['POP'] = pd.Categorical(Q_df['POP'], categories=[p for p in pop_order if p.upper() not in EXCLUDED_POPULATIONS], ordered=True)
        Q_df['max_anc'] = Q_df[[f'Anc_{i+1}' for i in range(K)]].max(axis=1)
        Q_df = Q_df.sort_values(['POP', 'max_anc'], ascending=[True, False])
        
        ind = np.arange(len(Q_df))
        bottom = np.zeros(len(Q_df))
        
        for k in range(K):
            values = Q_df[f'Anc_{k+1}'].values
            ax.bar(ind, values, bottom=bottom, width=1.0,
                   color=colors[k], linewidth=0, edgecolor='none')
            bottom += values
        
        amazonian_pops = set()
        if has_region:
            region_map = Q_df.groupby('POP')['REGION'].first()
            amazonian_pops = set(region_map[region_map == 'Amazon'].index)
        else:
            for pop in pop_counts.index:
                if pop.upper() in AMAZONIAN_FALLBACK:
                    amazonian_pops.add(pop)
        
        if idx == n_panels - 1:
            pop_groups = Q_df.groupby('POP', sort=False)
            pop_boundaries = []
            label_info = []
            micro_annotations = []
            
            current_pos = 0
            # IMPROVED: Wider vertical spacing and more stagger levels
            stagger_cycle = [1.05, 1.12, 1.19, 1.26]
            stagger_idx = 0
            small_pop_positions = []
            
            for pop_name, group in pop_groups:
                group_size = len(group)
                center = current_pos + group_size / 2
                count = pop_counts[pop_name]
                
                # Use abbreviated name for very small populations (n<=3)
                if count <= 3 and pop_name in ABBREVIATED_LABEL_MAP:
                    short_name = ABBREVIATED_LABEL_MAP[pop_name]
                else:
                    short_name = SHORT_LABEL_MAP.get(pop_name, pop_name)
                label = f'{short_name}\n(n={count})'
                
                width_pixels = get_block_width_pixels(ax, current_pos, current_pos + group_size,
                                                     fig_width, dpi)
                
                is_small_amazonian = (pop_name in amazonian_pops and count < 10)
                
                if is_small_amazonian:
                    # MODIFIED: Skip labels for small Amazonian populations (n<10)
                    # These populations are not included in the 28 analysed populations
                    # They still appear in the bars but without labels to avoid overlap
                    small_pop_positions.append(center)
                    stagger_idx += 1
                    
                    label_info.append({
                        'POP': pop_name,
                        'n': count,
                        'x_center': center,
                        'y_final': None,  # No label
                        'fontsize': None,
                        'stagger_level': 0
                    })
                else:
                    fontsize = 9 if width_pixels < 30 else (10 if width_pixels < 60 else 11)
                    ax.text(center, -0.02, label, ha='center', va='top',
                           rotation=40, fontsize=fontsize, transform=ax.get_xaxis_transform())
                    
                    label_info.append({
                        'POP': pop_name,
                        'n': count,
                        'x_center': center,
                        'y_final': -0.02,
                        'fontsize': fontsize,
                        'stagger_level': 0
                    })
                
                if current_pos > 0:
                    pop_boundaries.append(current_pos)
                current_pos += group_size
            
            if micro_annotations:
                micro_annotations = adjust_overlapping_labels(fig, micro_annotations)
                
                for i, ann in enumerate(micro_annotations):
                    for info in label_info:
                        if info['stagger_level'] > 0:
                            if abs(info['x_center'] - ann.get_position()[0]) < 0.1:
                                info['y_final'] = ann.get_position()[1]
                                info['fontsize'] = ann.get_fontsize()
                                break
            
            all_label_info[K] = label_info
        else:
            pop_groups = Q_df.groupby('POP', sort=False)
            pop_boundaries = []
            current_pos = 0
            for pop_name, group in pop_groups:
                if current_pos > 0:
                    pop_boundaries.append(current_pos)
                current_pos += len(group)
        
        for boundary in pop_boundaries:
            ax.axvline(x=boundary, color='black', alpha=0.4, linewidth=0.6, zorder=3)
        
        title = f'K={K}'
        if cv_error is not None:
            title += f' (CV={cv_error:.6f})'
        ax.set_title(title, fontsize=11, loc='left', pad=2, fontweight='bold')
        
        ax.set_ylim([0, 1])
        ax.set_xlim([0, len(Q_df)])
        ax.set_ylabel('', fontsize=9)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(['0', '0.5', '1'], fontsize=9)
        ax.set_xticks([])
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    fig.suptitle(f'ADMIXTURE (n={n_samples}; K={K_values[0]}..{K_values[-1]})',
                fontsize=14, fontweight='bold', y=0.98)
    
    # Create legend with meaningful labels for max K (use K=5 labels if max_K=5)
    max_K = max(K_values)
    legend_elements = [mpatches.Patch(facecolor=colors_dict[max_K][i], 
                                      edgecolor='black' if max_K == 5 else 'none',
                                      linewidth=0.5,
                                      label=get_ancestry_label(i, max_K))
                      for i in range(max_K)]
    fig.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(0.985, 0.5),
              frameon=True if max_K == 5 else False, fontsize=10,
              title='Ancestry' if max_K == 5 else None,
              edgecolor='0.8')
    
    plt.subplots_adjust(bottom=0.20, right=0.85, top=0.92, hspace=0.12)
    
    for ext in ['png', 'pdf']:
        output_file = f"{output_base}.{ext}"
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=ext)
        print(f"  Saved: {output_file}")
    
    plt.close()
    
    return all_label_info

def main():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    
    base_dir = Path('ANALYSIS/00-07-ADMIXTURE')
    pop_file = Path('ANALYSIS/00-08-PCA/pop_736.updated.tsv')
    order_file = base_dir / 'admixture_population_order.txt'
    output_dir = base_dir
    
    print("=" * 80)
    print("FIX OVERLAPPING LABELS FOR SMALL AMAZONIAN POPULATIONS")
    print("=" * 80)
    print(f"ADMIXTURE directory: {base_dir}")
    print(f"Population file: {pop_file}")
    print(f"Order file: {order_file}")
    
    if not order_file.exists():
        print("ERROR: Population order file not found!")
        return
    
    with open(order_file, 'r') as f:
        pop_order = [line.strip() for line in f if line.strip()]
    print(f"Loaded population order: {len(pop_order)} populations")
    
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
    
    print("\nLoading data...")
    fam_df = pd.read_csv(fam_file, sep=' ', header=None, usecols=[0, 1])
    fam_df.columns = ['FID', 'IID']
    print(f"  FAM: {len(fam_df)} samples")
    
    pop_df = parse_population_file(pop_file)
    has_region = 'REGION' in pop_df.columns
    
    merged_df = fam_df.merge(pop_df, on='IID', how='left')
    has_pop = merged_df['POP'].notna()
    n_with_pop = has_pop.sum()
    
    if n_with_pop == 0:
        print("ERROR: No samples with population assignments!")
        return
    
    merged_df = merged_df[has_pop].copy()
    print(f"  Samples with population: {n_with_pop}")
    
    pop_counts = merged_df.groupby('POP').size()
    
    amazonian_pops = set()
    if has_region:
        region_map = merged_df.groupby('POP')['REGION'].first()
        amazonian_pops = set(region_map[region_map == 'Amazon'].index)
        print(f"  Amazonian populations (from REGION): {amazonian_pops}")
    else:
        for pop in pop_counts.index:
            if pop.upper() in AMAZONIAN_FALLBACK:
                amazonian_pops.add(pop)
        print(f"  Amazonian populations (fallback): {amazonian_pops}")
    
    small_amazonian = [pop for pop in amazonian_pops if pop_counts[pop] < 10]
    print(f"  Small (n<10) Amazonian populations: {small_amazonian}")
    
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
    
    print(f"  Available K values: {available_K}")
    
    print("\nLoading ADMIXTURE results...")
    Q_matrices = {}
    Q_matrices_filtered = {}
    pop_df_dict = {}
    
    # Parse CV errors from logs
    cv_errors = parse_cv_from_logs(str(base_dir))
    
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
    
    print("Aligning ancestry components...")
    alignments = align_ancestry_components(Q_matrices_filtered, available_K)
    
    max_K = max(available_K)
    base_colors = get_okabe_ito_palette(max_K)
    
    colors_dict = {}
    for K in available_K:
        alignment = alignments[K]
        colors_dict[K] = [base_colors[alignment[i] if alignment[i] < len(base_colors) else i]
                         for i in range(K)]
    
    fig_width = max(26, n_with_pop / 24)
    print(f"  Figure width: {fig_width:.1f} inches")
    
    print("\nGenerating plots with fixed labels...")
    
    all_label_info = {}
    
    for K in available_K:
        output_base = output_dir / f'admixture_grouped_K{K}_ordered_smallfix'
        label_info = plot_admixture_with_fixed_labels(
            Q_matrices_filtered[K], pop_df_dict[K], K,
            colors_dict[K], cv_errors.get(K),
            output_base, pop_order, pop_counts, fig_width
        )
        all_label_info[f'K{K}'] = label_info
    
    panel_K = [K for K in available_K if 2 <= K <= 8]
    if panel_K:
        print(f"\nCreating multi-K panel (K={panel_K[0]}..{panel_K[-1]})...")
        output_base = output_dir / f'admixture_grouped_K{panel_K[0]}_K{panel_K[-1]}_ordered_panel_smallfix'
        panel_label_info = plot_multi_panel_fixed(
            {K: Q_matrices_filtered[K] for K in panel_K},
            {K: pop_df_dict[K] for K in panel_K},
            panel_K, colors_dict, cv_errors,
            output_base, pop_order, pop_counts, fig_width
        )
        all_label_info.update(panel_label_info)
    
    if all_label_info:
        label_df_list = []
        for k_name, labels in all_label_info.items():
            if labels:
                df = pd.DataFrame(labels)
                df['K'] = k_name
                label_df_list.append(df)
        
        if label_df_list:
            final_label_df = pd.concat(label_df_list, ignore_index=True)
            label_file = output_dir / 'labels_smallpops.csv'
            final_label_df.to_csv(label_file, index=False)
            print(f"\nSaved label information: {label_file}")
    
    # Plot CV curve
    print("\nGenerating CV error curve...")
    plot_and_save_cv(cv_errors, str(output_dir))
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Excluded populations (too small): {EXCLUDED_POPULATIONS}")
    print(f"Small Amazonian populations (n<10, unlabelled): {len([p for p in small_amazonian if p.upper() not in EXCLUDED_POPULATIONS])}")
    print(f"Total populations in data: {len(pop_order)}")
    print(f"Populations in figures: {len(pop_order) - len([p for p in pop_order if p.upper() in EXCLUDED_POPULATIONS])}")
    print(f"K values processed: {available_K}")
    if cv_errors:
        optimal_k = min(cv_errors, key=cv_errors.get)
        print(f"Optimal K: {optimal_k} (CV={cv_errors[optimal_k]:.6f})")
    print(f"Output directory: {output_dir}")
    print("\nAncestry labels used for K=5:")
    for i, label in ANCESTRY_LABELS_K5.items():
        print(f"  Component {i+1}: {label}")
    print("\nFigures generated successfully!")

if __name__ == "__main__":
    main()