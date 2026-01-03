#!/usr/bin/env python3
"""
Complete PCA analysis and Figure 3 generation - REVIEWER-REVISED VERSION
Changes from previous version:
- Panel D: SGDP samples now use triangles (distinct from Peru circles)
- Added clear legend distinguishing Peru vs SGDP marker styles
- Fixed sample count (736, not 746)
- Expanded regional/language classifications
- Improved population group assignments
All outputs go to ANALYSIS/00-08-PCA/
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse, Patch
from matplotlib.lines import Line2D
from adjustText import adjust_text
import re

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base directory
BASE_DIR = 'ANALYSIS/00-08-PCA'

# Input VCF files
PERU_VCF = os.path.join(BASE_DIR, 'peru.736.vcf.gz')
COMMON_VCF = os.path.join(BASE_DIR, 'common_variants_ibd_clean.vcf.gz')
SGDP_VCF = os.path.join(BASE_DIR, 'common_peru_sgdp.vcf.gz')

# Alternative VCF paths if above don't exist
ALT_PERU_VCF = 'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz'
ALT_SGDP_VCF = 'ANALYSIS/12-SGDP/common_variants_peru_sgdp.vcf.gz'
ALT_SGDP_VCF2 = 'ANALYSIS/12-SGDP/common_peru_sgdp.vcf.gz'

# Population files
POP_FILE = os.path.join(BASE_DIR, 'pop_736.updated.tsv')
ALT_POP_FILE = os.path.join(BASE_DIR, 'ii_28_populations.txt')
ALT_POP_FILE2 = 'ANALYSIS/10-PCA/ii_28_populations.txt'

# SGDP population file
SGDP_POP_FILE = 'ANALYSIS/00-08-PCA/sgdp-peru-sample-ids-population-ids-24-SDGP-24-close-pops.txt'
ALT_SGDP_POP_FILE = 'ANALYSIS/12-SGDP/sgdp-peru-sample-ids-population-ids-24-SDGP-24-close-pops.txt'

# Output PCA prefixes
PCA_PERU_PREFIX = os.path.join(BASE_DIR, 'pca_peru')
PCA_COMMON_PREFIX = os.path.join(BASE_DIR, 'pca_common')
PCA_SGDP_PREFIX = os.path.join(BASE_DIR, 'pca_sgdp')

# Alternative SGDP PCA path
ALT_SGDP_PCA = 'ANALYSIS/12-SGDP/pca_output'

# Peru populations (for filtering)
PERU_POPULATIONS = {
    'ANCASH', 'AREQUIPA', 'ASHANINKA_INS', 'AWAJUN', 'AYACUCHO',
    'CANDOSHI', 'CHACHAPOYAS', 'CHOPCCAS', 'CUSCO', 'IQUITOS',
    'JACARUS', 'JAQARUS', 'LAMAS', 'LAMBAYEQUE', 'LIMA', 'MATZES', 'MATSIGUENKAS',
    'MOCHES', 'MOQUEGUA', 'NAHUA', 'PUNO', 'QEROS', 'SHIPIBO_INS',
    'TACNA', 'TALLANES', 'TRUJILLO', 'TUMBES', 'UROS', 'HUARAZ',
    'AFRODESCENDENTS'
}

# ============================================================================
# REGIONAL AND LANGUAGE CLASSIFICATION (EXPANDED per reviewer feedback)
# ============================================================================

# Regional classification - EXPANDED
HIGHLAND_POPULATIONS = {
    'CHOPCCAS', 'CUSCO', 'AYACUCHO', 'QEROS', 'PUNO', 'UROS', 
    'ANCASH', 'CHACHAPOYAS', 'AREQUIPA', 'MOQUEGUA', 'TACNA', 
    'HUARAZ', 'JAQARUS', 'JACARUS'
}

AMAZONIAN_POPULATIONS = {
    'MATZES', 'AWAJUN', 'CANDOSHI', 'SHIPIBO_INS', 'ASHANINKA_INS',
    'MATSIGUENKAS', 'NAHUA', 'SHIPIBO', 'ASHANINKA', 'MACHIGUENGA',
    'IQUITOS'  # Added: Amazonian mestizo city
}

COASTAL_POPULATIONS = {
    'LIMA', 'TRUJILLO', 'MOCHES', 'LAMBAYEQUE', 'TUMBES', 
    'TALLANES', 'LAMAS', 'TALLAN', 'AFRODESCENDENTS'
}

# Language groups - EXPANDED (Jaqaru is an Aymaran language)
QUECHUA_POPULATIONS = {
    'CHOPCCAS', 'CUSCO', 'AYACUCHO', 'QEROS', 'ANCASH', 'CHACHAPOYAS',
    'HUARAZ'
}

AYMARA_POPULATIONS = {
    'UROS', 'PUNO', 'JAQARUS', 'JACARUS'  # Jaqaru is Aymaran family
}

# Population name standardization mapping
POPULATION_NAME_MAPPING = {
    'AFRODESCENDIENTES': 'AFRODESCENDENTS',
    'Afrodescendientes': 'AFRODESCENDENTS',
    'afrodescendientes': 'AFRODESCENDENTS',
}

def normalize_population_name(pop_name):
    """Normalize population names to standard format"""
    if pop_name is None or pd.isna(pop_name):
        return pop_name
    pop_str = str(pop_name).strip()
    # Check mapping first
    if pop_str in POPULATION_NAME_MAPPING:
        return POPULATION_NAME_MAPPING[pop_str]
    # Also check uppercase version
    if pop_str.upper() in POPULATION_NAME_MAPPING:
        return POPULATION_NAME_MAPPING[pop_str.upper()]
    return pop_str

def normalize_population_column(df, col='POP'):
    """Apply population name normalization to a DataFrame column"""
    if col in df.columns:
        df[col] = df[col].apply(normalize_population_name)
    return df

def get_region(population):
    """Get region for a population"""
    pop_upper = population.upper() if population else ''
    if pop_upper in {p.upper() for p in HIGHLAND_POPULATIONS}:
        return 'Highland'
    elif pop_upper in {p.upper() for p in AMAZONIAN_POPULATIONS}:
        return 'Amazonian'
    elif pop_upper in {p.upper() for p in COASTAL_POPULATIONS}:
        return 'Coastal'
    else:
        return 'Other'

def get_language_group(population):
    """Get language group for a population"""
    pop_upper = population.upper() if population else ''
    if pop_upper in {p.upper() for p in QUECHUA_POPULATIONS}:
        return 'Quechua'
    elif pop_upper in {p.upper() for p in AYMARA_POPULATIONS}:
        return 'Aymara'
    else:
        return 'Other'

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def check_plink():
    """Check if PLINK is available"""
    try:
        result = subprocess.run(['plink', '--version'], capture_output=True, text=True)
        print(f"✓ PLINK found")
        return True
    except:
        print("ERROR: PLINK not found in PATH. Please install PLINK.")
        return False

def find_file(primary, *alternatives):
    """Find first existing file from a list"""
    if os.path.exists(primary):
        return primary
    for alt in alternatives:
        if os.path.exists(alt):
            print(f"  Using alternate file: {alt}")
            return alt
    return None

def run_pca_if_needed(vcf_file, output_prefix):
    """Run PLINK PCA if eigenvec doesn't exist"""
    eigenvec_file = f"{output_prefix}.eigenvec"
    
    if os.path.exists(eigenvec_file):
        print(f"  ✓ Found existing: {os.path.basename(eigenvec_file)}")
        return True
    
    if not vcf_file or not os.path.exists(vcf_file):
        print(f"  ✗ VCF not found: {vcf_file}")
        return False
    
    print(f"  Running PCA on: {os.path.basename(vcf_file)}")
    cmd = f"plink --vcf {vcf_file} --pca 10 --double-id --out {output_prefix}"
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if os.path.exists(eigenvec_file):
            print(f"  ✓ PCA complete: {os.path.basename(eigenvec_file)}")
            return True
        else:
            print(f"  ✗ PCA failed - check PLINK output")
            return False
    except Exception as e:
        print(f"  ✗ Error running PCA: {e}")
        return False

def load_eigenvec(filepath):
    """Load eigenvec file from PLINK"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Eigenvec file not found: {filepath}")
    
    df = pd.read_csv(filepath, delim_whitespace=True, header=None)
    n_cols = df.shape[1]
    
    df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, n_cols-1)]
    df['IID'] = df['IID'].astype(str).str.strip()
    
    for i in range(n_cols-1, 11):
        df[f'PC{i}'] = 0
    
    return df

def load_eigenval(filepath):
    """Load eigenvalues and calculate variance explained percentages"""
    eigenval_file = filepath.replace('.eigenvec', '.eigenval')
    
    if not os.path.exists(eigenval_file):
        print(f"  Warning: Eigenvalue file not found: {eigenval_file}")
        return None
    
    try:
        eigenvalues = pd.read_csv(eigenval_file, header=None, names=['eigenval'])
        total_var = eigenvalues['eigenval'].sum()
        var_explained = (eigenvalues['eigenval'] / total_var * 100).values
        print(f"  Variance explained: PC1={var_explained[0]:.2f}%, PC2={var_explained[1]:.2f}%, PC3={var_explained[2]:.2f}%")
        return var_explained
    except Exception as e:
        print(f"  Warning: Could not load eigenvalues: {e}")
        return None

def load_population_file(filepath):
    """Load population assignments - HARDENED VERSION"""
    if not os.path.exists(filepath):
        print(f"Warning: Population file not found: {filepath}")
        return pd.DataFrame(columns=['FID', 'IID', 'POP'])
    
    try:
        df = pd.read_csv(filepath, delimiter='\t')
    except:
        try:
            df = pd.read_csv(filepath, sep=r'\s+')
        except:
            print(f"Error reading population file: {filepath}")
            return pd.DataFrame(columns=['FID', 'IID', 'POP'])
    
    print(f"    Detected columns in {os.path.basename(filepath)}: {df.columns.tolist()}")
    
    df.columns = [col.strip() for col in df.columns]
    
    pop_col = None
    for col_name in ['Population', 'population', 'POP', 'pop']:
        if col_name in df.columns:
            pop_col = col_name
            break
    
    if pop_col:
        df['POP'] = df[pop_col].astype(str).str.strip()
    else:
        print(f"    WARNING: No population column found in {filepath}")
        print(f"    Available columns: {df.columns.tolist()}")
    
    iid_col = None
    for col_name in ['IID', 'iid', 'Sample', 'sample', 'ID', 'id']:
        if col_name in df.columns:
            iid_col = col_name
            break
    
    if iid_col and iid_col != 'IID':
        df['IID'] = df[iid_col].astype(str).str.strip()
    elif 'IID' in df.columns:
        df['IID'] = df['IID'].astype(str).str.strip()
    else:
        print(f"    WARNING: No IID column found in {filepath}")
    
    if 'FID' not in df.columns:
        if 'fid' in df.columns:
            df['FID'] = df['fid'].astype(str).str.strip()
        elif 'IID' in df.columns:
            df['FID'] = df['IID']
        else:
            df['FID'] = ''
    else:
        df['FID'] = df['FID'].astype(str).str.strip()
    
    if 'IID' in df.columns and 'POP' in df.columns:
        print(f"    Sample IIDs: {df['IID'].head(3).tolist()}")
        print(f"    Sample POPs: {df['POP'].head(3).tolist()}")
    
    return df

def parse_population_from_iid(iid):
    """Parse population name from SGDP IID - COMPLETE VERSION"""
    patterns_with_groups = [
        (r'^(.+)_HGDP\d+$', 1),
        (r'^(.+)_NA\d+$', 1),
        (r'^(.+)_HG\d+$', 1),
        (r'^(.+)_[A-Z]{2,4}\d+$', 1),
        (r'^(.+)_[A-Z][a-z]+\d+$', 1),
        (r'^(.+)_SS\d+$', 1),
        (r'^(.+)_Ayodo_\d+[A-Z]?$', 1),
        (r'^(.+)_[a-z]+\d+$', 1),
        (r'^(.+)_[A-Z]\d+$', 1),
        (r'^(.+)_[A-Z]\-\d+$', 1),
        (r'^(.+)_[a-z]+\d+_[a-z]+$', 1),
    ]
    
    for pattern, group in patterns_with_groups:
        match = re.match(pattern, iid)
        if match:
            population = match.group(group)
            population = re.sub(r'_(South|North|Jew|Pandit|Ayodo)$', r'_\1', population)
            if population and len(population) > 2:
                return population
    
    if '_' in iid:
        parts = iid.split('_')
        first_part = parts[0]
        if first_part and len(first_part) >= 2 and not first_part.isdigit():
            return first_part
    
    simple_match = re.match(r'^([A-Za-z_]+?)_?[A-Z]*\d+[A-Z]?$', iid)
    if simple_match:
        pop = simple_match.group(1)
        if pop and len(pop) > 2 and not pop.isdigit():
            return pop
    
    return None

# ============================================================================
# COLOR PALETTE
# ============================================================================

def get_population_colors(n_pops):
    """Generate distinct colors for populations"""
    if n_pops <= 20:
        colors = list(plt.cm.tab20.colors)
    else:
        colors = (list(plt.cm.tab20.colors) + 
                 list(plt.cm.tab20b.colors) + 
                 list(plt.cm.tab20c.colors))
    
    return colors[:n_pops]

# ============================================================================
# ENHANCED SGDP MATCHING FUNCTION
# ============================================================================

def match_sgdp_populations(sgdp_data, sgdp_pop_file):
    """Enhanced SGDP population matching with robust fallback"""
    print("\n  === ENHANCED SGDP POPULATION MATCHING ===")
    
    diagnostics = {
        'total_sgdp': len(sgdp_data),
        'direct_matches': 0,
        'fallback_matches': 0,
        'unmatched': 0,
        'unique_pops': set(),
        'pop_counts': {},
        'match_details': []
    }
    
    if sgdp_pop_file and os.path.exists(sgdp_pop_file):
        print(f"  Loading SGDP population file: {sgdp_pop_file}")
        sgdp_pop = load_population_file(sgdp_pop_file)
        
        if not sgdp_pop.empty and 'IID' in sgdp_pop.columns and 'POP' in sgdp_pop.columns:
            sgdp_pop['IID'] = sgdp_pop['IID'].astype(str).str.strip()
            sgdp_pop['POP'] = sgdp_pop['POP'].astype(str).str.strip()
            
            sgdp_pop_only = sgdp_pop[~sgdp_pop['POP'].isin(PERU_POPULATIONS)].copy()
            print(f"    Filtered to {len(sgdp_pop_only)} SGDP entries (removed Peru)")
            
            sgdp_data = sgdp_data.merge(
                sgdp_pop_only[['IID', 'POP']], 
                on='IID', 
                how='left'
            )
            
            diagnostics['direct_matches'] = sgdp_data['POP'].notna().sum()
            print(f"    Direct matches: {diagnostics['direct_matches']}/{len(sgdp_data)}")
            
            for idx, row in sgdp_data[sgdp_data['POP'].notna()].iterrows():
                diagnostics['match_details'].append({
                    'IID': row['IID'],
                    'POP': row['POP'],
                    'source': 'mapping_file'
                })
        else:
            print("    WARNING: Population file missing required columns")
            sgdp_data['POP'] = None
    else:
        print("    WARNING: No SGDP population file found")
        sgdp_data['POP'] = None
    
    unmatched_mask = sgdp_data['POP'].isna()
    n_unmatched = unmatched_mask.sum()
    
    if n_unmatched > 0:
        print(f"\n  Applying fallback parsing for {n_unmatched} unmatched samples...")
        
        fallback_pops = []
        parse_examples = []
        
        for idx, row in sgdp_data[unmatched_mask].iterrows():
            iid = row['IID']
            parsed_pop = parse_population_from_iid(iid)
            
            if parsed_pop and parsed_pop not in PERU_POPULATIONS:
                sgdp_data.loc[idx, 'POP'] = parsed_pop
                fallback_pops.append(parsed_pop)
                
                diagnostics['match_details'].append({
                    'IID': iid,
                    'POP': parsed_pop,
                    'source': 'parsed_from_iid'
                })
                
                if len(parse_examples) < 5:
                    parse_examples.append(f"      {iid} → {parsed_pop}")
        
        diagnostics['fallback_matches'] = len(fallback_pops)
        print(f"    Fallback matches: {diagnostics['fallback_matches']}")
        if parse_examples:
            print("    Example parsings:")
            for example in parse_examples:
                print(example)
    
    diagnostics['unmatched'] = sgdp_data['POP'].isna().sum()
    diagnostics['unique_pops'] = set(sgdp_data['POP'].dropna().unique())
    
    pop_counts = sgdp_data['POP'].value_counts()
    diagnostics['pop_counts'] = pop_counts.to_dict()
    
    total_matched = diagnostics['direct_matches'] + diagnostics['fallback_matches']
    match_rate = 100 * total_matched / diagnostics['total_sgdp'] if diagnostics['total_sgdp'] > 0 else 0
    
    print(f"\n  === MATCHING SUMMARY ===")
    print(f"    Direct matches: {diagnostics['direct_matches']}")
    print(f"    Fallback matches: {diagnostics['fallback_matches']}")
    print(f"    Total matched: {total_matched}/{diagnostics['total_sgdp']} ({match_rate:.1f}%)")
    print(f"    Unique SGDP populations: {len(diagnostics['unique_pops'])}")
    
    if diagnostics['pop_counts']:
        print(f"\n    Top 10 SGDP populations by sample count:")
        for pop, count in list(diagnostics['pop_counts'].items())[:10]:
            print(f"      {pop}: {count} samples")
    
    write_diagnostic_files(diagnostics, sgdp_data)
    
    if diagnostics['unmatched'] > 0:
        print(f"\n    WARNING: {diagnostics['unmatched']} samples still unmatched")
        unmatched_iids = sgdp_data[sgdp_data['POP'].isna()]['IID'].head(50).tolist()
        print(f"    First unmatched IIDs (up to 50):")
        for iid in unmatched_iids[:10]:
            print(f"      {iid}")
        print(f"    (See sgdp_unmatched.tsv for full list)")
    
    return sgdp_data, diagnostics

def write_diagnostic_files(diagnostics, sgdp_data):
    """Write diagnostic TSV files"""
    output_dir = BASE_DIR
    
    matched_file = os.path.join(output_dir, 'sgdp_iid_to_pop_resolved.tsv')
    if diagnostics['match_details']:
        matched_df = pd.DataFrame(diagnostics['match_details'])
        matched_df.to_csv(matched_file, sep='\t', index=False)
        print(f"    ✓ Wrote: {matched_file}")
    
    if diagnostics['unmatched'] > 0:
        unmatched_file = os.path.join(output_dir, 'sgdp_unmatched.tsv')
        unmatched_df = sgdp_data[sgdp_data['POP'].isna()][['IID']].copy()
        unmatched_df.to_csv(unmatched_file, sep='\t', index=False)
        print(f"    ✓ Wrote: {unmatched_file}")

# ============================================================================
# CIRCLE DRAWING FUNCTIONS
# ============================================================================

def draw_confidence_ellipse(ax, data, color, label, alpha=0.15, edgecolor=None, linewidth=2):
    """Draw a confidence ellipse around a cluster of points"""
    if len(data) < 3:
        return
    
    x = data[:, 0]
    y = data[:, 1]
    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    
    width, height = 2.5 * np.sqrt(eigenvalues)
    
    if edgecolor is None:
        edgecolor = color
    
    ellipse = Ellipse((mean_x, mean_y), width, height, angle=angle,
                     facecolor=color, alpha=alpha, 
                     edgecolor=edgecolor, linewidth=linewidth,
                     label=label, zorder=1)
    ax.add_patch(ellipse)

def draw_confidence_ellipse_custom(ax, data, color, label, alpha=0.15, edgecolor=None, linewidth=2, std_multiplier=2.5):
    """Draw a confidence ellipse around a cluster of points with custom std multiplier"""
    if len(data) < 3:
        return
    
    x = data[:, 0]
    y = data[:, 1]
    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    
    width, height = std_multiplier * np.sqrt(eigenvalues)
    
    if edgecolor is None:
        edgecolor = color
    
    ellipse = Ellipse((mean_x, mean_y), width, height, angle=angle,
                     facecolor=color, alpha=alpha, 
                     edgecolor=edgecolor, linewidth=linewidth,
                     label=label, zorder=1)
    ax.add_patch(ellipse)

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_panel_a(ax, data_dict, var_explained, pop_colors):
    """Panel A: PC1 vs PC2 - Main population structure"""
    merged = data_dict['peru_data']
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel A: Data not available', 
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    pops = sorted(merged['POP'].dropna().unique())
    
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        color = pop_colors.get(pop, 'gray')
        ax.scatter(subset['PC1'], subset['PC2'],
                  label=pop,
                  color=color,
                  alpha=0.7, s=30)
    
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC1'].mean()
            cy = subset['PC2'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold'))
    
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    # Add variance explained to axis labels
    pc1_var = f"{var_explained[0]:.2f}%" if var_explained is not None else ""
    pc2_var = f"{var_explained[1]:.2f}%" if var_explained is not None else ""
    
    ax.set_xlabel(f'PC1 ({pc1_var})' if pc1_var else 'PC1', fontsize=11)
    ax.set_ylabel(f'PC2 ({pc2_var})' if pc2_var else 'PC2', fontsize=11)
    ax.set_title('PCA of ARRAY & WGS Samples by Population', fontsize=12, pad=15)
    
    if len(pops) <= 30:
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                 frameon=False, fontsize=7, ncol=1)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_panel_b(ax, data_dict, var_explained, pop_colors):
    """Panel B: PC3 vs PC2 with regional circles - Highland/Amazonian/Coastal"""
    merged = data_dict.get('common_data', data_dict['peru_data'])
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel B: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Add region classification
    merged['Region'] = merged['POP'].apply(get_region)
    
    # Debug: print population counts by region
    print(f"\n  Panel B Regional Classification:")
    for region in ['Highland', 'Amazonian', 'Coastal', 'Other']:
        pops = merged[merged['Region'] == region]['POP'].unique()
        if len(pops) > 0:
            print(f"    {region}: {sorted(pops)}")
    
    # Define region colors for circles
    region_colors = {
        'Highland': '#e74c3c',      # Red
        'Amazonian': '#27ae60',     # Green
        'Coastal': '#3498db'        # Blue
    }
    
    # Use tight ellipses to avoid overlap
    std_multipliers = {
        'Highland': 1.9,
        'Amazonian': 2.0,
        'Coastal': 1.8
    }
    
    for region, color in region_colors.items():
        region_data = merged[merged['Region'] == region][['PC3', 'PC2']].values
        if len(region_data) >= 3:
            std_mult = std_multipliers.get(region, 1.9)
            draw_confidence_ellipse_custom(ax, region_data, color, region, 
                                  alpha=0.18, linewidth=3.0, std_multiplier=std_mult)
    
    # Plot populations using GLOBAL colors
    pops = sorted(merged['POP'].dropna().unique())
    
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        color = pop_colors.get(pop, 'gray')
        ax.scatter(subset['PC3'], subset['PC2'],
                  label=pop,
                  color=color,
                  alpha=0.8, s=35, zorder=5)
    
    # Add population labels
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC3'].mean()
            cy = subset['PC2'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold', zorder=10))
    
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    # Add variance explained
    pc3_var = f"{var_explained[2]:.2f}%" if var_explained is not None else ""
    pc2_var = f"{var_explained[1]:.2f}%" if var_explained is not None else ""
    
    ax.set_xlabel(f'PC3 ({pc3_var})' if pc3_var else 'PC3', fontsize=11)
    ax.set_ylabel(f'PC2 ({pc2_var})' if pc2_var else 'PC2', fontsize=11)
    
    n_samples = len(merged)
    n_pops = len(pops)
    ax.set_title(f'PC3 vs PC2 with Regional Groupings ({n_samples} Samples)', 
                 fontsize=12, pad=15)
    
    # Create custom legend for regions
    legend_elements = [
        Patch(facecolor='#e74c3c', alpha=0.3, label='Highland'),
        Patch(facecolor='#27ae60', alpha=0.3, label='Amazonian'),
        Patch(facecolor='#3498db', alpha=0.3, label='Coastal')
    ]
    ax.legend(handles=legend_elements, loc='upper right', 
             frameon=True, fontsize=9, title='Region')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_panel_c(ax, data_dict, var_explained, pop_colors):
    """Panel C: PC3 vs PC1 with Quechua/Aymara language circles"""
    merged = data_dict.get('common_data', data_dict['peru_data'])
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel C: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Add language classification
    merged['Language'] = merged['POP'].apply(get_language_group)
    
    # Debug: print language classification
    print(f"\n  Panel C Language Classification:")
    for language in ['Quechua', 'Aymara', 'Other']:
        pops = merged[merged['Language'] == language]['POP'].unique()
        if len(pops) > 0:
            print(f"    {language}: {sorted(pops)}")
    
    # Define language colors for circles
    language_colors = {
        'Quechua': '#ff8c00',      # Dark Orange
        'Aymara': '#8b008b'        # Dark Magenta
    }
    
    # Use tight ellipses
    std_multipliers = {
        'Quechua': 2.0,
        'Aymara': 1.8
    }
    
    for language, color in language_colors.items():
        language_data = merged[merged['Language'] == language][['PC3', 'PC1']].values
        if len(language_data) >= 3:
            std_mult = std_multipliers.get(language, 2.0)
            draw_confidence_ellipse_custom(ax, language_data, color, language,
                                  alpha=0.22, linewidth=4.0, std_multiplier=std_mult)
            print(f"    Drew {language} ellipse with {len(language_data)} points")
    
    # Plot populations using GLOBAL colors
    pops = sorted(merged['POP'].dropna().unique())
    
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        color = pop_colors.get(pop, 'gray')
        ax.scatter(subset['PC3'], subset['PC1'],
                  label=pop,
                  color=color,
                  alpha=0.8, s=35, zorder=5)
    
    # Add all population labels
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC3'].mean()
            cy = subset['PC1'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold', zorder=10))
    
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    # Add variance explained
    pc3_var = f"{var_explained[2]:.2f}%" if var_explained is not None else ""
    pc1_var = f"{var_explained[0]:.2f}%" if var_explained is not None else ""
    
    ax.set_xlabel(f'PC3 ({pc3_var})' if pc3_var else 'PC3', fontsize=11)
    ax.set_ylabel(f'PC1 ({pc1_var})' if pc1_var else 'PC1', fontsize=11)
    
    n_samples = len(merged)
    n_pops = len(pops)
    ax.set_title(f'PC3 vs PC1 with Language Groupings ({n_samples} Samples)',
                 fontsize=12, pad=15)
    
    # Create custom legend for language groups
    legend_elements = [
        Patch(facecolor='#ff8c00', alpha=0.3, label='Quechua'),
        Patch(facecolor='#8b008b', alpha=0.3, label='Aymara')
    ]
    ax.legend(handles=legend_elements, loc='upper left', 
             frameon=True, fontsize=9, title='Language')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_panel_d(ax, data_dict, var_explained, pop_colors):
    """
    Panel D: Global context with SGDP
    REVIEWER CHANGE: Peru = filled circles, SGDP = triangles (different shapes)
    """
    peru_data = data_dict.get('sgdp_peru', data_dict['peru_data'])
    sgdp_data = data_dict.get('sgdp_other', pd.DataFrame())
    
    if peru_data.empty:
        ax.text(0.5, 0.5, 'Panel D: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Plot Peru samples (filled CIRCLES) using GLOBAL colors
    peru_pops = sorted(peru_data['POP'].dropna().unique())
    
    for pop in peru_pops:
        subset = peru_data[peru_data['POP'] == pop]
        color = pop_colors.get(pop, 'gray')
        ax.scatter(subset['PC1'], subset['PC2'],
                  color=color,
                  marker='o',  # Circles for Peru
                  alpha=0.7, s=30, zorder=3)
    
    # Process SGDP samples (TRIANGULAR markers - key reviewer change)
    n_sgdp_pops = 0
    sgdp_texts = []
    peru_texts = []
    
    if not sgdp_data.empty and 'POP' in sgdp_data.columns:
        sgdp_with_pop = sgdp_data.dropna(subset=['POP'])
        
        if not sgdp_with_pop.empty:
            centroids = sgdp_with_pop.groupby('POP')[['PC1', 'PC2']].mean().reset_index()
            pop_counts = sgdp_with_pop.groupby('POP').size()
            centroids['count'] = centroids['POP'].map(pop_counts)
            centroids = centroids.sort_values('POP')
            sgdp_pops = centroids['POP'].tolist()
            n_sgdp_pops = len(sgdp_pops)
            
            print(f"    Found {n_sgdp_pops} SGDP populations to plot")
            
            sgdp_cmap = plt.cm.Set3
            sgdp_colors = [sgdp_cmap(i / max(12, n_sgdp_pops)) for i in range(n_sgdp_pops)]
            
            # SGDP samples: use TRIANGLES (▲) instead of unfilled circles
            for i, pop in enumerate(sgdp_pops):
                subset = sgdp_with_pop[sgdp_with_pop['POP'] == pop]
                color = sgdp_colors[i % len(sgdp_colors)]
                
                ax.scatter(subset['PC1'], subset['PC2'],
                          color=color,
                          marker='^',  # TRIANGLES for SGDP
                          s=35,
                          alpha=0.7,
                          edgecolors='darkgray',
                          linewidth=0.5,
                          zorder=2)
            
            # Add SGDP population labels
            for idx, row in centroids.iterrows():
                label = row['POP']
                if len(label) > 15:
                    label = label.replace('_', ' ')
                    if len(label) > 15:
                        label = label.replace('Eskimo ', 'Esk. ')
                        label = label.replace('Yemenite ', 'Yem. ')
                        label = label.replace('Ju hoan North', 'Ju hoan N.')
                        if len(label) > 15:
                            label = label[:12] + '..'
                
                text = ax.text(row['PC1'], row['PC2'], label,
                             fontsize=6, style='italic', alpha=0.85,
                             ha='center', va='center', zorder=10)
                sgdp_texts.append(text)
    else:
        # If no population info, still plot SGDP with triangles
        if not sgdp_data.empty:
            ax.scatter(sgdp_data['PC1'], sgdp_data['PC2'],
                      color='lightblue',
                      marker='^',  # Triangles
                      s=25,
                      alpha=0.5,
                      zorder=2)
    
    # Add Peru population labels
    for pop in peru_pops:
        subset = peru_data[peru_data['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC1'].mean()
            cy = subset['PC2'].mean()
            text = ax.text(cx, cy, pop, fontsize=7, weight='bold',
                          ha='center', va='center', zorder=11)
            peru_texts.append(text)
    
    # Adjust SGDP labels
    if sgdp_texts:
        adjust_text(sgdp_texts, 
                   ax=ax,
                   force_text=(0.3, 0.3),
                   force_points=(0.2, 0.2),
                   expand_points=(1.2, 1.2),
                   expand_text=(1.1, 1.1),
                   only_move={'points':'', 'text':'xy'},
                   arrowprops=dict(arrowstyle='-', 
                                 connectionstyle='arc3,rad=0',
                                 color='lightgray', 
                                 lw=0.4, 
                                 alpha=0.4),
                   avoid_self=True,
                   lim=300)
    
    # Adjust Peru labels
    if peru_texts:
        adjust_text(peru_texts,
                   ax=ax,
                   add_objects=sgdp_texts,
                   force_text=(0.4, 0.4),
                   force_points=(0.2, 0.2),
                   expand_points=(1.3, 1.3),
                   expand_text=(1.2, 1.2),
                   only_move={'points':'', 'text':'xy'},
                   arrowprops=dict(arrowstyle='-',
                                 connectionstyle='arc3,rad=0.1',
                                 color='gray',
                                 lw=0.4,
                                 alpha=0.3),
                   avoid_self=True,
                   lim=300)
    
    # Handle overlaps
    all_texts = sgdp_texts + peru_texts
    if len(all_texts) > 1:
        try:
            renderer = ax.figure.canvas.get_renderer()
            for i, text1 in enumerate(all_texts):
                bbox1 = text1.get_window_extent(renderer=renderer)
                for j, text2 in enumerate(all_texts[i+1:], i+1):
                    bbox2 = text2.get_window_extent(renderer=renderer)
                    if bbox1.overlaps(bbox2):
                        if text1 in sgdp_texts and text2 in peru_texts:
                            text1.set_fontsize(5)
                        elif text2 in sgdp_texts and text1 in peru_texts:
                            text2.set_fontsize(5)
                        else:
                            text2.set_fontsize(5)
        except:
            pass
    
    # Add variance explained
    pc1_var = f"{var_explained[0]:.2f}%" if var_explained is not None else ""
    pc2_var = f"{var_explained[1]:.2f}%" if var_explained is not None else ""
    
    ax.set_xlabel(f'PC1 ({pc1_var})' if pc1_var else 'PC1', fontsize=11)
    ax.set_ylabel(f'PC2 ({pc2_var})' if pc2_var else 'PC2', fontsize=11)
    
    n_peru = len(peru_data)
    n_sgdp = len(sgdp_data) if not sgdp_data.empty else 0
    
    title = f'Global Context: Peru ({n_peru}) + SGDP ({n_sgdp} samples, {n_sgdp_pops} populations)'
    ax.set_title(title, fontsize=12, pad=15)
    
    # ADD LEGEND distinguishing Peru vs SGDP markers (REVIEWER REQUEST)
    # Place in lower left to avoid overlap with SGDP labels in upper right
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markersize=8, label=f'Peru (n={n_peru})'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='lightblue',
               markeredgecolor='darkgray', markersize=8, label=f'SGDP (n={n_sgdp})')
    ]
    ax.legend(handles=legend_elements, loc='lower left', frameon=True, 
              fontsize=8, fancybox=True, framealpha=0.9)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    
    print("="*70)
    print("PCA Analysis and Figure 3 Generation - REVIEWER-REVISED VERSION")
    print("Changes: Different shapes for SGDP (triangles) vs Peru (circles)")
    print("="*70)
    
    if not check_plink():
        sys.exit(1)
    
    os.makedirs(BASE_DIR, exist_ok=True)
    
    print("\nLocating input files...")
    
    peru_vcf = find_file(PERU_VCF, ALT_PERU_VCF, 
                        'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz')
    common_vcf = find_file(COMMON_VCF, peru_vcf)
    sgdp_vcf = find_file(SGDP_VCF, ALT_SGDP_VCF, ALT_SGDP_VCF2)
    
    pop_file = find_file(POP_FILE, ALT_POP_FILE, ALT_POP_FILE2)
    
    if not peru_vcf:
        print("ERROR: No Peru VCF file found")
        sys.exit(1)
    
    if not pop_file:
        print("ERROR: No population file found")
        sys.exit(1)
    
    print(f"  Peru VCF: {peru_vcf}")
    print(f"  Population file: {pop_file}")
    if sgdp_vcf:
        print(f"  SGDP VCF: {sgdp_vcf}")
    
    print("\nRunning PCA analyses...")
    
    run_pca_if_needed(peru_vcf, PCA_PERU_PREFIX)
    
    if common_vcf and common_vcf != peru_vcf:
        run_pca_if_needed(common_vcf, PCA_COMMON_PREFIX)
    
    if sgdp_vcf:
        run_pca_if_needed(sgdp_vcf, PCA_SGDP_PREFIX)
    
    print("\nLoading PCA results...")
    data_dict = {}
    
    pop_df = load_population_file(pop_file)
    pop_df = normalize_population_column(pop_df, 'POP')  # Normalize population names
    print(f"  Loaded {len(pop_df)} samples from Peru population file")
    if 'POP' in pop_df.columns:
        print(f"  Found {pop_df['POP'].nunique()} unique Peru populations")
    
    # Load variance explained
    var_explained_peru = None
    var_explained_common = None
    var_explained_sgdp = None
    
    if os.path.exists(f"{PCA_PERU_PREFIX}.eigenvec"):
        peru_pca = load_eigenvec(f"{PCA_PERU_PREFIX}.eigenvec")
        var_explained_peru = load_eigenval(f"{PCA_PERU_PREFIX}.eigenvec")
        
        data_dict['peru_data'] = peru_pca.merge(
            pop_df[['IID', 'POP']], on='IID', how='left'
        )
        data_dict['peru_data'] = normalize_population_column(data_dict['peru_data'], 'POP')
        print(f"  Peru PCA: {len(data_dict['peru_data'])} samples, "
              f"{data_dict['peru_data']['POP'].nunique()} populations")
    else:
        print("  Warning: Peru PCA not available")
        data_dict['peru_data'] = pd.DataFrame()
    
    if os.path.exists(f"{PCA_COMMON_PREFIX}.eigenvec"):
        common_pca = load_eigenvec(f"{PCA_COMMON_PREFIX}.eigenvec")
        var_explained_common = load_eigenval(f"{PCA_COMMON_PREFIX}.eigenvec")
        
        data_dict['common_data'] = common_pca.merge(
            pop_df[['IID', 'POP']], on='IID', how='left'
        )
        data_dict['common_data'] = normalize_population_column(data_dict['common_data'], 'POP')
        print(f"  Common PCA: {len(data_dict['common_data'])} samples")
    else:
        data_dict['common_data'] = data_dict.get('peru_data', pd.DataFrame()).copy()
        var_explained_common = var_explained_peru
    
    sgdp_eigenvec_file = f"{PCA_SGDP_PREFIX}.eigenvec"
    
    if not os.path.exists(sgdp_eigenvec_file) and os.path.exists(f"{ALT_SGDP_PCA}.eigenvec"):
        sgdp_eigenvec_file = f"{ALT_SGDP_PCA}.eigenvec"
        print(f"  Using alternative SGDP eigenvec: {sgdp_eigenvec_file}")
    
    if os.path.exists(sgdp_eigenvec_file):
        print(f"\n  Loading SGDP PCA from: {sgdp_eigenvec_file}")
        sgdp_pca = load_eigenvec(sgdp_eigenvec_file)
        var_explained_sgdp = load_eigenval(sgdp_eigenvec_file)
        
        print(f"    Loaded {len(sgdp_pca)} total samples from SGDP eigenvec")
        
        peru_ids = set(pop_df['IID'].values) if not pop_df.empty else set()
        
        if peru_ids:
            is_peru = sgdp_pca['IID'].isin(peru_ids)
            
            data_dict['sgdp_peru'] = sgdp_pca[is_peru].copy()
            data_dict['sgdp_peru'] = data_dict['sgdp_peru'].merge(
                pop_df[['IID', 'POP']], on='IID', how='left'
            )
            data_dict['sgdp_peru'] = normalize_population_column(data_dict['sgdp_peru'], 'POP')
            
            data_dict['sgdp_other'] = sgdp_pca[~is_peru].copy()
            print(f"    Split: {len(data_dict['sgdp_peru'])} Peru, {len(data_dict['sgdp_other'])} SGDP samples")
            
            sgdp_pop_file = find_file(SGDP_POP_FILE, ALT_SGDP_POP_FILE)
            
            data_dict['sgdp_other'], diagnostics = match_sgdp_populations(
                data_dict['sgdp_other'], 
                sgdp_pop_file
            )
            data_dict['sgdp_other'] = normalize_population_column(data_dict['sgdp_other'], 'POP')
            
            if diagnostics['total_sgdp'] > 0:
                final_match_rate = 100 * (diagnostics['direct_matches'] + diagnostics['fallback_matches']) / diagnostics['total_sgdp']
                if final_match_rate < 90:
                    print(f"\n  WARNING: Match rate {final_match_rate:.1f}% is below target 90%")
                    print("    Consider updating the SGDP population mapping file")
    else:
        print("  SGDP PCA not available")
    
    print("\nCreating Figure 3...")
    
    # Create GLOBAL color dictionary for consistent colors across all panels
    all_peru_pops = sorted(data_dict['peru_data']['POP'].dropna().unique())
    global_colors = get_population_colors(len(all_peru_pops))
    PERU_POP_COLORS = dict(zip(all_peru_pops, global_colors))
    print(f"  Created consistent color mapping for {len(PERU_POP_COLORS)} populations")
    
    fig = plt.figure(figsize=(20, 16), dpi=300, constrained_layout=False)
    
    gs = GridSpec(2, 2, figure=fig,
                  left=0.06, right=0.94,
                  bottom=0.06, top=0.88,
                  wspace=0.40,
                  hspace=0.35)
    
    print("  Plotting Panel A: PC1 vs PC2")
    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_a(ax_a, data_dict, var_explained_peru, PERU_POP_COLORS)
    ax_a.text(-0.12, 1.08, 'A', transform=ax_a.transAxes,
              fontsize=20, fontweight='bold')
    
    print("  Plotting Panel B: PC3 vs PC2 with Regional Circles")
    ax_b = fig.add_subplot(gs[0, 1])
    plot_panel_b(ax_b, data_dict, var_explained_common, PERU_POP_COLORS)
    ax_b.text(-0.12, 1.08, 'B', transform=ax_b.transAxes,
              fontsize=20, fontweight='bold')
    
    print("  Plotting Panel C: PC3 vs PC1 with Language Circles")
    ax_c = fig.add_subplot(gs[1, 0])
    plot_panel_c(ax_c, data_dict, var_explained_common, PERU_POP_COLORS)
    ax_c.text(-0.12, 1.08, 'C', transform=ax_c.transAxes,
              fontsize=20, fontweight='bold')
    
    print("  Plotting Panel D: Global Context with SGDP (triangles vs circles)")
    ax_d = fig.add_subplot(gs[1, 1])
    plot_panel_d(ax_d, data_dict, var_explained_sgdp, PERU_POP_COLORS)
    ax_d.text(-0.12, 1.08, 'D', transform=ax_d.transAxes,
              fontsize=20, fontweight='bold')
    
    # FIXED: Correct sample count (736, not 746)
    fig.suptitle('Figure 3. Principal Component Analysis (PCA) of Peruvian Populations\n' +
                 'after removing related individuals (736 samples)',
                 fontsize=16, fontweight='bold', y=0.965)
    
    print("\nSaving outputs...")
    
    output_png = os.path.join(BASE_DIR, 'Figure3_PCA_composite.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  ✓ PNG saved: {output_png}")
    
    output_pdf = os.path.join(BASE_DIR, 'Figure3_PCA_composite.pdf')
    plt.savefig(output_pdf, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  ✓ PDF saved: {output_pdf}")
    
    plt.close()
    
    print("\n" + "="*70)
    print("FINAL SUMMARY:")
    print("-" * 70)
    
    if 'peru_data' in data_dict and not data_dict['peru_data'].empty:
        n_peru = len(data_dict['peru_data'])
        n_peru_pops = data_dict['peru_data']['POP'].nunique()
        print(f"Peru samples: {n_peru} samples, {n_peru_pops} populations")
    
    if 'sgdp_other' in data_dict and not data_dict['sgdp_other'].empty:
        n_sgdp = len(data_dict['sgdp_other'])
        if 'POP' in data_dict['sgdp_other'].columns:
            n_sgdp_matched = data_dict['sgdp_other']['POP'].notna().sum()
            n_sgdp_pops = data_dict['sgdp_other']['POP'].nunique()
            print(f"SGDP samples: {n_sgdp} loaded, {n_sgdp_matched} matched to {n_sgdp_pops} populations")
        else:
            print(f"SGDP samples: {n_sgdp} loaded (no population info)")
    
    print("-" * 70)
    print("REVIEWER CHANGES IMPLEMENTED:")
    print("  ✓ Panel D: SGDP samples now use triangles (▲)")
    print("  ✓ Panel D: Peru samples use circles (●)")
    print("  ✓ Panel D: Legend distinguishes marker shapes")
    print("  ✓ Fixed sample count: 736 (was 746)")
    print("  ✓ Expanded Aymara populations (added JAQARUS)")
    print("  ✓ Added IQUITOS to Amazonian region")
    print("  ✓ Consistent population colors across all panels A-D")
    print("-" * 70)
    print("Figure 3 has been successfully created!")
    print(f"Output files saved to: {BASE_DIR}/")
    print("="*70)

if __name__ == "__main__":
    main()