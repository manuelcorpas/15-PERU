#!/usr/bin/env python3
"""
Complete PCA analysis and Figure 4 generation - ENHANCED SGDP MATCHING VERSION
Implements robust fallback matching for 345/345 SGDP samples
All outputs go to ANALYSIS/00-08-PCA/
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

# SGDP population file - NOW IN 00-08-PCA directory as requested
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
    'JACARUS', 'LAMAS', 'LAMBAYEQUE', 'LIMA', 'MATZES', 'MATSIGUENKAS',
    'MOCHES', 'MOQUEGUA', 'NAHUA', 'PUNO', 'QEROS', 'SHIPIBO_INS',
    'TACNA', 'TALLANES', 'TRUJILLO', 'TUMBES', 'UROS'
}

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
    
    # Load with flexible whitespace delimiter
    df = pd.read_csv(filepath, delim_whitespace=True, header=None)
    n_cols = df.shape[1]
    
    # Set column names
    df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, n_cols-1)]
    
    # Ensure IID is string and stripped
    df['IID'] = df['IID'].astype(str).str.strip()
    
    # Pad with zeros if fewer than 10 PCs
    for i in range(n_cols-1, 11):
        df[f'PC{i}'] = 0
    
    return df

def load_population_file(filepath):
    """Load population assignments - HARDENED VERSION"""
    if not os.path.exists(filepath):
        print(f"Warning: Population file not found: {filepath}")
        return pd.DataFrame(columns=['FID', 'IID', 'POP'])
    
    # Try tab-delimited first
    try:
        df = pd.read_csv(filepath, delimiter='\t')
    except:
        # Try space-delimited
        try:
            df = pd.read_csv(filepath, sep=r'\s+')
        except:
            print(f"Error reading population file: {filepath}")
            return pd.DataFrame(columns=['FID', 'IID', 'POP'])
    
    # Debug: print detected columns
    print(f"    Detected columns in {os.path.basename(filepath)}: {df.columns.tolist()}")
    
    # Normalize column names - be more flexible
    df.columns = [col.strip() for col in df.columns]
    
    # Handle Population/population/POP/pop column variations
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
    
    # Handle IID/Sample/ID variations
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
    
    # Handle FID
    if 'FID' not in df.columns:
        if 'fid' in df.columns:
            df['FID'] = df['fid'].astype(str).str.strip()
        elif 'IID' in df.columns:
            df['FID'] = df['IID']
        else:
            df['FID'] = ''
    else:
        df['FID'] = df['FID'].astype(str).str.strip()
    
    # Debug: show sample of data
    if 'IID' in df.columns and 'POP' in df.columns:
        print(f"    Sample IIDs: {df['IID'].head(3).tolist()}")
        print(f"    Sample POPs: {df['POP'].head(3).tolist()}")
    
    return df
def parse_population_from_iid(iid):
    """
    Parse population name from SGDP IID - COMPLETE VERSION
    Handles all SGDP naming patterns including simple prefix extraction
    """
    # First try specific known patterns with underscores
    patterns_with_groups = [
        (r'^(.+)_HGDP\d+$', 1),           # e.g., Piapoco_HGDP00702
        (r'^(.+)_NA\d+$', 1),              # e.g., Quechua_NA11200, Luhya_NA19044
        (r'^(.+)_HG\d+$', 1),              # e.g., Gambian_HG02464, Esan_HG03100
        (r'^(.+)_[A-Z]{2,4}\d+$', 1),     # e.g., Saharawi_SAH31, Mada_CAMD042
        (r'^(.+)_[A-Z][a-z]+\d+$', 1),    # e.g., Aleut_Ale22, Estonian_Est400
        (r'^(.+)_SS\d+$', 1),              # e.g., Kharia_SS6004481
        (r'^(.+)_Ayodo_\d+[A-Z]?$', 1),   # e.g., Somali_Ayodo_81S, Luo_Ayodo_502C
        (r'^(.+)_[a-z]+\d+$', 1),         # e.g., Zapotec_zapo0098, Armenian_armenia293
        (r'^(.+)_[A-Z]\d+$', 1),           # e.g., Kapu_K1, Irula_I3
        (r'^(.+)_[A-Z]\-\d+$', 1),        # e.g., Kusunda_K-2, Crete_TZ-11
        (r'^(.+)_[a-z]+\d+_[a-z]+$', 1),  # e.g., Tajik_tdj409_shugnan
    ]
    
    for pattern, group in patterns_with_groups:
        match = re.match(pattern, iid)
        if match:
            population = match.group(group)
            # Clean up known suffixes that might still be attached
            population = re.sub(r'_(South|North|Jew|Pandit|Ayodo)$', r'_\1', population)
            if population and len(population) > 2:
                return population
    
    # FALLBACK: If no specific pattern matches but there's an underscore,
    # take everything before the first underscore as the population name
    # This catches cases like "Tu_HGDP01355", "Yi_HGDP01179"
    if '_' in iid:
        parts = iid.split('_')
        first_part = parts[0]
        # Check if the first part looks like a population name (not just numbers)
        if first_part and len(first_part) >= 2 and not first_part.isdigit():
            return first_part
    
    # If no underscore pattern matches, check if it's a simple population_number format
    simple_match = re.match(r'^([A-Za-z_]+?)_?[A-Z]*\d+[A-Z]?$', iid)
    if simple_match:
        pop = simple_match.group(1)
        if pop and len(pop) > 2 and not pop.isdigit():
            return pop
    
    # Return None if we can't parse
    return None
# ============================================================================
# COLOR PALETTE
# ============================================================================

def get_population_colors(n_pops):
    """Generate distinct colors for populations"""
    if n_pops <= 20:
        colors = list(plt.cm.tab20.colors)
    else:
        # Combine multiple colormaps for more colors
        colors = (list(plt.cm.tab20.colors) + 
                 list(plt.cm.tab20b.colors) + 
                 list(plt.cm.tab20c.colors))
    
    return colors[:n_pops]

# ============================================================================
# ENHANCED SGDP MATCHING FUNCTION
# ============================================================================

def match_sgdp_populations(sgdp_data, sgdp_pop_file):
    """
    Enhanced SGDP population matching with robust fallback.
    Returns: (matched_data, diagnostics_dict)
    """
    print("\n  === ENHANCED SGDP POPULATION MATCHING ===")
    
    # Initialize diagnostics
    diagnostics = {
        'total_sgdp': len(sgdp_data),
        'direct_matches': 0,
        'fallback_matches': 0,
        'unmatched': 0,
        'unique_pops': set(),
        'pop_counts': {},
        'match_details': []  # For TSV output
    }
    
    # Step A: Load and clean SGDP population file
    if sgdp_pop_file and os.path.exists(sgdp_pop_file):
        print(f"  Loading SGDP population file: {sgdp_pop_file}")
        sgdp_pop = load_population_file(sgdp_pop_file)
        
        if not sgdp_pop.empty and 'IID' in sgdp_pop.columns and 'POP' in sgdp_pop.columns:
            # Clean up
            sgdp_pop['IID'] = sgdp_pop['IID'].astype(str).str.strip()
            sgdp_pop['POP'] = sgdp_pop['POP'].astype(str).str.strip()
            
            # Filter out Peru populations
            sgdp_pop_only = sgdp_pop[~sgdp_pop['POP'].isin(PERU_POPULATIONS)].copy()
            print(f"    Filtered to {len(sgdp_pop_only)} SGDP entries (removed Peru)")
            
            # Direct merge
            sgdp_data = sgdp_data.merge(
                sgdp_pop_only[['IID', 'POP']], 
                on='IID', 
                how='left'
            )
            
            # Count direct matches
            diagnostics['direct_matches'] = sgdp_data['POP'].notna().sum()
            print(f"    Direct matches: {diagnostics['direct_matches']}/{len(sgdp_data)}")
            
            # Record match details for matched samples
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
    
    # Step B: Fallback parsing for unmatched samples (more conservative)
    unmatched_mask = sgdp_data['POP'].isna()
    n_unmatched = unmatched_mask.sum()
    
    if n_unmatched > 0:
        print(f"\n  Applying fallback parsing for {n_unmatched} unmatched samples...")
        
        # Parse population from IID for unmatched samples
        fallback_pops = []
        parse_examples = []
        
        for idx, row in sgdp_data[unmatched_mask].iterrows():
            iid = row['IID']
            parsed_pop = parse_population_from_iid(iid)
            
            # Only use parsed population if it's valid and not a Peru population
            if parsed_pop and parsed_pop not in PERU_POPULATIONS:
                sgdp_data.loc[idx, 'POP'] = parsed_pop
                fallback_pops.append(parsed_pop)
                
                diagnostics['match_details'].append({
                    'IID': iid,
                    'POP': parsed_pop,
                    'source': 'parsed_from_iid'
                })
                
                # Collect examples for reporting
                if len(parse_examples) < 5:
                    parse_examples.append(f"      {iid} → {parsed_pop}")
        
        diagnostics['fallback_matches'] = len(fallback_pops)
        print(f"    Fallback matches: {diagnostics['fallback_matches']}")
        if parse_examples:
            print("    Example parsings:")
            for example in parse_examples:
                print(example)
    
    # Step C: Final statistics
    diagnostics['unmatched'] = sgdp_data['POP'].isna().sum()
    diagnostics['unique_pops'] = set(sgdp_data['POP'].dropna().unique())
    
    # Count samples per population
    pop_counts = sgdp_data['POP'].value_counts()
    diagnostics['pop_counts'] = pop_counts.to_dict()
    
    # Print summary
    total_matched = diagnostics['direct_matches'] + diagnostics['fallback_matches']
    match_rate = 100 * total_matched / diagnostics['total_sgdp'] if diagnostics['total_sgdp'] > 0 else 0
    
    print(f"\n  === MATCHING SUMMARY ===")
    print(f"    Direct matches: {diagnostics['direct_matches']}")
    print(f"    Fallback matches: {diagnostics['fallback_matches']}")
    print(f"    Total matched: {total_matched}/{diagnostics['total_sgdp']} ({match_rate:.1f}%)")
    print(f"    Unique SGDP populations: {len(diagnostics['unique_pops'])}")
    
    # Top 10 populations by count
    if diagnostics['pop_counts']:
        print(f"\n    Top 10 SGDP populations by sample count:")
        for pop, count in list(diagnostics['pop_counts'].items())[:10]:
            print(f"      {pop}: {count} samples")
    
    # Write diagnostic TSV files
    write_diagnostic_files(diagnostics, sgdp_data)
    
    # Report any remaining unmatched
    if diagnostics['unmatched'] > 0:
        print(f"\n    WARNING: {diagnostics['unmatched']} samples still unmatched")
        unmatched_iids = sgdp_data[sgdp_data['POP'].isna()]['IID'].head(50).tolist()
        print(f"    First unmatched IIDs (up to 50):")
        for iid in unmatched_iids[:10]:  # Show first 10 in console
            print(f"      {iid}")
        print(f"    (See sgdp_unmatched.tsv for full list)")
    
    return sgdp_data, diagnostics

def write_diagnostic_files(diagnostics, sgdp_data):
    """Write diagnostic TSV files"""
    output_dir = BASE_DIR
    
    # Write matched IID to POP mapping
    matched_file = os.path.join(output_dir, 'sgdp_iid_to_pop_resolved.tsv')
    if diagnostics['match_details']:
        matched_df = pd.DataFrame(diagnostics['match_details'])
        matched_df.to_csv(matched_file, sep='\t', index=False)
        print(f"    ✓ Wrote: {matched_file}")
    
    # Write unmatched IIDs if any
    if diagnostics['unmatched'] > 0:
        unmatched_file = os.path.join(output_dir, 'sgdp_unmatched.tsv')
        unmatched_df = sgdp_data[sgdp_data['POP'].isna()][['IID']].copy()
        unmatched_df.to_csv(unmatched_file, sep='\t', index=False)
        print(f"    ✓ Wrote: {unmatched_file}")

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_panel_a(ax, data_dict):
    """Panel A: PC1 vs PC2 - Main population structure"""
    merged = data_dict['peru_data']
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel A: Data not available', 
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Get unique populations
    pops = sorted(merged['POP'].dropna().unique())
    colors = get_population_colors(len(pops))
    pop_colors = dict(zip(pops, colors))
    
    # Plot each population
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        ax.scatter(subset['PC1'], subset['PC2'],
                  label=pop,
                  color=pop_colors[pop],
                  alpha=0.7, s=30)
    
    # Add population labels at centroids
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC1'].mean()
            cy = subset['PC2'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold'))
    
    # Adjust text positions
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    ax.set_xlabel('PC1', fontsize=11)
    ax.set_ylabel('PC2', fontsize=11)
    ax.set_title('PCA of ARRAY & WGS Samples by Population', fontsize=12, pad=15)
    
    # Legend
    if len(pops) <= 30:
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                 frameon=False, fontsize=7, ncol=1)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_panel_b(ax, data_dict):
    """Panel B: PC3 on X-axis, PC2 on Y-axis"""
    merged = data_dict.get('common_data', data_dict['peru_data'])
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel B: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    pops = sorted(merged['POP'].dropna().unique())
    colors = get_population_colors(len(pops))
    pop_colors = dict(zip(pops, colors))
    
    # Plot PC3 on X-axis, PC2 on Y-axis
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        ax.scatter(subset['PC3'], subset['PC2'],
                  label=pop,
                  color=pop_colors[pop],
                  alpha=0.7, s=30)
    
    # Add population labels
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC3'].mean()
            cy = subset['PC2'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold'))
    
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    ax.set_xlabel('PC3', fontsize=11)
    ax.set_ylabel('PC2', fontsize=11)
    
    n_samples = len(merged)
    n_pops = len(pops)
    ax.set_title(f'PCA Common Variants IBD Clean ({n_samples} Samples, {n_pops} Populations)', 
                 fontsize=12, pad=15)
    
    if len(pops) <= 30:
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                 frameon=False, fontsize=7, ncol=1)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_panel_c(ax, data_dict):
    """Panel C: PC3 vs PC4 - Fine-scale structure"""
    merged = data_dict.get('common_data', data_dict['peru_data'])
    
    if merged.empty:
        ax.text(0.5, 0.5, 'Panel C: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    pops = sorted(merged['POP'].dropna().unique())
    colors = get_population_colors(len(pops))
    pop_colors = dict(zip(pops, colors))
    
    # Plot PC3 vs PC4
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        ax.scatter(subset['PC3'], subset['PC4'],
                  label=pop,
                  color=pop_colors[pop],
                  alpha=0.7, s=30)
    
    # Add all population labels
    texts = []
    for pop in pops:
        subset = merged[merged['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC3'].mean()
            cy = subset['PC4'].mean()
            texts.append(ax.text(cx, cy, pop, fontsize=8, weight='bold'))
    
    if texts:
        adjust_text(texts, ax=ax,
                   expand_points=(1.2, 1.2),
                   arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    ax.set_xlabel('PC3', fontsize=11)
    ax.set_ylabel('PC4', fontsize=11)
    
    n_samples = len(merged)
    n_pops = len(pops)
    ax.set_title(f'PCA Common Variants IBD Clean ({n_samples} Samples, {n_pops} Populations)',
                 fontsize=12, pad=15)
    
    # Legend at bottom
    if len(pops) <= 40:
        ax.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center',
                 frameon=False, fontsize=6, ncol=6)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
def plot_panel_d(ax, data_dict):
    """Panel D: Global context with SGDP - ROBUST LABEL DE-OVERLAPPING"""
    peru_data = data_dict.get('sgdp_peru', data_dict['peru_data'])
    sgdp_data = data_dict.get('sgdp_other', pd.DataFrame())
    
    if peru_data.empty:
        ax.text(0.5, 0.5, 'Panel D: Data not available',
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        return
    
    # Plot Peru samples (filled markers)
    peru_pops = sorted(peru_data['POP'].dropna().unique())
    colors = get_population_colors(len(peru_pops))
    peru_colors = dict(zip(peru_pops, colors))
    
    for pop in peru_pops:
        subset = peru_data[peru_data['POP'] == pop]
        ax.scatter(subset['PC1'], subset['PC2'],
                  color=peru_colors[pop],
                  alpha=0.6, s=25, zorder=2)
    
    # Process SGDP samples
    n_sgdp_pops = 0
    sgdp_texts = []
    peru_texts = []
    
    if not sgdp_data.empty and 'POP' in sgdp_data.columns:
        sgdp_with_pop = sgdp_data.dropna(subset=['POP'])
        
        if not sgdp_with_pop.empty:
            # Get centroids and counts
            centroids = sgdp_with_pop.groupby('POP')[['PC1', 'PC2']].mean().reset_index()
            pop_counts = sgdp_with_pop.groupby('POP').size()
            centroids['count'] = centroids['POP'].map(pop_counts)
            centroids = centroids.sort_values('POP')
            sgdp_pops = centroids['POP'].tolist()
            n_sgdp_pops = len(sgdp_pops)
            
            print(f"    Found {n_sgdp_pops} SGDP populations to plot")
            
            # Colors for SGDP
            sgdp_cmap = plt.cm.Set3
            sgdp_colors = [sgdp_cmap(i / max(12, n_sgdp_pops)) for i in range(n_sgdp_pops)]
            
            # Plot SGDP populations with hollow markers
            for i, pop in enumerate(sgdp_pops):
                subset = sgdp_with_pop[sgdp_with_pop['POP'] == pop]
                color = sgdp_colors[i % len(sgdp_colors)]
                
                ax.scatter(subset['PC1'], subset['PC2'],
                          facecolors='none',
                          edgecolors=color,
                          linewidth=1.0,
                          s=30,
                          alpha=0.7,
                          zorder=3)
            
            # Add SGDP labels
            for idx, row in centroids.iterrows():
                # Simplify long names
                label = row['POP']
                if len(label) > 15:
                    label = label.replace('_', ' ')
                    if len(label) > 15:
                        # Abbreviate specific known long names
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
        if not sgdp_data.empty:
            ax.scatter(sgdp_data['PC1'], sgdp_data['PC2'],
                      facecolors='none',
                      edgecolors='lightblue',
                      linewidth=0.8,
                      s=25,
                      alpha=0.5,
                      zorder=3)
    
    # Add Peru population labels
    for pop in peru_pops:
        subset = peru_data[peru_data['POP'] == pop]
        if len(subset) > 0:
            cx = subset['PC1'].mean()
            cy = subset['PC2'].mean()
            text = ax.text(cx, cy, pop, fontsize=7, weight='bold',
                          ha='center', va='center', zorder=11)
            peru_texts.append(text)
    
    # ROBUST LABEL ADJUSTMENT
    # First adjust SGDP labels among themselves
    if sgdp_texts:
        adjust_text(sgdp_texts, 
                   ax=ax,
                   force_text=(0.3, 0.3),  # Gentle force
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
    
    # Then adjust Peru labels, avoiding SGDP labels
    if peru_texts:
        # Get bounding boxes of SGDP texts to avoid them
        sgdp_boxes = []
        if sgdp_texts:
            renderer = ax.figure.canvas.get_renderer()
            for text in sgdp_texts:
                bbox = text.get_window_extent(renderer=renderer)
                bbox = bbox.transformed(ax.transData.inverted())
                sgdp_boxes.append(bbox)
        
        adjust_text(peru_texts,
                   ax=ax,
                   add_objects=sgdp_texts,  # Avoid SGDP labels
                   force_text=(0.4, 0.4),  # Slightly stronger force for Peru
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
    
    # Fine-tune: if any labels still overlap significantly, make one smaller
    all_texts = sgdp_texts + peru_texts
    if len(all_texts) > 1:
        try:
            renderer = ax.figure.canvas.get_renderer()
            for i, text1 in enumerate(all_texts):
                bbox1 = text1.get_window_extent(renderer=renderer)
                for j, text2 in enumerate(all_texts[i+1:], i+1):
                    bbox2 = text2.get_window_extent(renderer=renderer)
                    if bbox1.overlaps(bbox2):
                        # Make the less important one smaller
                        if text1 in sgdp_texts and text2 in peru_texts:
                            text1.set_fontsize(5)
                        elif text2 in sgdp_texts and text1 in peru_texts:
                            text2.set_fontsize(5)
                        else:
                            # Both same type, make the one with fewer samples smaller
                            text2.set_fontsize(5)
        except:
            pass  # Skip if renderer not available
    
    ax.set_xlabel('PC1', fontsize=11)
    ax.set_ylabel('PC2', fontsize=11)
    
    n_peru = len(peru_data)
    n_sgdp = len(sgdp_data) if not sgdp_data.empty else 0
    
    title = f'PCA Peru ({n_peru} Samples, {len(peru_pops)} Populations) SGDP ({n_sgdp} Samples, {n_sgdp_pops} Populations)'
    ax.set_title(title, fontsize=12, pad=15)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function"""
    
    print("="*70)
    print("PCA Analysis and Figure 4 Generation - ENHANCED SGDP MATCHING VERSION")
    print("="*70)
    
    # Check PLINK
    if not check_plink():
        sys.exit(1)
    
    # Create output directory
    os.makedirs(BASE_DIR, exist_ok=True)
    
    # Find input files
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
    
    # Run PCAs if needed
    print("\nRunning PCA analyses...")
    
    # Main Peru PCA
    run_pca_if_needed(peru_vcf, PCA_PERU_PREFIX)
    
    # Common variants PCA
    if common_vcf and common_vcf != peru_vcf:
        run_pca_if_needed(common_vcf, PCA_COMMON_PREFIX)
    
    # SGDP PCA
    if sgdp_vcf:
        run_pca_if_needed(sgdp_vcf, PCA_SGDP_PREFIX)
    
    # Load data
    print("\nLoading PCA results...")
    data_dict = {}
    
    # Load Peru population file
    pop_df = load_population_file(pop_file)
    print(f"  Loaded {len(pop_df)} samples from Peru population file")
    if 'POP' in pop_df.columns:
        print(f"  Found {pop_df['POP'].nunique()} unique Peru populations")
    
    # Load Peru PCA
    if os.path.exists(f"{PCA_PERU_PREFIX}.eigenvec"):
        peru_pca = load_eigenvec(f"{PCA_PERU_PREFIX}.eigenvec")
        
        # Merge with population data
        data_dict['peru_data'] = peru_pca.merge(
            pop_df[['IID', 'POP']], on='IID', how='left'
        )
        print(f"  Peru PCA: {len(data_dict['peru_data'])} samples, "
              f"{data_dict['peru_data']['POP'].nunique()} populations")
    else:
        print("  Warning: Peru PCA not available")
        data_dict['peru_data'] = pd.DataFrame()
    
    # Load common variants PCA
    if os.path.exists(f"{PCA_COMMON_PREFIX}.eigenvec"):
        common_pca = load_eigenvec(f"{PCA_COMMON_PREFIX}.eigenvec")
        data_dict['common_data'] = common_pca.merge(
            pop_df[['IID', 'POP']], on='IID', how='left'
        )
        print(f"  Common PCA: {len(data_dict['common_data'])} samples")
    else:
        data_dict['common_data'] = data_dict.get('peru_data', pd.DataFrame()).copy()
    
    # Load SGDP PCA - WITH ENHANCED MATCHING
    sgdp_eigenvec_file = f"{PCA_SGDP_PREFIX}.eigenvec"
    
    # Check alternative SGDP PCA location
    if not os.path.exists(sgdp_eigenvec_file) and os.path.exists(f"{ALT_SGDP_PCA}.eigenvec"):
        sgdp_eigenvec_file = f"{ALT_SGDP_PCA}.eigenvec"
        print(f"  Using alternative SGDP eigenvec: {sgdp_eigenvec_file}")
    
    if os.path.exists(sgdp_eigenvec_file):
        print(f"\n  Loading SGDP PCA from: {sgdp_eigenvec_file}")
        sgdp_pca = load_eigenvec(sgdp_eigenvec_file)
        print(f"    Loaded {len(sgdp_pca)} total samples from SGDP eigenvec")
        
        # Get Peru IDs for separation
        peru_ids = set(pop_df['IID'].values) if not pop_df.empty else set()
        
        if peru_ids:
            # Separate Peru and SGDP samples
            is_peru = sgdp_pca['IID'].isin(peru_ids)
            
            # Peru samples in global context
            data_dict['sgdp_peru'] = sgdp_pca[is_peru].copy()
            data_dict['sgdp_peru'] = data_dict['sgdp_peru'].merge(
                pop_df[['IID', 'POP']], on='IID', how='left'
            )
            
            # SGDP samples (non-Peru)
            data_dict['sgdp_other'] = sgdp_pca[~is_peru].copy()
            print(f"    Split: {len(data_dict['sgdp_peru'])} Peru, {len(data_dict['sgdp_other'])} SGDP samples")
            
            # ==================================================================
            # ENHANCED SGDP POPULATION MATCHING
            # ==================================================================
            
            # Find SGDP population file
            sgdp_pop_file = find_file(SGDP_POP_FILE, ALT_SGDP_POP_FILE)
            
            # Apply enhanced matching
            data_dict['sgdp_other'], diagnostics = match_sgdp_populations(
                data_dict['sgdp_other'], 
                sgdp_pop_file
            )
            
            # Verify match rate
            if diagnostics['total_sgdp'] > 0:
                final_match_rate = 100 * (diagnostics['direct_matches'] + diagnostics['fallback_matches']) / diagnostics['total_sgdp']
                if final_match_rate < 90:
                    print(f"\n  WARNING: Match rate {final_match_rate:.1f}% is below target 90%")
                    print("    Consider updating the SGDP population mapping file")
    else:
        print("  SGDP PCA not available")
    
    # Create Figure 4 with fixed spacing
    print("\nCreating Figure 4...")
    
    # Use constrained_layout=False and manual spacing control
    fig = plt.figure(figsize=(20, 16), dpi=300, constrained_layout=False)
    
    # Adjust GridSpec to prevent overlap
    gs = GridSpec(2, 2, figure=fig,
                  left=0.06, right=0.94,
                  bottom=0.06, top=0.88,  # Reduced top to make room for suptitle
                  wspace=0.40,  # Horizontal spacing between panels
                  hspace=0.35)  # Vertical spacing between panels
    
    # Panel A
    print("  Plotting Panel A: PC1 vs PC2")
    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_a(ax_a, data_dict)
    ax_a.text(-0.12, 1.08, 'A', transform=ax_a.transAxes,
              fontsize=20, fontweight='bold')
    
    # Panel B
    print("  Plotting Panel B: PC3 vs PC2 (PC2 on Y-axis)")
    ax_b = fig.add_subplot(gs[0, 1])
    plot_panel_b(ax_b, data_dict)
    ax_b.text(-0.12, 1.08, 'B', transform=ax_b.transAxes,
              fontsize=20, fontweight='bold')
    
    # Panel C
    print("  Plotting Panel C: PC3 vs PC4")
    ax_c = fig.add_subplot(gs[1, 0])
    plot_panel_c(ax_c, data_dict)
    ax_c.text(-0.12, 1.08, 'C', transform=ax_c.transAxes,
              fontsize=20, fontweight='bold')
    
    # Panel D - simple version
    print("  Plotting Panel D: Global Context with SGDP")
    ax_d = fig.add_subplot(gs[1, 1])
    plot_panel_d(ax_d, data_dict)
    ax_d.text(-0.12, 1.08, 'D', transform=ax_d.transAxes,
              fontsize=20, fontweight='bold')
    
    # Main title - positioned to avoid overlap
    fig.suptitle('Figure 4. Principal Component Analysis (PCA) of Peruvian Populations\n' +
                 'after removing related individuals (746 samples)',
                 fontsize=16, fontweight='bold', y=0.965)  # Explicit y position
    
    # Save outputs
    print("\nSaving outputs...")
    
    output_png = os.path.join(BASE_DIR, 'Figure4_PCA_composite.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  ✓ PNG saved: {output_png}")
    
    output_pdf = os.path.join(BASE_DIR, 'Figure4_PCA_composite.pdf')
    plt.savefig(output_pdf, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  ✓ PDF saved: {output_pdf}")
    
    plt.close()
    
    # Final summary
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
    print("Figure 4 has been successfully created!")
    print(f"Output files saved to: {BASE_DIR}/")
    print("TSV diagnostics saved:")
    print(f"  - {BASE_DIR}/sgdp_iid_to_pop_resolved.tsv")
    if os.path.exists(os.path.join(BASE_DIR, 'sgdp_unmatched.tsv')):
        print(f"  - {BASE_DIR}/sgdp_unmatched.tsv")
    print("="*70)

if __name__ == "__main__":
    main()