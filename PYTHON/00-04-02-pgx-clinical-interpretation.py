#!/usr/bin/env python3
"""
00-04-02-pgx-clinical-interpretation.py

Clinical PGx figure builder for pharmacogenomics clinical interpretation.
Generates publication-ready figures in Nature Genetics production style.

Usage:
    python 00-04-02-pgx-clinical-interpretation.py

All outputs saved to ANALYSIS/00-05-PyPGx/clinical/
"""

import os
import sys
import warnings
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Set, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings('ignore')

# ============================================================================
# NATURE GENETICS PRODUCTION-QUALITY MATPLOTLIB CONFIGURATION
# ============================================================================

# Set Nature Genetics publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 16  # Nature standard for panel titles
plt.rcParams['axes.labelsize'] = 14   # Nature standard for axis labels
plt.rcParams['xtick.labelsize'] = 12  # Nature standard for tick labels
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['figure.titlesize'] = 18  # Main figure title
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.edgecolor'] = '#666666'
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['patch.linewidth'] = 0.8
plt.rcParams['savefig.dpi'] = 600  # Nature standard DPI

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE = Path("ANALYSIS/00-05-PyPGx")
RESULTS_DIR = BASE / "results"
CLINICAL_DIR = BASE / "clinical"
CLINICAL_DIR.mkdir(parents=True, exist_ok=True)

# Input files from 00-04-00
PGX_CLASSIFICATION = RESULTS_DIR / "PGx_classification.csv"
PGX_POPULATION_SUMMARY = RESULTS_DIR / "PGx_population_summary.csv"
FDA_PATIENT_MATCHES = RESULTS_DIR / "FDA_patient_matches.csv"
FDA_TRIGGER_COUNTS = RESULTS_DIR / "FDA_trigger_counts.csv"
WGS_SAMPLES = BASE / "wgs_samples.txt"
COHORT_RESULTS = RESULTS_DIR / "cohort_pypgx_results.tsv"

# Clinical output files
CLINICAL_FIGURE = CLINICAL_DIR / "clinical_pgx_figure.png"
CLINICAL_FIGURE_PDF = CLINICAL_DIR / "clinical_pgx_figure.pdf"
LOG_FILE = CLINICAL_DIR / "clinical_build.log"

# Panel outputs
PANEL_A_CSV = CLINICAL_DIR / "panelA_actionable_frequencies.csv"
PANEL_A_PNG = CLINICAL_DIR / "panelA_actionable_frequencies.png"
PANEL_B_CSV = CLINICAL_DIR / "panelB_fda_trigger_rates.csv"
PANEL_B_PNG = CLINICAL_DIR / "panelB_fda_trigger_rates.png"
PANEL_C_CSV = CLINICAL_DIR / "panelC_population_burden.csv"
PANEL_C_PNG = CLINICAL_DIR / "panelC_population_burden.png"

# FDA-linked genes of clinical importance
FDA_LINKED_GENES = ['CYP2C19', 'CYP2B6', 'CYP3A5', 'SLCO1B1', 'TPMT', 'NUDT15', 'DPYD', 'UGT1A1']

# Nature-style Okabe-Ito colorblind-safe palette (muted)
OKABE_ITO = {
    'vermillion': '#CC6677',  # Muted red - Poor/High risk
    'orange': '#DDAA33',       # Muted orange - Intermediate risk  
    'bluish_green': '#88CCAA', # Muted teal - Normal
    'sky_blue': '#88BBDD',    # Muted sky blue - Rapid
    'reddish_purple': '#AA4488', # Muted purple - Ultrarapid
    'navy': '#332288',         # Navy - WGS
    'teal': '#44AA99',         # Teal - Cohort
    'yellow': '#DDCC77',       # Muted yellow
    'grey': '#999999'          # Grey
}

# Phenotype colors using muted palette
PHENOTYPE_COLORS = {
    'Poor Metabolizer': OKABE_ITO['vermillion'],
    'Intermediate Metabolizer': OKABE_ITO['orange'],
    'Normal Metabolizer': OKABE_ITO['bluish_green'],
    'Rapid Metabolizer': OKABE_ITO['sky_blue'],
    'Ultrarapid Metabolizer': OKABE_ITO['reddish_purple'],
    'Poor Function': OKABE_ITO['vermillion'],
    'Decreased Function': OKABE_ITO['orange'],
    'Normal Function': OKABE_ITO['bluish_green'],
    'Increased Function': OKABE_ITO['sky_blue']
}

# Phenotype order
PHENOTYPE_ORDER = [
    'Poor Metabolizer',
    'Intermediate Metabolizer', 
    'Normal Metabolizer',
    'Rapid Metabolizer',
    'Ultrarapid Metabolizer',
    'Poor Function',
    'Decreased Function',
    'Normal Function',
    'Increased Function'
]

# Category abbreviations for Panel A
CATEGORY_ABBREVIATIONS = {
    'Amazon_Indigenous': 'Amazon-Ind',
    'Amazon_Mestizo': 'Amazon-Mes',
    'Andes_Indigenous': 'Andes-Ind',
    'Andes_Mestizo': 'Andes-Mes',
    'Coast_Indigenous': 'Coast-Ind',
    'Coast_Mestizo': 'Coast-Mes'
}

# Define which phenotypes are clinically actionable per gene
ACTIONABLE_PHENOTYPES = {
    'CYP2C19': {'Poor Metabolizer', 'Intermediate Metabolizer'},
    'CYP2B6': {'Poor Metabolizer', 'Intermediate Metabolizer'},
    'CYP3A5': {'Poor Metabolizer', 'Intermediate Metabolizer'},
    'SLCO1B1': {'Poor Function', 'Decreased Function'},
    'TPMT': {'Poor Metabolizer', 'Intermediate Metabolizer', 'Rapid Metabolizer'},
    'NUDT15': {'Poor Metabolizer', 'Intermediate Metabolizer', 'Rapid Metabolizer'},
    'DPYD': {'Poor Metabolizer', 'Intermediate Metabolizer'},
    'UGT1A1': {'Poor Metabolizer', 'Intermediate Metabolizer'}
}

# ============================================================================
# STATISTICAL UTILITIES
# ============================================================================

def wilson_ci(successes, n, alpha=0.05):
    """Calculate Wilson confidence interval for proportion."""
    if n == 0:
        return 0, 0
    lower, upper = proportion_confint(successes, n, alpha=alpha, method='wilson')
    return lower * 100, upper * 100  # Convert to percentage

def two_proportion_test(x1, n1, x2, n2):
    """Two-proportion z-test."""
    if n1 == 0 or n2 == 0:
        return 1.0
    
    p1 = x1 / n1
    p2 = x2 / n2
    p_pool = (x1 + x2) / (n1 + n2)
    
    if p_pool == 0 or p_pool == 1:
        return 1.0
    
    se = np.sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
    if se == 0:
        return 1.0
    
    z = (p1 - p2) / se
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    return p_value

# ============================================================================
# SETUP
# ============================================================================

def setup_logging():
    """Configure logging for clinical analysis."""
    if LOG_FILE.exists():
        LOG_FILE.unlink()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE, mode='w'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def check_required_inputs(logger):
    """Check that required input files exist."""
    required_files = {
        'PGx_classification.csv': PGX_CLASSIFICATION,
        'PGx_population_summary.csv': PGX_POPULATION_SUMMARY,
        'FDA_patient_matches.csv': FDA_PATIENT_MATCHES,
        'wgs_samples.txt': WGS_SAMPLES
    }
    
    missing = []
    for name, path in required_files.items():
        if not path.exists():
            missing.append(name)
            logger.critical(f"Missing required input: {path}")
    
    if missing:
        logger.critical(f"Cannot proceed - missing {len(missing)} required files")
        logger.critical("Please run 00-04-00-run-pypgx.py first")
        sys.exit(1)
    
    logger.info("All required input files present")
    return True

# ============================================================================
# PANEL A: ACTIONABLE PHENOTYPE FREQUENCIES
# ============================================================================

def process_panel_a_data(logger):
    """Process data for Panel A and return structured data."""
    
    if not PGX_POPULATION_SUMMARY.exists():
        logger.warning("Cannot process Panel A - population summary missing")
        return None
    
    # Load data
    pop_summary = pd.read_csv(PGX_POPULATION_SUMMARY)
    
    # Filter to FDA-linked genes only
    fda_genes = [g for g in FDA_LINKED_GENES if g in pop_summary.Gene.unique()]
    
    if not fda_genes:
        logger.warning("No FDA-linked genes found in population summary")
        return None
    
    logger.info(f"Panel A: Processing {len(fda_genes)} FDA-linked genes: {', '.join(fda_genes)}")
    
    panel_data = pop_summary[pop_summary.Gene.isin(fda_genes)].copy()
    
    # Standardize phenotype names
    phenotype_map = {
        'poor metabolizer': 'Poor Metabolizer',
        'intermediate metabolizer': 'Intermediate Metabolizer',
        'normal metabolizer': 'Normal Metabolizer',
        'extensive metabolizer': 'Normal Metabolizer',
        'rapid metabolizer': 'Rapid Metabolizer',
        'ultrarapid metabolizer': 'Ultrarapid Metabolizer',
        'ultra-rapid metabolizer': 'Ultrarapid Metabolizer',
        'poor function': 'Poor Function',
        'decreased function': 'Decreased Function',
        'normal function': 'Normal Function',
        'increased function': 'Increased Function'
    }
    
    if 'Phenotype' in panel_data.columns:
        panel_data['Phenotype'] = panel_data['Phenotype'].apply(
            lambda x: phenotype_map.get(str(x).lower(), str(x).title()) if pd.notna(x) else 'Unknown'
        )
    
    # Ensure Frequency column is numeric
    if 'Frequency' in panel_data.columns:
        panel_data['Frequency'] = pd.to_numeric(panel_data['Frequency'], errors='coerce').fillna(0)
    else:
        logger.warning("Frequency column not found in population summary")
        return None
    
    # Focus on Region×Group
    region_group_data = panel_data[panel_data['Level'] == 'Region_Group'].copy()
    
    # Prepare data for plotting
    plot_data = {}
    csv_data = []
    
    for gene in fda_genes:
        gene_data = region_group_data[region_group_data.Gene == gene].copy()
        
        if not gene_data.empty:
            # Pivot for stacked bar
            pivot = gene_data.pivot_table(
                index='Category', 
                columns='Phenotype', 
                values='Frequency',
                aggfunc='mean',
                fill_value=0
            )
            
            # Get all phenotypes present
            all_phenotypes = pivot.columns.tolist()
            if all_phenotypes:
                # Order phenotypes
                cols_ordered = [p for p in PHENOTYPE_ORDER if p in all_phenotypes]
                cols_unknown = [p for p in all_phenotypes if p not in PHENOTYPE_ORDER]
                cols_present = cols_ordered + cols_unknown
                
                pivot = pivot[cols_present]
                pivot_percent = pivot * 100
                
                # Get totals for each category
                totals = {}
                for category in pivot.index:
                    total = gene_data[gene_data['Category'] == category]['Total'].iloc[0]
                    totals[category] = total
                
                plot_data[gene] = {
                    'data': pivot_percent,
                    'totals': totals,
                    'phenotypes': cols_present
                }
                
                # Save data for CSV
                for _, row in gene_data.iterrows():
                    csv_data.append({
                        'Stratum': row['Category'],
                        'Gene': gene,
                        'Phenotype': row['Phenotype'],
                        'Percent': row['Frequency'] * 100,
                        'N': row['Total']
                    })
    
    # Save CSV once
    if csv_data:
        tidy_df = pd.DataFrame(csv_data)
        tidy_df.to_csv(PANEL_A_CSV, index=False)
        logger.info(f"Panel A: Saved tidy CSV to {PANEL_A_CSV}")
    
    return {'genes': fda_genes, 'plot_data': plot_data}

def plot_panel_a(ax, data, logger, standalone=False, fig=None):
    """Plot Panel A with Nature Genetics production styling."""
    
    if not data:
        if not standalone:
            ax.text(0.5, 0.5, 'Panel A: Data not available', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=12)
        return []
    
    fda_genes = data['genes']
    plot_data = data['plot_data']
    
    # Create subplots for each gene with shared axes
    n_genes = len(fda_genes)
    
    if standalone and fig:
        axes = []
        ax_ref = None
        for i in range(n_genes):
            if ax_ref is None:
                ax_gene = fig.add_subplot(n_genes, 1, i+1)
                ax_ref = ax_gene
            else:
                ax_gene = fig.add_subplot(n_genes, 1, i+1, sharex=ax_ref, sharey=ax_ref)
            axes.append(ax_gene)
    elif not standalone:
        from matplotlib.gridspec import GridSpecFromSubplotSpec
        # Increased vertical spacing for better separation between rows
        gs_a = GridSpecFromSubplotSpec(n_genes, 1, subplot_spec=ax.get_subplotspec(), 
                                       hspace=0.45)  # Increased from 0.25 for more space
        axes = []
        ax_ref = None
        for i in range(n_genes):
            if ax_ref is None:
                ax_gene = plt.subplot(gs_a[i])
                ax_ref = ax_gene
            else:
                ax_gene = plt.subplot(gs_a[i], sharex=ax_ref, sharey=ax_ref)
            axes.append(ax_gene)
    else:
        return []
    
    # Set consistent y-axis limits and ticks on reference axis
    if ax_ref:
        ax_ref.set_ylim(-5, 105)
        ax_ref.set_yticks([0, 25, 50, 75, 100])
    
    for i, gene in enumerate(fda_genes):
        ax_gene = axes[i]
        
        if gene in plot_data:
            gene_info = plot_data[gene]
            # Make a copy to avoid modifying original data
            pivot_percent = gene_info['data'].copy()
            totals = gene_info['totals']
            cols_present = gene_info['phenotypes']
            
            # Save original index before applying abbreviations
            original_index = pivot_percent.index.tolist()
            pivot_percent.index = [CATEGORY_ABBREVIATIONS.get(idx, idx) for idx in pivot_percent.index]
            
            if not pivot_percent.empty and pivot_percent.sum().sum() > 0:
                # Create stacked bar chart with muted colors
                colors = [PHENOTYPE_COLORS.get(p, OKABE_ITO['grey']) for p in cols_present]
                pivot_percent.plot(kind='bar', stacked=True, ax=ax_gene,
                                  color=colors,
                                  edgecolor='white', linewidth=0.5, width=0.6,
                                  legend=False)
                
                # Add percentage labels only for segments ≥15% with better contrast
                for container in ax_gene.containers:
                    labels = []
                    for v in container.datavalues:
                        if v >= 15:
                            labels.append(f'{v:.0f}')
                        else:
                            labels.append('')
                    ax_gene.bar_label(container, labels=labels, label_type='center', 
                                      fontsize=9, weight='normal', color='white')
                
                # Add N under each bar in grey italics - use original index for lookup
                for j, orig_cat in enumerate(original_index):
                    ax_gene.text(j, -8, f'n={int(totals[orig_cat])}',
                                ha='center', fontsize=9, color='#666666', 
                                style='italic')
                
                # Gene title only - reduced padding to bring it closer to plot
                ax_gene.set_title(f'{gene}', fontsize=13, fontweight='bold', pad=4)
                ax_gene.set_xlabel('')
                ax_gene.set_ylabel('')  # No individual y-labels
                
                # Hide y-axis labels/ticks on all subplots
                ax_gene.tick_params(axis='y', labelleft=False)
                
                # X-axis handling: only show labels on bottom subplot
                if i < n_genes - 1:
                    # Hide x tick labels for all but the bottom row
                    ax_gene.tick_params(axis='x', labelbottom=False)
                    ax_gene.set_xlabel('')
                else:
                    # Bottom row: show horizontal tick labels
                    ax_gene.tick_params(axis='x', labelbottom=True)
                    labels = ax_gene.get_xticklabels()
                    ax_gene.set_xticklabels(labels, rotation=0, ha='center', fontsize=11)
                
                # Subtle grid for Nature style
                ax_gene.grid(True, axis='y', alpha=0.1, linestyle='-', linewidth=0.5)
                ax_gene.set_axisbelow(True)
                ax_gene.spines['top'].set_visible(False)
                ax_gene.spines['right'].set_visible(False)
                ax_gene.spines['left'].set_linewidth(0.8)
                ax_gene.spines['bottom'].set_linewidth(0.8)
                
            else:
                ax_gene.text(0.5, 0.5, f'{gene}: No data',
                           ha='center', va='center', transform=ax_gene.transAxes, fontsize=10)
                ax_gene.set_title(f'{gene}', fontsize=13, fontweight='bold')
                ax_gene.axis('off')
        else:
            ax_gene.text(0.5, 0.5, f'{gene}: No data',
                       ha='center', va='center', transform=ax_gene.transAxes, fontsize=10)
            ax_gene.set_title(f'{gene}', fontsize=13, fontweight='bold')
            ax_gene.axis('off')
    
    # Get current figure for sup labels
    current_fig = plt.gcf()
    
    if not standalone:
        # Panel title in Nature style
        ax.set_title('A: Actionable phenotype frequencies (n=736)', 
                    fontweight='bold', fontsize=16, pad=30, loc='left')
        ax.axis('off')
        
        # Add single Y-axis label spanning all rows - closer to graphs
        # Position it closer to the left edge of the panel
        ax.text(-0.03, 0.5, 'Frequency (%)', transform=ax.transAxes, 
                rotation=90, va='center', ha='center', fontsize=12)
        
        # Single centered legend below panel
        if plot_data:
            legend_handles = []
            phenotypes_used = set()
            for gene_data in plot_data.values():
                phenotypes_used.update(gene_data['phenotypes'])
            
            for phenotype in PHENOTYPE_ORDER:
                if phenotype in phenotypes_used and phenotype in PHENOTYPE_COLORS:
                    legend_handles.append(mpatches.Patch(color=PHENOTYPE_COLORS[phenotype], 
                                                        label=phenotype))
            
            # Centered horizontal legend with proper spacing
            ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.08),
                     ncol=min(5, len(legend_handles)), frameon=False, fontsize=11, 
                     columnspacing=1.5, handletextpad=0.8)
    else:
        # For standalone figure, add the single axis labels
        current_fig.supylabel('Frequency (%)', x=0.02, fontsize=12)
        current_fig.supxlabel('Region × Group', y=0.02, fontsize=12)
    
    return fda_genes

# ============================================================================
# PANEL B: FDA TRIGGER PREVALENCE WITH STATISTICS
# ============================================================================

def process_panel_b_data(logger):
    """Process data for Panel B and return structured data."""
    
    if not FDA_PATIENT_MATCHES.exists() or not WGS_SAMPLES.exists():
        logger.warning("Required files missing for Panel B")
        return None
    
    # Constants
    WGS_N = 109
    COHORT_N = 736
    
    # Load data
    wgs_ids = set(open(WGS_SAMPLES).read().splitlines())
    matches = pd.read_csv(FDA_PATIENT_MATCHES)
    
    # Filter to FDA-linked genes
    fda_matches = matches[matches.Gene.isin(FDA_LINKED_GENES)].copy()
    
    if fda_matches.empty:
        logger.warning("No FDA triggers for linked genes")
        return None
    
    # Calculate trigger counts, rates, and statistics
    trigger_data = []
    
    for gene in FDA_LINKED_GENES:
        gene_matches = fda_matches[fda_matches.Gene == gene]
        
        if not gene_matches.empty:
            wgs_count = gene_matches[gene_matches.SampleID.isin(wgs_ids)].shape[0]
            cohort_count = gene_matches.shape[0]
            
            wgs_rate = (wgs_count / WGS_N) * 100
            cohort_rate = (cohort_count / COHORT_N) * 100
            
            # Calculate Wilson CIs
            wgs_lcl, wgs_ucl = wilson_ci(wgs_count, WGS_N)
            cohort_lcl, cohort_ucl = wilson_ci(cohort_count, COHORT_N)
            
            # Two-proportion test
            p_value = two_proportion_test(wgs_count, WGS_N, cohort_count, COHORT_N)
            
            trigger_data.append({
                'Gene': gene,
                'WGS_Count': wgs_count,
                'WGS_N': WGS_N,
                'WGS_Rate': wgs_rate,
                'WGS_Rate_LCL': wgs_lcl,
                'WGS_Rate_UCL': wgs_ucl,
                'Cohort_Count': cohort_count,
                'Cohort_N': COHORT_N,
                'Cohort_Rate': cohort_rate,
                'Cohort_Rate_LCL': cohort_lcl,
                'Cohort_Rate_UCL': cohort_ucl,
                'P_WGS_vs_Cohort': p_value
            })
    
    if not trigger_data:
        logger.warning("No trigger data to process")
        return None
    
    df = pd.DataFrame(trigger_data)
    df = df.sort_values('Cohort_Rate', ascending=False)  # Sort descending
    
    # Save tidy CSV
    tidy_df = []
    for _, row in df.iterrows():
        tidy_df.append({
            'Gene': row.Gene,
            'Sample_Type': 'WGS',
            'N': WGS_N,
            'Count': row.WGS_Count,
            'Rate_Per_100': row.WGS_Rate,
            'Rate_LCL': row.WGS_Rate_LCL,
            'Rate_UCL': row.WGS_Rate_UCL,
            'P_WGS_vs_Cohort': row.P_WGS_vs_Cohort
        })
        tidy_df.append({
            'Gene': row.Gene,
            'Sample_Type': 'Cohort',
            'N': COHORT_N,
            'Count': row.Cohort_Count,
            'Rate_Per_100': row.Cohort_Rate,
            'Rate_LCL': row.Cohort_Rate_LCL,
            'Rate_UCL': row.Cohort_Rate_UCL,
            'P_WGS_vs_Cohort': row.P_WGS_vs_Cohort
        })
    
    pd.DataFrame(tidy_df).to_csv(PANEL_B_CSV, index=False)
    logger.info(f"Panel B: Saved tidy CSV to {PANEL_B_CSV}")
    logger.info(f"Panel B: Processed {len(df)} genes with FDA triggers")
    
    return df

def plot_panel_b(ax, data, logger):
    """Plot Panel B with Nature Genetics production styling."""
    
    if data is None or data.empty:
        ax.text(0.5, 0.5, 'Panel B: Data not available', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        return False
    
    df = data.copy()  # Make a copy to avoid modifying original
    
    # Calculate x-axis limits and clip CIs to plotting window
    x_max = max(df.WGS_Rate_UCL.max(), df.Cohort_Rate_UCL.max())
    x_max = np.ceil(x_max + 5)  # Add margin
    
    # Clip CI bounds to prevent extreme whiskers
    df['WGS_Rate_LCL'] = np.clip(df.WGS_Rate_LCL, 0, x_max)
    df['WGS_Rate_UCL'] = np.clip(df.WGS_Rate_UCL, 0, x_max)
    df['Cohort_Rate_LCL'] = np.clip(df.Cohort_Rate_LCL, 0, x_max)
    df['Cohort_Rate_UCL'] = np.clip(df.Cohort_Rate_UCL, 0, x_max)
    
    # Compute tidy CI half-widths once
    wgs_err_low = df.WGS_Rate - df.WGS_Rate_LCL
    wgs_err_high = df.WGS_Rate_UCL - df.WGS_Rate
    coh_err_low = df.Cohort_Rate - df.Cohort_Rate_LCL
    coh_err_high = df.Cohort_Rate_UCL - df.Cohort_Rate
    
    # Create grouped bar plot with wider bars
    genes = df.Gene.tolist()
    y_pos = np.arange(len(genes))
    bar_width = 0.35  # Wider bars for visibility
    
    # Plot bars with Nature-style colors (lower z-order)
    bars1 = ax.barh(y_pos - bar_width/2, df.WGS_Rate, bar_width,
                    label='WGS (n=109)', color=OKABE_ITO['navy'], 
                    edgecolor='white', linewidth=0.5, alpha=0.9, zorder=2)
    bars2 = ax.barh(y_pos + bar_width/2, df.Cohort_Rate, bar_width,
                    label='Cohort (n=736)', color=OKABE_ITO['teal'], 
                    edgecolor='white', linewidth=0.5, alpha=0.9, zorder=2)
    
    # Draw error bars with elegant styling (higher z-order)
    for i in range(len(df)):
        # Determine if this is a tiny rate needing delicate CIs
        wgs_is_tiny = (df.iloc[i].WGS_Rate < 1.0) and (df.iloc[i].WGS_Rate_UCL <= 3.0)
        coh_is_tiny = (df.iloc[i].Cohort_Rate < 1.0) and (df.iloc[i].Cohort_Rate_UCL <= 3.0)
        
        # WGS error bars
        if wgs_is_tiny:
            ax.errorbar(df.iloc[i].WGS_Rate, y_pos[i] - bar_width/2,
                       xerr=[[wgs_err_low.iloc[i]], [wgs_err_high.iloc[i]]],
                       fmt='none', ecolor='#444444', elinewidth=1.2, capsize=3, capthick=1.2,
                       alpha=0.8, zorder=3)
        else:
            ax.errorbar(df.iloc[i].WGS_Rate, y_pos[i] - bar_width/2,
                       xerr=[[wgs_err_low.iloc[i]], [wgs_err_high.iloc[i]]],
                       fmt='none', ecolor='#444444', elinewidth=1.8, capsize=4, capthick=1.8,
                       alpha=0.8, zorder=3)
        
        # Cohort error bars
        if coh_is_tiny:
            ax.errorbar(df.iloc[i].Cohort_Rate, y_pos[i] + bar_width/2,
                       xerr=[[coh_err_low.iloc[i]], [coh_err_high.iloc[i]]],
                       fmt='none', ecolor='#444444', elinewidth=1.2, capsize=3, capthick=1.2,
                       alpha=0.8, zorder=3)
        else:
            ax.errorbar(df.iloc[i].Cohort_Rate, y_pos[i] + bar_width/2,
                       xerr=[[coh_err_low.iloc[i]], [coh_err_high.iloc[i]]],
                       fmt='none', ecolor='#444444', elinewidth=1.8, capsize=4, capthick=1.8,
                       alpha=0.8, zorder=3)
    
    # Add rates inside or at end of bars with consistent alignment
    for i, row in df.iterrows():
        # WGS rate
        if row.WGS_Rate >= 3:
            ax.text(row.WGS_Rate/2, y_pos[i] - bar_width/2, f'{row.WGS_Rate:.1f}',
                   va='center', ha='center', fontsize=10, color='white', weight='bold')
        else:
            label_x = min(row.WGS_Rate + 0.4, x_max - 1.0)
            ax.text(label_x, y_pos[i] - bar_width/2, f'{row.WGS_Rate:.1f}',
                   va='center', ha='left', fontsize=10, color='#333333')
        
        # Cohort rate
        if row.Cohort_Rate >= 3:
            ax.text(row.Cohort_Rate/2, y_pos[i] + bar_width/2, f'{row.Cohort_Rate:.1f}',
                   va='center', ha='center', fontsize=10, color='white', weight='bold')
        else:
            label_x = min(row.Cohort_Rate + 0.4, x_max - 1.0)
            ax.text(label_x, y_pos[i] + bar_width/2, f'{row.Cohort_Rate:.1f}',
                   va='center', ha='left', fontsize=10, color='#333333')
    
    # Add p-values close to bars, avoiding random floating
    for i, row in df.iterrows():
        if row.P_WGS_vs_Cohort < 0.001:
            p_text = 'p<0.001'
        elif row.P_WGS_vs_Cohort < 0.05:
            p_text = f'p={row.P_WGS_vs_Cohort:.3f}'
        else:
            p_text = f'p={row.P_WGS_vs_Cohort:.2f}'
        
        # Position just above cohort bar end
        p_x = min(row.Cohort_Rate + 0.5, x_max - 0.5)
        p_y = y_pos[i] + bar_width * 0.45
        ax.text(p_x, p_y, p_text, 
               va='bottom', ha='left', fontsize=9, style='italic', 
               color='#666666')
    
    # Nature-style formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=12)
    ax.set_xlabel('Trigger rate per 100 individuals', fontsize=14)
    ax.set_xlim(0, x_max)
    ax.set_title('B: FDA trigger rates per 100 individuals', 
                fontweight='bold', fontsize=16, pad=20, loc='left')
    
    # Elegant legend
    ax.legend(loc='lower right', frameon=False, fontsize=11)
    
    # Subtle grid with improved readability
    ax.grid(True, axis='x', alpha=0.15, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    
    return True

# ============================================================================
# PANEL C: POPULATION BURDEN WITH NATURE-STYLE HEATMAP
# ============================================================================

def process_panel_c_data(logger):
    """Process data for Panel C with fallback strategy."""
    
    if not PGX_POPULATION_SUMMARY.exists():
        logger.warning("Population summary file missing for Panel C")
        return None
    
    # Load data
    pop_summary = pd.read_csv(PGX_POPULATION_SUMMARY)
    
    # Standardize phenotype names
    phenotype_map = {
        'poor metabolizer': 'Poor Metabolizer',
        'intermediate metabolizer': 'Intermediate Metabolizer',
        'normal metabolizer': 'Normal Metabolizer',
        'extensive metabolizer': 'Normal Metabolizer',
        'rapid metabolizer': 'Rapid Metabolizer',
        'ultrarapid metabolizer': 'Ultrarapid Metabolizer',
        'ultra-rapid metabolizer': 'Ultrarapid Metabolizer',
        'poor function': 'Poor Function',
        'decreased function': 'Decreased Function',
        'normal function': 'Normal Function',
        'increased function': 'Increased Function'
    }
    
    pop_summary['Phenotype'] = pop_summary['Phenotype'].apply(
        lambda x: phenotype_map.get(str(x).lower(), str(x).title()) if pd.notna(x) else 'Unknown'
    )
    
    # Filter to FDA-linked genes
    pop_summary = pop_summary[pop_summary.Gene.isin(FDA_LINKED_GENES)].copy()
    
    # Try Population level first
    pop_level_data = pop_summary[pop_summary['Level'] == 'Population'].copy()
    
    # Determine which level to use and minimum N
    level_used = 'Population'
    min_n = 15
    data_to_use = pop_level_data
    
    if not pop_level_data.empty:
        # Check if we have populations with N >= 15
        pop_totals = pop_level_data.groupby('Category')['Total'].first()
        valid_pops = pop_totals[pop_totals >= min_n].index.tolist()
        
        if not valid_pops:
            # Try min_n = 10
            logger.info("No populations with N>=15, trying N>=10")
            min_n = 10
            valid_pops = pop_totals[pop_totals >= min_n].index.tolist()
            
            if not valid_pops:
                # Fall back to Region×Group
                logger.info("No populations with N>=10, falling back to Region×Group")
                level_used = 'Region_Group'
                data_to_use = pop_summary[pop_summary['Level'] == 'Region_Group'].copy()
                min_n = 0
    else:
        # No population data, use Region×Group
        logger.info("No Population level data, using Region×Group")
        level_used = 'Region_Group'
        data_to_use = pop_summary[pop_summary['Level'] == 'Region_Group'].copy()
        min_n = 0
    
    if data_to_use.empty:
        logger.warning("No data available for Panel C")
        return None
    
    # Filter by minimum N if applicable
    if min_n > 0:
        category_totals = data_to_use.groupby('Category')['Total'].first()
        valid_categories = category_totals[category_totals >= min_n].index.tolist()
        data_to_use = data_to_use[data_to_use.Category.isin(valid_categories)]
    
    if data_to_use.empty:
        logger.warning(f"No {level_used} with N>={min_n}")
        return None
    
    # Calculate burden percentages
    burden_data = []
    
    for category in data_to_use.Category.unique():
        category_data = data_to_use[data_to_use.Category == category]
        
        for gene in FDA_LINKED_GENES:
            gene_data = category_data[category_data.Gene == gene]
            
            if not gene_data.empty:
                # Get actionable phenotypes for this gene
                actionable = ACTIONABLE_PHENOTYPES.get(gene, set())
                
                # Sum frequencies for actionable phenotypes
                actionable_freq = gene_data[gene_data.Phenotype.isin(actionable)]['Frequency'].sum()
                
                # Get total N for this category-gene combination
                total_n = gene_data['Total'].iloc[0] if 'Total' in gene_data.columns else 0
                
                burden_pct = actionable_freq * 100  # Already a proportion, convert to %
                
                burden_data.append({
                    'Population': category,
                    'Gene': gene,
                    'BurdenPercent': burden_pct,
                    'N': total_n
                })
    
    if not burden_data:
        logger.warning("No burden data calculated")
        return None
    
    burden_df = pd.DataFrame(burden_data)
    
    # Create pivot table for heatmap
    pivot = burden_df.pivot_table(
        index='Population',
        columns='Gene', 
        values='BurdenPercent',
        fill_value=0
    )
    
    # Sort populations by total burden
    pivot['Total_Burden'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('Total_Burden', ascending=False)
    pivot = pivot.drop('Total_Burden', axis=1)
    
    # Save tidy CSV
    burden_df.to_csv(PANEL_C_CSV, index=False)
    logger.info(f"Panel C: Saved tidy CSV to {PANEL_C_CSV}")
    logger.info(f"Panel C: Used {level_used} level with min_n={min_n}")
    logger.info(f"Panel C: Analyzed {len(pivot.index)} {level_used}s")
    
    return {
        'pivot': pivot,
        'level_used': level_used,
        'min_n': min_n,
        'n_categories': len(pivot.index)
    }

def plot_panel_c(ax, data, logger):
    """Plot Panel C heatmap with Nature Genetics production styling."""
    
    if data is None:
        ax.text(0.5, 0.5, 'Panel C: Data not available', 
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        return 0, None
    
    pivot = data['pivot']
    level_used = data['level_used']
    min_n = data['min_n']
    
    # Create soft sequential colormap (yellow→orange→teal) capped at 80%
    colors = ['#FFFFEE', '#FFF7C7', '#FEDA77', '#FDAE61', '#F46D43', 
              '#D73027', '#A50026', '#762A83', '#5B4394']
    n_bins = 100
    cmap = sns.blend_palette(colors, n_colors=n_bins, as_cmap=True)
    
    # Prepare annotations with adaptive coloring
    annot = pivot.copy()
    annot_array = annot.values
    annot_text = np.zeros_like(annot_array, dtype=object)
    annot_colors = np.zeros_like(annot_array, dtype=object)
    
    for i in range(annot_array.shape[0]):
        for j in range(annot_array.shape[1]):
            val = annot_array[i, j]
            if val > 0:
                annot_text[i, j] = f'{val:.0f}'
                # Adaptive text color: white for dark cells, black for light
                annot_colors[i, j] = 'white' if val > 40 else 'black'
            else:
                annot_text[i, j] = ''
                annot_colors[i, j] = 'black'
    
    # Create heatmap with soft cap at 80% for better contrast
    vmax = min(80, pivot.max().max())
    
    sns.heatmap(pivot, annot=annot_text, fmt='', cmap=cmap,
                cbar_kws={'label': 'Burden (%)', 'shrink': 0.7, 'pad': 0.02}, 
                linewidths=0.5, linecolor='#E6E6E6',  # Light grey borders
                ax=ax, vmin=0, vmax=vmax, 
                square=False, annot_kws={'fontsize': 12, 'fontweight': 'normal'})
    
    # Apply adaptive text colors manually
    for i, (text, color_row) in enumerate(zip(ax.texts, annot_colors.flatten())):
        text.set_color(color_row)
    
    # Nature-style formatting
    ax.set_xlabel('Gene', fontsize=14, labelpad=12, fontweight='normal')
    ax.set_ylabel('Population', fontsize=14, labelpad=12, fontweight='bold')
    
    title = 'C: Population burden of actionable phenotypes (%)'
    if min_n > 0:
        title += f' (n≥{min_n})'
    ax.set_title(title, fontweight='bold', fontsize=16, pad=20, loc='left')
    
    # Rotate x-axis labels 30° with clear spacing
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=12)
    # Keep y-axis labels bold, left-aligned, readable
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=12, fontweight='bold')
    ax.tick_params(axis='both', which='major', length=0)
    
    return len(pivot.index), level_used

# ============================================================================
# MAIN CLINICAL FIGURE
# ============================================================================

def create_clinical_figure():
    """Create Nature Genetics production-quality clinical PGx figure."""
    
    logger = setup_logging()
    logger.info("="*70)
    logger.info("Starting clinical PGx figure generation")
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    logger.info("Style: Nature Genetics production quality")
    logger.info("="*70)
    
    # Check inputs
    check_required_inputs(logger)
    
    # Process all data first (CSVs are saved here)
    logger.info("Processing data for all panels...")
    panel_a_data = process_panel_a_data(logger)
    panel_b_data = process_panel_b_data(logger)
    panel_c_data = process_panel_c_data(logger)
    
    # Create integrated figure with ample whitespace
    logger.info("Creating Nature Genetics production-quality figure...")
    fig = plt.figure(figsize=(18, 22))  # Larger for better spacing
    gs = GridSpec(3, 1, figure=fig, height_ratios=[2.2, 1, 1.2], hspace=0.6)  # More vertical space
    
    # Panel A: Actionable phenotype frequencies
    logger.info("-" * 50)
    ax_a = fig.add_subplot(gs[0])
    genes_included = plot_panel_a(ax_a, panel_a_data, logger)
    
    # Panel B: FDA trigger prevalence
    logger.info("-" * 50)
    ax_b = fig.add_subplot(gs[1])
    panel_b_success = plot_panel_b(ax_b, panel_b_data, logger)
    
    # Panel C: Population burden
    logger.info("-" * 50)
    ax_c = fig.add_subplot(gs[2])
    if panel_c_data:
        n_categories, level_used = plot_panel_c(ax_c, panel_c_data, logger)
    else:
        n_categories, level_used = 0, None
    
    # Overall title - Nature Genetics style
    fig.suptitle('Clinical Pharmacogenomics in Peru',
                 fontsize=18, fontweight='bold', y=0.98, x=0.1, ha='left')
    
    # Adjust layout with generous margins
    plt.tight_layout(rect=[0.03, 0.03, 0.97, 0.95])
    
    # Save high-resolution figures (600 DPI for Nature)
    logger.info("Saving Nature Genetics production-quality figures...")
    plt.savefig(CLINICAL_FIGURE, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(CLINICAL_FIGURE_PDF, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none', format='pdf')
    plt.close()
    
    # Generate standalone panel figures
    logger.info("Generating standalone panel figures...")
    
    # Panel A standalone
    if panel_a_data:
        fig_a = plt.figure(figsize=(14, 12))
        plot_panel_a(None, panel_a_data, logger, standalone=True, fig=fig_a)
        plt.suptitle('Actionable Phenotype Frequencies by Region/Group', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(PANEL_A_PNG, dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()
    
    # Panel B standalone
    if panel_b_data is not None:
        fig_b = plt.figure(figsize=(12, 8))
        ax_b = fig_b.add_subplot(111)
        plot_panel_b(ax_b, panel_b_data, logger)
        plt.tight_layout()
        plt.savefig(PANEL_B_PNG, dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()
    
    # Panel C standalone
    if panel_c_data:
        fig_c = plt.figure(figsize=(12, 10))
        ax_c = fig_c.add_subplot(111)
        plot_panel_c(ax_c, panel_c_data, logger)
        plt.tight_layout()
        plt.savefig(PANEL_C_PNG, dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()
    
    # Final summary
    logger.info("="*70)
    logger.info("SUMMARY:")
    logger.info(f"  Figure styled for Nature Genetics production")
    logger.info(f"  Export DPI: 600")
    logger.info(f"  Genes included: {', '.join(genes_included) if genes_included else 'None'}")
    logger.info(f"  Panel C level: {level_used if level_used else 'None'}")
    logger.info(f"  Panel C categories: {n_categories}")
    logger.info(f"  Panel A: {'Generated' if genes_included else 'Skipped (no data)'}")
    logger.info(f"  Panel B: {'Generated' if panel_b_success else 'Skipped (no data)'}")
    logger.info(f"  Panel C: {'Generated' if n_categories > 0 else 'Skipped (no data)'}")
    logger.info("")
    logger.info("OUTPUT FILES:")
    logger.info(f"  Main figure (PNG): {CLINICAL_FIGURE}")
    logger.info(f"  Main figure (PDF): {CLINICAL_FIGURE_PDF}")
    logger.info(f"  Panel figures: {CLINICAL_DIR}/*.png")
    logger.info(f"  Tidy CSVs: {CLINICAL_DIR}/*.csv")
    logger.info(f"  Log: {LOG_FILE}")
    logger.info("="*70)
    logger.info("Nature-style production rendering applied")
    
    print(f"\n✔ Nature Genetics production-quality clinical PGx figure generated!")
    print(f"  Main outputs (600 DPI):")
    print(f"    • {CLINICAL_FIGURE}")
    print(f"    • {CLINICAL_FIGURE_PDF}")
    print(f"  Individual panels: {CLINICAL_DIR}/panel*.png")
    print(f"  Log: {LOG_FILE}")
    print(f"\n  Clinical insights generated:")
    
    if PANEL_A_CSV.exists():
        print(f"    ✔ {PANEL_A_CSV.name}")
    if PANEL_B_CSV.exists():
        print(f"    ✔ {PANEL_B_CSV.name}")
    if PANEL_C_CSV.exists():
        print(f"    ✔ {PANEL_C_CSV.name}")
    
    print(f"\n  Nature Genetics production style applied:")
    print(f"  • Font hierarchy: 18pt/16pt/14pt/12pt")
    print(f"  • Muted Okabe-Ito colorblind-safe palette")
    print(f"  • Ample whitespace between panels")
    print(f"  • 600 DPI export resolution")
    print(f"  Focus: {len(genes_included)} FDA-linked genes" if genes_included else "")
    print(f"  Coverage: {n_categories} {level_used}s analyzed" if level_used else "")
    print("\n  All functionality preserved with complete statistical analysis")

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    try:
        create_clinical_figure()
    except KeyboardInterrupt:
        print("\n\nFigure generation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error generating clinical figure: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)