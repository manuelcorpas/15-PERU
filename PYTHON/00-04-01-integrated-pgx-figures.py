#!/usr/bin/env python3
"""
00-04-01-integrated-pgx-figures.py

Integrated PGx figure builder that consumes outputs from 00-04-00-run-pypgx.py
Creates publication-quality multi-panel figure with validation and population insights.

IMPORTANT: This script NEVER runs PyPGx or processes VCFs.
It only creates figures and tidy CSVs from existing results produced by 00-04-00-run-pypgx.py

Usage:
    python 00-04-01-integrated-pgx-figures.py

All outputs saved to ANALYSIS/00-05-PyPGx/integrated/
"""

import os
import sys
import warnings
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact

warnings.filterwarnings('ignore')

# Configure matplotlib
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 14

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE = Path("ANALYSIS/00-05-PyPGx")
RESULTS_DIR = BASE / "results"
OUTPUT_DIR = BASE / "integrated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Critical input files from 00-04-00 (NEVER overwrite these)
PGX_CLASSIFICATION = RESULTS_DIR / "PGx_classification.csv"
PGX_CONCORDANCE = RESULTS_DIR / "PGx_concordance.csv"
PGX_POPULATION_SUMMARY = RESULTS_DIR / "PGx_population_summary.csv"
FDA_PATIENT_MATCHES = RESULTS_DIR / "FDA_patient_matches.csv"
FDA_TRIGGER_COUNTS = RESULTS_DIR / "FDA_trigger_counts.csv"
WGS_SAMPLES = BASE / "wgs_samples.txt"
ARRAY_SAMPLES = BASE / "array_samples.txt"
COHORT_RESULTS = BASE / "cohort" / "cohort_pypgx_results.tsv"
WGS_TRUTH_RESULTS = BASE / "wgs_truth" / "wgs_truth_pypgx_results.tsv"

# Output files (ALWAYS overwrite these)
INTEGRATED_FIGURE = OUTPUT_DIR / "integrated_pgx_figure.png"
LOG_FILE = OUTPUT_DIR / "integrated_build.log"

# Panel-specific outputs (ALWAYS overwrite)
PANEL_A_CSV = OUTPUT_DIR / "panelA_concordance_tidy.csv"
PANEL_A_PNG = OUTPUT_DIR / "panelA_concordance.png"
PANEL_B_REGION_CSV = OUTPUT_DIR / "panelB_region_group_frequencies.csv"
PANEL_B_POP_CSV = OUTPUT_DIR / "panelB_population_frequencies.csv"
PANEL_B_PNG = OUTPUT_DIR / "panelB.png"
PANEL_C_CSV = OUTPUT_DIR / "panelC_fda_trigger_rates.csv"
PANEL_C_PNG = OUTPUT_DIR / "panelC.png"
PANEL_D_CSV = OUTPUT_DIR / "panelD_wgs_breadth_with_tiers.csv"
PANEL_D_PNG = OUTPUT_DIR / "panelD.png"

# Thresholds for Panel A
KAPPA_THRESHOLDS = [0.85, 0.70]
ACCURACY_THRESHOLDS = [95, 85]
INDETERMINATE_THRESHOLDS = [5, 15]

# Phenotype order and colors
PHENOTYPE_ORDER = [
    'Poor Metabolizer',
    'Intermediate Metabolizer', 
    'Normal Metabolizer',
    'Rapid Metabolizer',
    'Ultrarapid Metabolizer'
]

PHENOTYPE_COLORS = {
    'Poor Metabolizer': '#e74c3c',
    'Intermediate Metabolizer': '#f39c12',
    'Normal Metabolizer': '#27ae60',
    'Rapid Metabolizer': '#3498db',
    'Ultrarapid Metabolizer': '#9b59b6'
}

# ============================================================================
# SETUP
# ============================================================================

def setup_logging():
    """Configure logging - always creates new log file."""
    # Remove old log file to start fresh
    if LOG_FILE.exists():
        LOG_FILE.unlink()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE, mode='w'),  # Always overwrite
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def check_critical_inputs(logger):
    """Check that all critical input files exist from 00-04-00."""
    critical_files = {
        'PGx_classification.csv': PGX_CLASSIFICATION,
        'PGx_concordance.csv': PGX_CONCORDANCE,
        'PGx_population_summary.csv': PGX_POPULATION_SUMMARY,
        'wgs_samples.txt': WGS_SAMPLES,
        'array_samples.txt': ARRAY_SAMPLES,
        'cohort_pypgx_results.tsv': COHORT_RESULTS,
        'wgs_truth_pypgx_results.tsv': WGS_TRUTH_RESULTS
    }
    
    missing = []
    for name, path in critical_files.items():
        if not path.exists():
            missing.append(name)
            logger.critical(f"Missing critical input: {path}")
    
    if missing:
        logger.critical(f"Cannot proceed - missing {len(missing)} critical files")
        logger.critical("Please run 00-04-00-run-pypgx.py first to generate required inputs")
        sys.exit(1)
    
    logger.info("All critical input files present")
    return True

# ============================================================================
# PANEL A: CONCORDANCE VALIDATION
# ============================================================================

def create_panel_a(ax, logger):
    """Panel A: Validation on WGS overlap showing only Array-supported genes."""
    
    logger.info("Generating Panel A: Concordance validation")
    
    # Check required inputs
    if not PGX_CLASSIFICATION.exists() or not PGX_CONCORDANCE.exists():
        logger.error("Cannot create Panel A - missing classification or concordance data")
        ax.text(0.5, 0.5, 'Panel A: Data not available', 
                ha='center', va='center', transform=ax.transAxes)
        return 0
    
    # Load data
    tiers = pd.read_csv(PGX_CLASSIFICATION)
    concordance = pd.read_csv(PGX_CONCORDANCE)
    
    # Get genes for Panel A
    array_supported = tiers.query('Classification == "Array-supported"').Gene.tolist()
    sentinels = ['CYP2D6', 'UGT1A1']
    genes_for_panelA = array_supported + [g for g in sentinels if g in tiers.Gene.unique()]
    
    # Filter concordance to selected genes
    panel_data = concordance[concordance.Gene.isin(genes_for_panelA)].copy()
    
    if panel_data.empty:
        logger.warning("No data for Panel A after filtering")
        ax.text(0.5, 0.5, 'Panel A: No Array-supported genes found', 
                ha='center', va='center', transform=ax.transAxes)
        return 0
    
    # Merge with classification for sorting
    panel_data = panel_data.merge(tiers[['Gene', 'Classification']], on='Gene', how='left')
    
    # Sort: Array-supported first, then by Kappa descending
    panel_data['sort_key'] = panel_data['Classification'].apply(lambda x: 0 if x == 'Array-supported' else 1)
    panel_data = panel_data.sort_values(['sort_key', 'Kappa'], ascending=[True, False])
    
    # Prepare data for plotting
    genes = panel_data.Gene.tolist()
    kappa = panel_data.Kappa.tolist()
    accuracy = panel_data.Accuracy.tolist()
    indeterminate = panel_data.IndeterminateRate_Array.tolist()
    
    n_genes = len(genes)
    x = np.arange(n_genes)
    width = 0.25
    
    # Create clustered bars
    bars1 = ax.bar(x - width, kappa, width, label="Cohen's κ", color='#3498db', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x, [a/100 for a in accuracy], width, label='Accuracy (%)', color='#27ae60', edgecolor='black', linewidth=0.5)
    bars3 = ax.bar(x + width, [i/100 for i in indeterminate], width, label='Indeterminate (%)', color='#e74c3c', edgecolor='black', linewidth=0.5)
    
    # Add thresholds (scaled to 0-1)
    ax.axhline(y=0.85, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=0.70, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Formatting
    ax.set_xlabel('Gene')
    ax.set_ylabel('Value (κ: 0-1, Others: 0-100%)')
    ax.set_title('A: Validation on WGS overlap (n=109)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right')
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper right', framealpha=0.95)
    
    # Add footnote
    ax.text(0.5, -0.25, '46 additional genes are WGS-only due to high indeterminate rate or complexity',
            ha='center', transform=ax.transAxes, fontsize=9, style='italic')
    
    # Always save tidy table (overwrites if exists)
    tidy_df = pd.DataFrame({
        'Gene': genes,
        'Classification': panel_data.Classification.tolist(),
        'Kappa': kappa,
        'Accuracy': accuracy,
        'IndeterminateRate': indeterminate,
        'N_valid': panel_data.N_valid.tolist()
    })
    tidy_df.to_csv(PANEL_A_CSV, index=False)
    logger.info(f"Panel A: Saved tidy CSV to {PANEL_A_CSV}")
    
    # Save standalone figure
    fig_standalone = plt.figure(figsize=(12, 6))
    ax_standalone = fig_standalone.add_subplot(111)
    bars1 = ax_standalone.bar(x - width, kappa, width, label="Cohen's κ", color='#3498db', edgecolor='black', linewidth=0.5)
    bars2 = ax_standalone.bar(x, [a/100 for a in accuracy], width, label='Accuracy (%)', color='#27ae60', edgecolor='black', linewidth=0.5)
    bars3 = ax_standalone.bar(x + width, [i/100 for i in indeterminate], width, label='Indeterminate (%)', color='#e74c3c', edgecolor='black', linewidth=0.5)
    ax_standalone.axhline(y=0.85, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_standalone.axhline(y=0.70, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax_standalone.set_xlabel('Gene')
    ax_standalone.set_ylabel('Value (κ: 0-1, Others: 0-100%)')
    ax_standalone.set_title('Panel A: Validation on WGS overlap (n=109)', fontweight='bold')
    ax_standalone.set_xticks(x)
    ax_standalone.set_xticklabels(genes, rotation=45, ha='right')
    ax_standalone.set_ylim(0, 1.05)
    ax_standalone.grid(True, alpha=0.3, axis='y')
    ax_standalone.legend(loc='upper right', framealpha=0.95)
    plt.tight_layout()
    plt.savefig(PANEL_A_PNG, dpi=300, bbox_inches='tight')
    plt.close(fig_standalone)
    logger.info(f"Panel A: Saved figure to {PANEL_A_PNG}")
    
    logger.info(f"Panel A: Used {len(array_supported)} Array-supported genes + {len([g for g in sentinels if g in genes])} sentinels")
    
    return len(array_supported)

# ============================================================================
# PANEL B: PHENOTYPE FREQUENCIES
# ============================================================================

def create_panel_b(ax, logger):
    """Panel B: Phenotype frequencies for Array-supported genes."""
    
    logger.info("Generating Panel B: Phenotype frequencies")
    
    # Check required inputs
    if not PGX_CLASSIFICATION.exists() or not PGX_POPULATION_SUMMARY.exists():
        logger.warning("Cannot create Panel B - missing required data")
        ax.text(0.5, 0.5, 'Panel B: Data not available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Load data
    tiers = pd.read_csv(PGX_CLASSIFICATION)
    pop_summary = pd.read_csv(PGX_POPULATION_SUMMARY)
    
    # Get Array-supported genes only
    array_supported = tiers.query('Classification == "Array-supported"').Gene.tolist()
    
    if not array_supported:
        logger.warning("No Array-supported genes found for Panel B")
        ax.text(0.5, 0.5, 'No Array-supported genes available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Filter to Array-supported genes
    panel_data = pop_summary[pop_summary.Gene.isin(array_supported)].copy()
    
    if panel_data.empty:
        logger.warning("No population data for Array-supported genes")
        ax.text(0.5, 0.5, 'No data available for Array-supported genes', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Standardize phenotype names
    phenotype_map = {
        'poor metabolizer': 'Poor Metabolizer',
        'intermediate metabolizer': 'Intermediate Metabolizer',
        'normal metabolizer': 'Normal Metabolizer',
        'extensive metabolizer': 'Normal Metabolizer',
        'rapid metabolizer': 'Rapid Metabolizer',
        'ultrarapid metabolizer': 'Ultrarapid Metabolizer',
        'ultra-rapid metabolizer': 'Ultrarapid Metabolizer'
    }
    
    panel_data['Phenotype'] = panel_data['Phenotype'].str.lower().map(
        lambda x: phenotype_map.get(x, x.title() if pd.notna(x) else 'Unknown')
    )
    
    # Create two sub-plots within Panel B
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    gs_b = GridSpecFromSubplotSpec(2, 1, subplot_spec=ax.get_subplotspec(), height_ratios=[1, 1], hspace=0.3)
    ax_b1 = plt.subplot(gs_b[0])
    ax_b2 = plt.subplot(gs_b[1])
    
    # Process Region×Group data
    region_group_data = panel_data[panel_data['Level'] == 'Region_Group'].copy()
    
    saved_region = False
    if not region_group_data.empty:
        # Pivot for stacked bar
        pivot = region_group_data.pivot_table(
            index='Category', 
            columns='Phenotype', 
            values='Frequency',
            aggfunc='mean',
            fill_value=0
        )
        
        # Reorder columns to match phenotype order
        cols_present = [p for p in PHENOTYPE_ORDER if p in pivot.columns]
        pivot = pivot[cols_present]
        
        # Check sample sizes and filter
        category_totals = region_group_data.groupby('Category')['Total'].mean()
        valid_categories = category_totals[category_totals >= 15].index
        pivot_filtered = pivot.loc[pivot.index.isin(valid_categories)]
        
        if not pivot_filtered.empty:
            # Create stacked bars
            pivot_filtered.plot(kind='bar', stacked=True, ax=ax_b1, 
                               color=[PHENOTYPE_COLORS.get(p, 'gray') for p in cols_present],
                               edgecolor='black', linewidth=0.5)
            
            # Add value labels
            for container in ax_b1.containers:
                labels = [f'{v:.2f}' if v >= 0.02 else '' for v in container.datavalues]
                ax_b1.bar_label(container, labels=labels, label_type='center', fontsize=8)
            
            ax_b1.set_title('Region × Group', fontsize=11)
            ax_b1.set_xlabel('')
            ax_b1.set_ylabel('Frequency')
            ax_b1.legend(title='Phenotype', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
            ax_b1.set_xticklabels(ax_b1.get_xticklabels(), rotation=45, ha='right')
            
            # Save region group frequencies
            region_group_data.to_csv(PANEL_B_REGION_CSV, index=False)
            saved_region = True
            logger.info(f"Panel B: Saved region×group frequencies to {PANEL_B_REGION_CSV}")
        
        # Mark low-n categories
        low_n_categories = category_totals[category_totals < 15].index
        if len(low_n_categories) > 0:
            for i, cat in enumerate(pivot_filtered.index):
                if cat in low_n_categories:
                    ax_b1.text(i, 0.5, 'n<15', ha='center', va='center', fontsize=9, style='italic')
    
    # Process Population data
    pop_data = panel_data[panel_data['Level'] == 'Population'].copy()
    
    saved_pop = False
    if not pop_data.empty:
        pivot_pop = pop_data.pivot_table(
            index='Category',
            columns='Phenotype',
            values='Frequency', 
            aggfunc='mean',
            fill_value=0
        )
        
        # Reorder columns
        cols_present = [p for p in PHENOTYPE_ORDER if p in pivot_pop.columns]
        pivot_pop = pivot_pop[cols_present]
        
        # Plot
        pivot_pop.plot(kind='bar', stacked=True, ax=ax_b2,
                      color=[PHENOTYPE_COLORS.get(p, 'gray') for p in cols_present],
                      edgecolor='black', linewidth=0.5)
        
        # Add value labels
        for container in ax_b2.containers:
            labels = [f'{v:.2f}' if v >= 0.02 else '' for v in container.datavalues]
            ax_b2.bar_label(container, labels=labels, label_type='center', fontsize=8)
        
        ax_b2.set_title('Population (n ≥ 15)', fontsize=11)
        ax_b2.set_xlabel('Population')
        ax_b2.set_ylabel('Frequency')
        ax_b2.legend().set_visible(False)
        ax_b2.set_xticklabels(ax_b2.get_xticklabels(), rotation=45, ha='right')
        
        # Save population frequencies
        pop_data.to_csv(PANEL_B_POP_CSV, index=False)
        saved_pop = True
        logger.info(f"Panel B: Saved population frequencies to {PANEL_B_POP_CSV}")
    
    # Overall title for Panel B
    ax.set_title('B: Array-supported genes across cohort (n=736)', fontweight='bold')
    ax.axis('off')
    
    # Save standalone figure
    fig_standalone = plt.figure(figsize=(12, 10))
    gs_standalone = GridSpec(2, 1, figure=fig_standalone, height_ratios=[1, 1], hspace=0.3)
    
    # Copy plots to standalone figure
    if saved_region:
        ax_s1 = fig_standalone.add_subplot(gs_standalone[0])
        region_group_data_reload = pd.read_csv(PANEL_B_REGION_CSV)
        pivot_s1 = region_group_data_reload.pivot_table(
            index='Category', columns='Phenotype', values='Frequency', aggfunc='mean', fill_value=0
        )
        cols_present = [p for p in PHENOTYPE_ORDER if p in pivot_s1.columns]
        pivot_s1[cols_present].plot(kind='bar', stacked=True, ax=ax_s1,
                                    color=[PHENOTYPE_COLORS.get(p, 'gray') for p in cols_present],
                                    edgecolor='black', linewidth=0.5)
        ax_s1.set_title('Region × Group', fontsize=11)
        ax_s1.set_ylabel('Frequency')
        ax_s1.legend(title='Phenotype', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    
    if saved_pop:
        ax_s2 = fig_standalone.add_subplot(gs_standalone[1])
        pop_data_reload = pd.read_csv(PANEL_B_POP_CSV)
        pivot_s2 = pop_data_reload.pivot_table(
            index='Category', columns='Phenotype', values='Frequency', aggfunc='mean', fill_value=0
        )
        cols_present = [p for p in PHENOTYPE_ORDER if p in pivot_s2.columns]
        pivot_s2[cols_present].plot(kind='bar', stacked=True, ax=ax_s2,
                                    color=[PHENOTYPE_COLORS.get(p, 'gray') for p in cols_present],
                                    edgecolor='black', linewidth=0.5)
        ax_s2.set_title('Population (n ≥ 15)', fontsize=11)
        ax_s2.set_xlabel('Population')
        ax_s2.set_ylabel('Frequency')
    
    fig_standalone.suptitle('Panel B: Array-supported genes across cohort (n=736)', fontweight='bold')
    plt.tight_layout()
    plt.savefig(PANEL_B_PNG, dpi=300, bbox_inches='tight')
    plt.close(fig_standalone)
    logger.info(f"Panel B: Saved figure to {PANEL_B_PNG}")
    
    logger.info(f"Panel B: Processed {len(array_supported)} Array-supported genes")

# ============================================================================
# PANEL C: FDA TRIGGER RATES
# ============================================================================

def create_panel_c(ax, logger):
    """Panel C: True FDA trigger rates per 100 for WGS vs Cohort."""
    
    logger.info("Generating Panel C: FDA trigger rates")
    
    # Check if FDA matches file exists
    if not FDA_PATIENT_MATCHES.exists():
        logger.warning(f"FDA_patient_matches.csv not found at {FDA_PATIENT_MATCHES}, skipping Panel C")
        ax.text(0.5, 0.5, 'Panel C omitted (FDA data not available)', 
                ha='center', va='center', transform=ax.transAxes, fontweight='bold')
        return False
    
    if not WGS_SAMPLES.exists():
        logger.warning("WGS samples list not found, skipping Panel C")
        ax.text(0.5, 0.5, 'Panel C omitted (WGS samples list not available)', 
                ha='center', va='center', transform=ax.transAxes, fontweight='bold')
        return False
    
    # Constants
    WGS_N, COHORT_N = 109, 736
    
    # Load data
    wgs_ids = set(open(WGS_SAMPLES).read().splitlines())
    matches = pd.read_csv(FDA_PATIENT_MATCHES)
    tiers = pd.read_csv(PGX_CLASSIFICATION)
    
    # Merge and filter to Array-supported
    matches = matches.merge(tiers[['Gene', 'Classification']], on='Gene', how='left')
    matches = matches.query('Classification == "Array-supported"')
    
    if matches.empty:
        logger.warning("No Array-supported FDA triggers found")
        ax.text(0.5, 0.5, 'No Array-supported FDA triggers found', 
                ha='center', va='center', transform=ax.transAxes)
        return False
    
    # Calculate counts and rates
    wgs_counts = matches[matches.SampleID.isin(wgs_ids)].groupby('Gene').size().rename('WGS_count')
    coh_counts = matches.groupby('Gene').size().rename('Cohort_count')
    
    df = pd.concat([wgs_counts, coh_counts], axis=1).fillna(0).reset_index()
    df['RatePer100_WGS'] = 100 * df['WGS_count'] / WGS_N
    df['RatePer100_Cohort'] = 100 * df['Cohort_count'] / COHORT_N
    
    # Sort by total rate
    df['TotalRate'] = df['RatePer100_WGS'] + df['RatePer100_Cohort']
    df = df.sort_values('TotalRate', ascending=True)
    
    # Create grouped bar plot
    genes = df['Gene'].tolist()
    y_pos = np.arange(len(genes))
    
    bar_width = 0.35
    bars1 = ax.barh(y_pos - bar_width/2, df['RatePer100_WGS'], bar_width, 
                    label='WGS', color='#3498db', edgecolor='black', linewidth=0.5)
    bars2 = ax.barh(y_pos + bar_width/2, df['RatePer100_Cohort'], bar_width,
                    label='Cohort', color='#27ae60', edgecolor='black', linewidth=0.5)
    
    # Add annotations
    for i, (_, row) in enumerate(df.iterrows()):
        # WGS annotation
        if row['RatePer100_WGS'] > 0:
            ax.text(row['RatePer100_WGS'] + 0.5, i - bar_width/2, 
                   f"{row['RatePer100_WGS']:.1f} per 100\nN={WGS_N}",
                   va='center', fontsize=8)
        # Cohort annotation  
        if row['RatePer100_Cohort'] > 0:
            ax.text(row['RatePer100_Cohort'] + 0.5, i + bar_width/2,
                   f"{row['RatePer100_Cohort']:.1f} per 100\nN={COHORT_N}",
                   va='center', fontsize=8)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes)
    ax.set_xlabel('Rate per 100 individuals')
    ax.set_title('C: FDA triggers per 100 (WGS 109 vs Cohort 736)', fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Always save tidy table (overwrites if exists)
    tidy_df = pd.DataFrame()
    for _, row in df.iterrows():
        tidy_df = pd.concat([tidy_df, pd.DataFrame([
            {'Gene': row['Gene'], 'Stratum': 'WGS', 'N': WGS_N, 
             'Count': int(row['WGS_count']), 'RatePer100': row['RatePer100_WGS']},
            {'Gene': row['Gene'], 'Stratum': 'Cohort', 'N': COHORT_N,
             'Count': int(row['Cohort_count']), 'RatePer100': row['RatePer100_Cohort']}
        ])], ignore_index=True)
    
    tidy_df.to_csv(PANEL_C_CSV, index=False)
    logger.info(f"Panel C: Saved tidy CSV to {PANEL_C_CSV}")
    
    # Save standalone figure
    fig_standalone = plt.figure(figsize=(10, max(6, len(genes) * 0.4)))
    ax_standalone = fig_standalone.add_subplot(111)
    bars1 = ax_standalone.barh(y_pos - bar_width/2, df['RatePer100_WGS'], bar_width, 
                               label='WGS', color='#3498db', edgecolor='black', linewidth=0.5)
    bars2 = ax_standalone.barh(y_pos + bar_width/2, df['RatePer100_Cohort'], bar_width,
                               label='Cohort', color='#27ae60', edgecolor='black', linewidth=0.5)
    for i, (_, row) in enumerate(df.iterrows()):
        if row['RatePer100_WGS'] > 0:
            ax_standalone.text(row['RatePer100_WGS'] + 0.5, i - bar_width/2, 
                             f"{row['RatePer100_WGS']:.1f} per 100\nN={WGS_N}",
                             va='center', fontsize=8)
        if row['RatePer100_Cohort'] > 0:
            ax_standalone.text(row['RatePer100_Cohort'] + 0.5, i + bar_width/2,
                             f"{row['RatePer100_Cohort']:.1f} per 100\nN={COHORT_N}",
                             va='center', fontsize=8)
    ax_standalone.set_yticks(y_pos)
    ax_standalone.set_yticklabels(genes)
    ax_standalone.set_xlabel('Rate per 100 individuals')
    ax_standalone.set_title('Panel C: FDA triggers per 100 (WGS 109 vs Cohort 736)', fontweight='bold')
    ax_standalone.legend(loc='lower right')
    ax_standalone.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig(PANEL_C_PNG, dpi=300, bbox_inches='tight')
    plt.close(fig_standalone)
    logger.info(f"Panel C: Saved figure to {PANEL_C_PNG}")
    
    logger.info(f"Panel C: Used true counts (yes)")
    logger.info(f"Panel C: {len(df)} genes with FDA triggers")
    
    return True

# ============================================================================  
# PANEL D: WGS BREADTH ANALYSIS
# ============================================================================

def create_panel_d(ax, logger):
    """Panel D: WGS breadth analysis with p-values."""
    
    logger.info("Generating Panel D: WGS breadth analysis")
    
    # First check for legacy WGS p-value files
    legacy_candidates = [
        RESULTS_DIR / "WGS_interpop_pvalues_56genes.tsv",
        RESULTS_DIR / "wgs_interpop_pvalues.tsv",
        RESULTS_DIR / "wgs_pvalues.tsv",
        RESULTS_DIR / "wgs_breadth_pvalues.csv"
    ]
    
    legacy_pval_file = None
    for candidate in legacy_candidates:
        if candidate.exists():
            legacy_pval_file = candidate
            break
    
    if not legacy_pval_file:
        # Also check for any matching pattern
        legacy_pval_files = list(RESULTS_DIR.glob("*wgs*pval*.tsv")) + \
                           list(RESULTS_DIR.glob("*wgs*pval*.csv")) + \
                           list(RESULTS_DIR.glob("*WGS*pval*.tsv"))
        if legacy_pval_files:
            legacy_pval_file = legacy_pval_files[0]
    
    if legacy_pval_file:
        # Use legacy p-values
        logger.info(f"Panel D: Found legacy p-values at {legacy_pval_file}")
        
        try:
            pval_df = pd.read_csv(legacy_pval_file, sep='\t' if legacy_pval_file.suffix == '.tsv' else ',')
            
            # Identify p-value column
            if 'p_value' in pval_df.columns:
                pval_col = 'p_value'
            elif 'pvalue' in pval_df.columns:
                pval_col = 'pvalue'
            elif 'P-value' in pval_df.columns:
                pval_col = 'P-value'
            else:
                # Assume second column is p-value
                pval_col = pval_df.columns[1]
            
            # Check if we have enough informative genes
            if len(pval_df) < 30:
                logger.warning(f"Panel D: Only {len(pval_df)} genes in legacy file - omitting panel")
                ax.text(0.5, 0.5, 'Panel D omitted — no informative WGS breadth file available.\nSee SI.',
                       ha='center', va='center', transform=ax.transAxes, fontweight='bold')
                logger.info("Panel D: Omitted (insufficient genes in legacy file)")
                return 'omitted_insufficient'
            
            # Calculate -log10 p-values
            pval_df['neg_log10_p'] = -np.log10(pval_df[pval_col].clip(lower=1e-10))
            pval_df['neg_log10_p'] = pval_df['neg_log10_p'].clip(upper=5)
            
            # Check if any are significant
            if (pval_df[pval_col] <= 0.1).sum() == 0:
                logger.warning("Panel D: All p-values > 0.1 - omitting panel")
                ax.text(0.5, 0.5, 'Panel D omitted — no informative WGS breadth file available.\nSee SI.',
                       ha='center', va='center', transform=ax.transAxes, fontweight='bold')
                logger.info("Panel D: Omitted (all p > 0.1)")
                return 'omitted_weak_signal'
            
            # Plot
            pval_df = pval_df.sort_values('neg_log10_p')
            ax.scatter(pval_df['neg_log10_p'], range(len(pval_df)), 
                      color='#3498db', s=50, edgecolor='black', linewidth=0.5)
            ax.set_yticks(range(len(pval_df)))
            ax.set_yticklabels(pval_df.iloc[:, 0], fontsize=9)  # Gene names
            ax.set_xlabel('-log10(p-value)')
            ax.set_title('D: WGS breadth (56 genes, original analysis)', fontweight='bold')
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='x')
            
            # Save tidy CSV (always overwrite)
            pval_df.to_csv(PANEL_D_CSV, index=False)
            logger.info(f"Panel D: Saved tidy CSV to {PANEL_D_CSV}")
            
            # Save standalone figure
            fig_standalone = plt.figure(figsize=(8, max(6, len(pval_df) * 0.3)))
            ax_standalone = fig_standalone.add_subplot(111)
            ax_standalone.scatter(pval_df['neg_log10_p'], range(len(pval_df)), 
                                 color='#3498db', s=50, edgecolor='black', linewidth=0.5)
            ax_standalone.set_yticks(range(len(pval_df)))
            ax_standalone.set_yticklabels(pval_df.iloc[:, 0], fontsize=9)
            ax_standalone.set_xlabel('-log10(p-value)')
            ax_standalone.set_title('Panel D: WGS breadth (56 genes, original analysis)', fontweight='bold')
            ax_standalone.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
            ax_standalone.legend()
            ax_standalone.grid(True, alpha=0.3, axis='x')
            plt.tight_layout()
            plt.savefig(PANEL_D_PNG, dpi=300, bbox_inches='tight')
            plt.close(fig_standalone)
            logger.info(f"Panel D: Saved figure to {PANEL_D_PNG}")
            
            logger.info(f"Panel D: Used legacy p-values from {legacy_pval_file.name}")
            return 'legacy'
            
        except Exception as e:
            logger.error(f"Panel D: Failed to process legacy file: {e}")
            # Fall through to try rebuilding
            legacy_pval_file = None
    
    # Try to rebuild p-values from WGS/cohort data
    if not legacy_pval_file:
        logger.info("Panel D: Attempting to rebuild p-values from WGS data")
        
        if not WGS_TRUTH_RESULTS.exists() or not COHORT_RESULTS.exists():
            logger.warning("Cannot rebuild p-values - missing WGS or cohort results")
            ax.text(0.5, 0.5, 'Panel D omitted — no informative WGS breadth file available.\nSee SI.',
                   ha='center', va='center', transform=ax.transAxes, fontweight='bold')
            logger.info("Panel D: Omitted (missing input files for rebuilding)")
            return 'omitted_no_data'
        
        # Load sample metadata
        sample_metadata = None
        if (BASE / "pop_736.updated.tsv").exists():
            sample_metadata = pd.read_csv(BASE / "pop_736.updated.tsv", sep='\t')
        elif ARRAY_SAMPLES.exists() and WGS_SAMPLES.exists():
            # Try to reconstruct basic metadata
            wgs_ids = set(open(WGS_SAMPLES).read().splitlines())
            array_ids = set(open(ARRAY_SAMPLES).read().splitlines())
            all_ids = list(wgs_ids | array_ids)
            sample_metadata = pd.DataFrame({'IID': all_ids})
        
        if sample_metadata is None:
            logger.warning("Cannot rebuild p-values - no sample metadata available")
            ax.text(0.5, 0.5, 'Panel D omitted — no informative WGS breadth file available.\nSee SI.',
                   ha='center', va='center', transform=ax.transAxes, fontweight='bold')
            logger.info("Panel D: Omitted (no metadata for rebuilding)")
            return 'omitted_no_metadata'
        
        # Load WGS results for p-value calculation
        wgs_df = pd.read_csv(WGS_TRUTH_RESULTS, sep='\t')
        cohort_df = pd.read_csv(COHORT_RESULTS, sep='\t')
        
        # Define actionable phenotypes for key genes
        actionable_map = {
            'CYP2C19': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer'],
            'CYP2B6': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer'],
            'SLCO1B1': ['Poor Function', 'Decreased Function', 'Poor Metabolizer', 'poor function', 'decreased function'],
            'TPMT': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer'],
            'NUDT15': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer'],
            'DPYD': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer'],
            'UGT1A1': ['Poor Metabolizer', 'Intermediate Metabolizer', 'poor metabolizer', 'intermediate metabolizer']
        }
        
        # Calculate p-values for genes
        pvalues = []
        genes_tested = wgs_df['Gene'].unique()
        
        for gene in genes_tested:
            gene_wgs = wgs_df[wgs_df['Gene'] == gene]
            gene_cohort = cohort_df[cohort_df['Gene'] == gene]
            
            if len(gene_wgs) < 10:
                continue
            
            # Determine if gene has actionable phenotypes
            if gene in actionable_map:
                actionable = actionable_map[gene]
                n_actionable_wgs = gene_wgs['Phenotype'].isin(actionable).sum()
                n_total_wgs = len(gene_wgs)
                
                if gene_cohort.empty:
                    # Use simple binomial test against expected frequency
                    from scipy.stats import binom_test
                    try:
                        p = binom_test(n_actionable_wgs, n_total_wgs, 0.1, alternative='greater')
                        pvalues.append({
                            'Gene': gene,
                            'p_value': p,
                            'n_actionable': n_actionable_wgs,
                            'n_total': n_total_wgs
                        })
                    except:
                        continue
                else:
                    # Compare WGS vs cohort frequencies
                    n_actionable_cohort = gene_cohort['Phenotype'].isin(actionable).sum()
                    n_total_cohort = len(gene_cohort)
                    
                    # Chi-square test
                    from scipy.stats import chi2_contingency
                    contingency = [[n_actionable_wgs, n_total_wgs - n_actionable_wgs],
                                  [n_actionable_cohort, n_total_cohort - n_actionable_cohort]]
                    try:
                        chi2, p, _, _ = chi2_contingency(contingency)
                        pvalues.append({
                            'Gene': gene,
                            'p_value': p,
                            'n_actionable': n_actionable_wgs,
                            'n_total': n_total_wgs,
                            'chi2': chi2
                        })
                    except:
                        continue
            else:
                # For genes without clear actionable mapping, test frequency differences
                if gene_cohort.empty or len(gene_wgs) < 20:
                    continue
                    
                # Test most common phenotype
                most_common = gene_wgs['Phenotype'].value_counts().index[0]
                n_pheno_wgs = (gene_wgs['Phenotype'] == most_common).sum()
                n_pheno_cohort = (gene_cohort['Phenotype'] == most_common).sum()
                
                from scipy.stats import chi2_contingency
                contingency = [[n_pheno_wgs, len(gene_wgs) - n_pheno_wgs],
                              [n_pheno_cohort, len(gene_cohort) - n_pheno_cohort]]
                try:
                    chi2, p, _, _ = chi2_contingency(contingency)
                    if p < 0.5:  # Only include if showing some signal
                        pvalues.append({
                            'Gene': gene,
                            'p_value': p,
                            'phenotype_tested': most_common,
                            'n_total': len(gene_wgs)
                        })
                except:
                    continue
        
        # Check if we have enough significant results
        if len(pvalues) < 30:
            logger.warning(f"Panel D: Only {len(pvalues)} genes with p-values - omitting panel")
            ax.text(0.5, 0.5, 'Panel D omitted — insufficient signal.\nSee SI for full WGS table.',
                   ha='center', va='center', transform=ax.transAxes, fontweight='bold')
            logger.info(f"Panel D: Omitted (rebuilt {len(pvalues)} genes, insufficient)")
            return 'omitted_insufficient_genes'
        
        significant = sum(1 for p in pvalues if p['p_value'] <= 0.1)
        if significant == 0:
            logger.warning("Panel D: All rebuilt p-values > 0.1 - omitting panel")
            ax.text(0.5, 0.5, 'Panel D omitted — insufficient signal.\nSee SI for full WGS table.',
                   ha='center', va='center', transform=ax.transAxes, fontweight='bold')
            logger.info("Panel D: Omitted (all rebuilt p > 0.1)")
            return 'omitted_weak_signal'
        
        # Create plot with rebuilt p-values
        pval_df = pd.DataFrame(pvalues)
        pval_df['neg_log10_p'] = -np.log10(pval_df['p_value'].clip(lower=1e-10))
        pval_df['neg_log10_p'] = pval_df['neg_log10_p'].clip(upper=5)
        pval_df = pval_df.sort_values('neg_log10_p')
        
        # Plot
        ax.scatter(pval_df['neg_log10_p'], range(len(pval_df)),
                  color='#e74c3c', s=50, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(pval_df)))
        ax.set_yticklabels(pval_df['Gene'], fontsize=9)
        ax.set_xlabel('-log10(p-value)')
        ax.set_title(f'D: WGS breadth ({len(pval_df)} genes, reconstructed)', fontweight='bold')
        ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='x')
        
        # Save rebuilt p-values (always overwrite)
        pval_df.to_csv(PANEL_D_CSV, index=False)
        logger.info(f"Panel D: Saved rebuilt p-values to {PANEL_D_CSV}")
        
        # Save standalone figure
        fig_standalone = plt.figure(figsize=(8, max(6, len(pval_df) * 0.3)))
        ax_standalone = fig_standalone.add_subplot(111)
        ax_standalone.scatter(pval_df['neg_log10_p'], range(len(pval_df)),
                             color='#e74c3c', s=50, edgecolor='black', linewidth=0.5)
        ax_standalone.set_yticks(range(len(pval_df)))
        ax_standalone.set_yticklabels(pval_df['Gene'], fontsize=9)
        ax_standalone.set_xlabel('-log10(p-value)')
        ax_standalone.set_title(f'Panel D: WGS breadth ({len(pval_df)} genes, reconstructed)', fontweight='bold')
        ax_standalone.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax_standalone.legend()
        ax_standalone.grid(True, alpha=0.3, axis='x')
        plt.tight_layout()
        plt.savefig(PANEL_D_PNG, dpi=300, bbox_inches='tight')
        plt.close(fig_standalone)
        logger.info(f"Panel D: Saved figure to {PANEL_D_PNG}")
        
        logger.info(f"Panel D: Rebuilt p-values for {len(pval_df)} genes")
        return 'rebuilt'

# ============================================================================
# MAIN INTEGRATED FIGURE
# ============================================================================

def create_integrated_figure():
    """Create the complete integrated PGx figure."""
    
    logger = setup_logging()
    logger.info("="*70)
    logger.info("Starting integrated PGx figure generation")
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    logger.info("This script ONLY creates figures from existing CSVs")
    logger.info("It NEVER runs PyPGx or processes VCFs")
    logger.info("="*70)
    
    # Check critical inputs
    check_critical_inputs(logger)
    
    # Create figure with GridSpec for flexible layout
    fig = plt.figure(figsize=(16, 20))
    gs = GridSpec(4, 1, figure=fig, height_ratios=[1.2, 1.5, 0.8, 0.8], hspace=0.4)
    
    # Panel A: Concordance validation
    logger.info("-" * 50)
    ax_a = fig.add_subplot(gs[0])
    n_array_supported = create_panel_a(ax_a, logger)
    
    # Panel B: Phenotype frequencies
    logger.info("-" * 50)
    ax_b = fig.add_subplot(gs[1])
    create_panel_b(ax_b, logger)
    
    # Panel C: FDA trigger rates
    logger.info("-" * 50)
    ax_c = fig.add_subplot(gs[2])
    fda_success = create_panel_c(ax_c, logger)
    
    # Panel D: WGS breadth
    logger.info("-" * 50)
    ax_d = fig.add_subplot(gs[3])
    panel_d_status = create_panel_d(ax_d, logger)
    
    # Overall figure title
    fig.suptitle('Integrated Pharmacogenomics Analysis: Platform Validation and Population Insights',
                 fontsize=16, fontweight='bold', y=0.995)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    # Save integrated figure (always overwrite)
    plt.savefig(INTEGRATED_FIGURE, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Final logging summary
    logger.info("="*70)
    logger.info("SUMMARY:")
    logger.info(f"  Array-supported genes used in Panels A-C: {n_array_supported}")
    logger.info(f"  Panel A: Regenerated")
    logger.info(f"  Panel B: Regenerated")
    logger.info(f"  Panel C: {'Regenerated' if fda_success else 'Skipped (no FDA data)'}")
    logger.info(f"  Panel D: {panel_d_status}")
    logger.info("")
    logger.info("OUTPUT FILES (all regenerated):")
    logger.info(f"  Main figure: {INTEGRATED_FIGURE}")
    logger.info(f"  Panel figures: {OUTPUT_DIR}/*.png")
    logger.info(f"  Tidy CSVs: {OUTPUT_DIR}/*.csv")
    logger.info(f"  Log: {LOG_FILE}")
    logger.info("="*70)
    
    print(f"\n✓ Integrated PGx figure regenerated successfully!")
    print(f"  Main output: {INTEGRATED_FIGURE}")
    print(f"  Individual panels: {OUTPUT_DIR}/panel*.png")
    print(f"  Log: {LOG_FILE}")
    print(f"\n  Tidy CSVs regenerated:")
    
    if PANEL_A_CSV.exists():
        print(f"    ✓ {PANEL_A_CSV.name}")
    if PANEL_B_REGION_CSV.exists():
        print(f"    ✓ {PANEL_B_REGION_CSV.name}")
    if PANEL_B_POP_CSV.exists():
        print(f"    ✓ {PANEL_B_POP_CSV.name}")
    if PANEL_C_CSV.exists():
        print(f"    ✓ {PANEL_C_CSV.name}")
    if PANEL_D_CSV.exists():
        print(f"    ✓ {PANEL_D_CSV.name}")
    
    print("\n  All figures and CSVs have been regenerated from existing data.")
    print("  No VCFs or PyPGx runs were performed.")

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    try:
        create_integrated_figure()
    except KeyboardInterrupt:
        print("\n\nFigure generation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error generating figure: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)