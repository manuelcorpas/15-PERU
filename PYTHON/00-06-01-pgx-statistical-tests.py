#!/usr/bin/env python3
"""
00-06-01-pgx-statistical-tests.py

Statistical comparison of pharmacogenomic phenotype frequencies between
Indigenous and Mestizo populations within each region.

Purpose:
    Address Reviewer 3's comment: "Compare metabolizer frequencies between
    Peruvian vs global, Indigenous vs mestizo, Andean vs Amazonian."

Usage:
    python 00-06-01-pgx-statistical-tests.py

Input:
    ANALYSIS/00-05-PyPGx/results/PGx_population_summary.csv
    ANALYSIS/00-05-PyPGx/results/cohort_pypgx_results.tsv

Output:
    ANALYSIS/00-12-PGX-STATS/indigenous_vs_mestizo_tests.csv
    ANALYSIS/00-12-PGX-STATS/region_comparison_tests.csv
    ANALYSIS/00-12-PGX-STATS/manuscript_text.txt
    ANALYSIS/00-12-PGX-STATS/supplementary_table.csv

Author: Generated for Nature Health genomics paper revision
Date: 2025
"""

import os
import sys
from pathlib import Path
from typing import Dict, List
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("ANALYSIS/00-05-PyPGx")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = Path("ANALYSIS/00-12-PGX-STATS")

# Input files - primary is population summary (already aggregated)
PGX_SUMMARY = RESULTS_DIR / "PGx_population_summary.csv"

# Population classification
POPULATION_CLASSIFICATION = {
    # Amazonian Indigenous
    'MATZES': ('Amazon', 'Indigenous'),
    'AWAJUN': ('Amazon', 'Indigenous'),
    'CANDOSHI': ('Amazon', 'Indigenous'),
    'SHIPIBO_INS': ('Amazon', 'Indigenous'),
    'SHIPIBO': ('Amazon', 'Indigenous'),
    'ASHANINKA_INS': ('Amazon', 'Indigenous'),
    'ASHANINKA': ('Amazon', 'Indigenous'),
    'MATSIGUENKAS': ('Amazon', 'Indigenous'),
    'NAHUA': ('Amazon', 'Indigenous'),
    
    # Amazonian Mestizo
    'IQUITOS': ('Amazon', 'Mestizo'),
    'LAMAS': ('Amazon', 'Mestizo'),
    
    # Andean Indigenous
    'CHOPCCAS': ('Andes', 'Indigenous'),
    'UROS': ('Andes', 'Indigenous'),
    'QEROS': ('Andes', 'Indigenous'),
    'JAQARUS': ('Andes', 'Indigenous'),
    'JACARUS': ('Andes', 'Indigenous'),
    
    # Andean Mestizo
    'CUSCO': ('Andes', 'Mestizo'),
    'AYACUCHO': ('Andes', 'Mestizo'),
    'PUNO': ('Andes', 'Mestizo'),
    'ANCASH': ('Andes', 'Mestizo'),
    'HUARAZ': ('Andes', 'Mestizo'),
    'AREQUIPA': ('Andes', 'Mestizo'),
    'MOQUEGUA': ('Andes', 'Mestizo'),
    'TACNA': ('Andes', 'Mestizo'),
    'CHACHAPOYAS': ('Andes', 'Mestizo'),
    
    # Coastal Indigenous
    'TALLANES': ('Coast', 'Indigenous'),
    'TALLAN': ('Coast', 'Indigenous'),
    'MOCHES': ('Coast', 'Indigenous'),
    
    # Coastal Mestizo
    'LIMA': ('Coast', 'Mestizo'),
    'TRUJILLO': ('Coast', 'Mestizo'),
    'LAMBAYEQUE': ('Coast', 'Mestizo'),
    'TUMBES': ('Coast', 'Mestizo'),
    
    # Afro-Peruvian (separate category)
    'AFRODESCENDIENTES': ('Coast', 'Afro-Peruvian'),
}

# Genes to analyze
FDA_LINKED_GENES = ['CYP2C19', 'CYP2B6', 'CYP3A5', 'SLCO1B1', 'TPMT', 'NUDT15']

# Actionable phenotypes by gene
ACTIONABLE_PHENOTYPES = {
    'CYP2C19': ['Poor Metabolizer', 'Intermediate Metabolizer'],
    'CYP2B6': ['Poor Metabolizer', 'Intermediate Metabolizer'],
    'CYP3A5': ['Poor Metabolizer', 'Intermediate Metabolizer'],
    'SLCO1B1': ['Poor Function', 'Decreased Function'],
    'TPMT': ['Poor Metabolizer', 'Intermediate Metabolizer'],
    'NUDT15': ['Poor Metabolizer', 'Intermediate Metabolizer'],
}

# ============================================================================
# SETUP
# ============================================================================

def setup_directories():
    """Create output directory."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Created output directory: {OUTPUT_DIR}")

def load_data() -> pd.DataFrame:
    """Load pharmacogenomics population summary data."""
    print("\n📖 Loading pharmacogenomics data...")
    
    if not PGX_SUMMARY.exists():
        print(f"   ❌ Input file not found: {PGX_SUMMARY}")
        sys.exit(1)
    
    df = pd.read_csv(PGX_SUMMARY)
    print(f"   Loaded {len(df)} records from population summary")
    
    # Show available columns
    print(f"   Columns: {list(df.columns)}")
    
    # Filter to Region_Group level (the most useful for our comparisons)
    df = df[df['Level'] == 'Region_Group'].copy()
    print(f"   Filtered to Region_Group level: {len(df)} records")
    
    # Show available genes
    genes = df['Gene'].unique()
    print(f"   Genes: {len(genes)} available")
    
    # Show regions and groups
    print(f"   Regions: {df['Region'].unique().tolist()}")
    print(f"   Groups: {df['Group'].unique().tolist()}")
    
    return df

# ============================================================================
# STATISTICAL TESTS
# ============================================================================

def calculate_actionable_counts(df: pd.DataFrame, gene: str, 
                                 group_col: str = 'Group') -> pd.DataFrame:
    """
    Calculate actionable vs non-actionable counts for a gene by group.
    """
    # Filter to gene
    gene_df = df[df['gene'] == gene].copy() if 'gene' in df.columns else df.copy()
    
    if gene_df.empty:
        return pd.DataFrame()
    
    # Identify actionable phenotypes for this gene
    actionable = ACTIONABLE_PHENOTYPES.get(gene, [])
    
    # Identify phenotype column
    pheno_col = None
    for col in ['phenotype', 'Phenotype', 'predicted_phenotype']:
        if col in gene_df.columns:
            pheno_col = col
            break
    
    if pheno_col is None:
        return pd.DataFrame()
    
    # Mark actionable
    gene_df['is_actionable'] = gene_df[pheno_col].isin(actionable)
    
    # Count by group
    counts = gene_df.groupby(group_col).agg(
        actionable=('is_actionable', 'sum'),
        total=('is_actionable', 'count')
    ).reset_index()
    
    counts['non_actionable'] = counts['total'] - counts['actionable']
    counts['pct_actionable'] = 100 * counts['actionable'] / counts['total']
    
    return counts

def chi_square_test(group1_action: int, group1_total: int,
                   group2_action: int, group2_total: int) -> Dict:
    """
    Perform chi-square test comparing two groups.
    Falls back to Fisher's exact test if expected counts < 5.
    """
    # Create contingency table
    # [Actionable, Non-Actionable]
    # [Group1_Act, Group1_NonAct]
    # [Group2_Act, Group2_NonAct]
    
    table = np.array([
        [group1_action, group1_total - group1_action],
        [group2_action, group2_total - group2_action]
    ])
    
    # Check minimum expected count
    row_totals = table.sum(axis=1)
    col_totals = table.sum(axis=0)
    total = table.sum()
    
    expected = np.outer(row_totals, col_totals) / total
    min_expected = expected.min()
    
    if min_expected < 5:
        # Use Fisher's exact test
        odds_ratio, p_value = fisher_exact(table)
        test_used = "Fisher's exact"
        statistic = odds_ratio
    else:
        # Use chi-square test
        chi2, p_value, dof, expected = chi2_contingency(table)
        test_used = "Chi-square"
        statistic = chi2
    
    return {
        'test': test_used,
        'statistic': statistic,
        'p_value': p_value,
        'min_expected': min_expected
    }

def compare_indigenous_vs_mestizo(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare actionable phenotype frequencies between Indigenous and Mestizo
    populations within each region using pre-aggregated summary data.
    """
    print("\n📊 Comparing Indigenous vs Mestizo within regions...")
    
    results = []
    
    for region in ['Amazon', 'Andes', 'Coast']:
        region_df = df[df['Region'] == region]
        
        for gene in FDA_LINKED_GENES:
            gene_df = region_df[region_df['Gene'] == gene]
            
            if gene_df.empty:
                continue
            
            # Get actionable phenotypes for this gene
            actionable = ACTIONABLE_PHENOTYPES.get(gene, [])
            
            # Get Indigenous data
            ind_df = gene_df[gene_df['Group'] == 'Indigenous']
            if ind_df.empty:
                continue
            
            # Sum actionable phenotype counts
            ind_action = ind_df[ind_df['Phenotype'].isin(actionable)]['Count'].sum()
            # Total is the same for all phenotypes (it's the group total)
            ind_total = ind_df['Total'].iloc[0] if len(ind_df) > 0 else 0
            
            # Get Mestizo data
            mes_df = gene_df[gene_df['Group'] == 'Mestizo']
            if mes_df.empty:
                continue
            
            mes_action = mes_df[mes_df['Phenotype'].isin(actionable)]['Count'].sum()
            mes_total = mes_df['Total'].iloc[0] if len(mes_df) > 0 else 0
            
            # Skip if insufficient sample size
            if ind_total < 5 or mes_total < 5:
                continue
            
            # Perform test
            test_result = chi_square_test(
                int(ind_action), int(ind_total),
                int(mes_action), int(mes_total)
            )
            
            results.append({
                'Region': region,
                'Gene': gene,
                'Indigenous_Actionable': int(ind_action),
                'Indigenous_Total': int(ind_total),
                'Indigenous_Pct': 100 * ind_action / ind_total if ind_total > 0 else 0,
                'Mestizo_Actionable': int(mes_action),
                'Mestizo_Total': int(mes_total),
                'Mestizo_Pct': 100 * mes_action / mes_total if mes_total > 0 else 0,
                'Test': test_result['test'],
                'Statistic': test_result['statistic'],
                'P_Value': test_result['p_value'],
                'Min_Expected': test_result['min_expected']
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Apply multiple testing correction
        results_df['P_Adjusted'] = multipletests(
            results_df['P_Value'], method='bonferroni'
        )[1]
        results_df['Significant'] = results_df['P_Adjusted'] < 0.05
    
    return results_df

def compare_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare actionable phenotype frequencies between regions
    using pre-aggregated summary data.
    """
    print("\n📊 Comparing between regions...")
    
    results = []
    region_pairs = [('Amazon', 'Andes'), ('Amazon', 'Coast'), ('Andes', 'Coast')]
    
    for gene in FDA_LINKED_GENES:
        gene_df = df[df['Gene'] == gene]
        
        if gene_df.empty:
            continue
        
        actionable = ACTIONABLE_PHENOTYPES.get(gene, [])
        
        for region1, region2 in region_pairs:
            # Get data for each region (combine Indigenous and Mestizo)
            r1_df = gene_df[gene_df['Region'] == region1]
            r2_df = gene_df[gene_df['Region'] == region2]
            
            if r1_df.empty or r2_df.empty:
                continue
            
            # Sum actionable counts across groups within each region
            r1_action = r1_df[r1_df['Phenotype'].isin(actionable)]['Count'].sum()
            # Get total (sum of totals for each group, but need to avoid double counting)
            # Since Total is per group, sum unique group totals
            r1_total = r1_df.groupby('Group')['Total'].first().sum()
            
            r2_action = r2_df[r2_df['Phenotype'].isin(actionable)]['Count'].sum()
            r2_total = r2_df.groupby('Group')['Total'].first().sum()
            
            if r1_total < 5 or r2_total < 5:
                continue
            
            test_result = chi_square_test(
                int(r1_action), int(r1_total),
                int(r2_action), int(r2_total)
            )
            
            results.append({
                'Gene': gene,
                'Comparison': f'{region1} vs {region2}',
                f'{region1}_Actionable': int(r1_action),
                f'{region1}_Total': int(r1_total),
                f'{region1}_Pct': 100 * r1_action / r1_total if r1_total > 0 else 0,
                f'{region2}_Actionable': int(r2_action),
                f'{region2}_Total': int(r2_total),
                f'{region2}_Pct': 100 * r2_action / r2_total if r2_total > 0 else 0,
                'Test': test_result['test'],
                'Statistic': test_result['statistic'],
                'P_Value': test_result['p_value']
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df['P_Adjusted'] = multipletests(
            results_df['P_Value'], method='bonferroni'
        )[1]
        results_df['Significant'] = results_df['P_Adjusted'] < 0.05
    
    return results_df

# ============================================================================
# REPORT GENERATION
# ============================================================================

def generate_manuscript_text(ind_vs_mes: pd.DataFrame, region_comp: pd.DataFrame):
    """
    Generate manuscript text for the Results section.
    """
    output_file = OUTPUT_DIR / "manuscript_text.txt"
    
    with open(output_file, 'w') as f:
        f.write("MANUSCRIPT TEXT FOR PHARMACOGENOMICS STATISTICAL COMPARISONS\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("INSERT INTO METHODS:\n")
        f.write("-" * 40 + "\n")
        f.write('"To quantify the significance of observed differences, we performed ')
        f.write('chi-square tests (or Fisher\'s exact tests when expected cell counts ')
        f.write('were <5) comparing actionable phenotype frequencies between Indigenous ')
        f.write('and mestizo groups within each region. P-values were adjusted for ')
        f.write('multiple testing using Bonferroni correction."\n\n')
        
        f.write("INSERT INTO RESULTS:\n")
        f.write("-" * 40 + "\n")
        
        # Find significant results
        if len(ind_vs_mes) > 0:
            sig_results = ind_vs_mes[ind_vs_mes['Significant'] == True]
            
            if len(sig_results) > 0:
                f.write('"Statistical comparison revealed significant differences in ')
                f.write('actionable phenotype frequencies between Indigenous and mestizo ')
                f.write('populations. ')
                
                # Highlight top findings
                for _, row in sig_results.head(3).iterrows():
                    p_str = "< 0.001" if row["P_Adjusted"] < 0.001 else f"= {row['P_Adjusted']:.3f}"
                    f.write(f'For {row["Gene"]} in {row["Region"]} populations, ')
                    f.write(f'Indigenous individuals showed {row["Indigenous_Pct"]:.0f}% ')
                    f.write(f'actionable phenotypes compared to {row["Mestizo_Pct"]:.0f}% ')
                    f.write(f'in mestizo groups (χ² = {row["Statistic"]:.1f}, ')
                    f.write(f'p {p_str}). ')
                
                f.write('These differences remained significant after Bonferroni correction ')
                f.write(f'for {len(ind_vs_mes)} comparisons."\n')
            else:
                f.write('"No comparisons reached statistical significance after ')
                f.write('Bonferroni correction."\n')
        
        f.write("\n\nSUPPLEMENTARY TABLE CAPTION:\n")
        f.write("-" * 40 + "\n")
        f.write('"Supplementary Table X. Statistical comparison of actionable ')
        f.write('pharmacogenomic phenotype frequencies between Indigenous and mestizo ')
        f.write('populations. Chi-square tests (or Fisher\'s exact tests for small ')
        f.write('sample sizes) were used to compare frequencies within each region. ')
        f.write('P-values were adjusted using Bonferroni correction for multiple ')
        f.write('comparisons. Significant comparisons (adjusted p < 0.05) are ')
        f.write('highlighted."\n')
    
    print(f"   ✓ Manuscript text: {output_file}")

def generate_supplementary_table(ind_vs_mes: pd.DataFrame, region_comp: pd.DataFrame):
    """
    Generate publication-ready supplementary table.
    """
    # Format Indigenous vs Mestizo table
    if len(ind_vs_mes) > 0:
        supp_df = ind_vs_mes.copy()
        supp_df['Indigenous (%)'] = supp_df.apply(
            lambda r: f"{r['Indigenous_Actionable']}/{r['Indigenous_Total']} ({r['Indigenous_Pct']:.1f}%)",
            axis=1
        )
        supp_df['Mestizo (%)'] = supp_df.apply(
            lambda r: f"{r['Mestizo_Actionable']}/{r['Mestizo_Total']} ({r['Mestizo_Pct']:.1f}%)",
            axis=1
        )
        supp_df['P-value'] = supp_df['P_Value'].apply(
            lambda p: '< 0.001' if p < 0.001 else f'{p:.3f}'
        )
        supp_df['Adjusted P'] = supp_df['P_Adjusted'].apply(
            lambda p: '< 0.001' if p < 0.001 else f'{p:.3f}'
        )
        supp_df['Sig.'] = supp_df['Significant'].apply(
            lambda s: '*' if s else ''
        )
        
        # Select columns for publication
        pub_cols = ['Region', 'Gene', 'Indigenous (%)', 'Mestizo (%)', 
                   'Test', 'P-value', 'Adjusted P', 'Sig.']
        pub_df = supp_df[pub_cols]
        
        output_file = OUTPUT_DIR / "supplementary_table_indigenous_mestizo.csv"
        pub_df.to_csv(output_file, index=False)
        print(f"   ✓ Supplementary table: {output_file}")
    
    # Format region comparison table
    if len(region_comp) > 0:
        output_file = OUTPUT_DIR / "supplementary_table_region_comparison.csv"
        region_comp.to_csv(output_file, index=False)
        print(f"   ✓ Region comparison table: {output_file}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    print("=" * 70)
    print("PHARMACOGENOMICS STATISTICAL COMPARISONS")
    print("=" * 70)
    print("Addressing Reviewer 3: Compare metabolizer frequencies between groups")
    print()
    
    # Setup
    setup_directories()
    
    # Check input files
    if not PGX_SUMMARY.exists():
        print(f"❌ ERROR: Input file not found: {PGX_SUMMARY}")
        print("   Please run the PyPGx pipeline first.")
        sys.exit(1)
    
    # Load data
    df = load_data()
    
    # Step 1: Indigenous vs Mestizo comparison
    print("\n" + "=" * 50)
    print("STEP 1: Indigenous vs Mestizo Comparison")
    print("=" * 50)
    
    ind_vs_mes = compare_indigenous_vs_mestizo(df)
    
    if len(ind_vs_mes) > 0:
        output_file = OUTPUT_DIR / "indigenous_vs_mestizo_tests.csv"
        ind_vs_mes.to_csv(output_file, index=False)
        print(f"\n   ✓ Results saved: {output_file}")
        
        sig_count = ind_vs_mes['Significant'].sum()
        print(f"\n   Comparisons performed: {len(ind_vs_mes)}")
        print(f"   Significant after Bonferroni: {sig_count}")
        
        if sig_count > 0:
            print("\n   Significant findings:")
            for _, row in ind_vs_mes[ind_vs_mes['Significant']].iterrows():
                print(f"   • {row['Gene']} ({row['Region']}): "
                      f"Indigenous {row['Indigenous_Pct']:.1f}% vs "
                      f"Mestizo {row['Mestizo_Pct']:.1f}% "
                      f"(p_adj = {row['P_Adjusted']:.4f})")
    else:
        print("   ⚠️ No comparisons could be performed (insufficient data)")
    
    # Step 2: Regional comparison
    print("\n" + "=" * 50)
    print("STEP 2: Regional Comparison")
    print("=" * 50)
    
    region_comp = compare_regions(df)
    
    if len(region_comp) > 0:
        output_file = OUTPUT_DIR / "region_comparison_tests.csv"
        region_comp.to_csv(output_file, index=False)
        print(f"\n   ✓ Results saved: {output_file}")
        
        sig_count = region_comp['Significant'].sum()
        print(f"\n   Comparisons performed: {len(region_comp)}")
        print(f"   Significant after Bonferroni: {sig_count}")
    
    # Step 3: Generate reports
    print("\n" + "=" * 50)
    print("STEP 3: Generating Reports")
    print("=" * 50)
    
    generate_manuscript_text(ind_vs_mes, region_comp)
    generate_supplementary_table(ind_vs_mes, region_comp)
    
    # Final summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\n📁 Output files in: {OUTPUT_DIR}/")
    print("\nKey outputs:")
    print("  • indigenous_vs_mestizo_tests.csv - Statistical comparisons")
    print("  • region_comparison_tests.csv - Regional comparisons")
    print("  • manuscript_text.txt - Ready-to-use text")
    print("  • supplementary_table_*.csv - Publication-ready tables")
    print("=" * 70)

if __name__ == "__main__":
    main()