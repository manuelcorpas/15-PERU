#!/usr/bin/env python3
"""
00-06-05-pca-variance-explained.py

Extract variance explained percentages from PCA eigenvalues and generate
updated axis labels for Figure 3.

Purpose:
    Address Reviewer 1 and 3's comments: "Add variance explained percentages
    to PCA axes."

Usage:
    python 00-06-05-pca-variance-explained.py

Input:
    - PLINK .eigenval files from PCA analysis

Output:
    ANALYSIS/00-16-PCA-VAR/variance_explained.csv
    ANALYSIS/00-16-PCA-VAR/figure_axis_labels.txt
    ANALYSIS/00-16-PCA-VAR/manuscript_text.txt

Author: Generated for Nature Health genomics paper revision
Date: 2025
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Possible eigenval file locations
EIGENVAL_FILES = [
    "ANALYSIS/00-08-PCA/pca_peru.eigenval",
    "ANALYSIS/00-08-PCA/pca_common.eigenval",
    "ANALYSIS/00-08-PCA/pca_sgdp.eigenval",
    "ANALYSIS/10-PCA/pca_output.eigenval",
    "ANALYSIS/12-SGDP/pca_output.eigenval"
]

# Output directory
OUTPUT_DIR = Path("ANALYSIS/00-16-PCA-VAR")

# ============================================================================
# SETUP
# ============================================================================

def setup_directories():
    """Create output directory."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Created output directory: {OUTPUT_DIR}")

# ============================================================================
# EIGENVALUE ANALYSIS
# ============================================================================

def load_eigenvalues(eigenval_path: str) -> Optional[np.ndarray]:
    """
    Load eigenvalues from PLINK .eigenval file.
    """
    if not os.path.exists(eigenval_path):
        return None
    
    try:
        eigenvals = pd.read_csv(eigenval_path, header=None, names=['eigenval'])
        return eigenvals['eigenval'].values
    except Exception as e:
        print(f"   ⚠️ Error reading {eigenval_path}: {e}")
        return None

def calculate_variance_explained(eigenvals: np.ndarray) -> pd.DataFrame:
    """
    Calculate variance explained by each PC.
    """
    total_var = eigenvals.sum()
    var_explained = (eigenvals / total_var) * 100
    cumulative = np.cumsum(var_explained)
    
    results = pd.DataFrame({
        'PC': [f'PC{i+1}' for i in range(len(eigenvals))],
        'Eigenvalue': eigenvals,
        'Variance_Explained_Pct': var_explained,
        'Cumulative_Pct': cumulative
    })
    
    return results

def find_eigenval_files() -> Dict[str, str]:
    """
    Find all available eigenval files.
    """
    found = {}
    
    for path in EIGENVAL_FILES:
        if os.path.exists(path):
            name = Path(path).stem.replace('_output', '').replace('.eigenval', '')
            found[name] = path
    
    # Also search for any .eigenval files in analysis directories
    for analysis_dir in ['ANALYSIS/00-08-PCA', 'ANALYSIS/10-PCA', 'ANALYSIS/12-SGDP']:
        if os.path.exists(analysis_dir):
            for f in os.listdir(analysis_dir):
                if f.endswith('.eigenval'):
                    name = f.replace('.eigenval', '')
                    if name not in found:
                        found[name] = os.path.join(analysis_dir, f)
    
    return found

# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_scree(var_df: pd.DataFrame, name: str, output_dir: Path):
    """
    Create scree plot showing variance explained by each PC.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Bar plot of variance explained
    ax = axes[0]
    ax.bar(range(1, len(var_df)+1), var_df['Variance_Explained_Pct'], 
           color='steelblue', edgecolor='black', alpha=0.7)
    ax.set_xlabel('Principal Component', fontsize=12)
    ax.set_ylabel('Variance Explained (%)', fontsize=12)
    ax.set_title(f'Variance Explained by Each PC\n({name})', fontsize=14)
    ax.set_xticks(range(1, len(var_df)+1))
    ax.grid(True, alpha=0.3)
    
    # Cumulative variance plot
    ax = axes[1]
    ax.plot(range(1, len(var_df)+1), var_df['Cumulative_Pct'], 
            'o-', color='steelblue', linewidth=2, markersize=8)
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50%')
    ax.axhline(y=80, color='green', linestyle='--', alpha=0.5, label='80%')
    ax.set_xlabel('Number of Principal Components', fontsize=12)
    ax.set_ylabel('Cumulative Variance Explained (%)', fontsize=12)
    ax.set_title(f'Cumulative Variance Explained\n({name})', fontsize=14)
    ax.set_xticks(range(1, len(var_df)+1))
    ax.set_ylim(0, 100)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    output_file = output_dir / f"scree_plot_{name}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   ✓ Scree plot: {output_file}")

# ============================================================================
# REPORT GENERATION
# ============================================================================

def generate_axis_labels(var_df: pd.DataFrame) -> Dict[str, str]:
    """
    Generate figure axis labels with variance explained percentages.
    """
    labels = {}
    
    for i in range(min(10, len(var_df))):
        pc = f'PC{i+1}'
        pct = var_df.iloc[i]['Variance_Explained_Pct']
        labels[pc] = f'{pc} ({pct:.1f}%)'
    
    return labels

def generate_reports(all_results: Dict[str, pd.DataFrame]):
    """
    Generate all output files.
    """
    print("\n📝 Generating reports...")
    
    # Combine all results
    combined = []
    for name, df in all_results.items():
        df_copy = df.copy()
        df_copy['Dataset'] = name
        combined.append(df_copy)
    
    if combined:
        combined_df = pd.concat(combined, ignore_index=True)
        
        # Save combined results
        output_file = OUTPUT_DIR / "variance_explained.csv"
        combined_df.to_csv(output_file, index=False)
        print(f"   ✓ Variance explained: {output_file}")
    
    # Generate axis labels for each dataset
    labels_file = OUTPUT_DIR / "figure_axis_labels.txt"
    with open(labels_file, 'w') as f:
        f.write("FIGURE AXIS LABELS WITH VARIANCE EXPLAINED\n")
        f.write("=" * 70 + "\n\n")
        
        for name, df in all_results.items():
            f.write(f"\n{name.upper()}\n")
            f.write("-" * 40 + "\n")
            
            labels = generate_axis_labels(df)
            for pc, label in labels.items():
                f.write(f"  {pc}: {label}\n")
            
            f.write("\n")
    
    print(f"   ✓ Axis labels: {labels_file}")
    
    # Generate manuscript text
    manuscript_file = OUTPUT_DIR / "manuscript_text.txt"
    with open(manuscript_file, 'w') as f:
        f.write("MANUSCRIPT TEXT FOR PCA VARIANCE EXPLAINED\n")
        f.write("=" * 70 + "\n\n")
        
        # Get Peru PCA results (main analysis)
        peru_df = all_results.get('pca_peru', all_results.get(list(all_results.keys())[0]))
        
        pc1_var = peru_df.iloc[0]['Variance_Explained_Pct']
        pc2_var = peru_df.iloc[1]['Variance_Explained_Pct']
        pc3_var = peru_df.iloc[2]['Variance_Explained_Pct']
        cum_3pc = peru_df.iloc[2]['Cumulative_Pct']
        
        f.write("INSERT INTO FIGURE 3 LEGEND:\n")
        f.write("-" * 40 + "\n")
        f.write(f'"Variance explained by each PC is indicated on axes. ')
        f.write(f'PC1 ({pc1_var:.1f}%) captures the primary axis of genetic ')
        f.write(f'variation separating Amazonian from Andean populations. ')
        f.write(f'PC2 ({pc2_var:.1f}%) distinguishes highland from coastal groups. ')
        f.write(f'The first three PCs together explain {cum_3pc:.1f}% of total variance."\n\n')
        
        f.write("UPDATED AXIS LABELS FOR FIGURE 3:\n")
        f.write("-" * 40 + "\n")
        f.write(f'Panel A x-axis: "PC1 ({pc1_var:.1f}%)"\n')
        f.write(f'Panel A y-axis: "PC2 ({pc2_var:.1f}%)"\n')
        f.write(f'Panel B x-axis: "PC3 ({pc3_var:.1f}%)"\n')
        f.write(f'Panel B y-axis: "PC2 ({pc2_var:.1f}%)"\n')
        f.write(f'Panel C x-axis: "PC3 ({pc3_var:.1f}%)"\n')
        f.write(f'Panel C y-axis: "PC1 ({pc1_var:.1f}%)"\n\n')
        
        f.write("FOR RESULTS TEXT (optional):\n")
        f.write("-" * 40 + "\n")
        f.write(f'"Principal component analysis revealed substantial genetic structure ')
        f.write(f'among Peruvian populations. The first principal component (PC1, ')
        f.write(f'{pc1_var:.1f}% variance explained) separated Amazonian from Andean ')
        f.write(f'populations, while PC2 ({pc2_var:.1f}%) distinguished highland from ')
        f.write(f'coastal groups (Fig. 3A)."\n')
    
    print(f"   ✓ Manuscript text: {manuscript_file}")
    
    return combined_df if combined else pd.DataFrame()

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    print("=" * 70)
    print("PCA VARIANCE EXPLAINED EXTRACTION")
    print("=" * 70)
    print("Addressing Reviewers 1 & 3: Add variance explained to PCA axes")
    print()
    
    # Setup
    setup_directories()
    
    # Find eigenval files
    print("\n📖 Searching for eigenvalue files...")
    eigenval_files = find_eigenval_files()
    
    if not eigenval_files:
        print("\n❌ No eigenvalue files found")
        print("   Please run PCA analysis first (plink --pca)")
        print("\n   Expected locations:")
        for path in EIGENVAL_FILES:
            print(f"     - {path}")
        sys.exit(1)
    
    print(f"   Found {len(eigenval_files)} eigenvalue file(s):")
    for name, path in eigenval_files.items():
        print(f"     - {name}: {path}")
    
    # Process each file
    print("\n" + "=" * 50)
    print("CALCULATING VARIANCE EXPLAINED")
    print("=" * 50)
    
    all_results = {}
    
    for name, path in eigenval_files.items():
        print(f"\n   Processing: {name}")
        
        eigenvals = load_eigenvalues(path)
        
        if eigenvals is not None:
            var_df = calculate_variance_explained(eigenvals)
            all_results[name] = var_df
            
            # Print summary
            print(f"   PC1: {var_df.iloc[0]['Variance_Explained_Pct']:.2f}%")
            print(f"   PC2: {var_df.iloc[1]['Variance_Explained_Pct']:.2f}%")
            print(f"   PC3: {var_df.iloc[2]['Variance_Explained_Pct']:.2f}%")
            print(f"   Cumulative (PC1-3): {var_df.iloc[2]['Cumulative_Pct']:.2f}%")
            
            # Create scree plot
            plot_scree(var_df, name, OUTPUT_DIR)
    
    if not all_results:
        print("\n❌ No eigenvalues could be loaded")
        sys.exit(1)
    
    # Generate reports
    print("\n" + "=" * 50)
    print("GENERATING REPORTS")
    print("=" * 50)
    
    combined_df = generate_reports(all_results)
    
    # Final summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    
    # Show main results
    main_name = 'pca_peru' if 'pca_peru' in all_results else list(all_results.keys())[0]
    main_df = all_results[main_name]
    
    print(f"\n📊 KEY VALUES FOR FIGURE 3 ({main_name}):")
    print(f"   PC1: {main_df.iloc[0]['Variance_Explained_Pct']:.1f}%")
    print(f"   PC2: {main_df.iloc[1]['Variance_Explained_Pct']:.1f}%")
    print(f"   PC3: {main_df.iloc[2]['Variance_Explained_Pct']:.1f}%")
    
    print(f"\n📁 Output files in: {OUTPUT_DIR}/")
    print("=" * 70)

if __name__ == "__main__":
    main()
