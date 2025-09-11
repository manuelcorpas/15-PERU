"""
Population Assignment Script for PCA Analysis
==============================================
This script assigns population labels to 747 genetic samples for PCA analysis.

Suggested filename: 00-03-02-population-assignment-pca.py

Input:
- VCF file with 747 samples from IBD analysis (ANALYSIS/00-06-IBD/ibd_clean.vcf)
- Demographic data with population codes (DATA/DemographicData_INS.xlsx)

Output:
- Population mapping file (ANALYSIS/00-07-PCA/ii_28_populations.txt)
- Sample IDs file (ANALYSIS/00-07-PCA/sample_ids.txt)
- Population summary plot (ANALYSIS/00-07-PCA/population_summary.png)

The script handles two types of samples:
1. Array samples (628): Matched using chip IDs with demographic data
2. WGS samples (119): Population extracted from sample name

Total: 747 samples across 28 populations
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def standardize_population_name(pop_name):
    """Standardize population names to consistent format"""
    if not pop_name:
        return pop_name
    
    # Convert to uppercase for consistency
    standardized = str(pop_name).upper()
    
    # Define mappings for inconsistent names
    mappings = {
        'AFRO_DES': 'AFRODESCENDIENTES',
        'ANCASH': 'HUARAZ',  # Ancash is the department, Huaraz is the capital
        'ASHANINKA_INS': 'ASHANINKA',
        'JACARUS': 'JAQARUS',
        'MATSIGUENKAS': 'MACHIGUENGA',
        'MOCHE': 'MOCHES',
        'MATSES': 'MATZES',  # Standardize to MATZES
        'QEROS': 'QUEROS',
        'SHIPIBO_INS': 'SHIPIBO',
        'TALLANES': 'TALLAN'
    }
    
    return mappings.get(standardized, standardized)

def extract_sample_ids_from_vcf(vcf_file, output_file):
    """Extract sample IDs from VCF file"""
    import subprocess
    
    # Use bcftools to extract sample IDs
    cmd = f"bcftools query -l {vcf_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    
    with open(output_file, 'w') as f:
        f.write(result.stdout)
    
    return result.stdout.strip().split('\n')

def create_population_mapping(sample_ids, demographic_file, output_file):
    """Create population mapping file for PCA"""
    
    # Read demographic data
    df_demo = pd.read_excel(demographic_file)
    
    # Create dictionary for array code to population mapping
    array_to_pop = {}
    for _, row in df_demo.iterrows():
        array_code = row['ARRAY_CODE']
        population = standardize_population_name(row['POPULATION_CODE'])
        array_to_pop[array_code] = population
    
    # Process samples and assign populations
    population_assignments = []
    chip_populations = {}  # Track populations by chip ID for fixing unknowns
    
    for sample_id in sample_ids:
        # Check if it's a WGS sample (contains hyphen)
        if '-' in sample_id:
            # Extract population from sample name
            population = sample_id.split('-')[0]
            # WGS populations are already in uppercase
            iid = sample_id
            fid = sample_id
        else:
            # Array sample - extract base pattern
            parts = sample_id.split('_')
            if len(parts) >= 4:
                # Extract base pattern (chip_position)
                base_pattern = f"{parts[0]}_{parts[1]}"
                chip_id = parts[0]  # Extract chip ID for tracking
                
                # Look up population
                if base_pattern in array_to_pop:
                    population = array_to_pop[base_pattern]
                    # Track this chip's population
                    if chip_id not in chip_populations:
                        chip_populations[chip_id] = population
                else:
                    population = 'UNKNOWN'
                    print(f"Warning: No population found for {base_pattern}")
                
                iid = sample_id
                fid = sample_id
            else:
                print(f"Warning: Unexpected sample format: {sample_id}")
                continue
        
        population_assignments.append({
            'IID': iid,
            'FID': fid,
            'Population': population,
            'chip_id': chip_id if '-' not in sample_id else None
        })
    
    # Fix UNKNOWN populations by looking at other samples from the same chip
    print("\nFixing UNKNOWN populations using chip-level inference...")
    fixed_count = 0
    
    for i, assignment in enumerate(population_assignments):
        if assignment['Population'] == 'UNKNOWN' and assignment['chip_id']:
            chip_id = assignment['chip_id']
            if chip_id in chip_populations:
                old_pop = assignment['Population']
                new_pop = chip_populations[chip_id]
                population_assignments[i]['Population'] = new_pop
                print(f"Fixed: {assignment['IID']} from {old_pop} to {new_pop} (chip {chip_id})")
                fixed_count += 1
    
    if fixed_count > 0:
        print(f"Fixed {fixed_count} samples using chip-level population inference")
    
    # Create output dataframe (remove the helper chip_id column)
    df_output = pd.DataFrame([{k: v for k, v in assignment.items() if k != 'chip_id'} 
                             for assignment in population_assignments])
    
    # Write to file
    df_output.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    pop_counts = df_output['Population'].value_counts().sort_index()
    print("\nPopulation Summary:")
    print("===================")
    for pop, count in pop_counts.items():
        print(f"{pop}: {count} samples")
    print(f"\nTotal populations: {len(pop_counts)}")
    print(f"Total samples: {len(df_output)}")
    
    # Check for unknown populations
    if 'UNKNOWN' in pop_counts:
        print(f"\nWarning: {pop_counts['UNKNOWN']} samples still with unknown population")
    else:
        print("\n✅ All samples successfully assigned to populations!")
    
    print(f"\nOutput saved to: {output_file}")
    
    return df_output

def create_population_summary_plot(df_populations, output_dir):
    """Create publication-ready plot of population sample counts"""
    
    # Calculate population counts and sort by count (descending)
    pop_counts = df_populations['Population'].value_counts().sort_values(ascending=True)
    
    # Set up the plot with publication-ready styling
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create color palette
    colors = plt.cm.Set3(np.linspace(0, 1, len(pop_counts)))
    
    # Create horizontal bar plot
    bars = ax.barh(range(len(pop_counts)), pop_counts.values, color=colors, 
                   edgecolor='black', linewidth=0.5, alpha=0.8)
    
    # Customize the plot
    ax.set_yticks(range(len(pop_counts)))
    ax.set_yticklabels(pop_counts.index, fontsize=10)
    ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
    ax.set_ylabel('Population', fontsize=12, fontweight='bold')
    ax.set_title('Sample Distribution Across Peruvian Populations\n(N = 747 samples)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add sample count labels on bars
    for i, (bar, count) in enumerate(zip(bars, pop_counts.values)):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2, 
                str(count), ha='left', va='center', fontsize=9, fontweight='bold')
    
    # Add grid for better readability
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Adjust layout and styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    
    # Set x-axis limits with some padding
    max_count = pop_counts.values.max()
    ax.set_xlim(0, max_count * 1.15)
    
    # Add summary statistics as text box
    total_samples = len(df_populations)
    total_populations = len(pop_counts)
    median_samples = pop_counts.median()
    
    stats_text = f'Total Populations: {total_populations}\n'
    stats_text += f'Total Samples: {total_samples}\n'
    stats_text += f'Median Samples/Pop: {median_samples:.0f}\n'
    stats_text += f'Range: {pop_counts.min()}-{pop_counts.max()}'
    
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, 
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8),
            verticalalignment='bottom', horizontalalignment='right',
            fontsize=9, fontfamily='monospace')
    
    # Tight layout
    plt.tight_layout()
    
    # Save the plot
    plot_file = os.path.join(output_dir, 'population_summary.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
    
    # Also save as PDF for publication
    plot_file_pdf = os.path.join(output_dir, 'population_summary.pdf')
    plt.savefig(plot_file_pdf, bbox_inches='tight', facecolor='white')
    
    print(f"Population summary plots saved:")
    print(f"  - PNG: {plot_file}")
    print(f"  - PDF: {plot_file_pdf}")
    
    # Show plot (comment out if running in headless environment)
    # plt.show()
    
    plt.close()
    
    return plot_file

# Main workflow
if __name__ == "__main__":
    # File paths
    vcf_file = 'ANALYSIS/00-06-IBD/ibd_clean.vcf'
    demographic_file = 'DATA/DemographicData_INS.xlsx'
    output_dir = 'ANALYSIS/00-07-PCA'
    output_file = os.path.join(output_dir, 'ii_28_populations.txt')
    sample_ids_file = os.path.join(output_dir, 'sample_ids.txt')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory created/verified: {output_dir}")
    
    print("\nStep 1: Extracting sample IDs from VCF...")
    sample_ids = extract_sample_ids_from_vcf(
        vcf_file, 
        sample_ids_file
    )
    print(f"Extracted {len(sample_ids)} sample IDs")
    print(f"Sample IDs saved to: {sample_ids_file}")
    
    print("\nStep 2: Creating population mapping...")
    df_populations = create_population_mapping(
        sample_ids,
        demographic_file,
        output_file
    )
    
    print("\nStep 3: Creating population summary plot...")
    plot_file = create_population_summary_plot(df_populations, output_dir)
    
    print("\n✅ Population mapping and visualization complete!")
    print(f"Ready for PCA analysis with {len(df_populations['Population'].unique())} populations")
    print(f"\nFiles created:")
    print(f"  - Population mapping: {output_file}")
    print(f"  - Sample IDs: {sample_ids_file}")
    print(f"  - Summary plot: {plot_file}")