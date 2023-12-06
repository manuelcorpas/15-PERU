import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import seaborn as sns

def run_plink_pca(vcf_file, output_prefix):
    """
    Run PLINK to perform PCA on a VCF file.
    """
    plink_cmd = f"plink --vcf {vcf_file} --pca --double-id --out {output_prefix}"
    subprocess.run(plink_cmd, shell=True, check=True)
'''
def plot_pca_results(pca_file, population_file):
    pca_data = pd.read_csv(pca_file, delim_whitespace=True, header=None)
    pca_data.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, pca_data.shape[1] - 1)]

    # Read and adjust the population data to match the PCA data format
    pop_data = pd.read_csv(population_file, delimiter='\t')
    
    # Check if the IID needs duplication based on '_' presence
    # Adjust the 'IID' by duplicating it only if it contains an underscore
    #pop_data['IID'] = pop_data['IID'].apply(lambda x: f"{x}_{x}" if '_' in x else x)

    # Merge the PCA and population data
    merged_data = pd.merge(pca_data, pop_data, on='IID')

    # Check if the merged data is empty
    if merged_data.empty:
        print("Merged data is empty. Possible mismatch in IDs.")
        return

    # Plotting with alphabetically sorted population labels
    plt.figure(figsize=(12, 8))
    unique_pops = sorted(merged_data['Population'].unique())
    texts = []

    for pop in unique_pops:
        subset = merged_data[merged_data['Population'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.7)

        # Calculate the centroid of each population group
        centroid = subset[['PC1', 'PC2']].mean()
        texts.append(plt.text(centroid['PC1'], centroid['PC2'], pop, fontsize=10, weight='bold'))

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA Peru (746 Samples, 28 Populations) SGDP (345 Samples, 164 Populations')
    plt.legend()

    # Adjust text positions
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.3))

    plt.show()
    # Print the count of populations
    print(f"Total number of populations plotted: {len(unique_pops)}")
'''

def plot_pca_results(pca_file, population_file):
    pca_data = pd.read_csv(pca_file, delim_whitespace=True, header=None)
    pca_data.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, pca_data.shape[1] - 1)] 

    pop_data = pd.read_csv(population_file, delimiter='\t')
    merged_data = pd.merge(pca_data, pop_data, on='IID')

    if merged_data.empty:
        print("Merged data is empty. Possible mismatch in IDs.")
        return

    plt.figure(figsize=(20, 15))  # Increased plot size
    unique_pops = sorted(merged_data['Population'].unique())
    texts = []

    for pop in unique_pops:
        subset = merged_data[merged_data['Population'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.7)

        centroid = subset[['PC1', 'PC2']].mean()
        texts.append(plt.text(centroid['PC1'], centroid['PC2'], pop, fontsize=9, weight='bold'))

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA Peru (746 Samples, 28 Populations) SGDP (46 Samples, 24 Close Populations)')
    plt.legend()
    # Remove or adjust the legend
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # Uncomment this line to move the legend outside

    #adjust_text(texts, expand_points=(1.2, 1.2), arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    #plt.show()
    # ... [previous code to set up the plot] ...

    adjust_text(texts, 
            expand_points=(1.2, 1.2)) 
            #expand_text=(1.2, 1.2), 
            #arrowprops=dict(arrowstyle="-", color='black', lw=0.5), 
            #force_text=0.5, 
            #force_points=0.5,
            #autoalign=True)

    plt.show()

    print(f"Total number of populations plotted: {len(unique_pops)}")

def main():
    vcf_file = 'ANALYSIS/12-SGDP/common_variants_peru_sgdp.vcf.gz'
    output_prefix = 'ANALYSIS/12-SGDP/pca_output'
    population_file = 'ANALYSIS/12-SGDP/sgdp-peru-sample-ids-population-ids-24-SDGP-24-close-pops.txt'  # Path to the file mapping samples to populations

    run_plink_pca(vcf_file, output_prefix)
    pca_file = f"{output_prefix}.eigenvec"

    if os.path.exists(pca_file):
        print("PCA completed successfully. Plotting results...")
        plot_pca_results(pca_file, population_file)
    else:
        print("PCA failed or output files not found.")

if __name__ == "__main__":
    main()

