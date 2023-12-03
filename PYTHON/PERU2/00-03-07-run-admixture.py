import pandas as pd
import subprocess
import os

# Path to your files
population_file_path = 'ANALYSIS/10-PCA/ii_28_populations.txt'
vcf_file_path = 'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz'
plink_executable_path = 'plink'  # if PLINK is used for conversion
admixture_executable_path = 'admixture'
plink_format_output = 'ANALYSIS/11-ADMIXTURE/common_variants_ibd_clean'

# Read population data
population_data = pd.read_csv(population_file_path, sep='\t')

# VCF file handling and conversion to PLINK format
# This step depends on the format of your VCF file and the tools you have.
# You might use pyvcf or call PLINK via subprocess for this.

# Example of calling PLINK to convert VCF to PLINK format
subprocess.run([plink_executable_path, '--vcf', vcf_file_path, '--make-bed', '--out', plink_format_output])

# Running ADMIXTURE
# Adjust the number of populations (K) as needed
K = 3  # Example number of ancestral populations
subprocess.run([admixture_executable_path, plink_format_output+'.bed', str(K)])

# Process ADMIXTURE output
# Depends on your analysis needs. The output files from ADMIXTURE will typically include ancestry proportions for each individual.

# Example: Read and process ADMIXTURE output
# admixture_output = pd.read_csv('plink_format_output.' + str(K) + '.Q', header=None, sep=' ')
# Do further processing as needed

print("ADMIXTURE analysis complete.")

