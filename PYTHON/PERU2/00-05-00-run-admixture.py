import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

# Path to your files
population_file_path = 'ANALYSIS/10-PCA/ii_28_populations.txt'
vcf_file_path = 'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz'
plink_executable_path = 'plink'  # if PLINK is used for conversion
admixture_executable_path = '/Users/apple/Software/admixture_macosx-1.3.0/admixture'
plink_format_output = 'ANALYSIS/11-ADMIXTURE/common_variants_ibd_clean'

'''
# Convert VCF to PLINK format
subprocess.run([plink_executable_path, '--vcf', vcf_file_path, '--make-bed', '--double-id', '--out', plink_format_output])

# Prune for LD
subprocess.run([plink_executable_path, '--bfile', plink_format_output, 
                '--indep-pairwise', '50', '5', '0.2', 
                '--out', plink_format_output+'_pruned'])

# Extract the pruned dataset
subprocess.run([plink_executable_path, '--bfile', plink_format_output, 
                '--extract', plink_format_output+'_pruned.prune.in', 
                '--make-bed', '--out', plink_format_output +'_pruned'])
'''
# Running ADMIXTURE with 8 threads
K = 6  # Adjust the number of populations as necessary
admixture_output_file = plink_format_output+'_pruned.' + str(K)
subprocess.run([admixture_executable_path, plink_format_output+'_pruned.bed', str(K), '-j8'])

print("ADMIXTURE analysis complete.")


