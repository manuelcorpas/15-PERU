import subprocess
import pandas as pd
import gzip
import os
import math

def vcf_to_plink(vcf_file, plink_prefix):
    # Convert VCF to PLINK format
    subprocess.run(f"plink --vcf {vcf_file} --make-bed --double-id --out {plink_prefix}", shell=True)

def ld_prune(plink_prefix, ld_prefix, window_size=50, step_size=5, threshold=0.2):
    """
    Perform LD pruning.
    :param plink_prefix: Prefix of the input PLINK files (before pruning).
    :param output_prefix: Prefix for the output PLINK files (after pruning).
    :param window_size: Window size in SNPs (default 50).
    :param step_size: Step size for the window (default 5).
    :param threshold: LD threshold for pruning (default 0.2).
    """
    # Perform LD pruning
    subprocess.run(
        f"plink --bfile {plink_prefix} --indep-pairwise {window_size} {step_size} {threshold} --out {ld_prefix}", 
        shell=True
    )
    # Create a new file with the pruned dataset
    subprocess.run(
        f"plink --bfile {plink_prefix} --extract {ld_prefix}.prune.in --make-bed --out {ld_prefix}_pruned",
        shell=True
    )
def create_population_subsets(population_file, ld_pruned, population_out):
    # Read the population file with headers
    pop_data = pd.read_csv(population_file, sep='\t')

    # Iterate over each population and create a subset PLINK file
    for population in pop_data['Population'].unique():
        subset_file = f"{population_out}/{population}.txt"
        # Filter for the current population and write FID and IID to the subset file
        pop_subset = pop_data[pop_data['Population'] == population]
        pop_subset[['FID', 'IID']].to_csv(subset_file, index=False, sep=' ', header=False)

        # Create a subset PLINK file for the population
        subprocess.run(f"plink --bfile {ld_pruned} --keep {subset_file} --make-bed --freq --out {population_out}/{population}", shell=True)

def plink_freq_to_treemix(population_out, treemix_in):
    """
    Convert PLINK frequency files to TreeMix input format.
    :param population_out: Directory where population .frq files are.
    :param treemix_in: Output file path for TreeMix input formatted data.
    """
    # Read and process each frequency file
   
    frq_files   = []
    populations = []
    pop_count   = 0
    for root, dirs, files in os.walk(population_out):
        for file in files:
            if file.endswith('.frq'):
                pop_count += 1
                file_path = os.path.join(root, file)
                frq_files.append({'name': file, 'path': file_path})
                populations.append(file.split('.')[0])
            if pop_count == 100:
                break

    snp_data = {}
    for frq_file in frq_files:
        freq_file = f"{frq_file['path']}"
        print(f"Reading file: {freq_file}")
        pop       = frq_file['name'].split('.')[0]        
        with open(freq_file, 'r') as file:
            df = pd.read_csv(freq_file, delim_whitespace=True)
            for index, row in df.iterrows():
                snp = row['SNP']
                if snp not in snp_data:
                    snp_data[snp] = {}
                # Extract allele counts
                minor_allele_freq  = row['MAF'] if not math.isnan(row['MAF']) else 0
                total_alleles      = row['NCHROBS'] if not math.isnan(row['NCHROBS']) else 0
                minor_allele_count = int(minor_allele_freq * total_alleles)
                ref_allele_count   = int(total_alleles - minor_allele_count)
                print(f"{snp},{pop},{ref_allele_count},{minor_allele_count}")
                snp_data[snp][pop] = f"{ref_allele_count},{minor_allele_count}"
    
    # Write to TreeMix format
    with gzip.open(treemix_in, 'wt') as f_out:
        # Write header
        f_out.write(' '.join(populations) + '\n')
        # Write SNP data
        for snp in snp_data:
            line = ' '.join(snp_data[snp].get(pop, '0,0') for pop in populations)
            f_out.write(line + '\n')

vcf_file         = 'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz'
plink_prefix     = 'ANALYSIS/13-TREEMIX/common_variants_peru'
ld_prefix        = 'ANALYSIS/13-TREEMIX/common_variants_peru_ld'
population_file  = 'ANALYSIS/10-PCA/ii_28_populations.txt'  # Path to your population file
ld_pruned        = 'ANALYSIS/13-TREEMIX/common_variants_peru_ld_pruned'  # Prefix for the pruned PLINK files
population_out   = 'ANALYSIS/13-TREEMIX/POPULATIONS'
treemix_in       = 'ANALYSIS/13-TREEMIX/treemix_in.gz'

vcf_to_plink(vcf_file, plink_prefix)
ld_prune(plink_prefix, ld_prefix)
create_population_subsets(population_file,ld_pruned,population_out)
plink_freq_to_treemix(population_out, treemix_in)




'''
# Example usage
vcf_file = "ANALYSIS/12-SGDP/common_variants_peru_sgdp.vcf.gz"
plink_output_prefix = "ANALYSIS/13-TREEMIX/common_variants_peru_sgdp"
population_file = "ANALYSIS/12-SGDP/sgdp-peru-sample-ids-population-ids-24-SDGP-24-close-pops.txt"
treemix_output_file = "ANALYSIS/13-TREEMIX/common_variants_peru_sgdp_treemix"
output_prefix = "ANALYSIS/13-TREEMIX/common_variants_peru_sgdp.pop"

#convert_vcf_to_plink(vcf_file, plink_output_prefix)
allele_counts_data = calculate_allele_frequencies(plink_output_prefix, population_file,output_prefix)
#format_for_treemix(allele_counts_data, treemix_output_file)
#population_data = pd.read_csv(population_file, sep='\t')
#population_map = pd.Series(population_data.Population.values, index=population_data.IID).to_dict()
#frq_file = 'ANALYSIS/13-TREEMIX/common_variants_peru_sgdp.frq'
#output_file = 'ANALYSIS/13-TREEMIX/treemix_input.txt.gz'

#process_frq_file(frq_file, population_map, output_file)
'''
