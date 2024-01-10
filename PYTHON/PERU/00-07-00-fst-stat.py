import subprocess

def run_command(command):
    """Runs a system command."""
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

def load_population_mapping(ii_population_file):
    """Load population mapping from ii_population file."""
    population_mapping = {}
    with open(ii_population_file, 'r') as file:
        for line in file:
            parts = line.strip().split()  # Adjust split method based on actual file format
            if len(parts) >= 2:
                sample_id, population = parts[0], parts[2]
                population_mapping[sample_id] = population
    return population_mapping

def convert_ped_map_to_geno_snp_ind(ped_file, map_file, geno_out, snp_out, ind_out, ii_population_file):
    default_allele1='A'
    default_allele2='G'
    # Load population mapping
    population_mapping = load_population_mapping(ii_population_file)

    with open(map_file, 'r') as f_map, open(snp_out, 'w') as f_snp:
        for line in f_map:
            chrom, snp_id, genetic_dist, position = line.strip().split()
            # Use default alleles if actual allele information is not available
            f_snp.write(f"{snp_id} {chrom} {genetic_dist} {position} {default_allele1} {default_allele2}\n")

    # Process PED file
    with open(ped_file, 'r') as f:
        ped_data = [line.strip().split() for line in f]

    with open(ind_out, 'w') as ind_file, open(geno_out, 'w') as geno_file:
        for person in ped_data:
            population = population_mapping.get(person[1], "Unknown")

            # IND file format: Sample_name Gender Population
            ind_file.write(f"{person[1]} U {population}\n")

            # Extract genotype data and write to GENO file
            genotypes = person[6:]  # Skip initial columns (family, individual, etc.)
            geno_line = ""
            for i in range(0, len(genotypes), 2):
                # Convert genotype pair to single character
                if genotypes[i] == '0' or genotypes[i + 1] == '0':
                    geno_line += '9'  # Missing data
                elif genotypes[i] == genotypes[i + 1]:
                    geno_line += '0'  # Homozygous for the first allele
                else:
                    geno_line += '1'  # Heterozygous (or homozygous for the second allele if necessary)
            geno_file.write(geno_line + "\n")

def main():
    # Define file names and paths
    working_dir         = "/Users/apple/CloudDocs/PERU/HEINER-GUIO/ANALYSIS/15-DI-STAT/"
    plink_fileset       = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean"  # PLINK fileset prefix
    converted_prefix    = "common_variants_ibd_clean_converted"            # Prefix for converted files
    output_dir          = "D-Stat-out"                                     # Output directory for admixtools
    ii_population_file  = "ANALYSIS/10-PCA/ii_28_populations.txt"

    # Step 1: Convert PLINK files to .ped and .map format
    #run_command(f"plink --bfile {plink_fileset} --recode --out {working_dir}/{converted_prefix}")

    # Step 2: Convert .ped and .map to .geno, .snp, and .ind (Example command, replace with actual conversion command)
    ped_file = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean_converted.ped"
    map_file = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean_converted.map"
    geno_out = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean_converted.geno"
    snp_out = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean_converted.snp"
    ind_out = "ANALYSIS/15-DI-STAT/common_variants_ibd_clean_converted.ind" 
    #convert_ped_map_to_geno_snp_ind(ped_file, map_file, geno_out, snp_out, ind_out, ii_population_file)
    
    # Step 3: Run admixtools2 in R
    r_script = f"""
    library(admixtools)
    setwd('{working_dir}')
    f2_stats <- extract_f2(data = '{converted_prefix}', minmaf = 0.05, maxmiss = 0.01, outdir = '{output_dir}', pref = '{converted_prefix}', overwrite = TRUE)
    """

    with open(f"{working_dir}/run_admixtools2.R", "w") as file:
        file.write(r_script)

    run_command(f"Rscript {working_dir}/run_admixtools2.R")

if __name__ == "__main__":
    main()



