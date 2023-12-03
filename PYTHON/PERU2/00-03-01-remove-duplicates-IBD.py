import subprocess
import pandas as pd

def run_ibd_analysis(plink_path, input_file, output_prefix):
    # Step 1: Convert VCF to PLINK binary format
    vcf_to_binary_cmd = [
        plink_path,
        '--vcf', input_file,
        '--double-id', 
        '--make-bed',
        '--out', output_prefix
    ]
    subprocess.run(vcf_to_binary_cmd, check=True)

    # Step 2: Run IBD analysis in PLINK
    ibd_cmd = [
        plink_path,
        '--bfile', output_prefix,  # Using the converted binary file
        '--genome',
        '--out', output_prefix
    ]
    subprocess.run(ibd_cmd, check=True)

    # Step 3: Parse the PLINK .genome output to identify duplicates
    genome_file = f'{output_prefix}.genome'
    df = pd.read_csv(genome_file, delim_whitespace=True)
    
    # Identify duplicates (e.g., PI_HAT close to 1)
    duplicates = df[df['PI_HAT'] > 0.95]  # Adjust this threshold as needed

    # Create a list of individuals to remove (consider both FID and IID)
    to_remove = set(duplicates['IID1'].tolist() + duplicates['IID2'].tolist())

    # Step 4: Write the list of individuals to a file with FID and IID being the same
    with open('09-IBD/to_remove.txt', 'w') as f:
        for iid in to_remove:
            f.write(f'{iid} {iid}\n')  # Writing IID as both FID and IID

    # Step 5: Remove duplicates using PLINK
    remove_cmd = [
        plink_path,
        '--bfile', output_prefix,
        '--remove', '09-IBD/to_remove.txt',
        '--make-bed',
        '--out', f'{output_prefix}_clean'
    ]
    subprocess.run(remove_cmd, check=True)

    # Step 6: Convert the cleaned dataset back to VCF
    final_vcf_cmd = [
        plink_path,
        '--bfile', f'{output_prefix}_clean',
        '--recode', 'vcf',
        '--out', f'{output_prefix}_clean'
    ]
    subprocess.run(final_vcf_cmd, check=True)


# Example usage
plink_path = 'plink'  # Path to the PLINK executable
input_file = 'INPUT/VCF/common_variants.vcf.gz'  # Path to input dataset (without .bed/.bim/.fam extension)
output_prefix = 'ANALYSIS/09-IBD/ibd'  # Output prefix for PLINK files

run_ibd_analysis(plink_path, input_file, output_prefix)

