import subprocess
import os

# I need to do this first because I start from a plink file
def convert_plink_to_vcf_gz(plink_prefix, output_vcf_gz):
    """
    Convert PLINK files to gzipped VCF format with an index.
    :param plink_prefix: Prefix for the PLINK files (without .bed, .bim, .fam).
    :param output_vcf_gz: Path for the output gzipped VCF file.
    """
    # Convert PLINK files to VCF
    vcf_temp = output_vcf_gz.replace('.vcf.gz', '.vcf')
    plink_command = f"plink --bfile {plink_prefix} --recode vcf --out {vcf_temp.replace('.vcf', '')}"
    subprocess.run(plink_command, shell=True)

    # Gzip the VCF file
    bgzip_command = f"bgzip -c {vcf_temp} > {output_vcf_gz}"
    subprocess.run(bgzip_command, shell=True)

    # Create an index for the gzipped VCF file
    tabix_command = f"tabix -p vcf {output_vcf_gz}"
    subprocess.run(tabix_command, shell=True)

    # Remove the temporary VCF file
    os.remove(vcf_temp)

# I want to rely on bcftools to do the intersection between my clean IBD files and SGDP
def run_bcftools_intersection_with_samples(vcf_file1, vcf_file2, output_file):
    try:
        # Create a temporary directory to store intermediate files
        temp_dir = "temp_bcftools_output"
        os.makedirs(temp_dir, exist_ok=True)

        # Step 1: Find common variants and output them into separate files
        command_isec = f"bcftools isec -n=2 {vcf_file1} {vcf_file2} -p {temp_dir}"
        subprocess.run(command_isec, shell=True, check=True)

        # Step 2: Compress the intermediate VCF files with bgzip
        common_variants_file1 = os.path.join(temp_dir, "0000.vcf")
        common_variants_file2 = os.path.join(temp_dir, "0001.vcf")
        command_compress1 = f"bgzip -c {common_variants_file1} > {common_variants_file1}.gz"
        command_compress2 = f"bgzip -c {common_variants_file2} > {common_variants_file2}.gz"
        subprocess.run(command_compress1, shell=True, check=True)
        subprocess.run(command_compress2, shell=True, check=True)

        # Step 3: Index the compressed VCF files
        command_index1 = f"bcftools index {common_variants_file1}.gz"
        command_index2 = f"bcftools index {common_variants_file2}.gz"
        subprocess.run(command_index1, shell=True, check=True)
        subprocess.run(command_index2, shell=True, check=True)

        # Step 4: Merge the samples from both VCFs for the common variants
        command_merge = f"bcftools merge {common_variants_file1}.gz {common_variants_file2}.gz -Oz -o {output_file}"
        subprocess.run(command_merge, shell=True, check=True)

        print(f"Common variants with samples from both VCFs are stored in {output_file}")

        # Optional: Clean up the temporary directory
        # shutil.rmtree(temp_dir)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

# Example usage
plink_prefix = "ANALYSIS/12-SGDP/cteam_extended.v4.maf0.1perc"  # Replace with your PLINK file prefix
output_vcf_gz = "INPUT/VCF/sgdp.cteam_extended.v4.maf0.1perc.vcf.gz"  # Replace with your desired output VCF gzipped file path


# BEWARE I'VE COMMENTED THIS OUT
#convert_plink_to_vcf_gz(plink_prefix, output_vcf_gz)


# Replace with your VCF file paths
vcf_file1 = 'INPUT/VCF/sgdp.cteam_extended.v4.maf0.1perc.vcf.gz'
vcf_file2 = 'INPUT/VCF/common_variants_ibd_clean.vcf.gz'

# Output file
output_file = 'INPUT/VCF/common_variants_peru_sgdp.vcf.gz'

# Run the function
run_bcftools_intersection_with_samples(vcf_file1, vcf_file2, output_file)


'''
The bcftools isec command is used with the -p option to specify an output directory where the result files will be stored. The result includes four files by default, where 0000.vcf and 0001.vcf contain variants unique to the first and second files, respectively, and 0002.vcf contains variants common to both.
The -Oz option compresses the output files in VCF.GZ format.

Uses bcftools isec to find common variants and outputs them into separate files in a temporary directory.
Uses bcftools merge to combine the samples from both VCF files for these common variants into a single file.
Optionally, you can uncomment the shutil.rmtree(temp_dir) line to clean up the temporary directory after the script runs.

In this script, after generating the intermediate VCF files (0000.vcf and 0001.vcf), we compress them using bgzip before merging them with bcftools merge. This should resolve the issue you encountered.
'''

