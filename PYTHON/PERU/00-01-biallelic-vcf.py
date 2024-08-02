import subprocess

def run_bcftools(input_vcf, output_vcf):
    try:
        # Command to filter biallelic variants
        filter_cmd = [
            "bcftools", "view",
            "-m2", "-M2", "-v", "snps",
            "-Oz", "-o", output_vcf,
            input_vcf
        ]
        
        # Run the filter command
        subprocess.run(filter_cmd, check=True)
        print(f"Biallelic variant filtering completed. Output saved to {output_vcf}")

        # Command to index the output VCF file
        index_cmd = ["bcftools", "index", output_vcf]

        # Run the index command
        subprocess.run(index_cmd, check=True)
        print(f"Indexing of {output_vcf} completed.")

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bcftools: {e}")
    except FileNotFoundError as e:
        print(f"bcftools not found: {e}")

if __name__ == "__main__":
    input_vcf = "ANALYSIS/02-PGEN-PSAM-PVAR-2-VCF/Peru.joint.vcf.gz"
    output_vcf = "ANALYSIS/02-PGEN-PSAM-PVAR-2-VCF/Peru.biallelic.vcf.gz"

    run_bcftools(input_vcf, output_vcf)

