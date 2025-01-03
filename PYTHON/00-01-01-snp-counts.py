import pysam
from collections import defaultdict
import csv

# Define sample groups
sample_groups = {
    "CHOPCCAS": [],
    "CUSCO": [],
    "IQUITOS": [],
    "MATZES": [],
    "MOCHES": [],
    "TRUJILLO": [],
    "UROS": []
}

# Path to the VCF file
vcf_file = "ANALYSIS/00-01-GEN-DIV/Peru.joint.biallelic_snps.vcf.gz"
output_file = "ANALYSIS/00-01-GEN-DIV/individual_snp_counts.csv"

def process_vcf(vcf_file, sample_groups):
    """
    Process the VCF file to count SNPs for each individual within each group.
    """
    with pysam.VariantFile(vcf_file, "r") as vcf:
        # Convert header samples to a list for easier access
        samples = list(vcf.header.samples)

        # Map sample names to groups
        for sample in samples:
            for group in sample_groups:
                if sample.startswith(group):
                    sample_groups[group].append(sample)

        # Initialize a dictionary to hold SNP counts for each individual
        individual_snp_counts = {sample: 0 for group in sample_groups for sample in sample_groups[group]}

        # Process each SNP in the VCF
        for record in vcf.fetch():
            for sample in individual_snp_counts.keys():
                genotype = record.samples[sample]["GT"]
                if genotype:
                    # Count SNPs as presence of any non-reference allele
                    if 1 in genotype:  # 1 indicates an alternate allele
                        individual_snp_counts[sample] += 1

    return individual_snp_counts

def write_results(individual_snp_counts, sample_groups, output_file):
    """
    Write the SNP counts to a CSV file, grouped by sample group.
    """
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Group", "Individual", "SNP Count"])

        for group, samples in sample_groups.items():
            for sample in samples:
                writer.writerow([group, sample, individual_snp_counts.get(sample, 0)])

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    # Process the VCF file
    individual_snp_counts = process_vcf(vcf_file, sample_groups)

    # Write results to a CSV file
    write_results(individual_snp_counts, sample_groups, output_file)

