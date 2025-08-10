#!/usr/bin/env python

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

# Output file for counts of private variants
output_file = "ANALYSIS/00-01-GEN-DIV/private_variant_counts.csv"


def process_vcf_for_private_variants(vcf_path, sample_groups):
    """
    Process the VCF to count private variants in each group.

    A variant is 'private' to a group if:
       - That group has at least one sample carrying an alternate allele
       - No other group carries that variant (all reference or missing in other groups)
    """
    # Map each sample to its group
    sample_to_group = {}

    # Open the VCF
    with pysam.VariantFile(vcf_path, "r") as vcf:
        all_samples = list(vcf.header.samples)

        # Assign each sample to the correct group based on prefix
        for sample in all_samples:
            for group_prefix in sample_groups:
                if sample.startswith(group_prefix):
                    sample_groups[group_prefix].append(sample)
                    sample_to_group[sample] = group_prefix
                    break

        # Initialize a private-variant counter for each group
        private_variant_counts = {group: 0 for group in sample_groups}

        # Iterate over each variant in the VCF
        for record in vcf.fetch():
            groups_with_alt = set()

            # Check which groups carry an alternate allele
            for sample in all_samples:
                gt = record.samples[sample]["GT"]
                if gt is None:
                    continue
                # If there's a '1' in the genotype, that means an alt allele is present
                if 1 in gt:
                    grp = sample_to_group[sample]
                    groups_with_alt.add(grp)

            # If exactly one group has alt alleles, it is 'private' to that group
            if len(groups_with_alt) == 1:
                only_group = next(iter(groups_with_alt))
                private_variant_counts[only_group] += 1

    return private_variant_counts


def write_private_variant_counts(private_variant_counts, output_path):
    """
    Write the private variant counts per group to a CSV file.
    """
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Group", "PrivateVariantCount"])
        for group, count in private_variant_counts.items():
            writer.writerow([group, count])

    print(f"Private variant counts written to {output_path}")


if __name__ == "__main__":
    # Process the VCF to find private variants
    private_variant_counts = process_vcf_for_private_variants(vcf_file, sample_groups)

    # Write the results to a CSV
    write_private_variant_counts(private_variant_counts, output_file)

