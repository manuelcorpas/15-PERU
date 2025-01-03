import pysam
from collections import defaultdict
import matplotlib.pyplot as plt

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

def process_vcf(vcf_file, sample_groups):
    """
    Process the VCF file to compute total and normalized alternative allele counts.
    """
    with pysam.VariantFile(vcf_file, "r") as vcf:
        # Convert header samples to a list for easier indexing
        samples = list(vcf.header.samples)

        # Map sample names to groups
        for sample in samples:
            for group in sample_groups:
                if sample.startswith(group):
                    sample_groups[group].append(sample)

        # Precompute group indices
        group_indices = {
            group: [samples.index(sample) for sample in sample_groups[group]]
            for group in sample_groups
        }

        # Initialize counters for alternative alleles
        alt_allele_counts = {group: 0 for group in sample_groups}

        # Process each record in the VCF
        for record in vcf.fetch():
            for group, indices in group_indices.items():
                for idx in indices:
                    genotype = record.samples[idx]["GT"]
                    if genotype:
                        # Count occurrences of the alternative allele (1)
                        alt_allele_counts[group] += genotype.count(1)

    # Calculate normalized counts (per sample in the group)
    normalized_counts = {
        group: alt_allele_counts[group] / len(sample_groups[group]) if sample_groups[group] else 0
        for group in sample_groups
    }

    return alt_allele_counts, normalized_counts

def plot_results(normalized_counts):
    """
    Plot the normalized alternative allele counts.
    """
    groups = list(normalized_counts.keys())
    normalized_values = list(normalized_counts.values())

    # Create a bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(groups, normalized_values, color='skyblue')
    plt.xlabel("Groups")
    plt.ylabel("Average Alternative Alleles per Individual")
    plt.title("Comparison of Alternative Alleles by Group")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Process the VCF and compute allele counts
    alt_allele_counts, normalized_counts = process_vcf(vcf_file, sample_groups)

    # Print the results
    print("Total and Normalized Alternative Alleles by Group:")
    for group in sample_groups:
        total = alt_allele_counts[group]
        normalized = normalized_counts[group]
        print(f"{group}: Total = {total}, Normalized = {normalized:.2f}")

    # Plot the normalized counts
    plot_results(normalized_counts)

