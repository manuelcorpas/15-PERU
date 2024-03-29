import sys
import mariadb
import glob
import os
from pathlib import Path


# Function to parse VCF file and extract sample IDs
def parse_vcf_sample_ids(vcf_path):
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("#CHROM"):
                parts = line.strip().split('\t')
                return parts[9:]  # Sample IDs start from the 10th column in VCF header
    return []

# Function to analyze zygosity from ZYG string, given sample IDs
def analyze_zygosity(zyg_str, sample_ids):
    zyg_values = zyg_str.split(':')
    allele_frequencies = {sample_id: 0 for sample_id in sample_ids}  # Initialize frequencies
    
    for sample_id, zyg in zip(sample_ids, zyg_values):
        # Simplified analysis: count '1' alleles (homozygous or heterozygous variants)
        allele_frequencies[sample_id] += zyg.count('1')
    
    return allele_frequencies

# Main script
if __name__ == "__main__":
    vcf_path = "INPUT/VCF/Peru.joint.vcf"
    sample_ids = parse_vcf_sample_ids(vcf_path)


    try:
        # Connect to MariaDB
        db = mariadb.connect(
            host='localhost',
            user='admin',
            password='root',  # Your actual password
            database='peru'
        )
    except mariadb.Error as e:
        print(f"Error connecting to MariaDB: {e}")
        sys.exit(1)

    cur = db.cursor()

    # Fetch ZYG column for high impact variants
    cur.execute("SELECT ZYG FROM 00_04_VEP_WGS_HIGH_IMPACT")

    # Initialize a dictionary to hold overall allele frequencies across all variants
    overall_allele_frequencies = {sample_id: 0 for sample_id in sample_ids}

    for (zyg,) in cur:
        # Analyze zygosity for the current variant
        allele_frequencies = analyze_zygosity(zyg, sample_ids)
        
        # Update overall allele frequencies
        for sample_id, freq in allele_frequencies.items():
            overall_allele_frequencies[sample_id] += freq

    # Close database connection
    db.close()

    # Print overall allele frequencies for each sample
    for sample_id, freq in overall_allele_frequencies.items():
        print(f"Sample ID: {sample_id}, Allele Frequency: {freq}")

