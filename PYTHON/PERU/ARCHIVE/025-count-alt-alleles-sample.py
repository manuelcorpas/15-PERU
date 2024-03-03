import subprocess

def count_alternative_alleles(vcf_gz_file):
    # Index the VCF file (if not already indexed)
    #subprocess.run(['bcftools', 'index', vcf_gz_file], check=True)

    # Use bcftools to extract genotype information
    cmd = ['bcftools', 'query', '-f', '%SAMPLE\t%GT\n', vcf_gz_file]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error executing bcftools:")
        print(result.stderr)
        exit(1)

    alt_allele_counts = {}

    for line in result.stdout.strip().split('\n'):
        sample, genotype = line.split('\t')
        if sample not in alt_allele_counts:
            alt_allele_counts[sample] = 0
        # Count non-reference and non-missing alleles
        alt_allele_counts[sample] += sum(allele not in ['0', '.'] for allele in genotype.replace('|', '/').split('/'))

    return alt_allele_counts

# Example usage
vcf_file = 'INPUT/VCF/Peru.joint150WG.vcf.gz'
alternative_allele_counts = count_alternative_alleles(vcf_file)
for sample, count in alternative_allele_counts.items():
    print(f"{sample}: {count}")

