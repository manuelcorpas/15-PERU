# The script reads the VCF file line by line.
# It identifies the sample names based on the #CHROM line.
# For each variant line, it parses the genotype (GT) field from each sample column.
# It identifies and counts alternative alleles in the genotype, ignoring reference alleles.
# The script then prints out the count of alternative alleles for each sample and variant 
# position where at least one alternative allele is present

import re

def count_alternative_alleles(vcf_filename):
    with open(vcf_filename, 'r') as vcf_file:
        sample_names = []
        for line in vcf_file:
            line = line.strip()

            # Skip comment lines
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    columns = line.split("\t")
                    sample_names = columns[9:]
                    print(f"Found samples: {', '.join(sample_names)}")
                continue

            columns = line.split("\t")
            info_field = columns[7]
            sample_fields = columns[9:]

            # Parse each sample field to count alternative alleles
            for i, sample_field in enumerate(sample_fields):
                genotype_field = sample_field.split(":")[0]
                alleles = re.findall(r'[0-9]+', genotype_field)

                # Count the number of alternative alleles (ignoring the reference allele '0')
                alt_allele_count = sum(1 for allele in alleles if allele != '0')
                
                if alt_allele_count > 0:
                    print(f"Sample '{sample_names[i]}' has {alt_allele_count} alternative alleles at position {columns[1]}.")

# Replace 'your_file.vcf' with the name of your input VCF file
count_alternative_alleles('your_file.vcf')

