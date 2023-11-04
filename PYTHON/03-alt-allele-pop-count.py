# The script reads the VCF file line by line.
# It identifies the sample names based on the #CHROM line.
# For each variant line, it parses the genotype (GT) field from each sample column.
# It identifies and counts alternative alleles in the genotype, ignoring reference alleles.
# The script then prints out the count of alternative alleles for each sample and variant 
# position where at least one alternative allele is present
# The script counts alternative alleles for those coordinates that are queried from the database


import mariadb
import glob, os
from pathlib import Path
import configparser
import re
import pysam
import sys

config   = configparser.ConfigParser()
config.read_file(open(r'CONF/mariadb.conf'))

try:
    db = mariadb.connect(
        host   = config.get('peru','host'),
        user   = config.get('peru','user'     ),
        passwd = config.get('peru','password' ),
        db     = config.get('peru','database'))

except mariadb.Error as e:
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)


db.autocommit = True
cursor = db.cursor(dictionary=True)

def count_alt_alleles(vcf_path, region):
    """
    Counts alternative alleles per sample in a specified region of a VCF file.
    
    Parameters:
    - vcf_path: str, path to the VCF file
    - region: str, specific region in format 'chr:start-end'
    """
    try:
        # Open the VCF file
        vcf_file = pysam.VariantFile(vcf_path)
        
        # Extract the samples
        samples = vcf_file.header.samples
        
        # Dictionary to hold the counts of alternative alleles per sample
        alt_allele_counts = {sample: 0 for sample in samples}
        
        # Fetch the variants in the specified region
        for variant in vcf_file.fetch(region=region):
            for sample, sample_data in variant.samples.items():
                # Count alternative alleles
                if sample_data.alleles:
                    alt_alleles = [allele for allele in sample_data.alleles if allele != variant.ref]
                    alt_allele_counts[sample] += len(alt_alleles)
        
        results = {}
 
        # Output the results
        for sample, count in alt_allele_counts.items():
            parts = sample.split('-')
            Population = parts[0]
            Sample_no  = parts[1]
            if not Population in results:
                results[Population] = 0
            results[Population] += count
            #print(region,end='\t')       
            #print(f"{Population}\t{Sample_no}: {count}")
            
        #for pop in results:
        #    print(f"{region}\t{pop}\t{results[pop]}")
        return results;
        # Close the VCF file
        vcf_file.close()
        
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

# Peru.joint150WG.vcf.gz
vcf_path = "INPUT/VCF/Peru.joint150WG.vcf.gz"

# Unknown clinical significance of high impact in protein coding regions:
sql1 = "SELECT Chromosome,Chr_Start,Chr_End FROM 02_UP_VEP_EX_COMMON_VAR_1_CONSEQ WHERE IMPACT ='HIGH' AND CLIN_SIG = '-' AND BIOTYPE ='protein_coding'"

cursor.execute(sql1)

results = cursor.fetchall()
for result in results:
    Chromosome = result['Chromosome']
    Chr_Start  = result['Chr_Start']
    Chr_End    = result['Chr_End']
    region     = str(Chromosome) + ':' + str(Chr_Start) + "-" + str(Chr_End)
    #print(region,sep='\t')
    variant = count_alt_alleles(vcf_path, region)
    print (region,variant)
    
    
'''
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
count_alternative_alleles('INPUT/VCF/split_34.vcf')
'''
