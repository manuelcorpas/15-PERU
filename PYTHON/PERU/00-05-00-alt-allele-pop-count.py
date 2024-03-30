# The script reads the VCF file line by line.
# It identifies the sample names based on the #CHROM line.
# For each variant line, it parses the genotype (GT) field from each sample column.
# It identifies and counts alternative alleles in the genotype, ignoring reference alleles.
# The script then prints out the count of alternative alleles for each sample and variant 
# position where at least one alternative allele is present
# The script counts alternative alleles for those coordinates that are queried from the database


import pymysql  # Replace 'import mariadb' with this
import glob, os
from pathlib import Path
import configparser
import re
import sys
import pysam

config = configparser.ConfigParser()
config.read_file(open(r'CONF/mariadb.conf'))

try:
    db = pymysql.connect(
        host=config.get('peru', 'host'),
        user=config.get('peru', 'user'),
        password=config.get('peru', 'password'),  # 'passwd' is changed to 'password'
        database=config.get('peru', 'database')  # 'db' is changed to 'database'
    )

except pymysql.Error as e:  # 'mariadb.Error' is changed to 'pymysql.Error'
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)

db.autocommit(True)  # 'db.autocommit = True' is changed to 'db.autocommit(True)'
cursor = db.cursor(pymysql.cursors.DictCursor)  # 'dictionary=True' is changed to 'pymysql.cursors.DictCursor'

def count_alt_alleles(vcf_path, region):
    """
    Counts alternative alleles per sample in a specified region of a VCF file using tabix.

    Parameters:
    - vcf_path: str, path to the VCF file
    - region: str, specific region in format 'chr:start-end'
    """
    try:
        # Initialize a dictionary to hold the counts of alternative alleles per sample
        alt_allele_counts = {}

        # Open the VCF file with Tabix
        vcf = pysam.TabixFile(vcf_path)

        # Extract the sample names from the VCF file's header
        header = vcf.header
        samples = [line for line in header if line.startswith('#CHROM')][0].strip().split('\t')[9:]

        # Initialize counts for each sample
        for sample in samples:
            alt_allele_counts[sample] = 0

        # Fetch the specified region
        for record in vcf.fetch(region=region):
            parts = record.strip().split('\t')

            # Extract the genotype information for each sample
            genotypes = parts[9:]

            # Iterate through the genotypes and count alternative alleles
            for i, genotype in enumerate(genotypes):
                # Extract the genotype information
                genotype_parts = re.split(r'[/|]', genotype.split(':')[0])

                # Count the alternative alleles
                alt_allele_count = genotype_parts.count('1')

                # Update the counts dictionary
                alt_allele_counts[samples[i]] += alt_allele_count

        vcf.close()

        # Output the results
        results = {}
        for sample, count in alt_allele_counts.items():
            parts = sample.split('-')
            Population = parts[0]
            Sample_no = parts[1]
            if Population not in results:
                results[Population] = 0
            results[Population] += count

        # Return the aggregated counts by population
        return results

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)


# Peru.joint150WG.vcf.gz
vcf_path = "INPUT/VCF/Peru.joint150WG.vcf.gz"

# Unknown clinical significance of high impact in protein coding regions:
sql1 = "SELECT Chromosome,Chr_position,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG FROM 00_04_VEP_WGS_HIGH_IMPACT WHERE CLIN_SIG NOT like '%benign%'"
sql2 = "INSERT INTO 00_05_ALT_ALLELE_POP_COUNT(LOCATION,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG,CHOPCCAS,CUSCO,IQUITOS,MATZES,MOCHES,TRUJILLO,UROS) VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}')"
cursor.execute(sql1)

results = cursor.fetchall()
regions = {}
for result in results:
    Chromosome = result['Chromosome']
    Chr_Start  = result['Chr_position']
    Chr_End    = result['Chr_position']
    SYMBOL      = result['SYMBOL']
    CADD_PHRED  = result['CADD_PHRED']
    Consequence = result['Consequence']
    CLIN_SIG    = result['CLIN_SIG']
    region     = str(Chromosome) + ':' + str(Chr_Start) + "-" + str(Chr_End)
    print(region,sep='\t')
    variant = count_alt_alleles(vcf_path, region)
    if not region in regions:
        regions[region] = 0
        cursor.execute(sql2.format(*[region,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG,variant['CHOPCCAS'],variant['CUSCO'],variant['IQUITOS'],variant['MATZES'],variant['MOCHES'],variant['TRUJILLO'],variant['UROS']]))
        print (region,variant)
    regions[region] +=1

for loc in regions:
    if regions[loc] > 1:    
        print('Location with more than 1 VEP consequences: ',loc,regions[loc],sep='\t')    
