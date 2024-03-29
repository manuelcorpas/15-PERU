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
sql1 = "SELECT Chromosome,Chr_position,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG FROM 00_04_VEP_WGS_HIGH_IMPACT WHERE CLIN_SIG NOT like '%benign%'"
sql2 = "INSERT INTO 00_05_ALT_ALLELE_POP_COUNT(LOCATION,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG,CHOPCCAS,CUSCO,IQUITOS,MATZES,MOCHES,TRUJILLO,UROS) VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}')"
cursor.execute(sql1)

results = cursor.fetchall()
regions = {}
for result in results:
    Chromosome = result['Chromosome']
    Chr_Start  = result['Chr_Start']
    Chr_End    = result['Chr_End']
    SYMBOL      = result['SYMBOL']
    CADD_PHRED  = result['CADD_PHRED']
    Consequence = result['Consequence']
    CLIN_SIG    = result['CLIN_SIG']
    region     = str(Chromosome) + ':' + str(Chr_Start) + "-" + str(Chr_End)
    #print(region,sep='\t')
    variant = count_alt_alleles(vcf_path, region)
    if not region in regions:
        regions[region] = 0
        #cursor.execute(sql2.format(*[region,SYMBOL,CADD_PHRED,Consequence,CLIN_SIG,variant['CHOPCCAS'],variant['CUSCO'],variant['IQUITOS'],variant['MATZES'],variant['MOCHES'],variant['TRUJILLO'],variant['UROS']]))
        print (region,variant)
    regions[region] +=1

for loc in regions:
    if regions[loc] > 1:    
        print('Location with more than 1 VEP consequences: ',loc,regions[loc],sep='\t')    
