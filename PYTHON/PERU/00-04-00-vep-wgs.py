# Assigns VEP annotation to biallelic SNPs


import mysql.connector
import configparser
import sys
import os
import glob
from pathlib import Path 
# Reading database configuration
config = configparser.ConfigParser()
config.read_file(open(r'CONF/mariadb.conf'))

# Connecting to the database with explicit charset and collation
try:
    db = mysql.connector.connect(
        host=config.get('peru', 'host'),
        user=config.get('peru', 'user'),
        password=config.get('peru', 'password'),
        database=config.get('peru', 'database'),
        charset='utf8mb4',
        collation='utf8mb4_general_ci'
    )
except mysql.connector.Error as e:
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)

db.autocommit = True
cursor = db.cursor()

# Explicitly set collation for the connection
try:
    cursor.execute("SET collation_connection = 'utf8mb4_general_ci'")
except mysql.connector.Error as e:
    print(f"Error setting collation: {e}")
    sys.exit(1)

# Function to process multiallelic sites
def process_multiallelic_site(line):
    parts = line.split('\t')
    chromosome, position, _, ref, alts = parts[:5]
    genotypes = parts[9:]

    biallelic_lines = []

    for i, alt in enumerate(alts.split(','), start=1):
        new_genotypes = []
        for gt in genotypes:
            if gt == '0/0':
                new_gt = '0/0'
            elif gt == f'0/{i}':
                new_gt = '0/1'
            elif gt == f'{i}/{i}':
                new_gt = '1/1'
            else:
                new_gt = './.'  # or handle as appropriate
            new_genotypes.append(new_gt)

        new_line = '\t'.join([chromosome, position, '.', ref, alt] + parts[5:9] + new_genotypes)
        biallelic_lines.append(new_line)

    return biallelic_lines


db.autocommit = True
cursor = db.cursor(dictionary=True)
sql1="SELECT Location,Allele,Consequence,IMPACT,SYMBOL,BIOTYPE,Existing_variation,SIFT,PolyPhen,AF,CLIN_SIG,CADD_PHRED,LOEUF " \
    +"FROM 00_02_VEP_WGS WHERE Location = '{0}' AND Allele = '{1}' "
insert_statement = ("INSERT INTO `peru`.`00_04_VEP_WGS` (`Chromosome`, `Chr_position`, `REF`, `ALT`, `Consequence`, `IMPACT`, `SYMBOL`, `BIOTYPE`, `Existing_variation`, `SIFT`, `PolyPhen`, `AF`, `CLIN_SIG`, `CADD_PHRED`, `LOEUF`, `ZYG`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)")

with open('ANALYSIS/VEP/TXT/failed_vcf_lines_v2.txt', 'w') as failed_file:
    failed_file.write('')  # This will clear the file or create it if it doesn't exist


file = "INPUT/VCF/Peru.joint.vcf"
count = 0

severity_ranking = {
    # High Impact
    'transcript_ablation': 1,
    'splice_acceptor_variant': 2,
    'splice_donor_variant': 3,
    'stop_gained': 4,
    'frameshift_variant': 5,
    'stop_lost': 6,
    'start_lost': 7,
    'transcript_amplification': 8,
    'inframe_insertion': 9,
    'inframe_deletion': 10,
    'missense_variant': 11,
    'protein_altering_variant': 12,
    'splice_region_variant': 13,
    'incomplete_terminal_codon_variant': 14,
    'start_retained_variant': 15,
    'stop_retained_variant': 16,
    'synonymous_variant': 17,
    'coding_sequence_variant': 18,

    # Moderate Impact
    'regulatory_region_ablation': 19,
    'regulatory_region_amplification': 20,
    'TF_binding_site_ablation': 21,
    'TF_binding_site_amplification': 22,
    'regulatory_region_variant': 23,
    'TF_binding_site_variant': 24,
    'mature_miRNA_variant': 25,
    '5_prime_UTR_variant': 26,
    '3_prime_UTR_variant': 27,

    # Low Impact
    'non_coding_transcript_exon_variant': 28,
    'intron_variant': 29,
    'NMD_transcript_variant': 30,
    'non_coding_transcript_variant': 31,
    'upstream_gene_variant': 32,
    'downstream_gene_variant': 33,
    'TFBS_ablation': 34,
    'TFBS_amplification': 35,
    'TF_binding_site_variant': 36,
    'regulatory_region_ablation': 37,
    'regulatory_region_amplification': 38,
    'feature_elongation': 39,
    'regulatory_region_variant': 40,
    'feature_truncation': 41,
    'intergenic_variant': 42,

    # Modifier Impact
    'incomplete_terminal_codon_variant': 43,
    'start_retained_variant': 44,
    'stop_retained_variant': 45,
    'synonymous_variant': 46,
    'coding_sequence_variant': 47,
    'mature_miRNA_variant': 48,
    '5_prime_UTR_variant': 49,
    '3_prime_UTR_variant': 50,
    'non_coding_transcript_exon_variant': 51,
    'intron_variant': 52,
    'NMD_transcript_variant': 53,
    'non_coding_transcript_variant': 54,
    'upstream_gene_variant': 55,
    'downstream_gene_variant': 56,
    'TFBS_ablation': 57,
    'TFBS_amplification': 58,
    'TF_binding_site_variant': 59,
    'regulatory_region_ablation': 60,
    'regulatory_region_amplification': 61,
    'feature_elongation': 62,
    'regulatory_region_variant': 63,
    'feature_truncation': 64,
    'intergenic_variant': 65
}


import traceback

def update_or_insert_record(data, failed_file_path='ANALYSIS/VEP/TXT/failed_vcf_lines_v2.txt'):
    Chromosome, Chr_position, REF, ALT, consequence, *_ = data
    new_consequence_score = severity_ranking.get(consequence, float('inf'))

    try:
        # Use 'ALT' instead of 'Allele' based on your table schema
        cursor.execute("SELECT Consequence FROM `00_04_VEP_WGS` WHERE Chromosome = %s AND Chr_position = %s AND REF = %s AND ALT = %s", (Chromosome, Chr_position, REF, ALT))
        existing_record = cursor.fetchone()

        if existing_record:
            existing_consequence_score = severity_ranking.get(existing_record['Consequence'], float('inf'))
            if new_consequence_score < existing_consequence_score:
                # Update the existing record; ensure the SQL matches your table schema
                cursor.execute("UPDATE `00_04_VEP_WGS` SET Consequence=%s, IMPACT=%s, SYMBOL=%s, BIOTYPE=%s, Existing_variation=%s, SIFT=%s, PolyPhen=%s, AF=%s, CLIN_SIG=%s, CADD_PHRED=%s, LOEUF=%s, ZYG=%s WHERE Chromosome = %s AND Chr_position = %s AND REF = %s AND ALT = %s", 
                               (consequence, data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], Chromosome, Chr_position, REF, ALT))
        else:
            # Insert new record
            cursor.execute(insert_statement, data)
    except Exception as e:
        # Log the detailed error and traceback
        error_details = traceback.format_exc()
        print(f"Error: {e}\nDetails: {error_details}")
        with open(failed_file_path, 'a') as failed_file:
            failed_file.write(f"Failed to process record for Chromosome: {Chromosome}, Position: {Chr_position}, REF: {REF}, ALT: {ALT}, Error: {e}\nDetails: {error_details}\n")

def process_line(parts):
    Chromosome   = parts[0]
    Chr_position = parts[1]
    REF          = parts[3]
    ALT          = parts[4]
    zyg          = ':'.join(parts[9:])
    Location     = Chromosome + ':' + Chr_position + '-' + Chr_position
    Allele       = ALT
    query = sql1.format(Location, Allele)
    cursor.execute(query)
    results = cursor.fetchall()

    for result in results:
        try:
            consequence = result['Consequence']
            consequence_severities = [(severity_ranking.get(c.strip(), float('inf')), c.strip()) for c in consequence.split(',')]
            top_consequence = min(consequence_severities, key=lambda x: x[0])[1]  # Select the top-ranked consequence
            result['Consequence'] = top_consequence  # Update result with the top-ranked consequence

            # Prepare the data for insertion or update
            data = (Chromosome, Chr_position, REF, ALT, top_consequence, result['IMPACT'], result['SYMBOL'], result['BIOTYPE'], result['Existing_variation'], result['SIFT'], result['PolyPhen'], result['AF'], result['CLIN_SIG'], result['CADD_PHRED'], result['LOEUF'], zyg)

            # Use the new function to decide whether to insert or update
            update_or_insert_record(data)
        except KeyError as e:
            # Handle cases where a KeyError occurs
            print(f"Error processing line: {e}")

with open(file, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue  # Skip header lines
        parts = line.strip().split('\t')
        if len(parts) < 10: 
            continue  # Skip incomplete lines

        # Determine if the line is multiallelic and split if necessary
        if ',' in parts[4]:  # Check if ALT field is multiallelic
            biallelic_lines = process_multiallelic_site(line)
            for biallelic_line in biallelic_lines:
                # Process each biallelic line as usual
                b_parts = biallelic_line.strip().split('\t')
                process_line(b_parts)
        else:
            # Process the line as usual
            process_line(parts)
            
# Close the database connection
db.close()

