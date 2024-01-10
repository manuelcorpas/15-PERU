import sys
import mariadb
import glob
import os
from pathlib import Path

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

# Establish a connection to the MariaDB database
try:
    db = mariadb.connect(
        host='localhost',
        user='admin',
        password='root',  # Fill in your password here
        database='peru')
except mariadb.Error as e:
    print(f"Error connecting to MariaDB: {e}")
    sys.exit(1)

db.autocommit = True
cursor = db.cursor(dictionary=True)
sql1="SELECT Location,Allele,Consequence,IMPACT,SYMBOL,BIOTYPE,Existing_variation,SIFT,PolyPhen,AF,CLIN_SIG,CADD_PHRED,LOEUF " \
    +"FROM 02_VEP_UPLOAD WHERE Location = '{0}' AND Allele = '{1}' "
insert_statement = ("INSERT INTO `peru`.`03_00_VEP_WGS` (`Chromosome`, `Chr_position`, `REF`, `ALT`, `Consequence`, `IMPACT`, `SYMBOL`, `BIOTYPE`, `Existing_variation`, `SIFT`, `PolyPhen`, `AF`, `CLIN_SIG`, `CADD_PHRED`, `LOEUF`, `ZYG`) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)")


file = "ANALYSIS/02-PGEN-PSAM-PVAR-2-VCF/Peru.joint.vcf"
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

    # Process results and insert into the database
    top_result = None
    min_severity_rank = float('inf')
    for result in results:
        consequence = result['Consequence']
        consequence_severities = [severity_ranking.get(c.strip(), float('inf')) for c in consequence.split(',')]
        highest_severity = min(consequence_severities)
        if highest_severity < min_severity_rank:
            min_severity_rank = highest_severity
            top_result = result

    if top_result:
        location = top_result['Location']
        allele = top_result['Allele']
        consequence = top_result['Consequence']
        impact = top_result['IMPACT']
        symbol = top_result['SYMBOL']
        biotype = top_result['BIOTYPE']
        existing_variation = top_result['Existing_variation']
        sift = top_result['SIFT']
        polyphen = top_result['PolyPhen']
        af = top_result['AF']
        clin_sig = top_result['CLIN_SIG']
        cadd_phred = top_result['CADD_PHRED']
        loeuf = top_result['LOEUF']

        data = (Chromosome, Chr_position, REF, ALT, consequence, impact, symbol, biotype, existing_variation, sift, polyphen, af, clin_sig, cadd_phred, loeuf, zyg)
        cursor.execute(insert_statement, data)
        print(data)
    else:
        with open('ANALYSIS/VEP/TXT/failed_vcf_lines_v2.txt', 'a') as failed_file:
            failed_file.write('\t'.join(parts) + '\n')

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

