import csv
import mariadb
import glob
import os
from pathlib import Path
import configparser
import sys

# Reading database configuration
config = configparser.ConfigParser()
config.read_file(open(r'CONF/mariadb.conf'))

# Connecting to the database
try:
    db = mariadb.connect(
        host=config.get('peru', 'host'),
        user=config.get('peru', 'user'),
        passwd=config.get('peru', 'password'),
        db=config.get('peru', 'database'))
except mariadb.Error as e:
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)

db.autocommit = True
cursor = db.cursor()

def truncate_value(value, max_length):
    """Truncate the value to a specified max length."""
    return value[:max_length]

def process_file(file_path, cursor):
    with open(file_path) as f:
        reader = csv.reader(f, delimiter="\t")
        headers = next(reader, None)  # Read the file headers

        headers = [header.lstrip('#') for header in headers if header != 'PUBMED']

        sql_columns = ', '.join(headers)
        sql_insert_base = f"INSERT INTO 02_VEP_UPLOAD ({sql_columns}) VALUES "

        for row_index, row in enumerate(reader, start=1):
            # Exclude the 'PUBMED' value from each row
            row = [value for header, value in zip(headers, row) if header != 'PUBMED']

            # Truncate MOTIF_NAME if present
            if 'MOTIF_NAME' in headers:
                motif_index = headers.index('MOTIF_NAME')
                row[motif_index] = truncate_value(row[motif_index], 255)  # Assuming 255 is the max length

            if len(row) != len(headers):
                print(f"Row {row_index} skipped: column count does not match value count.")
                print(f"Problematic row content: {row}")
                continue

            sql_values = ', '.join(["'{}'".format(value.replace("'", "''")) for value in row])
            sql_query = sql_insert_base + f"({sql_values});"

            try:
                cursor.execute(sql_query)
            except mariadb.Error as e:
                print(f"Error inserting data at row {row_index}: {e}")

os.chdir("ANALYSIS/VEP/TXT/")

for file in glob.glob("*.txt"):
    print(Path(file).stem)
    process_file(file, cursor)

db.close()



'''


sql1 = "INSERT INTO 02_UP_VEP_EX_COMMON_VAR_1_CONSEQ(Uploaded_variation,Chromosome,Chr_Start,Chr_End,Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,SIFT,PolyPhen,AF,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,CLIN_SIG,SOMATIC,PHENO,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,CADD_PHRED,CADD_RAW,LOEUF) VALUES('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}','{8}','{9}','{10}','{11}','{12}','{13}','{14}','{15}','{16}','{17}','{18}','{19}','{20}','{21}','{22}','{23}','{24}','{25}','{26}','{27}','{28}','{29}','{30}','{31}','{32}','{33}','{34}','{35}','{36}','{37}','{38}','{39}','{40}','{41}','{42}','{43}','{44}','{45}','{46}','{47}','{48}','{49}','{50}','{51}')"


os.chdir("ANALYSIS/VEP/TXT/")

for file in glob.glob("*.txt"):
    print (Path(file).stem)
    with open(file) as f:
        reader = csv.reader(f,delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) == 51:
                Location = row[1]
                arr1 = Location.split(':')
                arr2 = arr1[1].split('-')
                Chromosome = arr1[0]
                Chr_Start  = arr2[0]
                Chr_End    = arr2[1]
                new_row    = [row[0],Chromosome,Chr_Start,Chr_End]+row[2:]
                cursor.execute(sql1.format(*new_row))
'''
