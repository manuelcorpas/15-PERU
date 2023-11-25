import pandas as pd
from sqlalchemy import create_engine
import os

# Database credentials
username = 'admin'
password = 'root'
host = 'localhost'
database = 'peru2'
port = 3306  # default MariaDB port

# SQLAlchemy engine for MariaDB
engine = create_engine(f'mysql+pymysql://{username}:{password}@{host}:{port}/{database}')

# Directory containing the files
directory = 'INPUT/VEP/WGS_DATA/'

# Define the column names (modify these to match your SQL table's columns)
columns = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "SIFT",
    "PolyPhen",
    "AF",
    "gnomADg_AF",
    "gnomADg_AFR_AF",
    "gnomADg_AMI_AF",
    "gnomADg_AMR_AF",
    "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF",
    "gnomADg_FIN_AF",
    "gnomADg_MID_AF",
    "gnomADg_NFE_AF",
    "gnomADg_OTH_AF",
    "gnomADg_SAS_AF",
    "FREQS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "CADD_PHRED",
    "CADD_RAW",
    "LOEUF"
]

# Iterate over files and upload data
for filename in os.listdir(directory):
    if filename.endswith(".txt"):  # Assuming all files are .txt format
        file_path = os.path.join(directory, filename)
        # Read the file without headers and assign column names
        data = pd.read_csv(file_path, sep='\t', header=None, names=columns, skiprows=1)
        # Upload data to the database
        data.to_sql('00_02_VEP_WGS', con=engine, if_exists='append', index=False)
        print(f"Uploaded {filename}")






'''
# Directory containing the files
directory = 'INPUT/VEP/ARRAY_DATA/'

# Iterate over files and upload data
for filename in os.listdir(directory):
    if filename.endswith(".txt"):  # Assuming all files are .txt format
        file_path = os.path.join(directory, filename)
        data = pd.read_csv(file_path, sep='\t', header=None, skiprows=1)  # Change sep if needed
        data.to_sql('00_01_VEP_ARRAY', con=engine, if_exists='append', index=False)  # Appends data to the table
        print(f"Uploaded {filename}")
'''
