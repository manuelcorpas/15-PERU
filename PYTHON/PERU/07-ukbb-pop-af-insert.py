import os
import zipfile
import pandas as pd
import csv
import mariadb
import os
import configparser
import glob

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

cursor = db.cursor(dictionary=True)


sql1 = "INSERT INTO 07_UKBB_POP_AF_INSERT(GENE,POPULATION,ALLELE,ALLELES_OBSERVED,ALLELES_TOTAL,FREQUENCY) VALUES('{0}','{1}','{2}','{3}','{4}','{5}')"

def extract_and_print_data(root_directory):

    pattern = root_directory + '*_UKBB_frequencies.tsv'
    tsv_files = glob.glob(pattern)
    
    for file in tsv_files:
        gene_name = os.path.basename(file).split('_UKBB_frequencies.tsv')[0]
        df = pd.read_csv(file, sep='\t')
        for index, row in df.iterrows():
            print(gene_name,row.iloc[0],row.iloc[1],row.iloc[2],row.iloc[3],row.iloc[4],sep='\t')
            input_data = [gene_name,row.iloc[0],row.iloc[1].replace("'", ""),row.iloc[2],row.iloc[3],row.iloc[4]]
            cursor.execute(sql1.format(*input_data))

    db.commit()

# Example usage
extract_and_print_data('INPUT/PHARMGKB/')
#print(df)

# Replace '/path/to/search' with the actual root directory path where the search should start.


#process_directories('ANALYSIS/06-PyPGX/')



