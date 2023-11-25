import os
import zipfile
import pandas as pd
import csv
import mariadb
import os
import configparser

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

sql1 = "INSERT INTO 05_PGX_INSERT(GENE,SAMPLE,ALLELE,EFFECT) VALUES('{0}','{1}','{2}','{3}')"


def extract_and_print_data(root_directory):
    # Loop through the directories and subdirectories
    for root, dirs, files in os.walk(root_directory):
        # Extract the part of the directory name after 'pipeline_'
        dir_name = os.path.basename(root)
        if dir_name.startswith('pipeline_'):
            pipeline_id = dir_name[len('pipeline_'):]

            # If 'results.zip' exists in this directory
            if 'results.zip' in files:
                # Construct the full path to the zip file
                zip_path = os.path.join(root, 'results.zip')

                # Unzip the file
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    # Extract to the same directory
                    zip_ref.extractall(root)

                    # Assuming the first item in the zip is the directory
                    unzipped_dir = zip_ref.namelist()[0].split('/')[0]
                    unzipped_path = os.path.join(root, unzipped_dir)

                    # Check for 'data.tsv' in the unzipped directory
                    if 'data.tsv' in os.listdir(unzipped_path):
                        data_file_path = os.path.join(unzipped_path, 'data.tsv')

                        # Print the pipeline ID
                        #print(f"Pipeline ID: {pipeline_id}")

                        # Read and print the contents of 'data.tsv'
                        df = pd.read_csv(data_file_path,sep='\t')
                        for index, row in df.iterrows():
                            print(pipeline_id,row.iloc[0],row.iloc[1],row.iloc[2], sep='\t')
                            input_data = [pipeline_id,row.iloc[0],row.iloc[1],row.iloc[2]]
                            cursor.execute(sql1.format(*input_data))
                            db.commit()
                        #df.to_excel(root_directory, index=False)

# Example usage
extract_and_print_data('ANALYSIS/06-PyPGX/')
#print(df)

# Replace '/path/to/search' with the actual root directory path where the search should start.


#process_directories('ANALYSIS/06-PyPGX/')



