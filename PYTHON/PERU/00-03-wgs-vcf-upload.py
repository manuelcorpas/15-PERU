import mariadb
import glob
import os
from pathlib import Path

def parse_vcf_file(file_path, db):
    cursor = db.cursor()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue  # Skip incomplete lines

            chromosome = parts[0]
            position = parts[1]
            ref = parts[3]
            alt = parts[4]
            zyg = ':'.join(parts[9:])

            # Insert into database
            try:
                cursor.execute("INSERT INTO 00_03_WGS_VCF_UPLOAD (Chromosome, Chr_position, REF, ALT, ZYG) VALUES (?, ?, ?, ?, ?)", 
                               (chromosome, position, ref, alt, zyg))
            except mariadb.Error as e:
                print(f"Error inserting data: {e}")
            db.commit()

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

file = "ANALYSIS/02-PGEN-PSAM-PVAR-2-VCF/Peru.joint.vcf"

parse_vcf_file(file, db)

# Close the database connection
db.close()

