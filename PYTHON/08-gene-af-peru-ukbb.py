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

sql1 = "SELECT DISTINCT GENE FROM 07_UKBB_POP_AF_INSERT ORDER BY GENE ASC"
sql2 = "SELECT POPULATION,ALLELE,COUNTER FROM 06_SUMMARISE_AF_BY_POP WHERE GENE = '{0}'"
sql3 = "SELECT p.FROM 07_UKBB_POP_AF_INSERT as uk, 06_SUMMARISE_AF_BY_POP as peru WHERE ALLELE LIKE '%Reference%';"

cursor.execute(sql1)
genes = cursor.fetchall()

for row in genes:
    GENE = row['GENE']
    cursor.execute(sql2.format(GENE))
    populations = cursor.fetchall()
    for pop in populations:
        POPULATION = pop['POPULATION']
        ALLELE     = pop['ALLELE']
        COUNTER    = pop['COUNTER']
        print(GENE,POPULATION,ALLELE,COUNTER,sep='\t')

db.commit()




