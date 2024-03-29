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

sql1 = "INSERT INTO 00_00_SAMPLE(ARRAY_CODE,INTERNAL_CODE,POPULATION,DEPARTMENT,STATUS) VALUES('{0}','{1}','{2}','{3}','{4}')"

os.chdir("INPUT/SAMPLES")
file = "DemographicData_INS.csv"

with open(file) as f:
    reader = csv.reader(f,delimiter=",")
    next(reader, None)
    for row in reader:
        print(row)
        cursor.execute(sql1.format(*row))

db.commit()


