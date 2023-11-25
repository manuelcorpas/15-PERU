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

sql1 = "INSERT INTO 00_PHARMGKB(drug,fda,ema,swissmedic,hcsc,pmda) VALUES('{0}','{1}','{2}','{3}','{4}','{5}')"

os.chdir("INPUT/PHARMGKB")
file = "all-data.tsv"

with open(file) as f:
    reader = csv.reader(f,delimiter="\t")
    next(reader, None)
    for row in reader:
        print(row)
        cursor.execute(sql1.format(*row))

db.commit()


