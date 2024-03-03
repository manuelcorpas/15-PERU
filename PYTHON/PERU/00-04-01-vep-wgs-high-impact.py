import sys
import mariadb
import glob
import os
from pathlib import Path


# Establish a connection to MariaDB database
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
cur = db.cursor(dictionary=True)

cur.execute('''
CREATE TABLE IF NOT EXISTS 00_04_VEP_WGS_HIGH_IMPACT AS
SELECT *
FROM 00_04_VEP_WGS
WHERE IMPACT = 'HIGH';
''')

db.close()

