import pandas as pd
from sqlalchemy import create_engine

# MariaDB connection (modify with your actual credentials and database details)
username = 'admin'
password = 'root'
host = 'localhost'
database = 'peru2'
port = 3306  # default MariaDB port

# SQLAlchemy engine for MariaDB
engine = create_engine(f'mysql+pymysql://{username}:{password}@{host}:{port}/{database}')

# Load the CSV file
file_path = 'INPUT/SAMPLES/peru-biobank-pruned-columns.csv'
data = pd.read_csv(file_path)
print(data.head())

# Filter out rows where 'cod microarray' or 'cod genoma' is blank
#filtered_data = data[data['cod microarray'].notna() | data['cod genoma'].notna()]
#print(filtered_data.head())

# Upload filtered data to the MariaDB database
#filtered_data.to_sql('00_00_PERU_BIOBANK_PRUNED', con=engine, if_exists='replace', index=False)


data.to_sql('00_00_PERU_BIOBANK_PRUNED', con=engine, if_exists='replace', index=False)
