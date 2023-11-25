import pandas as pd
from sqlalchemy import create_engine

# MariaDB connection (modify with your actual credentials and database details)
username = 'admin'
password = 'root'
host = 'localhost'
database = 'peru'
port = 3306  # default MariaDB port

# SQLAlchemy engine for MariaDB
engine = create_engine(f'mysql+pymysql://{username}:{password}@{host}:{port}/{database}')

# Load the CSV file
file_path = 'INPUT/SAMPLES/genoma_peruano_db.csv'
data = pd.read_csv(file_path)

# Optional: Convert data to the appropriate format if needed
# For example, converting date strings to datetime objects, handling NaNs, etc.

# Upload data to the MariaDB database
data.to_sql('00_GENOMA_PERUANO', con=engine, if_exists='replace', index=False)

