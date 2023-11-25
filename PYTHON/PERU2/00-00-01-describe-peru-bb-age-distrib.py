from sqlalchemy import create_engine
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# MariaDB connection (modify with your actual credentials and database details)
username = 'admin'
password = 'root'
host = 'localhost'
database = 'peru2'
port = 3306  # default MariaDB port

# SQLAlchemy engine for MariaDB
engine = create_engine(f'mysql+pymysql://{username}:{password}@{host}:{port}/{database}')

# Query to select all data from the table
query = 'SELECT * FROM 00_00_PERU_BIOBANK_PRUNED'

# Reading the data into a pandas DataFrame
df = pd.read_sql(query, engine)

# Replace invalid date entries with 'NaT' (Not a Time)
df['110 cuando nacio'] = pd.to_datetime(df['110 cuando nacio'], errors='coerce', format='%d/%m/%Y')


# Filter the DataFrame to include only records with at least one of the two columns populated
df = df[df['cod microarray'].notna() | df['cod genoma'].notna()]

# Drop rows where '110 cuando nacio' is NaT (if you want to exclude them from the analysis)
df = df.dropna(subset=['110 cuando nacio'])

# Calculate age for the rows with valid dates
current_year = datetime.now().year
df['age'] = current_year - df['110 cuando nacio'].dt.year

# Remove rows with ages less than 0 or greater than 120
df = df[(df['age'] >= 0) & (df['age'] <= 120)]

# Calculate the average (mean) and standard deviation
age_mean = df['age'].mean()
age_std = df['age'].std()

# Plotting
plt.figure(figsize=(10, 6))
count, bins, ignored = plt.hist(df['age'], bins=range(int(df['age'].min()), int(df['age'].max()) + 1, 1), alpha=0.7, color='blue', edgecolor='black')

# Plot vertical lines for the average and standard deviation
plt.axvline(age_mean, color='red', linestyle='dashed', linewidth=2)
plt.axvline(age_mean - age_std, color='green', linestyle='dashed', linewidth=2)
plt.axvline(age_mean + age_std, color='green', linestyle='dashed', linewidth=2)

# Annotate the mean, standard deviation, and total number of records on the plot
plt.text(age_mean, plt.ylim()[1]*0.9, f'Mean: {age_mean:.2f}', ha='center', color='red')
plt.text(age_mean - age_std, plt.ylim()[1]*0.95, f'-1 STD: {age_mean - age_std:.2f}', ha='center', color='green')
plt.text(age_mean + age_std, plt.ylim()[1]*0.95, f'+1 STD: {age_mean + age_std:.2f}', ha='center', color='green')

# Annotate the total number of records
plt.text(bins[0], plt.ylim()[1]*0.8, f'Total records: {len(df)}', ha='left', color='black')

plt.title('Age Distribution with Mean, Standard Deviation, and Total Records (Filtered Data)')
plt.xlabel('Age')
plt.ylabel('Frequency')

plt.show()



'''
# For 'codigo' column
print('--- codigo ---')
print('Number of unique values:', df['codigo'].nunique())
print('Number of missing values:', df['codigo'].isnull().sum())
print('Data type:', df['codigo'].dtype)
print('Sample values:', df['codigo'].dropna().unique()[:5])

# For 'cod microarray' column
print('--- cod microarray ---')
print('Number of unique values:', df['cod microarray'].nunique())
print('Number of missing values:', df['cod microarray'].isnull().sum())
print('Data type:', df['cod microarray'].dtype)
print('Sample values:', df['cod microarray'].dropna().unique()[:5])

# ... Repeat for each column ...

# For 'b11 bilirr indirecta' column
print('--- b11 bilirr indirecta ---')
print('Number of unique values:', df['b11 bilirr indirecta'].nunique())
print('Number of missing values:', df['b11 bilirr indirecta'].isnull().sum())
print('Data type:', df['b11 bilirr indirecta'].dtype)
print('Sample values:', df['b11 bilirr indirecta'].dropna().unique()[:5])
'''
