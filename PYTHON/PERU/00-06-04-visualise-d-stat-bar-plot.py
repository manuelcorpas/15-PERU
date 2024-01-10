import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Replace with the path to your D-statistics results file
file_path = 'ANALYSIS/15-DI-STAT/d_statistics_results.txt'

# Load the data into a DataFrame
# Ensure that the delimiter matches the format of your file
df = pd.read_csv(file_path, delim_whitespace=True)



# Unique combinations of pop1 and pop4
unique_combinations = df[['pop1', 'pop4']].drop_duplicates()

# Plotting multiple bar plots
for index, row in unique_combinations.iterrows():
    subset_df = df[(df['pop1'] == row['pop1']) & (df['pop4'] == row['pop4'])]
    plt.figure(figsize=(12, 6))
    sns.barplot(data=subset_df, x='pop2', y='est', hue='pop3')
    plt.title(f'D-statistics for Pop1: {row["pop1"]} and Pop4: {row["pop4"]}')
    plt.xlabel('Population 2')
    plt.ylabel('Estimated D-statistic (est)')
    plt.xticks(rotation=45)
    plt.legend(title='Population 3', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()
