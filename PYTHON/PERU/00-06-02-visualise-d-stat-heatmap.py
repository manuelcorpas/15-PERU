import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Replace with the path to your D-statistics results file
file_path = 'ANALYSIS/15-DI-STAT/d_statistics_results.txt'

# Load the data into a DataFrame
# Ensure that the delimiter matches the format of your file
df = pd.read_csv(file_path, delim_whitespace=True)

# Pivot the DataFrame to create a matrix suitable for a heatmap
# The matrix will have pop2 as rows, pop3 as columns, and est values as cell values
pivot_df = df.pivot_table(index='pop2', columns='pop3', values='est', aggfunc=np.mean)

# Create the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(pivot_df, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('D-statistics Heatmap')
plt.xlabel('Population 3')
plt.ylabel('Population 2')
plt.show()





