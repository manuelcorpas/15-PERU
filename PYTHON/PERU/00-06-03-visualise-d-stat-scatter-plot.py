import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Replace with the path to your D-statistics results file
file_path = 'ANALYSIS/15-DI-STAT/d_statistics_results.txt'

# Load the data into a DataFrame
# Ensure that the delimiter matches the format of your file
df = pd.read_csv(file_path, delim_whitespace=True)


# Scatter plot of D-statistic vs. Z-score
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='est', y='z', hue='pop1', style='pop4', s=100)
plt.title('Scatter Plot of D-statistics vs. Z-scores')
plt.xlabel('Estimated D-statistic (est)')
plt.ylabel('Z-score (z)')
plt.legend(title='Populations', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.show()



