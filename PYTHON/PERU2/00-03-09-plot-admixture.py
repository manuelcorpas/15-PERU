import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Path to your files
population_file_path = 'ANALYSIS/10-PCA/ii_28_populations.txt'
vcf_file_path = 'ANALYSIS/10-PCA/common_variants_ibd_clean.vcf.gz'
plink_executable_path = 'plink'  # if PLINK is used for conversion
admixture_executable_path = 'admixture'
plink_format_output = 'ANALYSIS/11-ADMIXTURE/common_variants_ibd_clean'
K=3

# Read ADMIXTURE output (assuming K=3)
admixture_output_file = plink_format_output+'_pruned.' + str(K) + '.Q'
admixture_data = pd.read_csv(admixture_output_file, header=None, sep='\s+')

# Read population data
population_data = pd.read_csv(population_file_path, sep='\t')
population_data = population_data[['IID', 'Population']]  # Keep only necessary columns

# The columns in admixture_data are the ancestral populations, and rows are individuals
# Make sure the individual IDs in both DataFrames align correctly
# Assuming the order and number of individuals are the same in both DataFrames

# Correctly merging ADMIXTURE output with population data
merged_data = pd.concat([population_data, admixture_data], axis=1)

# Convert the ancestral proportion columns to numeric types
numeric_columns = merged_data.columns[2:]
merged_data[numeric_columns] = merged_data[numeric_columns].apply(pd.to_numeric, errors='coerce')

# Drop any rows with NaN values resulting from the conversion
merged_data.dropna(subset=numeric_columns, inplace=True)

# Group data by population and calculate mean ancestral proportions for each population
grouped_data = merged_data.groupby('Population')[numeric_columns].mean()

# Specify the ancestry column to sort by (e.g., the first ancestry column)
sort_ancestry = numeric_columns[0]  # Adjust this if necessary

# Sort the populations based on the chosen ancestry
grouped_data = grouped_data.sort_values(by=sort_ancestry, ascending=False)

# Prepare data for stacked bar plot
barWidth = 0.85
indices = np.arange(grouped_data.shape[0])

# Plotting
plt.figure(figsize=(12, 8))

bottom = np.zeros(grouped_data.shape[0])

for ancestry in numeric_columns:
    ancestry_means = grouped_data[ancestry]
    plt.bar(indices, ancestry_means, bottom=bottom, edgecolor='white', width=barWidth, label=f'Ancestry {ancestry}')
    bottom += ancestry_means.values

plt.xlabel('Population', fontweight='bold')
plt.xticks(indices, grouped_data.index, rotation='vertical')  # Set x-tick labels to be vertical
plt.ylabel('Average Ancestral Proportion')
plt.title('Average Ancestral Proportions by Population (Sorted by ' + str(sort_ancestry) + ')')
plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
plt.tight_layout()
plt.show()
'''
### Explanation of the Visualization Code:

1. **Read ADMIXTURE Output**: The script reads the ADMIXTURE output file, which contains the ancestry proportions for each individual.

2. **Read Population Data**: It then reads the population data file and merges it with the ADMIXTURE results. This step assumes that the order of individuals in the ADMIXTURE output and the population file is the same. If not, you'll need to adjust the script to correctly match individuals between the two files.

3. **Data Transformation**: The data is transformed (or "melted") for easier plotting. This restructures the dataframe so that each row represents a single ancestry proportion for an individual.

4. **Plotting**: The script uses seaborn to create a bar plot showing the average admixture proportions for each ancestral population across the different populations in your dataset.

Please ensure that you have `pandas`, `matplotlib`, and `seaborn` installed in your Python environment to run this script. Adjust the plot settings as necessary for your specific needs or preferences.
'''
