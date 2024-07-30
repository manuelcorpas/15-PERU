import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Load the dataset
file_path = 'ANALYSIS/04-ALLELE-COUNT/total_ac_vs_gnomad_amr.csv'
data = pd.read_csv(file_path)

# Initialize a list to store p-values
p_values = []

# Adding logging for debugging
log = []

# Iterate over each row in the DataFrame
for index, row in data.iterrows():
    # Extract the necessary counts for the contingency table
    total_allele_count = row['TOTAL Allele Count']
    total_allele_number = row['TOTAL Allele Number']
    gnomad_allele_count = row['GNOMAD_AMR Allele Count']
    gnomad_allele_number = row['GNOMAD_AMR Allele Number']
    
    # Construct the contingency table
    table = [
        [total_allele_count if total_allele_count > 0 else 0.1, total_allele_number - total_allele_count if (total_allele_number - total_allele_count) > 0 else 0.1],
        [gnomad_allele_count if gnomad_allele_count > 0 else 0.1, gnomad_allele_number - gnomad_allele_count if (gnomad_allele_number - gnomad_allele_count) > 0 else 0.1]
    ]

    # Log the constructed contingency table
    log.append(f"Row {index}: Contingency Table: {table}")

    # Perform Fisher's exact test
    try:
        _, p_value = fisher_exact(table, alternative='two-sided')
        log.append(f"Row {index}: Fisher's Exact Test computed successfully.")
    except ValueError as e:
        log.append(f"Row {index}: Fisher's Exact Test ValueError: {e}")
        p_value = pd.NA

    p_values.append(p_value)

# Applying the Benjamini-Hochberg procedure to control the false discovery rate
rejected, p_values_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# Add the list of corrected p-values as a new column to the DataFrame
data['p_value_corrected'] = p_values_corrected

# Save the updated DataFrame to a new CSV file
output_file_path = 'ANALYSIS/04-ALLELE-COUNT/updated_total_ac_vs_gnomad_amr.csv'
data.to_csv(output_file_path, index=False)

# Optionally, save or print the log for debugging
log_file_path = 'ANALYSIS/04-ALLELE-COUNT/fisher_test_log.txt'
with open(log_file_path, 'w') as log_file:
    log_file.write('\n'.join(log))

print("Fisher's Exact Test with multiple testing correction completed. Results saved to:", output_file_path)
print(f"Debug log saved to: {log_file_path}")

