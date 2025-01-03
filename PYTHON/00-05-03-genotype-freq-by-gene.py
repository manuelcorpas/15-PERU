import pandas as pd

# File paths
input_file = "ANALYSIS/00-05-PyPGx/all_genes_results.tsv"
output_file = "ANALYSIS/00-05-PyPGx/genotype_frequency_by_gene.tsv"

# Load the data
data = pd.read_csv(input_file, sep="\t")

# Ensure necessary columns are present
required_columns = ["Gene", "Genotype", "Phenotype"]
if not all(col in data.columns for col in required_columns):
    raise ValueError(f"Input file must contain {', '.join(required_columns)} columns.")

# Group data by Gene and Genotype
grouped = data.groupby(["Gene", "Genotype", "Phenotype"]).size().reset_index(name="No of Subjects with Genotype")

# Calculate frequencies
total_counts = grouped.groupby("Gene")["No of Subjects with Genotype"].transform("sum")
grouped["Genotype Frequency"] = (grouped["No of Subjects with Genotype"] / total_counts).round(4)
grouped["Genotype Frequency (%)"] = (grouped["Genotype Frequency"] * 100).round(4)

# Reorder columns for better readability
grouped = grouped[["Gene", "Genotype", "No of Subjects with Genotype", "Genotype Frequency", "Genotype Frequency (%)", "Phenotype"]]

# Save to file
grouped.to_csv(output_file, sep="\t", index=False)

print(f"Genotype frequencies and associated phenotypes by gene saved to {output_file}")

