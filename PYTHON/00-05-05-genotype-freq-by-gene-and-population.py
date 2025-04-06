import pandas as pd

# File paths
input_file = "ANALYSIS/00-05-PyPGx/all_genes_results.tsv"
output_file = "ANALYSIS/00-05-PyPGx/00-05-06-genotype-freq-by-gene-and-population.csv"

# Load the data
data = pd.read_csv(input_file, sep="\t")

# Debug: Check input file
print("Input file loaded. Columns:", data.columns)

# Ensure necessary columns are present
required_columns = ["Unnamed: 0", "Gene", "Genotype", "Phenotype"]
if not all(col in data.columns for col in required_columns):
    raise ValueError(f"Input file must contain {', '.join(required_columns)} columns.")

# Extract sample group from the first column
if "Unnamed: 0" in data.columns:
    data[["Sample Group", "Sample ID"]] = data["Unnamed: 0"].str.split("-", n=1, expand=True)
    print("Sample groups extracted. Example:")
    print(data[["Sample Group", "Sample ID", "Gene"]].head())
else:
    raise ValueError("Expected 'Unnamed: 0' column for sample group information.")

# Group data by Sample Group, Gene, and Genotype
grouped = data.groupby(["Sample Group", "Gene", "Genotype", "Phenotype"]).size().reset_index(name="No of Subjects with Genotype")

# Debug: Check grouping results
if grouped.empty:
    raise ValueError("No data after grouping. Verify input structure.")
print("Grouping results:")
print(grouped.head())

# Calculate frequencies for each Sample Group and Gene
total_counts = grouped.groupby(["Sample Group", "Gene"])["No of Subjects with Genotype"].transform("sum")
grouped["Genotype Frequency"] = (grouped["No of Subjects with Genotype"] / total_counts).round(4)
grouped["Genotype Frequency (%)"] = (grouped["Genotype Frequency"] * 100).round(4)

# Debug: Check frequency calculations
if grouped.empty:
    raise ValueError("No data after frequency calculation. Verify grouping logic.")
print("Frequency calculations:")
print(grouped.head())

# Reorder columns for better readability
grouped = grouped[["Sample Group", "Gene", "Genotype", "Phenotype", "No of Subjects with Genotype", "Genotype Frequency", "Genotype Frequency (%)"]]

# Sort by Gene
grouped = grouped.sort_values(by=["Gene", "Sample Group"])

# Save to file
grouped.to_csv(output_file, index=False)
print(f"Genotype frequencies by population and gene saved to {output_file}")

