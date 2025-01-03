import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

# Load the input file
path = "ANALYSIS/00-05-PyPGx/"
input_file = path + "genotype_frequency_by_gene.tsv"
df = pd.read_csv(input_file, sep="\t")

# Ensure correct column names
df.columns = df.columns.str.strip()

# Consolidate data by summing frequencies for each Gene-Phenotype combination
consolidated_df = df.groupby(["Gene", "Phenotype"], as_index=False).agg({
    "Genotype Frequency (%)": "sum"
})

# Sort genes alphabetically
consolidated_df = consolidated_df.sort_values(by="Gene")

# Generate a horizontal stacked bar chart
plt.figure(figsize=(14, 10))

genes = consolidated_df["Gene"].unique()
phenotypes = sorted(consolidated_df["Phenotype"].unique())
y_positions = range(len(genes))

# Use a larger colormap for distinct phenotype colours
cmap = get_cmap("tab20c")
phenotype_colors = {phenotype: cmap(i / len(phenotypes)) for i, phenotype in enumerate(phenotypes)}

for i, gene in enumerate(genes):
    gene_data = consolidated_df[consolidated_df["Gene"] == gene]
    cumulative = 0

    for _, row in gene_data.iterrows():
        phenotype = row["Phenotype"]
        freq = row["Genotype Frequency (%)"]
        plt.barh(i, freq, left=cumulative, color=phenotype_colors[phenotype], edgecolor="black", height=0.8)
        cumulative += freq

# Add a legend for phenotypes
handles = [plt.Rectangle((0, 0), 1, 1, color=phenotype_colors[phenotype]) for phenotype in phenotypes]
plt.legend(handles, phenotypes, title="Phenotypes", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.yticks(y_positions, genes, fontsize=10)
plt.xticks(fontsize=10)
plt.xlabel("Phenotype Frequency (%)", fontsize=12)
plt.ylabel("Genes", fontsize=12)
plt.title("Predicted Phenotype Breakdown Per Pharmacogene", fontsize=14)  # Updated title
plt.grid(axis="x", linestyle="--", alpha=0.6)
plt.tight_layout()

# Save the plot
output_plot = path + "genotype_frequency_by_gene_alphabetical.png"
plt.savefig(output_plot, dpi=300)
plt.show()

print(f"Consolidated horizontal stacked bar chart saved as {output_plot}")

