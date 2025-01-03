import os
import zipfile
import pandas as pd

# Directory containing PyPGx results
output_dir = "ANALYSIS/00-05-PyPGx/results"
output_file = "ANALYSIS/00-05-PyPGx/all_genes_results.tsv"

def extract_zip_to_dataframe(zip_path, target_file="data.tsv"):
    """Extract a specified file from a zip archive and load it into a pandas DataFrame."""
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            # Find the target file within any subdirectory
            matching_files = [name for name in zf.namelist() if name.endswith(target_file)]
            if matching_files:
                with zf.open(matching_files[0]) as file:
                    df = pd.read_csv(file, sep="\t")
                return df
            else:
                print(f"File {target_file} not found in {zip_path}. Archive contents: {zf.namelist()}")
                return None
    except Exception as e:
        print(f"Error reading {zip_path}: {e}")
        return None

def collect_results(output_dir):
    """Collect all data.tsv files from results.zip in each gene directory."""
    all_results = []

    for gene in os.listdir(output_dir):
        gene_dir = os.path.join(output_dir, gene)
        if os.path.isdir(gene_dir):
            results_zip = os.path.join(gene_dir, "results.zip")
            if os.path.exists(results_zip):
                gene_name = gene.split("pipeline_")[-1]
                df = extract_zip_to_dataframe(results_zip)
                if df is not None:
                    df["Gene"] = gene_name  # Add a column for the gene name
                    all_results.append(df)

    return pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()

# Process all results and save to a single file
all_results_df = collect_results(output_dir)
if not all_results_df.empty:
    all_results_df.to_csv(output_file, sep="\t", index=False)
    print(f"Combined results saved to {output_file}")
else:
    print("No valid data found in results.")

