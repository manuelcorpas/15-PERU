import os
import subprocess
import logging
import shutil

# Configure logging
logging.basicConfig(filename='pypgx_analysis.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Input files and directories
vcf_file = "ANALYSIS/00-05-PyPGx/Peru.joint.vcf.gz"  # Updated VCF file
gene_list_file = "INPUT/PHARMACOGENES/PHARMACOGENES.txt"
output_dir = "ANALYSIS/00-05-PyPGx/results"
os.makedirs(output_dir, exist_ok=True)

def run_ngs_pipeline_for_gene(gene_name):
    """Run pypgx NGS pipeline for a specific gene."""
    gene_output_dir = f"{output_dir}/pipeline_{gene_name}"

    # Remove the output directory if it already exists
    if os.path.exists(gene_output_dir):
        logging.warning(f"Directory already exists for gene {gene_name}, removing: {gene_output_dir}")
        shutil.rmtree(gene_output_dir)

    try:
        command = (
            f"pypgx run-ngs-pipeline {gene_name} "
            f"{gene_output_dir} "
            f"--variants {vcf_file} "
            f"--assembly GRCh37"
        )
        logging.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running pipeline for gene {gene_name}: {e}")
        raise

def process_genes(file_path):
    """Process all genes listed in the input file."""
    try:
        with open(file_path, 'r') as file:
            for line in file:
                gene_name = line.strip()
                logging.info(f"Processing gene: {gene_name}")
                run_ngs_pipeline_for_gene(gene_name)
    except Exception as e:
        logging.error(f"Error processing genes: {e}")
        raise

# Process genes from the list
if not os.path.exists(gene_list_file):
    logging.error(f"Gene list file not found: {gene_list_file}")
    raise FileNotFoundError(f"Gene list file not found: {gene_list_file}")

logging.info("Starting processing of gene list...")
process_genes(gene_list_file)
logging.info("Processing complete.")
print("Pipeline completed successfully.")

