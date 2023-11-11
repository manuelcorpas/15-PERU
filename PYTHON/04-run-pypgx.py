import subprocess

def run_command_for_gene(gene_name):
    # Replace this with the actual command you want to run
    command = f"pypgx run-chip-pipeline {gene_name} ANALYSIS/06-PyPGX/pipeline_{gene_name} ANALYSIS/06-PyPGX/Full_INS_hg37_Autosomic.vcf.gz"
    print(command)
    subprocess.run(command, shell=True)

def process_genes(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            gene_name = line.strip()
            run_command_for_gene(gene_name)

# Replace 'your_file.txt' with the path to your file
process_genes('INPUT/PHARMACOGENES/PHARMACOGENES.txt')

