import subprocess
import gzip
from ete3 import Tree

def extract_and_read_tree(gz_file_path):
    # Assuming the tree is in Newick format
    with gzip.open(gz_file_path, 'rt') as f:
        return Tree(f.read().strip())


# Step 1: Prepare Input Data
# This step depends on your data format and might require external tools or scripts.
# Replace this with commands or Python code to prepare your TreeMix input data.

input_file = 'ANALYSIS/13-TREEMIX/treemix_in.gz'  # TreeMix input file
treemix_output = 'ANALYSIS/13-TREEMIX/treemix_output'  # Output file prefix

# Step 2: Run TreeMix
# Root: AFRO_DES
# Grouped together SNPs to account for LD (-k)
# Build the ML G graph with 10 migration events
treemix_command = f'treemix -i {input_file} -root AFRO_DES -k 1000 -o {treemix_output}'
subprocess.run(treemix_command, shell=True)

tree = extract_and_read_tree('ANALYSIS/13-TREEMIX/treemix_output.treeout.gz')

# Render the tree
tree.show()
