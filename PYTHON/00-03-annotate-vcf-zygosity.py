import subprocess
import os

def run_command(command, output_file=None):
    """Runs a shell command and optionally writes output to a file."""
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        if output_file:
            with open(output_file, 'w') as f:
                f.write(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e.stderr}")
        raise

def extract_zygosities(vcf_file, output_tsv):
    """Extract zygosity information from a VCF file."""
    command = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' {vcf_file}"
    run_command(command, output_tsv)


def label_zygosities(input_tsv, output_tsv):
    """Label zygosities based on genotype."""
    with open(input_tsv, 'r') as infile, open(output_tsv, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            genotype = fields[5]
            if genotype == "0/0":
                zygosity = "Homozygous Reference"
            elif genotype in ["0/1", "1/0"]:
                zygosity = "Heterozygous"
            elif genotype == "1/1":
                zygosity = "Homozygous Alternate"
            else:
                zygosity = "Unknown"
            outfile.write('\t'.join(fields + [zygosity]) + '\n')


import tempfile

def sort_and_index_tsv(tsv_file):
    """Sort, compress, and index the annotation file."""
    # Sort the annotation file
    sorted_file = f"{tsv_file}.sorted"
    run_command(f"sort -k1,1 -k2,2n {tsv_file} > {sorted_file}")

    # Compress the sorted file
    compressed_file = f"{sorted_file}.gz"
    run_command(f"bgzip -c {sorted_file} > {compressed_file}")

    # Index the compressed file
    run_command(f"tabix -s 1 -b 2 -e 2 {compressed_file}")

    # Return the path to the compressed and indexed file
    return compressed_file

def annotate_vcf(vcf_file, tsv_file, output_vcf):
    """Annotate the VCF file with zygosity information."""
    # Sort, compress, and index the annotation file
    indexed_tsv = sort_and_index_tsv(tsv_file)

    # Create a temporary file for the header
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".txt") as header_file:
        header_file.write(
            '##INFO=<ID=ZYGOSITY,Number=1,Type=String,Description="Zygosity of the variant (Homozygous Reference, Heterozygous, Homozygous Alternate)">\n'
        )
        header_path = header_file.name

    # Annotate the VCF file
    command = (
        f"bcftools annotate -a {indexed_tsv} -h {header_path} "
        f"-c CHROM,POS,ID,REF,ALT,ZYGOSITY {vcf_file} -Oz -o {output_vcf}"
    )
    run_command(command)

    # Clean up temporary header file
    os.remove(header_path)

def validate_vcf(vcf_file):
    """Validate the VCF file by indexing and viewing the header."""
    try:
        run_command(f"bcftools index {vcf_file}")
        header = run_command(f"bcftools view {vcf_file} | head")
        print(header)
    except Exception as e:
        print("Validation failed.", str(e))


def main():
    # File paths (adjust as needed)
    original_vcf = "ANALYSIS/00-03-LIFTOVER/Peru.joint.vcf.gz"
    extracted_tsv = "ANALYSIS/00-03-LIFTOVER/zygosities.tsv"
    labeled_tsv = "ANALYSIS/00-03-LIFTOVER/zygosities_with_labels.tsv"
    annotated_vcf = "ANALYSIS/00-03-LIFTOVER/Peru.joint.GRCh38.annotated.vcf.gz"

    # Steps
    print("Step 1: Extracting zygosities from VCF file...")
    extract_zygosities(original_vcf, extracted_tsv)

    print("Step 2: Labelling zygosities...")
    label_zygosities(extracted_tsv, labeled_tsv)

    print("Step 3: Annotating VCF with zygosities...")
    annotate_vcf(original_vcf, labeled_tsv, annotated_vcf)

    print("Step 4: Validating annotated VCF file...")
    validate_vcf(annotated_vcf)

    print("Annotation completed successfully.")

if __name__ == "__main__":
    main()
