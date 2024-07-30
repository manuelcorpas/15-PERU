mport subprocess
import os
import argparse

def annotate_f5_gene(input_vcf, output_vcf):
    """
    Annotate the F5 gene in the input VCF file using panno.

    :param input_vcf: Path to the input VCF file.
    :param output_vcf: Path to save the annotated output VCF file.
    """
    # Construct the command to run panno
    command = ["panno", "annotate", input_vcf, output_vcf, "--gene", "F5"]

    # Run the panno command
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"Annotation completed successfully. Output saved to {output_vcf}")
    else:
        print(f"Annotation failed. Error: {result.stderr}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate the F5 gene in a VCF file using panno.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input VCF file.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the annotated output VCF file.")

    args = parser.parse_args()

    input_vcf = args.input
    output_vcf = args.output

    # Annotate the VCF file with the F5 gene
    annotate_f5_gene(input_vcf, output_vcf)

