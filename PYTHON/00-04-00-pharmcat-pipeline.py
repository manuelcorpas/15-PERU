import os
import subprocess
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command, description):
    """Run a shell command and log its progress."""
    logging.info(f"{description}...")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        logging.error(f"Failed: {description}")
        raise RuntimeError(f"Command failed: {command}")
    logging.info(f"Completed: {description}")

def preprocess_vcf(input_vcf, bed_file, ref_fasta, ref_vcf, output_dir):
    """Preprocess the VCF file for PharmCAT."""
    sorted_vcf = os.path.join(output_dir, "Peru.joint.pharmacogenes.sorted.vcf.gz")
    preprocessed_vcf = os.path.join(output_dir, "Peru.joint.pharmacogenes.sorted.preprocessed.vcf.bgz")

    # Extract pharmacogene regions
    run_command(
        f"bcftools view -R {bed_file} {input_vcf} -Oz -o {output_dir}/Peru.joint.pharmacogenes.vcf.gz",
        "Extract pharmacogene regions"
    )

    # Sort and index the VCF
    run_command(
        f"bcftools sort {output_dir}/Peru.joint.pharmacogenes.vcf.gz -Oz -o {sorted_vcf}",
        "Sort pharmacogene VCF"
    )
    run_command(f"bcftools index {sorted_vcf}", "Index pharmacogene VCF")

    # Preprocess VCF for PharmCAT
    run_command(
        f"python3 ~/CloudDocs/SOFTWARE/PharmCAT/pharmcat-preprocessor/pharmcat_vcf_preprocessor.py \
        -vcf {sorted_vcf} \
        -refFna {ref_fasta} \
        -refVcf {ref_vcf} \
        -o {output_dir}/",
        "Preprocess VCF for PharmCAT"
    )

    return preprocessed_vcf

def run_pharmcat(preprocessed_vcf, output_dir, pharmcat_jar):
    """Run PharmCAT on the preprocessed VCF."""
    run_command(
        f"java -jar {pharmcat_jar} --matcher-vcf {preprocessed_vcf} --output-dir {output_dir}",
        "Run PharmCAT"
    )

def main():
    # Define file paths
    input_vcf = "ANALYSIS/00-04-PHARMCAT/Peru.joint.GRCh38.annotated.vcf.gz"
    bed_file = "ANALYSIS/00-04-PHARMCAT/pharmacogenes_GRCh38_variant_annotation_filtered_corrected.bed"
    ref_fasta = "ANALYSIS/00-04-PHARMCAT/Homo_sapiens.UCSC.GRCh38.fa"
    ref_vcf = "ANALYSIS/00-04-PHARMCAT/pharmcat_positions.vcf.bgz"
    pharmcat_jar = "ANALYSIS/00-04-PHARMCAT/pharmcat-2.15.3-all.jar"
    output_dir = "ANALYSIS/00-04-PHARMCAT/results"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    try:
        logging.info("Starting PharmCAT pipeline...")

        # Step 1: Preprocess the VCF
        preprocessed_vcf = preprocess_vcf(input_vcf, bed_file, ref_fasta, ref_vcf, output_dir)

        # Step 2: Run PharmCAT
        run_pharmcat(preprocessed_vcf, output_dir, pharmcat_jar)

        logging.info("PharmCAT preprocessing and report generation completed. Please run the result parser script for further analysis.")

    except Exception as e:
        logging.error(f"Pipeline failed: {e}")

if __name__ == "__main__":
    main()

