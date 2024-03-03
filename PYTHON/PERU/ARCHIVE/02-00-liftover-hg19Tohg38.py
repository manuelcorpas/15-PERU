import subprocess
import os

def parse_original_vcf(vcf_file):
    """
    Parse the original VCF file to extract and store the header, along with necessary information for each variant.
    """
    vcf_header = []
    vcf_data = {}
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                vcf_header.append(line.strip())
                continue
            parts = line.strip().split('\t')
            chrom, pos, id, ref, alt, qual, filter, info, format, *genotypes = parts
            vcf_data[(chrom, pos)] = (id, ref, alt, qual, filter, info, format, genotypes)
    return vcf_header, vcf_data

def vcf_to_bed(vcf_file, bed_file):
    """
    Convert VCF to UCSC BED format. Only include chromosomal positions.
    """
    with open(vcf_file, 'r') as vcf, open(bed_file, 'w') as bed:
        for line in vcf:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            chrom, pos = parts[0], parts[1]
            bed.write(f"{chrom}\t{int(pos)-1}\t{pos}\n")

def run_liftover(input_bed, output_bed, unmapped_bed, chain_file, liftover_path):
    """
    Run the UCSC LiftOver tool to convert coordinates from one genome build to another.
    """
    subprocess.run([liftover_path, input_bed, chain_file, output_bed, unmapped_bed], check=True)

def convert_bed_to_vcf(bed_file, original_vcf_data, original_vcf_header, output_vcf_file):
    """
    Convert the BED file back to VCF format using the data and header from the original VCF.
    """
    with open(bed_file, 'r') as bed, open(output_vcf_file, 'w') as vcf:
        # Write the header
        for header_line in original_vcf_header:
            vcf.write(header_line + "\n")

        # Write the VCF content
        for line in bed:
            chrom, start, end = line.strip().split('\t')
            pos = str(int(start) + 1)  # Convert back to 1-based position
            if (chrom, pos) in original_vcf_data:
                id, ref, alt, qual, filter, info, format, genotypes = original_vcf_data[(chrom, pos)]
                vcf_entry = [chrom, pos, id, ref, alt, qual, filter, info, format] + genotypes
                vcf.write('\t'.join(vcf_entry) + "\n")

def process_vcf_files(input_dir, output_dir, chain_file_path, liftover_executable_path):
    """
    Process all VCF files in the given input directory and save the lifted-over VCFs in the output directory.
    """
    for filename in os.listdir(input_dir):
        if filename.endswith(".vcf"):
            original_vcf = os.path.join(input_dir, filename)
            intermediate_bed = os.path.join(output_dir, filename.replace('.vcf', '_intermediate.bed'))
            lifted_bed = os.path.join(output_dir, filename.replace('.vcf', '_lifted.bed'))
            unmapped_bed_file = os.path.join(output_dir, filename.replace('.vcf', '_unmapped.bed'))
            new_vcf = os.path.join(output_dir, filename.replace('.vcf', '_hg38.vcf'))

            print(f"Processing {filename}")

            # Convert original VCF to BED format
            vcf_to_bed(original_vcf, intermediate_bed)

            # Perform the liftover
            run_liftover(intermediate_bed, lifted_bed, unmapped_bed_file, chain_file_path, liftover_executable_path)

            # Parse the original VCF file
            vcf_header, vcf_data = parse_original_vcf(original_vcf)

            # Convert lifted-over BED back to VCF
            convert_bed_to_vcf(lifted_bed, vcf_data, vcf_header, new_vcf)

            print(f"Completed processing {filename}")

# File paths and parameters
input_directory = "ANALYSIS/17-LIFTOVER/SPLIT"
output_directory = "ANALYSIS/17-LIFTOVER/SPLIT"
chain_file_path = "ANALYSIS/17-LIFTOVER/hg19ToHg38.over.chain.gz"
liftover_executable_path = 'ANALYSIS/17-LIFTOVER/liftover'

# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)

# Process all VCF files in the input directory
process_vcf_files(input_directory, output_directory, chain_file_path, liftover_executable_path)
