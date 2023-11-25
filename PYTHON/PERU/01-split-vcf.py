import os

os.chdir("ANALYSIS/02-PGEN-PSAM-PVAR-2-VCF/")

def split_vcf(input_filename, max_size_mb=500):
    max_size = max_size_mb * 1024 * 1024  # Convert MB to bytes

    with open(input_filename, 'r') as file:
        headers = []
        current_chromosome = None
        output_file = None
        output_file_number = 0

        for line in file:
            # Collecting headers
            if line.startswith("#"):
                headers.append(line)
                continue

            # Getting chromosome from the line
            chromosome = line.split('\t')[0]

            # Check whether to start a new file
            if (output_file is None) or \
               (chromosome != current_chromosome) or \
               ((os.path.getsize(output_file.name) + len(line)) > max_size):
               
                if output_file:
                    output_file.close()

                output_file_number += 1
                output_file = open(f"split_{output_file_number}.vcf", 'w')
                output_file.writelines(headers)

            # Writing the line to the file
            output_file.write(line)

            # Updating the current chromosome
            current_chromosome = chromosome

        if output_file:
            output_file.close()

    print(f"Splitting completed: {output_file_number} files created.")

# Replace 'your_file.vcf' with the name of your input VCF file
split_vcf('Peru.joint.vcf')

