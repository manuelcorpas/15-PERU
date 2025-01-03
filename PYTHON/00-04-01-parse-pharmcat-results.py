import os
import json
import pandas as pd
import logging
from bs4 import BeautifulSoup

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_pharmcat_json_html(directory, output_file):
    """Parse PharmCAT match.json, phenotype.json, and HTML files to extract metabolizer statuses grouped by population, drug, gene, and haplotype."""
    data = []

    # Iterate over files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".phenotype.json"):
            base_filename = filename.replace(".phenotype.json", "")
            phenotype_filepath = os.path.join(directory, filename)
            match_filepath = os.path.join(directory, f"{base_filename}.match.json")
            html_filepath = os.path.join(directory, f"{base_filename}.report.html")

            logging.info(f"Processing files: {filename}, {base_filename}.match.json, and {base_filename}.report.html")

            # Determine population based on filename (assumes population is part of the filename)
            population = "Unknown"
            for pop in ["MOCHES", "MATZES", "CUSCO", "TRUJILLO", "CHOPCCAS", "UROS", "IQUITOS"]:
                if pop in filename.upper():
                    population = pop
                    break

            if population == "Unknown":
                logging.warning(f"Population could not be determined for file: {filename}")

            # Load phenotype JSON content
            try:
                with open(phenotype_filepath, 'r') as phenotype_file:
                    phenotype_data = json.load(phenotype_file)

                    # Extract drug-level data from phenotype JSON
                    drug_reports = phenotype_data.get("drugReports", {})
                    for drug, report in drug_reports.items():
                        phenotypes = report.get("phenotypes", [])
                        gene = report.get("gene", "Unknown")
                        haplotypes = report.get("haplotypes", [])

                        # Ensure each phenotype is paired with its corresponding haplotype
                        for i, phenotype in enumerate(phenotypes):
                            haplotype = haplotypes[i] if i < len(haplotypes) else "Unknown"
                            data.append({
                                "Population": population,
                                "Drug": drug,
                                "Gene": gene,
                                "Haplotype": haplotype,
                                "Metabolizer Status": phenotype
                            })
            except json.JSONDecodeError as e:
                logging.error(f"Error decoding JSON in phenotype file {phenotype_filepath}: {e}")

            # Parse HTML file for additional drug and gene information
            if os.path.exists(html_filepath):
                try:
                    with open(html_filepath, 'r') as html_file:
                        soup = BeautifulSoup(html_file, 'html.parser')

                        # Locate drug and gene information in HTML
                        for guideline in soup.find_all('section', class_='guideline'):
                            drug_name = guideline.find('h3').text.strip()
                            gene_element = guideline.find('p', class_='gene-info')  # Adjust selector as needed
                            gene_name = gene_element.text.strip() if gene_element else "Unknown"

                            phenotype_elements = guideline.select('div.hint + p')
                            metabolizer_status = phenotype_elements[0].text.strip() if phenotype_elements else "Unknown"

                            genotype_elements = guideline.select('table.genotype-table td')  # Adjust selector as needed
                            haplotypes = [elem.text.strip() for elem in genotype_elements] if genotype_elements else ["Unknown"]

                            for haplotype in haplotypes:
                                data.append({
                                    "Population": population,
                                    "Drug": drug_name,
                                    "Gene": gene_name,
                                    "Haplotype": haplotype,
                                    "Metabolizer Status": metabolizer_status
                                })
                except Exception as e:
                    logging.error(f"Error parsing HTML file {html_filepath}: {e}")

    # Convert to DataFrame
    if not data:
        logging.warning("No data was collected from the provided files.")
        return

    df = pd.DataFrame(data)

    # Group by population, drug, gene, haplotype, and metabolizer status, and calculate frequencies
    grouped = df.groupby(["Population", "Drug", "Gene", "Haplotype", "Metabolizer Status"]).size().reset_index(name="Frequency")

    # Save to CSV
    grouped.to_csv(output_file, index=False)
    logging.info(f"Results grouped by frequency saved to {output_file}")

def main():
    input_directory = "ANALYSIS/00-04-PHARMCAT/results"
    output_file = os.path.join("ANALYSIS/00-04-PHARMCAT/results", "metabolizer_status_by_population_and_drug.csv")
    parse_pharmcat_json_html(input_directory, output_file)

if __name__ == "__main__":
    main()

