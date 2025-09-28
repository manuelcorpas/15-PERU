#!/usr/bin/env python3
"""
00-04-run-pypgx.py

Complete end-to-end PyPGx workflow pipeline.
Runs automatically with no command-line arguments.
All analyses performed using GRCh37 reference genome.

Usage:
    python 00-04-run-pypgx.py

Outputs saved to ANALYSIS/00-05-PyPGx/
"""

import os
import sys
import subprocess
import shutil
import zipfile
import logging
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, cohen_kappa_score
from scipy import stats
import requests
from bs4 import BeautifulSoup
from rapidfuzz import fuzz

# ============================================================================
# CONFIGURATION - HARD-CODED CONSTANTS
# ============================================================================

# Input files
JOINT_VCF = Path("ANALYSIS/00-05-PyPGx/Peru.joint.vcf.gz")
SAMPLE_FILE = Path("ANALYSIS/00-05-PyPGx/pop_736.updated.tsv")  # Contains all 736 samples
GENE_LIST = Path("INPUT/PHARMACOGENES/PHARMACOGENES.txt")
FDA_CSV = Path("ANALYSIS/00-05-PyPGx/fda_pharmacogenomics_associations.csv")

# Additional paths for compatibility
IBD_DIR = Path("ANALYSIS/00-06-IBD/paper936k")
COHORT_PREFIX = IBD_DIR / "merged_936k_final"

# Output directories
BASE_DIR = Path("ANALYSIS/00-05-PyPGx")
COHORT_DIR = BASE_DIR / "cohort"
COHORT_RAW_DIR = COHORT_DIR / "raw_zips"
WGS_TRUTH_DIR = BASE_DIR / "wgs_truth"
WGS_RAW_DIR = WGS_TRUTH_DIR / "raw_zips"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
LOGS_DIR = BASE_DIR / "logs"

# Assembly
ASSEMBLY = "GRCh37"

# Classification thresholds
ARRAY_SUPPORTED_THRESHOLDS = {
    "kappa_min": 0.85,
    "accuracy_min": 95.0,
    "indeterminate_max": 5.0
}

PARTIAL_SUPPORTED_THRESHOLDS = {
    "kappa_range": (0.70, 0.85),
    "accuracy_range": (85.0, 95.0),
    "indeterminate_range": (5.0, 15.0)
}

# Force these genes to WGS-only (complex loci)
WGS_ONLY_OVERRIDE = {"CYP2D6", "HLA-A", "HLA-B"}

# Fuzzy matching threshold for FDA mapping
FUZZY_THRESHOLD = 85

# Logging setup
LOG_FILE = LOGS_DIR / "pypgx_pipeline.log"

# ============================================================================
# SETUP FUNCTIONS
# ============================================================================

def setup_logging():
    """Configure logging."""
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def ensure_directories():
    """Create necessary output directories."""
    dirs = [
        COHORT_DIR, COHORT_RAW_DIR, WGS_TRUTH_DIR, WGS_RAW_DIR,
        RESULTS_DIR, FIGURES_DIR, LOGS_DIR
    ]
    for dir_path in dirs:
        dir_path.mkdir(parents=True, exist_ok=True)

def check_prerequisites():
    """Check that all required input files exist."""
    required_files = [JOINT_VCF, SAMPLE_FILE, GENE_LIST]
    
    missing = []
    for file_path in required_files:
        if not file_path.exists():
            missing.append(str(file_path))
    
    if missing:
        logger.error(f"Missing required files: {missing}")
        sys.exit(1)
    
    # Check PyPGx is installed
    try:
        result = subprocess.run(["pypgx", "--version"], capture_output=True, text=True)
        logger.info(f"PyPGx version: {result.stdout.strip()}")
    except:
        logger.error("PyPGx not found. Please install with: pip install pypgx")
        sys.exit(1)

def load_sample_lists():
    """Load and identify WGS vs array samples from TSV file."""
    # Read the sample file
    df_samples = pd.read_csv(SAMPLE_FILE, sep='\t')
    
    wgs_samples = []
    array_samples = []
    
    for _, row in df_samples.iterrows():
        # IID is the sample identifier
        sample_id = row['IID']
        
        # Identify sample type by ID pattern
        # Array samples have long IDs with underscores (e.g., 9512704011_R01C01_9512704011_R01C01)
        # WGS samples have short IDs with dashes (e.g., CHOPCCAS-250)
        if '_' in sample_id and len(sample_id) > 30:
            array_samples.append(sample_id)
        else:
            wgs_samples.append(sample_id)
    
    # Check for reasonable split
    if len(wgs_samples) == 0 or len(array_samples) == 0:
        logger.error("Failed to identify WGS vs array samples")
        logger.error(f"Found {len(wgs_samples)} WGS and {len(array_samples)} array samples")
        sys.exit(1)
    
    # Check total
    total = len(wgs_samples) + len(array_samples)
    logger.info(f"Loaded {len(wgs_samples)} WGS and {len(array_samples)} array samples (total: {total})")
    
    if total != 736:
        logger.warning(f"Total samples = {total}, expected 736")
    
    return wgs_samples, array_samples

def load_gene_list():
    """Load pharmacogenes to analyze."""
    genes = []
    with open(GENE_LIST, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.append(gene)
    
    logger.info(f"Loaded {len(genes)} genes to analyze")
    return genes

def load_metadata():
    """Generate sample metadata from the sample file."""
    # Read the sample file with population info
    df = pd.read_csv(SAMPLE_FILE, sep='\t')
    
    # Map populations to regions
    region_map = {
        # Amazon
        "IQUITOS": "Amazon", "MATZES": "Amazon", "SHIPIBO": "Amazon",
        "ASHANINKA": "Amazon", "CANDOSHI": "Amazon", "MACHIGUENGA": "Amazon",
        "NAHUA": "Amazon", "AWAJUN": "Amazon",
        # Andes  
        "CUSCO": "Andes", "CHOPCCAS": "Andes", "UROS": "Andes", 
        "PUNO": "Andes", "QUEROS": "Andes", "AYACUCHO": "Andes", 
        "JAQARUS": "Andes", "HUARAZ": "Andes", "CHACHAPOYAS": "Andes",
        # Coast
        "TRUJILLO": "Coast", "MOCHES": "Coast", "LIMA": "Coast",
        "LAMBAYEQUE": "Coast", "TUMBES": "Coast", "TALLAN": "Coast",
        "TACNA": "Coast", "MOQUEGUA": "Coast", "AREQUIPA": "Coast",
        # Other
        "AFRODESCENDIENTES": "Other", "LAMAS": "Other"
    }
    
    # Map populations to groups (simplified - adjust as needed)
    group_map = {
        # Indigenous groups
        "IQUITOS": "Mestizo", "MATZES": "Indigenous", "SHIPIBO": "Indigenous",
        "ASHANINKA": "Indigenous", "CANDOSHI": "Indigenous", 
        "MACHIGUENGA": "Indigenous", "NAHUA": "Indigenous", "AWAJUN": "Indigenous",
        "CUSCO": "Mestizo", "CHOPCCAS": "Indigenous", "UROS": "Indigenous",
        "PUNO": "Mestizo", "QUEROS": "Indigenous", "AYACUCHO": "Mestizo",
        "JAQARUS": "Indigenous", "HUARAZ": "Mestizo", "CHACHAPOYAS": "Mestizo",
        "TRUJILLO": "Mestizo", "MOCHES": "Indigenous", "LIMA": "Mestizo",
        "LAMBAYEQUE": "Mestizo", "TUMBES": "Mestizo", "TALLAN": "Indigenous",
        "TACNA": "Mestizo", "MOQUEGUA": "Mestizo", "AREQUIPA": "Mestizo",
        "AFRODESCENDIENTES": "Other", "LAMAS": "Indigenous"
    }
    
    # Create metadata dataframe
    metadata_df = pd.DataFrame({
        'SampleID': df['IID'],
        'Population': df['POP'],
        'Region': df['POP'].map(region_map).fillna('Unknown'),
        'Group': df['POP'].map(group_map).fillna('Unknown')
    })
    
    logger.info(f"Generated metadata for {len(metadata_df)} samples")
    logger.info(f"Regions: {metadata_df['Region'].value_counts().to_dict()}")
    logger.info(f"Groups: {metadata_df['Group'].value_counts().to_dict()}")
    
    return metadata_df

# ============================================================================
# PYPGX FUNCTIONS
# ============================================================================

def extract_samples_from_vcf(vcf_path: Path, sample_list: List[str], output_vcf: Path):
    """Extract specific samples from VCF using bcftools."""
    sample_file = output_vcf.parent / f"{output_vcf.stem}_samples.txt"
    
    with open(sample_file, 'w') as f:
        for sample in sample_list:
            f.write(f"{sample}\n")
    
    cmd = [
        "bcftools", "view",
        "-S", str(sample_file),
        "-o", str(output_vcf),
        "-O", "z",
        str(vcf_path)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        subprocess.run(["bcftools", "index", "-t", str(output_vcf)], check=True)
        logger.info(f"Created subset VCF: {output_vcf}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to extract samples: {e}")
        return False

def run_pypgx_for_gene(vcf_path: Path, gene: str, output_dir: Path):
    """Run PyPGx for a single gene."""
    gene_dir = output_dir / f"pypgx_{gene}"
    
    # Clean up if exists
    if gene_dir.exists():
        shutil.rmtree(gene_dir)
    
    cmd = [
        "pypgx", "run-ngs-pipeline", gene,
        str(gene_dir),
        "--variants", str(vcf_path),
        "--assembly", ASSEMBLY
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Move results.zip to raw_zips
        results_zip = gene_dir / "results.zip"
        if results_zip.exists():
            raw_dir = output_dir.parent / "raw_zips"
            raw_dir.mkdir(exist_ok=True)
            dest = raw_dir / f"{gene}_results.zip"
            shutil.move(str(results_zip), str(dest))
            
            # Clean up gene directory
            shutil.rmtree(gene_dir)
            
            return dest
        return None
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"PyPGx failed for {gene}: {e.stderr}")
        if gene_dir.exists():
            shutil.rmtree(gene_dir)
        return None

def extract_pypgx_results(zip_path: Path, gene: str):
    """Extract results from PyPGx output zip."""
    results = []
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            # Find data.tsv file
            data_files = [n for n in zf.namelist() if n.endswith("data.tsv")]
            if not data_files:
                return results
            
            with zf.open(data_files[0]) as f:
                df = pd.read_csv(f, sep="\t")
                
                if df.empty:
                    return results
                
                # Process each sample
                for idx, row in df.iterrows():
                    # Get sample ID from first column
                    sample_id = row.iloc[0] if not pd.isna(row.iloc[0]) else ""
                    
                    phenotype = str(row.get("Phenotype", ""))
                    
                    result = {
                        "SampleID": sample_id,
                        "Gene": gene,
                        "Phenotype": phenotype,
                        "Haplotype1": str(row.get("Haplotype1", "")),
                        "Haplotype2": str(row.get("Haplotype2", "")),
                        "Genotype": str(row.get("Genotype", "")),
                        "Score": row.get("Score", None),
                        "IndeterminateFlag": "Indeterminate" in phenotype
                    }
                    results.append(result)
    
    except Exception as e:
        logger.error(f"Error extracting results from {zip_path}: {e}")
    
    return results

def run_pypgx_batch(vcf_path: Path, genes: List[str], output_tag: str, output_dir: Path):
    """Run PyPGx for multiple genes on a VCF."""
    all_results = []
    successful_genes = []
    failed_genes = []
    
    for i, gene in enumerate(genes, 1):
        logger.info(f"[{output_tag}] Processing {gene} ({i}/{len(genes)})")
        
        zip_path = run_pypgx_for_gene(vcf_path, gene, output_dir)
        
        if zip_path:
            results = extract_pypgx_results(zip_path, gene)
            if results:
                all_results.extend(results)
                successful_genes.append(gene)
            else:
                failed_genes.append(gene)
        else:
            failed_genes.append(gene)
    
    # Create DataFrame
    if all_results:
        df = pd.DataFrame(all_results)
    else:
        df = pd.DataFrame()
    
    logger.info(f"[{output_tag}] Successful: {len(successful_genes)}, Failed: {len(failed_genes)}")
    if failed_genes:
        logger.warning(f"[{output_tag}] Failed genes: {failed_genes}")
    
    return df

# ============================================================================
# CONCORDANCE ANALYSIS
# ============================================================================

def compute_concordance(cohort_df: pd.DataFrame, wgs_df: pd.DataFrame):
    """Compute concordance metrics between cohort and WGS-truth results."""
    concordance_results = []
    confusion_matrices = {}
    discordant_calls = {}
    
    # Get common genes
    genes = sorted(set(cohort_df["Gene"]) & set(wgs_df["Gene"]))
    
    for gene in genes:
        gene_cohort = cohort_df[cohort_df["Gene"] == gene].copy()
        gene_wgs = wgs_df[wgs_df["Gene"] == gene].copy()
        
        # Merge on SampleID
        merged = pd.merge(
            gene_cohort[["SampleID", "Phenotype", "IndeterminateFlag"]],
            gene_wgs[["SampleID", "Phenotype", "IndeterminateFlag"]],
            on="SampleID",
            suffixes=("_array", "_wgs")
        )
        
        if merged.empty:
            continue
        
        # Filter out indeterminate calls
        valid_mask = ~merged["IndeterminateFlag_array"] & ~merged["IndeterminateFlag_wgs"]
        valid_data = merged[valid_mask].copy()
        
        # Calculate metrics
        n_valid = len(valid_data)
        indeterminate_rate = (merged["IndeterminateFlag_array"].sum() / len(merged)) * 100
        
        if n_valid > 0:
            array_pheno = valid_data["Phenotype_array"]
            wgs_pheno = valid_data["Phenotype_wgs"]
            
            # Accuracy
            accuracy = accuracy_score(wgs_pheno, array_pheno) * 100
            
            # Cohen's Kappa
            try:
                kappa = cohen_kappa_score(wgs_pheno, array_pheno)
            except:
                kappa = 0
            
            # Macro F1
            try:
                macro_f1 = f1_score(wgs_pheno, array_pheno, average='macro')
            except:
                macro_f1 = 0
            
            # Confusion matrix
            labels = sorted(set(array_pheno) | set(wgs_pheno))
            cm = confusion_matrix(wgs_pheno, array_pheno, labels=labels)
            confusion_matrices[gene] = (cm, labels)
            
            # Discordant calls
            discordant = valid_data[valid_data["Phenotype_array"] != valid_data["Phenotype_wgs"]]
            if not discordant.empty:
                discordant_calls[gene] = discordant[["SampleID", "Phenotype_array", "Phenotype_wgs"]]
        else:
            accuracy = kappa = macro_f1 = 0
        
        concordance_results.append({
            "Gene": gene,
            "Accuracy": accuracy,
            "MacroF1": macro_f1,
            "Kappa": kappa,
            "IndeterminateRate_Array": indeterminate_rate,
            "Npairs": len(merged),
            "N_valid": n_valid
        })
    
    concordance_df = pd.DataFrame(concordance_results)
    
    # Save concordance table
    concordance_df.to_csv(RESULTS_DIR / "PGx_concordance.csv", index=False)
    logger.info(f"Saved concordance results for {len(concordance_df)} genes")
    
    # Save confusion matrices
    for gene, (cm, labels) in confusion_matrices.items():
        save_confusion_matrix(gene, cm, labels)
    
    # Save discordant calls
    for gene, discordant in discordant_calls.items():
        discordant.to_csv(RESULTS_DIR / f"discordant_calls_{gene}.csv", index=False)
    
    return concordance_df

def save_confusion_matrix(gene: str, cm: np.ndarray, labels: List[str]):
    """Save confusion matrix as image."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=labels, yticklabels=labels,
                cbar_kws={'label': 'Count'}, ax=ax)
    
    ax.set_xlabel('Array Phenotype')
    ax.set_ylabel('WGS Phenotype')
    ax.set_title(f'Confusion Matrix: {gene}')
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / f"confmat_{gene}.png", dpi=150)
    plt.close()

# ============================================================================
# GENE CLASSIFICATION
# ============================================================================

def classify_genes(concordance_df: pd.DataFrame):
    """Classify genes based on concordance metrics."""
    classifications = []
    
    for _, row in concordance_df.iterrows():
        gene = row["Gene"]
        kappa = row["Kappa"]
        accuracy = row["Accuracy"]
        indeterminate = row["IndeterminateRate_Array"]
        n_pairs = row["Npairs"]
        
        # Force override for known complex genes
        if gene in WGS_ONLY_OVERRIDE:
            classification = "WGS-only"
            notes = "Known complex locus"
        # Check sample size
        elif n_pairs < 10:
            classification = "WGS-only"
            notes = f"Insufficient data (n={n_pairs})"
        # Array-supported
        elif (kappa >= ARRAY_SUPPORTED_THRESHOLDS["kappa_min"] and
              accuracy >= ARRAY_SUPPORTED_THRESHOLDS["accuracy_min"] and
              indeterminate <= ARRAY_SUPPORTED_THRESHOLDS["indeterminate_max"]):
            classification = "Array-supported"
            notes = "High concordance"
        # Partially supported
        elif ((PARTIAL_SUPPORTED_THRESHOLDS["kappa_range"][0] <= kappa < PARTIAL_SUPPORTED_THRESHOLDS["kappa_range"][1]) or
              (PARTIAL_SUPPORTED_THRESHOLDS["accuracy_range"][0] <= accuracy < PARTIAL_SUPPORTED_THRESHOLDS["accuracy_range"][1]) or
              (PARTIAL_SUPPORTED_THRESHOLDS["indeterminate_range"][0] <= indeterminate < PARTIAL_SUPPORTED_THRESHOLDS["indeterminate_range"][1])):
            classification = "Partially supported"
            notes = "Moderate concordance"
        # WGS-only
        else:
            classification = "WGS-only"
            notes = "Low concordance"
        
        classifications.append({
            "Gene": gene,
            "Accuracy": round(accuracy, 1),
            "MacroF1": round(row["MacroF1"], 3),
            "Kappa": round(kappa, 3),
            "IndeterminateRate": round(indeterminate, 1),
            "Npairs": n_pairs,
            "Classification": classification,
            "Notes": notes
        })
    
    classification_df = pd.DataFrame(classifications)
    classification_df.to_csv(RESULTS_DIR / "PGx_classification.csv", index=False)
    
    # Log summary
    class_counts = classification_df["Classification"].value_counts()
    logger.info(f"Gene classifications:\n{class_counts.to_string()}")
    
    return classification_df

# ============================================================================
# POPULATION SUMMARIES
# ============================================================================

def create_population_summaries(cohort_df: pd.DataFrame, classification_df: pd.DataFrame, metadata_df: pd.DataFrame):
    """Create population-level summaries for array-supported genes."""
    # Get array-supported genes
    array_genes = classification_df[
        classification_df["Classification"] == "Array-supported"
    ]["Gene"].tolist()
    
    if not array_genes:
        logger.warning("No array-supported genes found for population summaries")
        return pd.DataFrame()
    
    # Merge with metadata
    merged = pd.merge(cohort_df, metadata_df, on="SampleID", how="left")
    
    # Filter to array-supported genes and exclude indeterminate
    filtered = merged[
        (merged["Gene"].isin(array_genes)) & 
        (~merged["IndeterminateFlag"])
    ].copy()
    
    if filtered.empty:
        logger.warning("No valid data for population summaries")
        return pd.DataFrame()
    
    summaries = []
    
    for gene in array_genes:
        gene_data = filtered[filtered["Gene"] == gene]
        
        if gene_data.empty:
            continue
        
        # By Region × Group
        for region in gene_data["Region"].unique():
            for group in gene_data["Group"].unique():
                subset = gene_data[(gene_data["Region"] == region) & (gene_data["Group"] == group)]
                if len(subset) > 0:
                    pheno_counts = subset["Phenotype"].value_counts()
                    total = len(subset)
                    
                    for phenotype, count in pheno_counts.items():
                        summaries.append({
                            "Gene": gene,
                            "Level": "Region_Group",
                            "Category": f"{region}_{group}",
                            "Region": region,
                            "Group": group,
                            "Phenotype": phenotype,
                            "Count": count,
                            "Total": total,
                            "Frequency": count / total
                        })
        
        # By Population (if n >= 15)
        for pop in gene_data["Population"].unique():
            pop_data = gene_data[gene_data["Population"] == pop]
            if len(pop_data) >= 15:
                pheno_counts = pop_data["Phenotype"].value_counts()
                total = len(pop_data)
                
                for phenotype, count in pheno_counts.items():
                    summaries.append({
                        "Gene": gene,
                        "Level": "Population",
                        "Category": pop,
                        "Region": pop_data["Region"].iloc[0],
                        "Group": pop_data["Group"].iloc[0] if "Group" in pop_data.columns else "",
                        "Phenotype": phenotype,
                        "Count": count,
                        "Total": total,
                        "Frequency": count / total
                    })
    
    if summaries:
        summary_df = pd.DataFrame(summaries)
        summary_df.to_csv(RESULTS_DIR / "PGx_population_summary.csv", index=False)
        logger.info(f"Created {len(summary_df)} population summaries")
        
        # Create phenotype frequency plots
        create_phenotype_frequency_plots(summary_df, array_genes)
        
        return summary_df
    else:
        logger.warning("No population summaries generated")
        return pd.DataFrame()

def create_phenotype_frequency_plots(summary_df: pd.DataFrame, genes: List[str]):
    """Create stacked bar plots of phenotype frequencies."""
    for gene in genes:
        gene_data = summary_df[summary_df["Gene"] == gene]
        
        if gene_data.empty:
            continue
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # Plot 1: Region × Group
        ax = axes[0]
        region_group_data = gene_data[gene_data["Level"] == "Region_Group"].copy()
        if not region_group_data.empty:
            pivot = region_group_data.pivot_table(
                index="Category", columns="Phenotype", values="Frequency", fill_value=0
            )
            pivot.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
            ax.set_title(f'{gene}: Phenotype Frequencies by Region × Group')
            ax.set_xlabel('Region × Group')
            ax.set_ylabel('Frequency')
            ax.legend(title='Phenotype', bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        # Plot 2: Population (n >= 15)
        ax = axes[1]
        pop_data = gene_data[gene_data["Level"] == "Population"].copy()
        if not pop_data.empty:
            pivot = pop_data.pivot_table(
                index="Category", columns="Phenotype", values="Frequency", fill_value=0
            )
            pivot.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
            ax.set_title(f'{gene}: Phenotype Frequencies by Population (n ≥ 15)')
            ax.set_xlabel('Population')
            ax.set_ylabel('Frequency')
            ax.legend(title='Phenotype', bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / f"phenotype_frequencies_{gene}.png", dpi=300, bbox_inches='tight')
        plt.close()

# ============================================================================
# FDA MAPPING
# ============================================================================

def parse_fda_data():
    """Parse FDA pharmacogenomics associations from website if CSV doesn't exist."""
    if FDA_CSV.exists():
        logger.info(f"FDA data already exists: {FDA_CSV}")
        return
    
    logger.info("Fetching FDA pharmacogenomics associations...")
    
    url = "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations"
    
    try:
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception(f"Failed to fetch URL: {url}")
        
        soup = BeautifulSoup(response.text, 'html.parser')
        table = soup.find('table')
        
        if not table:
            raise Exception("Could not find table on FDA page")
        
        headers = [header.text.strip() for header in table.find_all('th')]
        
        rows = []
        for row in table.find_all('tr')[1:]:
            cells = [cell.text.strip() for cell in row.find_all('td')]
            if cells:
                rows.append(cells)
        
        df = pd.DataFrame(rows, columns=headers)
        df.to_csv(FDA_CSV, index=False)
        logger.info(f"FDA data saved to {FDA_CSV}")
        
    except Exception as e:
        logger.error(f"Failed to parse FDA data: {e}")
        logger.error("Please manually download FDA associations to: {FDA_CSV}")
        sys.exit(1)

def extract_phenotype_category(phenotype_text: str):
    """Extract standardized category from phenotype text."""
    text_lower = phenotype_text.lower()
    
    patterns = [
        (r"\bpoor metabolizer", "PM"),
        (r"\bintermediate metabolizer", "IM"),
        (r"\bnormal metabolizer", "NM"),
        (r"\bextensive metabolizer", "NM"),
        (r"\bultra[-\s]?rapid metabolizer", "UM"),
        (r"\brapid metabolizer", "RM"),
    ]
    
    for pattern, category in patterns:
        if re.search(pattern, text_lower):
            return category
    
    return None

def map_to_fda_recommendations(cohort_df: pd.DataFrame, classification_df: pd.DataFrame, metadata_df: pd.DataFrame):
    """Map phenotypes to FDA drug recommendations using deterministic matching."""
    # Parse FDA data if needed
    parse_fda_data()
    
    # Load FDA data
    fda_df = pd.read_csv(FDA_CSV)
    
    # Standardize gene names
    fda_df["Gene"] = fda_df["Gene"].str.upper()
    cohort_df["Gene"] = cohort_df["Gene"].str.upper()
    
    # Merge with metadata
    merged_df = pd.merge(cohort_df, metadata_df, on="SampleID", how="left")
    
    # Add classification
    merged_df = pd.merge(merged_df, classification_df[["Gene", "Classification"]], 
                         on="Gene", how="left")
    
    # Extract phenotype categories
    merged_df["PhenotypeCategory"] = merged_df["Phenotype"].apply(extract_phenotype_category)
    fda_df["AffectedCategories"] = fda_df["Affected Subgroups+"].apply(
        lambda x: extract_phenotype_category(x) if pd.notna(x) else None
    )
    
    # Match patients to FDA recommendations
    patient_matches = []
    fuzzy_matches = []
    
    for _, patient in merged_df.iterrows():
        if patient["IndeterminateFlag"] or pd.isna(patient["PhenotypeCategory"]):
            continue
        
        gene = patient["Gene"]
        category = patient["PhenotypeCategory"]
        
        # Find FDA entries for this gene
        gene_fda = fda_df[fda_df["Gene"] == gene]
        
        for _, fda_row in gene_fda.iterrows():
            fda_category = fda_row["AffectedCategories"]
            
            # Deterministic match on category
            if category and fda_category and category == fda_category:
                patient_matches.append({
                    "SampleID": patient["SampleID"],
                    "Region": patient.get("Region", "Unknown"),
                    "Group": patient.get("Group", "Unknown"),
                    "Population": patient.get("Population", "Unknown"),
                    "Gene": gene,
                    "GeneClassification": patient.get("Classification", "Unknown"),
                    "Phenotype": patient["Phenotype"],
                    "PhenotypeCategory": category,
                    "Drug": fda_row["Drug"],
                    "Recommendation": fda_row.get("Description of Gene-Drug Interaction", ""),
                    "MatchType": "Deterministic"
                })
            # Fallback to fuzzy if no category info
            elif not fda_category and patient["Phenotype"] and fda_row.get("Affected Subgroups+"):
                score = fuzz.partial_ratio(
                    patient["Phenotype"].lower(),
                    str(fda_row["Affected Subgroups+"]).lower()
                )
                if score >= FUZZY_THRESHOLD:
                    fuzzy_matches.append({
                        "SampleID": patient["SampleID"],
                        "Region": patient.get("Region", "Unknown"),
                        "Group": patient.get("Group", "Unknown"),
                        "Population": patient.get("Population", "Unknown"),
                        "Gene": gene,
                        "GeneClassification": patient.get("Classification", "Unknown"),
                        "Phenotype": patient["Phenotype"],
                        "PhenotypeCategory": category,
                        "Drug": fda_row["Drug"],
                        "Recommendation": fda_row.get("Description of Gene-Drug Interaction", ""),
                        "MatchType": f"Fuzzy (score={score})"
                    })
    
    # Log fuzzy matches
    if fuzzy_matches:
        logger.info(f"Used fuzzy matching for {len(fuzzy_matches)} FDA matches")
    
    # Combine matches
    all_matches = patient_matches + fuzzy_matches
    
    if not all_matches:
        logger.warning("No FDA matches found")
        return pd.DataFrame(), pd.DataFrame()
    
    matches_df = pd.DataFrame(all_matches)
    
    # Save patient-level matches
    matches_df.to_csv(RESULTS_DIR / "FDA_patient_matches.csv", index=False)
    logger.info(f"Saved {len(matches_df)} FDA patient matches")
    
    # Aggregate trigger counts
    trigger_counts = []
    
    # By Region × Group
    region_group_counts = matches_df.groupby(
        ["Region", "Group", "Gene", "GeneClassification"]
    ).size().reset_index(name="Count")
    region_group_counts["Level"] = "Region_Group"
    region_group_counts["Category"] = region_group_counts["Region"] + "_" + region_group_counts["Group"]
    trigger_counts.append(region_group_counts)
    
    # By Population (if n >= 15)
    pop_sample_counts = merged_df.groupby("Population")["SampleID"].nunique()
    valid_pops = pop_sample_counts[pop_sample_counts >= 15].index
    
    pop_matches = matches_df[matches_df["Population"].isin(valid_pops)]
    if not pop_matches.empty:
        pop_counts = pop_matches.groupby(
            ["Population", "Gene", "GeneClassification"]
        ).size().reset_index(name="Count")
        pop_counts["Level"] = "Population"
        pop_counts["Category"] = pop_counts["Population"]
        pop_counts["Region"] = pop_counts["Population"].map(
            merged_df.groupby("Population")["Region"].first()
        )
        pop_counts["Group"] = pop_counts["Population"].map(
            merged_df.groupby("Population")["Group"].first()
        )
        trigger_counts.append(pop_counts)
    
    if trigger_counts:
        trigger_df = pd.concat(trigger_counts, ignore_index=True)
        trigger_df.to_csv(RESULTS_DIR / "FDA_trigger_counts.csv", index=False)
        logger.info(f"Saved FDA trigger counts: {len(trigger_df)} entries")
        return matches_df, trigger_df
    else:
        return matches_df, pd.DataFrame()

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def create_concordance_figure(concordance_df: pd.DataFrame):
    """Create comprehensive concordance visualization."""
    # Sort by kappa
    plot_df = concordance_df[concordance_df["N_valid"] > 0].copy()
    plot_df = plot_df.sort_values("Kappa", ascending=False)
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    # Kappa subplot
    ax = axes[0]
    colors = ['green' if k >= 0.85 else 'orange' if k >= 0.70 else 'red' 
              for k in plot_df["Kappa"]]
    bars = ax.bar(range(len(plot_df)), plot_df["Kappa"], color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(y=0.85, color='green', linestyle='--', alpha=0.5, label='Array-supported')
    ax.axhline(y=0.70, color='orange', linestyle='--', alpha=0.5, label='Partially supported')
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df["Gene"], rotation=45, ha='right')
    ax.set_ylabel("Cohen's Kappa", fontsize=12)
    ax.set_title("Concordance Metrics by Gene", fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Accuracy subplot
    ax = axes[1]
    colors = ['green' if a >= 95 else 'orange' if a >= 85 else 'red' 
              for a in plot_df["Accuracy"]]
    ax.bar(range(len(plot_df)), plot_df["Accuracy"], color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(y=95, color='green', linestyle='--', alpha=0.5)
    ax.axhline(y=85, color='orange', linestyle='--', alpha=0.5)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df["Gene"], rotation=45, ha='right')
    ax.set_ylabel("Accuracy (%)", fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Indeterminate Rate subplot
    ax = axes[2]
    colors = ['green' if i <= 5 else 'orange' if i <= 15 else 'red' 
              for i in plot_df["IndeterminateRate_Array"]]
    ax.bar(range(len(plot_df)), plot_df["IndeterminateRate_Array"], color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(y=5, color='green', linestyle='--', alpha=0.5)
    ax.axhline(y=15, color='orange', linestyle='--', alpha=0.5)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df["Gene"], rotation=45, ha='right')
    ax.set_ylabel("Array Indeterminate Rate (%)", fontsize=12)
    ax.set_xlabel("Gene", fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "concordance_by_gene.png", dpi=300, bbox_inches='tight')
    plt.close()

def create_fda_trigger_figure(trigger_df: pd.DataFrame):
    """Create FDA trigger visualization."""
    if trigger_df.empty:
        logger.warning("No FDA triggers to visualize")
        return
    
    # Aggregate by gene
    gene_totals = trigger_df.groupby(["Gene", "GeneClassification"])["Count"].sum().reset_index()
    gene_totals = gene_totals.sort_values("Count", ascending=True)
    
    # Color map
    color_map = {
        "Array-supported": "#2ecc71",
        "Partially supported": "#f39c12",
        "WGS-only": "#e74c3c"
    }
    colors = [color_map.get(c, "#95a5a6") for c in gene_totals["GeneClassification"]]
    
    # Create lollipop plot
    fig, ax = plt.subplots(figsize=(12, max(8, len(gene_totals) * 0.4)))
    
    y_pos = np.arange(len(gene_totals))
    
    # Draw lines
    for i, (_, row) in enumerate(gene_totals.iterrows()):
        ax.plot([0, row["Count"]], [i, i], 'k-', alpha=0.3, linewidth=1)
    
    # Draw bars and dots
    ax.barh(y_pos, gene_totals["Count"], color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.scatter(gene_totals["Count"], y_pos, color=colors, s=100, zorder=5, edgecolors='black', linewidth=1)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(gene_totals["Gene"], fontsize=11)
    ax.set_xlabel("Number of FDA Recommendation Triggers", fontsize=12)
    ax.set_title("FDA Drug Recommendations by Gene and Platform Support", fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=label, edgecolor='black') 
                      for label, color in color_map.items()]
    ax.legend(handles=legend_elements, title="Gene Classification", 
              loc='lower right', framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "FDA_triggers.png", dpi=300, bbox_inches='tight')
    plt.close()

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Execute the complete PyPGx pipeline."""
    
    print("=" * 70)
    print("    PyPGx PHARMACOGENOMICS ANALYSIS PIPELINE")
    print("    No arguments needed - all inputs are pre-configured")
    print("=" * 70)
    
    # Initialize
    global logger
    logger = setup_logging()
    
    logger.info("Starting PyPGx pipeline")
    logger.info(f"Assembly: {ASSEMBLY}")
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    
    # Setup
    ensure_directories()
    check_prerequisites()
    
    # Load inputs
    wgs_samples, array_samples = load_sample_lists()
    genes = load_gene_list()
    metadata_df = load_metadata()
    
    # Save sample lists for reference
    wgs_list_file = BASE_DIR / "wgs_samples.txt"
    array_list_file = BASE_DIR / "array_samples.txt"
    
    with open(wgs_list_file, 'w') as f:
        for sample in wgs_samples:
            f.write(f"{sample}\n")
    logger.info(f"Saved WGS sample list to {wgs_list_file}")
    
    with open(array_list_file, 'w') as f:
        for sample in array_samples:
            f.write(f"{sample}\n")
    logger.info(f"Saved array sample list to {array_list_file}")
    
    # Step 1: Cohort PyPGx run
    print("\n" + "=" * 70)
    print("STEP 1: Running PyPGx on full cohort (736 samples)")
    print("=" * 70)
    
    cohort_results_file = COHORT_DIR / "cohort_pypgx_results.tsv"
    
    if cohort_results_file.exists():
        logger.info(f"Loading existing cohort results from {cohort_results_file}")
        cohort_df = pd.read_csv(cohort_results_file, sep='\t')
    else:
        cohort_df = run_pypgx_batch(JOINT_VCF, genes, "cohort", COHORT_DIR)
        if not cohort_df.empty:
            cohort_df.to_csv(cohort_results_file, sep='\t', index=False)
            logger.info(f"Saved cohort results: {len(cohort_df)} records")
    
    # Step 2: WGS-truth PyPGx run
    print("\n" + "=" * 70)
    print("STEP 2: Running PyPGx on WGS truth set (109 samples)")
    print("=" * 70)
    
    wgs_vcf = BASE_DIR / "wgs_subset.vcf.gz"
    wgs_results_file = WGS_TRUTH_DIR / "wgs_truth_pypgx_results.tsv"
    
    if wgs_results_file.exists():
        logger.info(f"Loading existing WGS results from {wgs_results_file}")
        wgs_df = pd.read_csv(wgs_results_file, sep='\t')
    else:
        # Extract WGS samples
        if not wgs_vcf.exists():
            logger.info("Creating WGS subset VCF...")
            if not extract_samples_from_vcf(JOINT_VCF, wgs_samples, wgs_vcf):
                logger.error("Failed to create WGS subset VCF")
                sys.exit(1)
        
        wgs_df = run_pypgx_batch(wgs_vcf, genes, "wgs_truth", WGS_TRUTH_DIR)
        if not wgs_df.empty:
            wgs_df.to_csv(wgs_results_file, sep='\t', index=False)
            logger.info(f"Saved WGS results: {len(wgs_df)} records")
    
    # Step 3: Concordance analysis
    print("\n" + "=" * 70)
    print("STEP 3: Computing concordance between array and WGS")
    print("=" * 70)
    
    # Filter cohort to WGS samples only
    cohort_wgs = cohort_df[cohort_df["SampleID"].isin(wgs_samples)]
    concordance_df = compute_concordance(cohort_wgs, wgs_df)
    
    # Step 4: Gene tiering
    print("\n" + "=" * 70)
    print("STEP 4: Classifying genes by platform support")
    print("=" * 70)
    
    classification_df = classify_genes(concordance_df)
    
    # Step 5: Population summaries
    print("\n" + "=" * 70)
    print("STEP 5: Creating population-level summaries")
    print("=" * 70)
    
    summary_df = create_population_summaries(cohort_df, classification_df, metadata_df)
    
    # Step 6: FDA mapping
    print("\n" + "=" * 70)
    print("STEP 6: Mapping to FDA recommendations")
    print("=" * 70)
    
    matches_df, trigger_df = map_to_fda_recommendations(cohort_df, classification_df, metadata_df)
    
    # Step 7: Create figures
    print("\n" + "=" * 70)
    print("STEP 7: Creating publication-quality figures")
    print("=" * 70)
    
    create_concordance_figure(concordance_df)
    logger.info("Created concordance figure")
    
    if not trigger_df.empty:
        create_fda_trigger_figure(trigger_df)
        logger.info("Created FDA triggers figure")
    
    # Final summary
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE!")
    print("=" * 70)
    
    print("\nOutputs created in ANALYSIS/00-05-PyPGx/:")
    print("\n  Results:")
    print("    - cohort/cohort_pypgx_results.tsv")
    print("    - wgs_truth/wgs_truth_pypgx_results.tsv")
    print("    - results/PGx_concordance.csv")
    print("    - results/PGx_classification.csv")
    print("    - results/PGx_population_summary.csv")
    print("    - results/FDA_patient_matches.csv")
    print("    - results/FDA_trigger_counts.csv")
    print("    - results/confmat_<GENE>.png")
    print("    - results/discordant_calls_<GENE>.csv")
    
    print("\n  Figures:")
    print("    - figures/concordance_by_gene.png")
    print("    - figures/phenotype_frequencies_<GENE>.png")
    print("    - figures/FDA_triggers.png")
    
    print(f"\nLog file: {LOG_FILE}")
    print("\nAnalysis complete and ready for manuscript!")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nPipeline failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
        