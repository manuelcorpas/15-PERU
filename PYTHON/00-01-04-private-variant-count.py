#!/usr/bin/env python
"""
Figure 2B: Private variant counts per population
OPTIMIZED VERSION for high-core systems (M3 Ultra: 32 cores, 256GB RAM)

Optimizations:
1. Uses cyvcf2 instead of pysam (faster C-based VCF parsing)
2. Multiprocessing across chromosomes
3. Numpy vectorized genotype operations
4. Pre-computed sample-to-group index arrays
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
import time

# Try cyvcf2 first (faster), fall back to pysam
try:
    from cyvcf2 import VCF
    USE_CYVCF2 = True
    print("Using cyvcf2 (optimized)")
except ImportError:
    import pysam
    USE_CYVCF2 = False
    print("cyvcf2 not found, using pysam (slower)")

# Define sample groups
GROUP_NAMES = ["CHOPCCAS", "CUSCO", "IQUITOS", "MATZES", "MOCHES", "TRUJILLO", "UROS"]

# Path to the VCF file
vcf_file = "ANALYSIS/00-01-GEN-DIV/Peru.joint.biallelic_snps.vcf.gz"

# Output files
output_csv = "ANALYSIS/00-01-GEN-DIV/private_variant_counts.csv"
output_png = "ANALYSIS/00-01-GEN-DIV/Figure2B_private_variants.png"
output_pdf = "ANALYSIS/00-01-GEN-DIV/Figure2B_private_variants.pdf"

# Number of cores to use (leave 2 for system)
N_CORES = max(1, cpu_count() - 2)


def get_sample_group_indices(samples):
    """
    Pre-compute sample-to-group mapping as numpy array for vectorized operations.
    Returns: group_indices array where group_indices[i] = group_id for sample i
    """
    group_indices = np.full(len(samples), -1, dtype=np.int8)
    sample_counts = {g: 0 for g in GROUP_NAMES}
    
    for i, sample in enumerate(samples):
        for g_idx, group in enumerate(GROUP_NAMES):
            if sample.startswith(group):
                group_indices[i] = g_idx
                sample_counts[group] += 1
                break
    
    return group_indices, sample_counts


def process_chunk_cyvcf2(args):
    """
    Process a chromosome/region using cyvcf2.
    Returns dict of private variant counts per group.
    """
    vcf_path, region, group_indices, n_groups = args
    
    vcf = VCF(vcf_path)
    private_counts = np.zeros(n_groups, dtype=np.int64)
    variant_count = 0
    
    # Iterate over variants in this region
    for variant in vcf(region):
        variant_count += 1
        
        # Get genotypes as numpy array: shape (n_samples, 3) -> [allele1, allele2, phased]
        gt = variant.gt_types  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
        
        # Sample has alt if gt_type is 1 (HET) or 2 (HOM_ALT)
        has_alt = (gt == 1) | (gt == 2)
        
        # Find which groups have at least one sample with alt
        groups_with_alt = set()
        for g_idx in range(n_groups):
            mask = group_indices == g_idx
            if np.any(has_alt & mask):
                groups_with_alt.add(g_idx)
        
        # If exactly one group has alt, it's private
        if len(groups_with_alt) == 1:
            private_counts[list(groups_with_alt)[0]] += 1
    
    vcf.close()
    return private_counts, variant_count


def process_chunk_pysam(args):
    """
    Process a chromosome/region using pysam (fallback).
    """
    vcf_path, region, group_indices, n_groups = args
    
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    private_counts = np.zeros(n_groups, dtype=np.int64)
    variant_count = 0
    
    try:
        iterator = vcf.fetch(region=region) if region else vcf.fetch()
    except:
        iterator = vcf.fetch()
    
    for record in iterator:
        variant_count += 1
        
        # Check which groups have alt alleles
        groups_with_alt = set()
        for i, sample in enumerate(samples):
            gt = record.samples[sample]["GT"]
            if gt is not None and 1 in gt:
                g_idx = group_indices[i]
                if g_idx >= 0:
                    groups_with_alt.add(g_idx)
        
        if len(groups_with_alt) == 1:
            private_counts[list(groups_with_alt)[0]] += 1
    
    vcf.close()
    return private_counts, variant_count


def get_chromosomes(vcf_path):
    """Get list of chromosomes/contigs from VCF."""
    if USE_CYVCF2:
        vcf = VCF(vcf_path)
        chroms = vcf.seqnames
        vcf.close()
    else:
        vcf = pysam.VariantFile(vcf_path)
        chroms = list(vcf.header.contigs)
        vcf.close()
    return chroms


def process_vcf_parallel(vcf_path, n_cores=None):
    """
    Process VCF in parallel across chromosomes.
    """
    if n_cores is None:
        n_cores = N_CORES
    
    print(f"Using {n_cores} cores for parallel processing")
    
    # Get sample info
    if USE_CYVCF2:
        vcf = VCF(vcf_path)
        samples = vcf.samples
        vcf.close()
    else:
        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
        vcf.close()
    
    # Pre-compute group indices
    group_indices, sample_counts = get_sample_group_indices(samples)
    n_groups = len(GROUP_NAMES)
    
    print(f"Found {len(samples)} samples across {n_groups} groups")
    for g, c in sample_counts.items():
        print(f"  {g}: {c} samples")
    
    # Get chromosomes for parallel processing
    chromosomes = get_chromosomes(vcf_path)
    print(f"Processing {len(chromosomes)} chromosomes/contigs...")
    
    # Prepare arguments for each chunk
    process_func = process_chunk_cyvcf2 if USE_CYVCF2 else process_chunk_pysam
    args_list = [(vcf_path, chrom, group_indices, n_groups) for chrom in chromosomes]
    
    # Process in parallel
    start_time = time.time()
    total_counts = np.zeros(n_groups, dtype=np.int64)
    total_variants = 0
    
    with Pool(n_cores) as pool:
        results = pool.map(process_func, args_list)
    
    # Aggregate results
    for counts, n_variants in results:
        total_counts += counts
        total_variants += n_variants
    
    elapsed = time.time() - start_time
    print(f"\nProcessed {total_variants:,} variants in {elapsed:.1f} seconds")
    print(f"Rate: {total_variants/elapsed:,.0f} variants/second")
    
    # Convert to dict
    private_variant_counts = {GROUP_NAMES[i]: int(total_counts[i]) for i in range(n_groups)}
    
    return private_variant_counts, sample_counts


def process_vcf_single_pass(vcf_path):
    """
    Single-threaded optimized version (for comparison or if VCF isn't indexed).
    Uses numpy vectorization for speed.
    """
    print("Running single-threaded optimized version...")
    
    if USE_CYVCF2:
        vcf = VCF(vcf_path)
        samples = vcf.samples
    else:
        vcf = pysam.VariantFile(vcf_path)
        samples = list(vcf.header.samples)
    
    group_indices, sample_counts = get_sample_group_indices(samples)
    n_groups = len(GROUP_NAMES)
    
    private_counts = np.zeros(n_groups, dtype=np.int64)
    variant_count = 0
    
    start_time = time.time()
    report_interval = 1_000_000
    
    if USE_CYVCF2:
        for variant in vcf:
            variant_count += 1
            
            gt = variant.gt_types
            has_alt = (gt == 1) | (gt == 2)
            
            groups_with_alt = set()
            for g_idx in range(n_groups):
                if np.any(has_alt & (group_indices == g_idx)):
                    groups_with_alt.add(g_idx)
            
            if len(groups_with_alt) == 1:
                private_counts[list(groups_with_alt)[0]] += 1
            
            if variant_count % report_interval == 0:
                elapsed = time.time() - start_time
                rate = variant_count / elapsed
                print(f"  Processed {variant_count:,} variants ({rate:,.0f}/sec)...")
    else:
        for record in vcf.fetch():
            variant_count += 1
            
            groups_with_alt = set()
            for i, sample in enumerate(samples):
                gt = record.samples[sample]["GT"]
                if gt is not None and 1 in gt:
                    g_idx = group_indices[i]
                    if g_idx >= 0:
                        groups_with_alt.add(g_idx)
            
            if len(groups_with_alt) == 1:
                private_counts[list(groups_with_alt)[0]] += 1
            
            if variant_count % report_interval == 0:
                elapsed = time.time() - start_time
                rate = variant_count / elapsed
                print(f"  Processed {variant_count:,} variants ({rate:,.0f}/sec)...")
    
    vcf.close()
    
    elapsed = time.time() - start_time
    print(f"\nProcessed {variant_count:,} variants in {elapsed:.1f} seconds")
    print(f"Rate: {variant_count/elapsed:,.0f} variants/second")
    
    private_variant_counts = {GROUP_NAMES[i]: int(private_counts[i]) for i in range(n_groups)}
    
    return private_variant_counts, sample_counts


def write_private_variant_counts(private_variant_counts, output_path):
    """Write the private variant counts per group to a CSV file."""
    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Group", "PrivateVariantCount"])
        for group, count in private_variant_counts.items():
            writer.writerow([group, count])
    print(f"Private variant counts written to {output_path}")


def plot_private_variants(private_variant_counts, sample_counts, save_png=None, save_pdf=None):
    """
    Plot bar chart of private variant counts per population.
    Style matched to Figure 2A for visual consistency.
    """
    df = pd.DataFrame([
        {"Group": group, "Private Variants": count, "n": sample_counts.get(group, 0)}
        for group, count in private_variant_counts.items()
    ])
    
    df = df.sort_values("Private Variants", ascending=True).reset_index(drop=True)
    
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sns.barplot(
        data=df,
        x="Group",
        y="Private Variants",
        palette="Set2",
        edgecolor="black",
        linewidth=1.5,
        ax=ax
    )
    
    for i, (idx, row) in enumerate(df.iterrows()):
        ax.annotate(
            f'{row["Private Variants"]:,}',
            xy=(i, row["Private Variants"]),
            xytext=(0, 5),
            textcoords='offset points',
            ha='center',
            va='bottom',
            fontsize=10,
            fontweight='bold'
        )
        ax.annotate(
            f'(n={row["n"]})',
            xy=(i, row["Private Variants"]),
            xytext=(0, 18),
            textcoords='offset points',
            ha='center',
            va='bottom',
            fontsize=9,
            color='dimgray'
        )
    
    ax.set_title("Private Variant Counts by Population (Ordered by Count)", fontsize=14)
    ax.set_xlabel("Population", fontsize=12)
    ax.set_ylabel("Number of Private Variants", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    ymax = df["Private Variants"].max()
    ax.set_ylim(0, ymax * 1.15)
    plt.tight_layout()
    
    if save_png:
        plt.savefig(save_png, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_png}")
    if save_pdf:
        plt.savefig(save_pdf, format='pdf', bbox_inches='tight')
        print(f"Saved: {save_pdf}")
    
    plt.show()
    return fig, ax


if __name__ == "__main__":
    print("=" * 60)
    print("Figure 2B: Private Variant Analysis (Optimized)")
    print(f"System: {cpu_count()} cores available, using {N_CORES}")
    print("=" * 60)
    
    # Check if VCF is indexed (required for parallel processing)
    vcf_path = Path(vcf_file)
    index_exists = vcf_path.with_suffix('.gz.tbi').exists() or vcf_path.with_suffix('.gz.csi').exists()
    
    if index_exists:
        print(f"\nVCF index found - using parallel processing")
        private_variant_counts, sample_counts = process_vcf_parallel(vcf_file, n_cores=N_CORES)
    else:
        print(f"\nNo VCF index found - using single-threaded mode")
        print("TIP: Run 'tabix -p vcf {vcf_file}' to enable parallel processing")
        private_variant_counts, sample_counts = process_vcf_single_pass(vcf_file)
    
    # Print results
    print("\n" + "=" * 60)
    print("Private Variant Counts:")
    for group, count in sorted(private_variant_counts.items(), key=lambda x: x[1]):
        print(f"  {group}: {count:,}")
    print("=" * 60)
    
    # Write CSV
    write_private_variant_counts(private_variant_counts, output_csv)
    
    # Plot
    print("\nGenerating Figure 2B...")
    plot_private_variants(
        private_variant_counts,
        sample_counts,
        save_png=output_png,
        save_pdf=output_pdf
    )
    
    print("\nDone!")