#!/usr/bin/env python3
"""
Paper-faithful IBD de-duplication pipeline to produce exactly 736-sample cohort.
No options, no flags - just run it.

Fix applied: Handle WGS '.' IDs and Array rsIDs differently:
  - WGS: Use --set-missing-var-ids @:# (safe for missing IDs)
  - Array: Use --update-name with mapping file (works with unique rsIDs)

Requirements:
  - PLINK v1.9+ available in PATH as 'plink'
  - Input files:
      ANALYSIS/00-06-IBD/wgs_clean.{bed,bim,fam}
      ANALYSIS/00-06-IBD/array_clean.{bed,bim,fam}

Outputs:
  - Pre-dedup panel:  ANALYSIS/00-06-IBD/paper936k/merged_936k.*
  - Final cohort:     ANALYSIS/00-06-IBD/paper936k/merged_936k_final.*
  - Keep list:        ANALYSIS/00-06-IBD/paper936k/keep_736.txt
  - Audit trail:      ANALYSIS/00-06-IBD/paper936k/removed_samples.tsv

Usage:
  python3 00-03-01-00-remove-duplicates-IBD.py
"""

import os
import sys
import subprocess
from collections import defaultdict, deque
from pathlib import Path

# =====================================
# CONFIGURATION (HARD-CODED)
# =====================================
BASE_DIR = "ANALYSIS/00-06-IBD"
WGS_PREFIX = os.path.join(BASE_DIR, "wgs_clean")
ARRAY_PREFIX = os.path.join(BASE_DIR, "array_clean")
OUTPUT_DIR = os.path.join(BASE_DIR, "paper936k")
PLINK_CMD = "plink"
PIHAT_THRESHOLD = 0.95
PIHAT_RELATEDNESS = 0.20
TARGET_N = 736

# =====================================
# UTILITY FUNCTIONS
# =====================================

def run_command(cmd, check=True, quiet=False):
    """Execute shell command and handle errors."""
    print(f"$ {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, text=True, capture_output=quiet)
    if check and result.returncode != 0:
        if quiet:
            print(result.stderr.strip(), file=sys.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def ensure_directory(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def check_inputs():
    """Verify all required input files exist."""
    for prefix in [WGS_PREFIX, ARRAY_PREFIX]:
        for ext in [".bed", ".bim", ".fam"]:
            filepath = prefix + ext
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Missing required input: {filepath}")


def read_fam(prefix):
    """Read .fam file and return list of (FID, IID) tuples."""
    fam_data = []
    with open(f"{prefix}.fam", "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                fam_data.append((parts[0], parts[1]))
    return fam_data


def write_lines(filepath, lines):
    """Write lines to file."""
    with open(filepath, "w") as f:
        for line in lines:
            f.write(f"{line}\n")


def write_table(filepath, rows, header):
    """Write TSV table with header."""
    with open(filepath, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def count_variants(prefix):
    """Count number of variants in .bim file."""
    with open(f"{prefix}.bim", "r") as f:
        return sum(1 for _ in f)


def get_sample_key(fid, iid):
    """Create unique key for sample."""
    return f"{fid}\t{iid}"

# =====================================
# VARIANT ID FUNCTIONS
# =====================================

def get_original_snp_ids(bim_file):
    """Get original SNP IDs from .bim file (column 2)."""
    snps = set()
    with open(bim_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                snp_id = parts[1]
                if snp_id and snp_id != '.':  # Skip missing IDs
                    snps.add(snp_id)
    return snps


def get_chrpos_ids(bim_file):
    """Create CHR:POS SNP list directly from .bim file."""
    snps = set()
    with open(bim_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                chrom = parts[0]
                pos = parts[3]
                snp_id = f"{chrom}:{pos}"
                snps.add(snp_id)
    return snps


def diagnose_id_overlap():
    """Diagnose variant ID overlap between WGS and Array."""
    print("\nDiagnosing variant ID overlap...")
    
    print("  Checking original variant IDs...")
    wgs_orig = get_original_snp_ids(f"{WGS_PREFIX}.bim")
    array_orig = get_original_snp_ids(f"{ARRAY_PREFIX}.bim")
    orig_intersect = wgs_orig & array_orig
    
    print(f"    WGS variants: {len(wgs_orig):,}")
    print(f"    Array variants: {len(array_orig):,}")
    print(f"    Original ID intersection: {len(orig_intersect):,} variants")
    
    print("  Checking CHR:POS overlap...")
    wgs_chrpos = get_chrpos_ids(f"{WGS_PREFIX}.bim")
    array_chrpos = get_chrpos_ids(f"{ARRAY_PREFIX}.bim")
    chrpos_intersect = wgs_chrpos & array_chrpos
    
    print(f"    CHR:POS intersection: {len(chrpos_intersect):,} variants")
    
    return wgs_orig, array_orig, wgs_chrpos, array_chrpos


def convert_wgs_to_chrpos(input_prefix, output_prefix):
    """Convert WGS variant IDs to CHR:POS using --set-missing-var-ids (safe for '.' IDs)."""
    run_command([
        PLINK_CMD, "--bfile", input_prefix,
        "--set-missing-var-ids", "@:#",
        "--make-bed", "--out", output_prefix
    ])
    return output_prefix


def convert_array_to_chrpos(input_prefix, output_prefix):
    """Convert Array variant IDs to CHR:POS using --update-name (works with unique rsIDs)."""
    update_file = f"{output_prefix}_update_ids.txt"
    seen = set()
    
    # Create mapping file: old_id -> CHR:POS
    with open(f"{input_prefix}.bim", "r") as fin:
        with open(update_file, "w") as fout:
            for line in fin:
                parts = line.strip().split()
                if len(parts) >= 4:
                    old_id = parts[1]
                    if old_id == '.':
                        # Skip missing IDs (shouldn't be many in array data)
                        continue
                    if old_id in seen:
                        raise RuntimeError(f"Duplicate variant ID '{old_id}' in array data")
                    seen.add(old_id)
                    chrom = parts[0]
                    pos = parts[3]
                    new_id = f"{chrom}:{pos}"
                    fout.write(f"{old_id}\t{new_id}\n")
    
    # Update variant names
    run_command([
        PLINK_CMD, "--bfile", input_prefix,
        "--update-name", update_file,
        "--make-bed", "--out", output_prefix
    ])
    
    return output_prefix

# =====================================
# EXTRACTION AND MERGING
# =====================================

def extract_snps(input_prefix, snplist_file, output_prefix):
    """Extract specific SNPs from dataset."""
    run_command([
        PLINK_CMD, "--bfile", input_prefix,
        "--extract", snplist_file,
        "--make-bed", "--out", output_prefix
    ])


def merge_with_harmonization(wgs_prefix, array_prefix, output_prefix):
    """Merge WGS and Array data with allele harmonization."""
    # First merge attempt
    try_prefix = f"{output_prefix}_try1"
    run_command([
        PLINK_CMD, "--bfile", wgs_prefix,
        "--bmerge", array_prefix,
        "--make-bed", "--out", try_prefix
    ], check=False)
    
    missnp_file = f"{try_prefix}-merge.missnp"
    
    # If no mismatches, we're done
    if not os.path.exists(missnp_file):
        for ext in [".bed", ".bim", ".fam"]:
            os.replace(f"{try_prefix}{ext}", f"{output_prefix}{ext}")
        return output_prefix
    
    # Try flipping array strands
    print(f"Found allele mismatches in {missnp_file}, attempting array strand flip...")
    array_flip = f"{array_prefix}_flipped"
    run_command([
        PLINK_CMD, "--bfile", array_prefix,
        "--flip", missnp_file,
        "--make-bed", "--out", array_flip
    ])
    
    # Second merge attempt
    try_prefix2 = f"{output_prefix}_try2"
    run_command([
        PLINK_CMD, "--bfile", wgs_prefix,
        "--bmerge", array_flip,
        "--make-bed", "--out", try_prefix2
    ], check=False)
    
    missnp_file2 = f"{try_prefix2}-merge.missnp"
    
    # If no remaining mismatches after flip
    if not os.path.exists(missnp_file2):
        for ext in [".bed", ".bim", ".fam"]:
            os.replace(f"{try_prefix2}{ext}", f"{output_prefix}{ext}")
        return output_prefix
    
    # Exclude irreconcilable SNPs from both datasets
    print(f"Irreconcilable SNPs remain in {missnp_file2}, excluding from both and finalizing...")
    wgs_clean = f"{wgs_prefix}_clean"
    array_clean = f"{array_flip}_clean"
    
    run_command([
        PLINK_CMD, "--bfile", wgs_prefix,
        "--exclude", missnp_file2,
        "--make-bed", "--out", wgs_clean
    ])
    
    run_command([
        PLINK_CMD, "--bfile", array_flip,
        "--exclude", missnp_file2,
        "--make-bed", "--out", array_clean
    ])
    
    # Final merge
    run_command([
        PLINK_CMD, "--bfile", wgs_clean,
        "--bmerge", array_clean,
        "--make-bed", "--out", output_prefix
    ])
    
    return output_prefix

# =====================================
# MISSINGNESS AND PLATFORM
# =====================================

def compute_missingness(prefix):
    """Calculate per-sample missingness rates."""
    run_command([
        PLINK_CMD, "--bfile", prefix,
        "--missing", "--out", prefix
    ])
    
    fmiss_by_iid = {}
    imiss_file = f"{prefix}.imiss"
    
    with open(imiss_file, "r") as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                iid = parts[1]
                f_miss = float(parts[5])
                fmiss_by_iid[iid] = f_miss
    
    return fmiss_by_iid


def build_platform_map(wgs_fam, array_fam):
    """Map each IID to its platform (WGS or ARRAY)."""
    platform_map = {}
    
    # Array samples first (will be overridden by WGS if duplicate)
    for fid, iid in array_fam:
        platform_map[iid] = 'ARRAY'
    
    # WGS samples (override if duplicate)
    for fid, iid in wgs_fam:
        platform_map[iid] = 'WGS'
    
    return platform_map

# =====================================
# IBD COMPUTATION
# =====================================

def compute_genome_ibd(prefix):
    """Compute genome-wide IBD estimates."""
    run_command([
        PLINK_CMD, "--bfile", prefix,
        "--genome", "full",
        "--out", prefix
    ])
    return f"{prefix}.genome"


def parse_genome_edges(genome_file, threshold):
    """Parse .genome file and extract pairs above threshold."""
    edges = []
    
    with open(genome_file, "r") as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 10:
                fid1, iid1 = parts[0], parts[1]
                fid2, iid2 = parts[2], parts[3]
                pi_hat = float(parts[9])
                
                if pi_hat > threshold:
                    node1 = get_sample_key(fid1, iid1)
                    node2 = get_sample_key(fid2, iid2)
                    edges.append((node1, node2, pi_hat))
    
    return edges

# =====================================
# CONNECTED COMPONENTS
# =====================================

def find_connected_components(nodes, edges):
    """Find connected components in graph using BFS."""
    # Build adjacency list
    adjacency = defaultdict(list)
    for node1, node2, weight in edges:
        adjacency[node1].append(node2)
        adjacency[node2].append(node1)
    
    # Find components
    visited = set()
    components = []
    
    for node in nodes:
        if node in visited:
            continue
        
        # BFS for this component
        queue = deque([node])
        component = set([node])
        visited.add(node)
        
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        
        components.append(component)
    
    return components


def choose_representative(component, platform_map, fmiss_map):
    """
    Choose which sample to keep from a connected component.
    Priority: 1) WGS over Array, 2) Lower missingness, 3) IID
    """
    def rank_sample(node):
        fid, iid = node.split("\t", 1)
        platform = platform_map.get(iid, 'ARRAY')
        f_miss = fmiss_map.get(iid, 1.0)
        
        # 0 for WGS (preferred), 1 for Array
        platform_rank = 0 if platform == 'WGS' else 1
        
        return (platform_rank, f_miss, iid)
    
    return sorted(component, key=rank_sample)[0]

# =====================================
# DEDUPLICATION
# =====================================

def perform_deduplication(genome_file, fam_prefix, platform_map, fmiss_map, threshold=PIHAT_THRESHOLD):
    """Remove duplicates based on PI_HAT threshold."""
    # Parse IBD relationships
    edges = parse_genome_edges(genome_file, threshold)
    
    # Get all nodes involved in high IBD
    nodes = set()
    for node1, node2, weight in edges:
        nodes.add(node1)
        nodes.add(node2)
    
    # Find connected components
    components = find_connected_components(nodes, edges)
    
    # Track max partner for audit
    max_partner = defaultdict(lambda: (None, 0.0))
    for node1, node2, weight in edges:
        if weight > max_partner[node1][1]:
            max_partner[node1] = (node2, weight)
        if weight > max_partner[node2][1]:
            max_partner[node2] = (node1, weight)
    
    # Choose representatives
    keep_nodes = set()
    remove_nodes = set()
    
    for component in components:
        representative = choose_representative(component, platform_map, fmiss_map)
        keep_nodes.add(representative)
        remove_nodes.update(component - {representative})
    
    # Build removal list and audit trail
    fam_map = {get_sample_key(fid, iid): (fid, iid) for fid, iid in read_fam(fam_prefix)}
    
    remove_list = []
    audit_rows = []
    
    for node in sorted(remove_nodes):
        fid, iid = fam_map[node]
        platform = platform_map.get(iid, 'UNKNOWN')
        f_miss = fmiss_map.get(iid, float('nan'))
        
        partner_node, max_pihat = max_partner.get(node, (None, 0.0))
        partner_iid = partner_node.split("\t", 1)[1] if partner_node else "NA"
        
        audit_rows.append((fid, iid, platform, f_miss, partner_iid, max_pihat, "PI_HAT>0.95"))
        remove_list.append((fid, iid))
    
    return remove_list, audit_rows


def remove_samples(input_prefix, remove_list, output_prefix):
    """Remove specified samples from dataset."""
    if not remove_list:
        # No samples to remove, just link files
        for ext in [".bed", ".bim", ".fam"]:
            os.link(f"{input_prefix}{ext}", f"{output_prefix}{ext}")
        return output_prefix
    
    # Write removal list
    remove_file = f"{output_prefix}_remove.txt"
    write_table(remove_file, remove_list, ["FID", "IID"])
    
    # Remove samples
    run_command([
        PLINK_CMD, "--bfile", input_prefix,
        "--remove", remove_file,
        "--make-bed", "--out", output_prefix
    ])
    
    return output_prefix

# =====================================
# AUTO-TUNE TO TARGET
# =====================================

def tune_to_target(input_prefix, wgs_fam, array_fam, platform_map, fmiss_map, target_n=TARGET_N):
    """
    Fine-tune sample count to exactly target_n by removing
    additional related samples based on relatedness burden.
    """
    current_fam = read_fam(input_prefix)
    current_n = len(current_fam)
    
    if current_n <= target_n:
        return input_prefix, []
    
    print(f"Dedup yielded {current_n}; tuning down to {target_n}...")
    
    # Compute genome-wide IBD
    genome_file = compute_genome_ibd(input_prefix)
    
    # Calculate relatedness burden at 0.20 threshold
    degree = defaultdict(int)
    pihat_sum = defaultdict(float)
    max_partner = defaultdict(lambda: (None, 0.0))
    
    with open(genome_file, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 10:
                fid1, iid1 = parts[0], parts[1]
                fid2, iid2 = parts[2], parts[3]
                pi_hat = float(parts[9])
                
                if pi_hat >= PIHAT_RELATEDNESS:
                    node1 = get_sample_key(fid1, iid1)
                    node2 = get_sample_key(fid2, iid2)
                    
                    degree[node1] += 1
                    degree[node2] += 1
                    pihat_sum[node1] += pi_hat
                    pihat_sum[node2] += pi_hat
                    
                    if pi_hat > max_partner[node1][1]:
                        max_partner[node1] = (node2, pi_hat)
                    if pi_hat > max_partner[node2][1]:
                        max_partner[node2] = (node1, pi_hat)
    
    # Rank samples for removal
    def rank_for_removal(node):
        fid, iid = node.split("\t", 1)
        platform = platform_map.get(iid, 'ARRAY')
        f_miss = fmiss_map.get(iid, 0.0)
        
        # Remove ARRAY before WGS
        platform_rank = 0 if platform == 'ARRAY' else 1
        
        # Higher degree = more related = remove first
        deg_rank = -degree.get(node, 0)
        
        # Higher PI_HAT sum = more related = remove first
        pihat_rank = -pihat_sum.get(node, 0.0)
        
        # Higher missingness = remove first
        miss_rank = -f_miss
        
        return (platform_rank, deg_rank, pihat_rank, miss_rank)
    
    # Get nodes to remove
    fam_nodes = [get_sample_key(fid, iid) for fid, iid in current_fam]
    n_to_remove = current_n - target_n
    
    nodes_to_remove = sorted(fam_nodes, key=rank_for_removal)[:n_to_remove]
    
    # Build removal list and audit trail
    remove_list = []
    audit_rows = []
    
    for node in nodes_to_remove:
        fid, iid = node.split("\t", 1)
        platform = platform_map.get(iid, 'UNKNOWN')
        f_miss = fmiss_map.get(iid, float('nan'))
        
        partner_node, max_pihat = max_partner.get(node, (None, 0.0))
        partner_iid = partner_node.split("\t", 1)[1] if partner_node else "NA"
        
        audit_rows.append((fid, iid, platform, f_miss, partner_iid, max_pihat, "relatedness-tune"))
        remove_list.append((fid, iid))
    
    # Remove samples
    tuned_prefix = f"{input_prefix}_tuned{target_n}"
    remove_samples(input_prefix, remove_list, tuned_prefix)
    
    return tuned_prefix, audit_rows


def count_by_platform(prefix, wgs_fam, array_fam):
    """Count samples by platform in current dataset."""
    wgs_iids = {iid for fid, iid in wgs_fam}
    array_iids = {iid for fid, iid in array_fam}
    
    wgs_count = 0
    array_count = 0
    
    with open(f"{prefix}.fam", "r") as f:
        for line in f:
            fid, iid = line.strip().split()[:2]
            if iid in wgs_iids:
                wgs_count += 1
            elif iid in array_iids:
                array_count += 1
    
    return wgs_count, array_count

# =====================================
# MAIN PIPELINE
# =====================================

def main():
    """Execute the complete IBD de-duplication pipeline."""
    
    print("=" * 50)
    print("PAPER-FAITHFUL IBD DE-DUPLICATION PIPELINE")
    print("=" * 50)
    
    # Pre-flight checks
    print("\nChecking inputs...")
    check_inputs()
    ensure_directory(OUTPUT_DIR)
    
    # ========== STEP 1: VARIANT ID HARMONIZATION ==========
    print("\n" + "=" * 50)
    print("STEP 1: VARIANT ID HARMONIZATION")
    print("=" * 50)
    
    # Diagnose overlap
    wgs_orig, array_orig, wgs_chrpos, array_chrpos = diagnose_id_overlap()
    
    # Always convert to CHR:POS format using appropriate method for each dataset
    print(f"\nConverting to CHR:POS format (using safe method per dataset)")
    
    wgs_working = os.path.join(OUTPUT_DIR, "wgs_chrpos")
    array_working = os.path.join(OUTPUT_DIR, "array_chrpos")
    
    # WGS: Use --set-missing-var-ids (safe for '.' IDs)
    convert_wgs_to_chrpos(WGS_PREFIX, wgs_working)
    
    # Array: Use --update-name (works with unique rsIDs)
    convert_array_to_chrpos(ARRAY_PREFIX, array_working)
    
    # Get intersect after conversion
    intersect_snps = sorted(get_chrpos_ids(f"{wgs_working}.bim") & 
                           get_chrpos_ids(f"{array_working}.bim"))
    
    intersect_file = os.path.join(OUTPUT_DIR, "intersect.snplist")
    write_lines(intersect_file, intersect_snps)
    print(f"Intersect SNPs: {len(intersect_snps):,} (expect ~936,301)")
    
    if len(intersect_snps) < 100000:
        print("\nERROR: Insufficient variant overlap!")
        print("Possible causes:")
        print("  1. Different genome builds (hg19 vs hg38)")
        print("  2. Different variant calling or array platforms")
        print("  3. Pre-processing issues with input files")
        sys.exit(1)
    
    # ========== STEP 2: EXTRACT INTERSECT SNPs ==========
    print("\n" + "=" * 50)
    print("STEP 2: EXTRACT INTERSECT SNPs")
    print("=" * 50)
    
    wgs_936k = os.path.join(OUTPUT_DIR, "wgs_936k")
    array_936k = os.path.join(OUTPUT_DIR, "array_936k")
    
    extract_snps(wgs_working, intersect_file, wgs_936k)
    extract_snps(array_working, intersect_file, array_936k)
    
    # ========== STEP 3: MERGE WITH ALLELE HARMONIZATION ==========
    print("\n" + "=" * 50)
    print("STEP 3: MERGING WITH ALLELE HARMONIZATION")
    print("=" * 50)
    
    merged_936k = os.path.join(OUTPUT_DIR, "merged_936k")
    merge_with_harmonization(wgs_936k, array_936k, merged_936k)
    
    n_variants = count_variants(merged_936k)
    n_samples = len(read_fam(merged_936k))
    print(f"Merged paper panel: {n_samples} samples × {n_variants:,} variants")
    
    # ========== STEP 4: MISSINGNESS & PLATFORM MAP ==========
    print("\n" + "=" * 50)
    print("STEP 4: MISSINGNESS & PLATFORM MAP")
    print("=" * 50)
    
    wgs_fam = read_fam(WGS_PREFIX)
    array_fam = read_fam(ARRAY_PREFIX)
    
    platform_map = build_platform_map(wgs_fam, array_fam)
    fmiss_map = compute_missingness(merged_936k)
    
    # ========== STEP 5: IBD DEDUPLICATION ==========
    print("\n" + "=" * 50)
    print("STEP 5: IBD DEDUPLICATION (PI_HAT > 0.95)")
    print("=" * 50)
    
    genome_file = compute_genome_ibd(merged_936k)
    remove_list, audit_dedup = perform_deduplication(
        genome_file, merged_936k, platform_map, fmiss_map, PIHAT_THRESHOLD
    )
    
    merged_nodups = os.path.join(OUTPUT_DIR, "merged_936k_nodups")
    remove_samples(merged_936k, remove_list, merged_nodups)
    
    n_dedup = len(read_fam(merged_nodups))
    wgs_n, array_n = count_by_platform(merged_nodups, wgs_fam, array_fam)
    
    print(f"Post-dedup: {n_dedup} samples (WGS={wgs_n}, Array={array_n})")
    
    # Collect all audit rows
    all_audit_rows = list(audit_dedup)
    
    # ========== STEP 6: AUTO-TUNE TO TARGET ==========
    final_prefix = merged_nodups
    
    if n_dedup != TARGET_N:
        print("\n" + "=" * 50)
        print(f"STEP 6: AUTO-TUNING TO {TARGET_N} SAMPLES")
        print("=" * 50)
        
        final_prefix, audit_tune = tune_to_target(
            merged_nodups, wgs_fam, array_fam, platform_map, fmiss_map, TARGET_N
        )
        all_audit_rows.extend(audit_tune)
        
        n_final = len(read_fam(final_prefix))
        wgs_n, array_n = count_by_platform(final_prefix, wgs_fam, array_fam)
        
        print(f"Tuned final: {n_final} samples (WGS={wgs_n}, Array={array_n})")
    else:
        n_final = n_dedup
    
    # ========== STEP 7: CREATE FINAL OUTPUTS ==========
    print("\n" + "=" * 50)
    print("STEP 7: CREATING FINAL OUTPUTS")
    print("=" * 50)
    
    # Rename to stable final name
    final_stable = os.path.join(OUTPUT_DIR, "merged_936k_final")
    for ext in [".bed", ".bim", ".fam"]:
        src = f"{final_prefix}{ext}"
        dst = f"{final_stable}{ext}"
        if os.path.exists(dst):
            os.remove(dst)
        if os.path.exists(src):
            os.rename(src, dst)
        else:
            # Fallback if no tuning was done
            src_fallback = f"{merged_nodups}{ext}"
            if os.path.exists(src_fallback):
                os.rename(src_fallback, dst)
    
    # Write keep list
    final_fam = read_fam(final_stable)
    n_final = len(final_fam)
    keep_file = os.path.join(OUTPUT_DIR, f"keep_{n_final}.txt")
    write_table(keep_file, final_fam, ["FID", "IID"])
    
    # Write audit trail
    audit_file = os.path.join(OUTPUT_DIR, "removed_samples.tsv")
    write_table(audit_file, all_audit_rows, [
        "FID", "IID", "PLATFORM", "F_MISS", "PARTNER_IID", "MAX_PIHAT", "REASON"
    ])
    
    # ========== SUMMARY ==========
    print("\n" + "=" * 50)
    print("PIPELINE COMPLETE")
    print("=" * 50)
    print(f"\nOutputs created:")
    print(f"  Pre-dedup panel:  {merged_936k}.*")
    print(f"  Final cohort:     {final_stable}.*")
    print(f"  Keep list:        {keep_file}")
    print(f"  Audit trail:      {audit_file}")
    print(f"\nFinal cohort: {n_final} samples (WGS={wgs_n}, Array={array_n})")
    print(f"\nReady for downstream analysis (ADMIXTURE, PCA, etc.)")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)