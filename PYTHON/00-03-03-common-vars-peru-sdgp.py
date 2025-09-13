#!/usr/bin/env python3
"""
Peru-SGDP Common Variants Pipeline
Converts Peru PLINK data, ensures 736 samples, normalizes contigs, and finds common variants with SGDP.
Run without arguments - uses default paths.
"""

import subprocess
import sys
import os
import shutil
import tempfile
from pathlib import Path
import multiprocessing
import re


# DEFAULT PATHS
FAM_PATH = 'ANALYSIS/00-06-IBD/paper936k/merged_936k_final.fam'
PERU_PLINK_PREFIX = 'ANALYSIS/00-06-IBD/paper936k/merged_936k_final'  # Peru genotype data (same prefix as FAM)
SGDP_VCF_PATH = 'INPUT/VCF/sgdp.cteam_extended.v4.maf0.1perc.vcf.gz'
OUTPUT_DIR = 'ANALYSIS/00-08-PCA'  # All outputs go here


def check_tools():
    """Check that required tools are available on PATH."""
    required_tools = ['bcftools', 'plink', 'bgzip', 'tabix']
    missing = []
    
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    
    if missing:
        print(f"ERROR: Missing required tools: {', '.join(missing)}", file=sys.stderr)
        print("Please ensure these tools are installed and on your PATH", file=sys.stderr)
        sys.exit(1)
    
    print("✓ All required tools found on PATH")


def run_command(cmd, error_msg, capture_output=False):
    """Run a shell command with error handling."""
    try:
        if capture_output:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            return result.stdout.strip()
        else:
            subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {error_msg}", file=sys.stderr)
        print(f"Command failed: {cmd}", file=sys.stderr)
        if hasattr(e, 'stderr') and e.stderr:
            print(f"Error details: {e.stderr}", file=sys.stderr)
        sys.exit(2)


def load_fam_samples(fam_path):
    """Load sample IIDs from FAM file in order."""
    if not os.path.exists(fam_path):
        print(f"ERROR: FAM file not found: {fam_path}", file=sys.stderr)
        sys.exit(2)
    
    samples = []
    with open(fam_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 2:
                print(f"ERROR: Invalid FAM file format at line {line_num}", file=sys.stderr)
                sys.exit(2)
            iid = parts[1]  # Second column is IID
            samples.append(iid)
    
    if len(samples) != 736:
        print(f"ERROR: Expected 736 samples in FAM file, found {len(samples)}", file=sys.stderr)
        sys.exit(2)
    
    print(f"✓ Loaded {len(samples)} samples from FAM file")
    return samples


def write_sample_list(samples, output_file):
    """Write sample list to file (one per line)."""
    with open(output_file, 'w') as f:
        for sample in samples:
            f.write(f"{sample}\n")


def detect_contig_style(vcf_path):
    """Detect contig naming style of a VCF file."""
    cmd = f"bcftools index -s {vcf_path}"
    output = run_command(cmd, f"Failed to get contig stats from {vcf_path}", capture_output=True)
    
    contigs = set()
    for line in output.split('\n'):
        if line.strip():
            parts = line.split('\t')
            if len(parts) >= 1:
                contigs.add(parts[0])
    
    # Check for chr-prefixed style
    has_chr = any(c.startswith('chr') for c in contigs)
    has_numeric = any(c.isdigit() or c in ['X', 'Y', 'MT', 'M'] for c in contigs)
    has_numeric_codes = any(c in ['23', '24', '25'] for c in contigs)
    
    if has_chr:
        style = 'chr'
    elif has_numeric or has_numeric_codes:
        style = 'numeric'
    else:
        style = 'unknown'
    
    print(f"  Detected contig style: {style}")
    print(f"  Sample contigs: {', '.join(sorted(list(contigs)[:5]))}")
    
    return style, contigs


def make_rename_map(src_style, src_contigs, dst_style, mapping_file):
    """Create a contig renaming map file."""
    mappings = []
    
    if src_style == dst_style:
        return False  # No renaming needed
    
    if src_style == 'numeric' and dst_style == 'chr':
        # Numeric to chr-prefixed
        for i in range(1, 23):
            if str(i) in src_contigs:
                mappings.append(f"{i}\tchr{i}")
        
        # Handle sex chromosomes and mitochondria
        if '23' in src_contigs:
            mappings.append("23\tchrX")
        if '24' in src_contigs:
            mappings.append("24\tchrY")
        if '25' in src_contigs:
            mappings.append("25\tchrM")
        
        if 'X' in src_contigs:
            mappings.append("X\tchrX")
        if 'Y' in src_contigs:
            mappings.append("Y\tchrY")
        if 'MT' in src_contigs:
            mappings.append("MT\tchrM")
        if 'M' in src_contigs:
            mappings.append("M\tchrM")
    
    elif src_style == 'chr' and dst_style == 'numeric':
        # Chr-prefixed to numeric
        for i in range(1, 23):
            if f"chr{i}" in src_contigs:
                mappings.append(f"chr{i}\t{i}")
        
        if 'chrX' in src_contigs:
            mappings.append("chrX\tX")
        if 'chrY' in src_contigs:
            mappings.append("chrY\tY")
        if 'chrM' in src_contigs:
            mappings.append("chrM\tMT")
        if 'chrMT' in src_contigs:
            mappings.append("chrMT\tMT")
    
    if mappings:
        with open(mapping_file, 'w') as f:
            for mapping in mappings:
                f.write(f"{mapping}\n")
        print(f"  Created contig mapping with {len(mappings)} entries")
        return True
    
    return False


def normalize_contigs(peru_vcf, sgdp_vcf, work_dir):
    """Normalize contig names in Peru VCF to match SGDP."""
    print("\n3. Normalizing contig names...")
    
    # Detect contig styles
    print("  Checking SGDP contigs...")
    sgdp_style, sgdp_contigs = detect_contig_style(sgdp_vcf)
    
    print("  Checking Peru contigs...")
    peru_style, peru_contigs = detect_contig_style(peru_vcf)
    
    if peru_style == sgdp_style:
        print(f"✓ Contig styles already match ({sgdp_style})")
        return peru_vcf, sgdp_style, peru_style
    
    # Need to rename Peru contigs to match SGDP
    print(f"  Peru uses '{peru_style}', SGDP uses '{sgdp_style}' → will rename Peru contigs")
    
    mapping_file = os.path.join(work_dir, 'rename_chrs.map')
    if make_rename_map(peru_style, peru_contigs, sgdp_style, mapping_file):
        # Apply renaming
        peru_normalized = peru_vcf.replace('.vcf.gz', '.norm.vcf.gz')
        cmd = f"bcftools annotate --rename-chrs {mapping_file} -Oz -o {peru_normalized} {peru_vcf}"
        run_command(cmd, "Failed to rename contigs")
        
        # Index the normalized VCF
        cmd = f"tabix -p vcf {peru_normalized}"
        run_command(cmd, "Failed to index normalized VCF")
        
        # Copy mapping file to output directory
        output_mapping = os.path.join(OUTPUT_DIR, 'rename_chrs.map')
        shutil.copy(mapping_file, output_mapping)
        print(f"✓ Renamed contigs and saved mapping to {output_mapping}")
        
        # Verify renaming worked
        print("  Verifying normalized contigs...")
        new_style, new_contigs = detect_contig_style(peru_normalized)
        if new_style != sgdp_style:
            print(f"ERROR: Contig normalization failed. Expected {sgdp_style}, got {new_style}", file=sys.stderr)
            sys.exit(2)
        
        return peru_normalized, sgdp_style, new_style
    else:
        print("ERROR: Could not create contig mapping", file=sys.stderr)
        sys.exit(2)


def convert_plink_to_vcf(plink_prefix, output_vcf, sample_file, threads, fam_iids):
    """Convert PLINK files to VCF format with IID sample names."""
    print(f"\nConverting PLINK to VCF format...")
    
    # Check if PLINK files exist
    for ext in ['.bed', '.bim', '.fam']:
        if not os.path.exists(plink_prefix + ext):
            print(f"ERROR: PLINK file not found: {plink_prefix}{ext}", file=sys.stderr)
            sys.exit(2)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    
    # Create temporary VCF
    temp_vcf = output_vcf.replace('.vcf.gz', '.temp.vcf')
    temp_prefix = temp_vcf.replace('.vcf', '')
    
    # Convert PLINK to VCF using vcf-iid to get IID as sample names
    cmd = f"plink --bfile {plink_prefix} --keep-fam {sample_file} --recode vcf-iid --out {temp_prefix}"
    run_command(cmd, "Failed to convert PLINK to VCF")
    
    # Compress with bgzip
    cmd = f"bgzip -c {temp_vcf} > {output_vcf}"
    run_command(cmd, "Failed to compress VCF")
    
    # Index
    cmd = f"tabix -p vcf {output_vcf}"
    run_command(cmd, "Failed to index VCF")
    
    # Clean up temp file
    if os.path.exists(temp_vcf):
        os.remove(temp_vcf)
    
    # Verify sample names match FAM IIDs
    vcf_samples = get_vcf_samples(output_vcf)
    
    if len(vcf_samples) != 736:
        print(f"ERROR: Expected 736 samples in converted VCF, got {len(vcf_samples)}", file=sys.stderr)
        sys.exit(2)
    
    # Check if sample names match exactly
    if vcf_samples != fam_iids:
        print("  Sample names don't match FAM IIDs, reheadering...")
        
        # Create reheader map
        reheader_map = output_vcf.replace('.vcf.gz', '.reheader.txt')
        with open(reheader_map, 'w') as f:
            for old_name, new_name in zip(vcf_samples, fam_iids):
                f.write(f"{old_name}\t{new_name}\n")
        
        # Apply reheader
        reheadered_vcf = output_vcf.replace('.vcf.gz', '.rehead.vcf.gz')
        cmd = f"bcftools reheader -s {reheader_map} -o {reheadered_vcf} {output_vcf}"
        run_command(cmd, "Failed to reheader VCF")
        
        # Index reheadered VCF
        cmd = f"tabix -p vcf {reheadered_vcf}"
        run_command(cmd, "Failed to index reheadered VCF")
        
        # Replace original with reheadered
        os.remove(output_vcf)
        os.rename(reheadered_vcf, output_vcf)
        
        # Verify again
        vcf_samples = get_vcf_samples(output_vcf)
        if vcf_samples != fam_iids:
            print("ERROR: Sample names still don't match after reheadering", file=sys.stderr)
            sys.exit(2)
        
        print("✓ Reheadered VCF to match FAM IIDs")
    
    print(f"✓ Converted PLINK to VCF with {len(vcf_samples)} samples (IID format)")
    return output_vcf


def get_vcf_samples(vcf_path):
    """Get list of samples from VCF file."""
    cmd = f"bcftools query -l {vcf_path}"
    output = run_command(cmd, f"Failed to query samples from {vcf_path}", capture_output=True)
    return output.split('\n') if output else []


def filter_noncanonical_contigs(vcf_path, work_dir):
    """Remove non-canonical contigs from VCF."""
    print("\n  Checking for non-canonical contigs...")
    
    # Get all contigs
    cmd = f"bcftools index -s {vcf_path}"
    output = run_command(cmd, "Failed to get contig stats", capture_output=True)
    
    all_contigs = set()
    for line in output.split('\n'):
        if line.strip():
            parts = line.split('\t')
            if len(parts) >= 1:
                all_contigs.add(parts[0])
    
    # Define canonical contigs
    canonical = set()
    for i in range(1, 23):
        canonical.update([str(i), f"chr{i}"])
    canonical.update(['X', 'Y', 'MT', 'M', 'chrX', 'chrY', 'chrM', 'chrMT'])
    canonical.update(['23', '24', '25'])  # Numeric codes for sex/mito
    
    # Find non-canonical contigs
    noncanonical = all_contigs - canonical
    
    if noncanonical:
        print(f"  Found {len(noncanonical)} non-canonical contigs: {', '.join(sorted(noncanonical))}")
        
        # Create filtered VCF
        filtered_vcf = vcf_path.replace('.vcf.gz', '.filtered.vcf.gz')
        
        # Build exclusion list
        exclude_list = '^' + ','.join(noncanonical)
        cmd = f"bcftools view -t {exclude_list} {vcf_path} -Oz -o {filtered_vcf}"
        run_command(cmd, "Failed to filter non-canonical contigs")
        
        # Index filtered VCF
        cmd = f"tabix -p vcf {filtered_vcf}"
        run_command(cmd, "Failed to index filtered VCF")
        
        print(f"✓ Removed {len(noncanonical)} non-canonical contigs")
        return filtered_vcf
    else:
        print("  No non-canonical contigs found")
        return vcf_path


def intersect_vcfs(sgdp_vcf, peru_vcf, output_dir, threads):
    """Perform correct intersection of SGDP and Peru VCFs."""
    print("\nIntersecting SGDP and Peru VCFs for common variants...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Method A: Two-pass isec with write-which option
    common_sgdp = os.path.join(output_dir, "common.sgdp.vcf.gz")
    common_peru = os.path.join(output_dir, "common.peru.vcf.gz")
    
    # Get common sites as they appear in SGDP (n=2 means sites in both, w1 means write first file's version)
    cmd1 = f"bcftools isec -n=2 -w1 --threads {threads} -Oz -o {common_sgdp} {sgdp_vcf} {peru_vcf}"
    run_command(cmd1, "Failed to extract common variants from SGDP")
    
    # Get common sites as they appear in Peru (w2 means write second file's version)
    cmd2 = f"bcftools isec -n=2 -w2 --threads {threads} -Oz -o {common_peru} {sgdp_vcf} {peru_vcf}"
    run_command(cmd2, "Failed to extract common variants from Peru")
    
    # Index both files
    run_command(f"bcftools index --threads {threads} {common_sgdp}", "Failed to index common SGDP VCF")
    run_command(f"bcftools index --threads {threads} {common_peru}", "Failed to index common Peru VCF")
    
    # Count variants in each
    sgdp_count = run_command(f"bcftools view -H {common_sgdp} | wc -l", 
                             "Failed to count SGDP variants", capture_output=True)
    
    print(f"✓ Common variants extracted: {sgdp_count.strip()} sites")
    
    return common_sgdp, common_peru, int(sgdp_count.strip())


def merge_vcfs(sgdp_vcf, peru_vcf, output_vcf, threads):
    """Merge the common variant VCFs."""
    print("\nMerging SGDP and Peru samples on common variants...")
    
    cmd = f"bcftools merge --threads {threads} {sgdp_vcf} {peru_vcf} -Oz -o {output_vcf}"
    run_command(cmd, "Failed to merge VCFs")
    
    # Index the merged VCF
    run_command(f"bcftools index --threads {threads} {output_vcf}", "Failed to index merged VCF")
    
    # Verify merged samples
    merged_samples = get_vcf_samples(output_vcf)
    print(f"✓ Merged VCF contains {len(merged_samples)} total samples")
    
    return merged_samples


def write_sample_check(output_dir, peru_samples, sgdp_samples, merged_samples,
                      sgdp_style, peru_style_before, peru_style_after, common_sites):
    """Write sample check TSV file with contig style information."""
    check_file = os.path.join(output_dir, "sample_check.tsv")
    
    # Get first 5 samples from each source
    peru_first5 = ','.join(peru_samples[:5])
    sgdp_first5 = ','.join(sgdp_samples[:5]) if len(sgdp_samples) >= 5 else ','.join(sgdp_samples)
    merged_first5 = ','.join(merged_samples[:5])
    
    with open(check_file, 'w') as f:
        f.write("source\tn_samples\tfirst5_ids\n")
        f.write(f"peru_fam\t{len(peru_samples)}\t{peru_first5}\n")
        f.write(f"peru_vcf\t{len(peru_samples)}\t{peru_first5}\n")
        f.write(f"sgdp\t{len(sgdp_samples)}\t{sgdp_first5}\n")
        f.write(f"merged\t{len(merged_samples)}\t{merged_first5}\n")
        f.write("\n")
        f.write("contig_audit\tvalue\n")
        f.write(f"contig_style_sgdp\t{sgdp_style}\n")
        f.write(f"contig_style_peru_before\t{peru_style_before}\n")
        f.write(f"contig_style_peru_after\t{peru_style_after}\n")
        f.write(f"n_common_sites\t{common_sites}\n")
    
    print(f"✓ Wrote sample check to {check_file}")


def verify_peru_in_merged(merged_samples, peru_samples, output_dir):
    """Verify all Peru samples are in the merged VCF."""
    peru_set = set(peru_samples)
    merged_set = set(merged_samples)
    
    peru_in_merged = peru_set & merged_set
    missing_peru = peru_set - merged_set
    
    if len(peru_in_merged) != 736:
        print(f"ERROR: Expected 736 Peru samples in merged VCF, found {len(peru_in_merged)}", file=sys.stderr)
        
        if missing_peru:
            print(f"  Missing {len(missing_peru)} Peru samples:", file=sys.stderr)
            for i, sample in enumerate(sorted(missing_peru)[:10], 1):
                print(f"    {i}. {sample}", file=sys.stderr)
            if len(missing_peru) > 10:
                print(f"    ... and {len(missing_peru)-10} more", file=sys.stderr)
            
            # Write full list of missing samples
            missing_file = os.path.join(output_dir, 'peru_missing_in_merged.txt')
            with open(missing_file, 'w') as f:
                for sample in sorted(missing_peru):
                    f.write(f"{sample}\n")
            print(f"  Full list written to {missing_file}", file=sys.stderr)
        
        sys.exit(2)
    
    print(f"✓ All 736 Peru samples present in merged VCF")
    return True


def print_summary(peru_count, sgdp_count, common_site_count, merged_total, output_dir):
    """Print final summary."""
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Peru samples (from FAM):  {peru_count}")
    print(f"SGDP samples:             {sgdp_count}")
    print(f"Common variant sites:     {common_site_count}")
    print(f"Merged VCF total samples: {merged_total}")
    print("\nOutput files:")
    print(f"  - {os.path.join(output_dir, 'peru.736.vcf.gz')}")
    print(f"  - {os.path.join(output_dir, 'common.sgdp.vcf.gz')}")
    print(f"  - {os.path.join(output_dir, 'common.peru.vcf.gz')}")
    print(f"  - {os.path.join(output_dir, 'common_peru_sgdp.vcf.gz')}")
    print(f"  - {os.path.join(output_dir, 'sample_check.tsv')}")
    print(f"  - {os.path.join(output_dir, 'rename_chrs.map')} (if contigs were renamed)")
    print("="*60)


def main():
    print("Peru-SGDP Common Variants Pipeline")
    print("="*60)
    
    # Use all available CPU threads
    threads = multiprocessing.cpu_count()
    print(f"Using {threads} threads")
    
    # Check tools
    check_tools()
    
    # Check that input files exist
    print("\nChecking input files...")
    
    if not os.path.exists(FAM_PATH):
        print(f"ERROR: FAM file not found: {FAM_PATH}", file=sys.stderr)
        sys.exit(2)
    print(f"✓ Found FAM file: {FAM_PATH}")
    
    # Check for Peru PLINK files
    for ext in ['.bed', '.bim', '.fam']:
        if not os.path.exists(PERU_PLINK_PREFIX + ext):
            print(f"ERROR: Peru PLINK file not found: {PERU_PLINK_PREFIX}{ext}", file=sys.stderr)
            sys.exit(2)
    print(f"✓ Found Peru PLINK files: {PERU_PLINK_PREFIX}.*")
    
    if not os.path.exists(SGDP_VCF_PATH):
        print(f"ERROR: SGDP VCF file not found: {SGDP_VCF_PATH}", file=sys.stderr)
        sys.exit(2)
    print(f"✓ Found SGDP VCF: {SGDP_VCF_PATH}")
    
    # Create work directory
    work_dir = tempfile.mkdtemp(prefix='peru_sgdp_work_')
    print(f"\nWorking directory: {work_dir}")
    
    try:
        # Load FAM samples
        print("\n1. Loading authoritative sample list...")
        peru_samples = load_fam_samples(FAM_PATH)
        
        # Create sample file for PLINK --keep-fam (needs FID and IID)
        fid_iid_file = os.path.join(work_dir, 'peru.keep.txt')
        with open(FAM_PATH, 'r') as f_in, open(fid_iid_file, 'w') as f_out:
            for line in f_in:
                parts = line.strip().split()
                if len(parts) >= 2:
                    f_out.write(f"{parts[0]}\t{parts[1]}\n")  # FID IID
        
        # Convert PLINK to VCF with IID sample names
        print("\n2. Converting Peru PLINK data to VCF...")
        peru_output = os.path.join(OUTPUT_DIR, 'peru.736.vcf.gz')
        peru_vcf = convert_plink_to_vcf(PERU_PLINK_PREFIX, peru_output, fid_iid_file, threads, peru_samples)
        
        # Normalize contig names to match SGDP
        peru_vcf_final, sgdp_style, peru_style_after = normalize_contigs(peru_vcf, SGDP_VCF_PATH, work_dir)
        peru_style_before = peru_style_after if peru_vcf == peru_vcf_final else 'numeric'  # Assume original was numeric if changed
        
        # Filter non-canonical contigs if present
        peru_vcf_final = filter_noncanonical_contigs(peru_vcf_final, work_dir)
        
        # Get SGDP sample info
        print("\n4. Checking SGDP VCF...")
        sgdp_samples = get_vcf_samples(SGDP_VCF_PATH)
        print(f"✓ SGDP VCF contains {len(sgdp_samples)} samples")
        if len(sgdp_samples) > 0:
            print(f"  First 5 samples: {', '.join(sgdp_samples[:5])}")
        
        # Intersect VCFs using the normalized Peru VCF
        print("\n5. Finding common variants...")
        common_sgdp, common_peru, common_sites = intersect_vcfs(SGDP_VCF_PATH, peru_vcf_final, 
                                                                OUTPUT_DIR, threads)
        
        # Merge VCFs
        print("\n6. Merging samples on common variants...")
        merged_vcf = os.path.join(OUTPUT_DIR, 'common_peru_sgdp.vcf.gz')
        merged_samples = merge_vcfs(common_sgdp, common_peru, merged_vcf, threads)
        
        # Verify Peru samples in merged VCF (exact string match)
        verify_peru_in_merged(merged_samples, peru_samples, OUTPUT_DIR)
        
        # Write sample check file with contig audit
        print("\n7. Writing validation files...")
        write_sample_check(OUTPUT_DIR, peru_samples, sgdp_samples, merged_samples,
                          sgdp_style, peru_style_before, peru_style_after, common_sites)
        
        # Write merged sample list for transparency
        merged_sample_file = os.path.join(OUTPUT_DIR, 'merged.samples.txt')
        write_sample_list(merged_samples, merged_sample_file)
        print(f"✓ Wrote merged sample list to {merged_sample_file}")
        
        # Print summary
        print_summary(736, len(sgdp_samples), common_sites, len(merged_samples), OUTPUT_DIR)
        
    finally:
        # Clean up work directory
        shutil.rmtree(work_dir)
        print(f"\n✓ Cleaned up temporary directory")
    
    print("\n✓ Pipeline completed successfully")
    return 0


if __name__ == '__main__':
    sys.exit(main())