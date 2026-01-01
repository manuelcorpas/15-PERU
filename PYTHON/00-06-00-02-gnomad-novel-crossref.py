#!/usr/bin/env python3
"""
00-06-00-03-gnomad-novel-crossref.py

Cross-reference VEP-defined novel variants against gnomAD annotation.

This script properly:
1. Reads VEP SPLIT.ARCHIVE files (matching the 1.6M novel in your summary)
2. Deduplicates by position (VEP outputs multiple lines per variant for each transcript)
3. Cross-references against gnomAD-annotated VCF

Usage:
    python 00-06-00-03-gnomad-novel-crossref.py

Author: For Nature Health genomics paper revision
Date: December 2025
"""

import os
import sys
import glob
import subprocess
import tempfile
from datetime import datetime
from collections import defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================

# Use SPLIT.ARCHIVE - this matches your 1.6M novel summary
VEP_DIR = "ANALYSIS/VEP/TXT/SPLIT.ARCHIVE"
VEP_PATTERN = "split_*.txt"

GNOMAD_VCF = "ANALYSIS/00-11-GNOMAD-NOVEL/Peru.joint.gnomad_annotated.vcf.gz"
OUTPUT_DIR = "ANALYSIS/00-11-GNOMAD-NOVEL"


# ============================================================================
# STEP 1: Parse VEP files and deduplicate
# ============================================================================

def parse_vep_files(vep_dir, pattern):
    """
    Parse VEP output files and deduplicate by position.
    
    VEP outputs multiple lines per variant (one per transcript).
    We deduplicate by (chrom, pos) to get unique variants.
    
    Returns: (novel_positions, existing_positions) as sets of (chrom, pos) tuples
    """
    
    novel_positions = set()
    existing_positions = set()
    
    vep_files = sorted(glob.glob(os.path.join(vep_dir, pattern)))
    
    if not vep_files:
        print(f"❌ No VEP files found: {os.path.join(vep_dir, pattern)}")
        return None, None
    
    print(f"   Found {len(vep_files)} VEP files")
    
    for vep_file in vep_files:
        filename = os.path.basename(vep_file)
        print(f"   {filename}...", end=' ', flush=True)
        
        col_idx = None
        lines_read = 0
        
        with open(vep_file, 'r') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('##'):
                    continue
                
                # Parse header to get column indices
                if line.startswith('#'):
                    cols = line.lstrip('#').strip().split('\t')
                    col_idx = {col: i for i, col in enumerate(cols)}
                    continue
                
                if col_idx is None:
                    continue
                
                lines_read += 1
                parts = line.strip().split('\t')
                
                try:
                    # Get Location column (format: "1:10492-10492" or "1:10492")
                    location_idx = col_idx.get('Location', 1)
                    location = parts[location_idx]
                    
                    if ':' not in location:
                        continue
                    
                    chrom = location.split(':')[0]
                    pos_part = location.split(':')[1]
                    pos = pos_part.split('-')[0]  # Take start position
                    
                    # Get Existing_variation column
                    existing_idx = col_idx.get('Existing_variation', 19)
                    if existing_idx >= len(parts):
                        existing_var = "-"
                    else:
                        existing_var = parts[existing_idx].strip()
                    
                    # Create position key for deduplication
                    pos_key = (chrom, pos)
                    
                    # Classify as novel or existing
                    # Novel = no rsID (Existing_variation is "-" or empty)
                    if existing_var == "-" or existing_var == "" or existing_var == ".":
                        novel_positions.add(pos_key)
                    else:
                        existing_positions.add(pos_key)
                        
                except (IndexError, ValueError):
                    continue
        
        print(f"{lines_read:,} lines")
    
    # Handle positions that appear as both novel and existing
    # (different annotations for same position) - prioritize existing
    novel_only = novel_positions - existing_positions
    
    return novel_only, existing_positions


# ============================================================================
# STEP 2: Query gnomAD VCF for positions
# ============================================================================

def query_gnomad_vcf(gnomad_vcf, positions):
    """
    Query gnomAD-annotated VCF to check which positions have gnomAD_AF.
    
    Returns: (in_gnomad, not_in_gnomad, not_found) as sets
    """
    
    print(f"   Checking {len(positions):,} positions against gnomAD VCF...")
    
    # Group positions by chromosome
    by_chrom = defaultdict(set)
    for chrom, pos in positions:
        by_chrom[chrom].add(int(pos))
    
    # Create regions file for bcftools query
    regions_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    
    # Sort chromosomes: 1-22, then X, Y
    def chrom_sort_key(c):
        if c.isdigit():
            return (0, int(c))
        elif c == 'X':
            return (1, 0)
        elif c == 'Y':
            return (1, 1)
        else:
            return (2, c)
    
    for chrom in sorted(by_chrom.keys(), key=chrom_sort_key):
        for pos in sorted(by_chrom[chrom]):
            regions_file.write(f"{chrom}\t{pos}\n")
    regions_file.close()
    
    # Query VCF
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%INFO/gnomAD_AF\\n' -R {regions_file.name} {gnomad_vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    os.unlink(regions_file.name)
    
    if result.returncode != 0:
        print(f"   ⚠️  bcftools warning: {result.stderr[:100]}")
    
    # Parse results
    in_gnomad = set()
    not_in_gnomad = set()
    
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        
        chrom, pos, gnomad_af = parts[0], parts[1], parts[2]
        pos_key = (chrom, pos)
        
        if pos_key in positions:
            if gnomad_af and gnomad_af != "." and gnomad_af != "":
                in_gnomad.add(pos_key)
            else:
                not_in_gnomad.add(pos_key)
    
    # Positions not found in VCF
    found = in_gnomad | not_in_gnomad
    not_found = positions - found
    
    return in_gnomad, not_in_gnomad, not_found


# ============================================================================
# STEP 3: Generate reports
# ============================================================================

def generate_reports(output_dir, stats):
    """Generate summary and manuscript text files."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate truly novel (not in gnomAD + not found in VCF)
    truly_novel = stats['not_in_gnomad'] + stats['not_found']
    pct_truly_novel = 100 * truly_novel / stats['novel'] if stats['novel'] > 0 else 0
    pct_in_gnomad = 100 - pct_truly_novel
    
    # Summary file
    summary_file = os.path.join(output_dir, 'gnomad_novel_verification_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("gnomAD VERIFICATION OF VEP-DEFINED NOVEL VARIANTS\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("VEP VARIANT COUNTS (DEDUPLICATED)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total unique positions:            {stats['total']:>12,}\n")
        f.write(f"  With rsID (existing):            {stats['existing']:>12,}\n")
        f.write(f"  Without rsID (novel):            {stats['novel']:>12,}\n\n")
        
        f.write("gnomAD STATUS OF NOVEL VARIANTS\n")
        f.write("-" * 50 + "\n")
        f.write(f"Novel variants checked:            {stats['novel']:>12,}\n")
        f.write(f"  Present in gnomAD:               {stats['in_gnomad']:>12,} ({pct_in_gnomad:.1f}%)\n")
        f.write(f"  Absent from gnomAD:              {truly_novel:>12,} ({pct_truly_novel:.1f}%)\n")
        f.write(f"    - In VCF, no gnomAD_AF:        {stats['not_in_gnomad']:>12,}\n")
        f.write(f"    - Not in VCF:                  {stats['not_found']:>12,}\n\n")
        
        f.write("=" * 70 + "\n")
    
    # Manuscript text
    manuscript_file = os.path.join(output_dir, 'manuscript_text_gnomad_novel.txt')
    with open(manuscript_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("MANUSCRIPT TEXT\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("RESULTS SECTION:\n")
        f.write("-" * 50 + "\n")
        f.write(f'"To further characterize the {stats["novel"]:,} variants lacking dbSNP ')
        f.write(f'identifiers, we cross-referenced them against gnomAD v2.1.1. ')
        f.write(f'Of these, {pct_truly_novel:.0f}% ({truly_novel:,}) were also absent from gnomAD, ')
        f.write(f'representing truly novel variation not captured in major population databases. ')
        f.write(f'The remaining {pct_in_gnomad:.0f}% ({stats["in_gnomad"]:,}) were present in gnomAD ')
        f.write(f'but had not been assigned dbSNP rsIDs at the time of annotation."\n\n')
        
        f.write("REVIEWER RESPONSE:\n")
        f.write("-" * 50 + "\n")
        f.write(f'"As requested, we verified our novel variants against gnomAD v2.1.1. ')
        f.write(f'Of the {stats["novel"]:,} variants identified as novel by VEP (lacking dbSNP rsIDs), ')
        f.write(f'{pct_truly_novel:.0f}% ({truly_novel:,}) were also absent from gnomAD, ')
        f.write(f'confirming their novelty. This analysis has been added to the revised manuscript."\n\n')
        
        f.write("=" * 70 + "\n")
    
    return summary_file, manuscript_file, truly_novel, pct_truly_novel


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("gnomAD VERIFICATION OF VEP-DEFINED NOVEL VARIANTS")
    print("=" * 70)
    print()
    
    # Check paths
    if not os.path.exists(VEP_DIR):
        print(f"❌ VEP directory not found: {VEP_DIR}")
        sys.exit(1)
    print(f"✓ VEP directory: {VEP_DIR}")
    
    if not os.path.exists(GNOMAD_VCF):
        print(f"❌ gnomAD VCF not found: {GNOMAD_VCF}")
        sys.exit(1)
    print(f"✓ gnomAD VCF: {GNOMAD_VCF}")
    
    start_time = datetime.now()
    
    # Step 1: Parse VEP files
    print()
    print("🔄 Step 1: Parsing VEP files (deduplicating by position)...")
    print()
    
    novel_positions, existing_positions = parse_vep_files(VEP_DIR, VEP_PATTERN)
    
    if novel_positions is None:
        sys.exit(1)
    
    total_unique = len(novel_positions) + len(existing_positions)
    
    print()
    print(f"   UNIQUE VARIANTS (deduplicated):")
    print(f"   Total:     {total_unique:,}")
    print(f"   Existing:  {len(existing_positions):,} ({100*len(existing_positions)/total_unique:.1f}%)")
    print(f"   Novel:     {len(novel_positions):,} ({100*len(novel_positions)/total_unique:.1f}%)")
    
    # Step 2: Query gnomAD VCF
    print()
    print("🔄 Step 2: Cross-referencing with gnomAD...")
    print()
    
    in_gnomad, not_in_gnomad, not_found = query_gnomad_vcf(GNOMAD_VCF, novel_positions)
    
    # Step 3: Generate reports
    print()
    print("🔄 Step 3: Generating reports...")
    
    stats = {
        'total': total_unique,
        'existing': len(existing_positions),
        'novel': len(novel_positions),
        'in_gnomad': len(in_gnomad),
        'not_in_gnomad': len(not_in_gnomad),
        'not_found': len(not_found)
    }
    
    summary_file, manuscript_file, truly_novel, pct_truly_novel = generate_reports(OUTPUT_DIR, stats)
    
    elapsed = datetime.now() - start_time
    
    # Print results
    print()
    print("=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"   Time: {elapsed}")
    print()
    print(f"📊 KEY FINDINGS:")
    print(f"   Novel variants (VEP):      {len(novel_positions):,}")
    print(f"   • In gnomAD:               {len(in_gnomad):,} ({100-pct_truly_novel:.1f}%)")
    print(f"   • TRULY NOVEL:             {truly_novel:,} ({pct_truly_novel:.1f}%)")
    print(f"     (absent from both dbSNP AND gnomAD)")
    print()
    print(f"📁 Output files:")
    print(f"   • {summary_file}")
    print(f"   • {manuscript_file}")
    print()
    print("=" * 70)
    
    # Print manuscript text
    pct_in_gnomad = 100 - pct_truly_novel
    print()
    print("📝 MANUSCRIPT TEXT:")
    print("-" * 70)
    print(f'"To further characterize the {len(novel_positions):,} variants lacking')
    print(f'dbSNP identifiers, we cross-referenced them against gnomAD v2.1.1.')
    print(f'Of these, {pct_truly_novel:.0f}% ({truly_novel:,}) were also absent from gnomAD,')
    print(f'representing truly novel variation. The remaining {pct_in_gnomad:.0f}%')
    print(f'({len(in_gnomad):,}) were present in gnomAD but lacked dbSNP identifiers."')
    print("-" * 70)


if __name__ == "__main__":
    main()