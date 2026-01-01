#!/usr/bin/env python3
"""
00-06-00-01-gnomad-annotate.py

Annotate VCF with gnomAD v2.1.1 allele frequencies using bcftools.

This adds gnomAD_AF, gnomAD_AC, gnomAD_AN fields to your VCF.

Features:
    - Parallel chromosome processing (--parallel N)
    - Extracts chromosomes first for speed
    - Auto-detects chromosome naming (chr1 vs 1)
    - Handles chr prefix mismatch with --add-chr

Prerequisites:
    - bcftools and tabix installed
    - gnomAD VCF files downloaded (run 00-06-00-00-gnomad-download.py first)

Usage:
    python 00-06-00-01-gnomad-annotate.py --chr 22 --input-vcf DATA/Peru.joint.vcf.gz
    python 00-06-00-01-gnomad-annotate.py --input-vcf DATA/Peru.joint.vcf.gz --parallel 12

Author: For Nature Health genomics paper revision
Date: December 2025
"""

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_INPUT_VCF = "DATA/Peru.joint.vcf.gz"
DEFAULT_GNOMAD_DIR = "GNOMAD"
DEFAULT_OUTPUT_DIR = "ANALYSIS/00-11-GNOMAD-NOVEL"

ALL_CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X']

# Default parallel workers
DEFAULT_WORKERS = max(1, multiprocessing.cpu_count() // 2)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def run_cmd(cmd, shell=False, check=False):
    """Run command and return result."""
    if shell:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
    return result


def check_tools():
    """Check if required tools are installed."""
    tools = {'bcftools': False, 'tabix': False}
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], capture_output=True, text=True)
            tools[tool] = result.returncode == 0
        except FileNotFoundError:
            tools[tool] = False
    return tools


def detect_chrom_naming(vcf_path):
    """Detect if VCF uses 'chr' prefix."""
    try:
        result = run_cmd(['bcftools', 'view', '-H', vcf_path])
        for line in result.stdout.split('\n')[:10]:
            if line:
                chrom = line.split('\t')[0]
                return chrom.startswith('chr')
        return False
    except:
        return False


def index_vcf(vcf_path, force=False):
    """Index VCF file."""
    tbi_path = f"{vcf_path}.tbi"
    if not force and os.path.exists(tbi_path):
        return True
    try:
        run_cmd(['tabix', '-f', '-p', 'vcf', vcf_path], check=True)
        return True
    except:
        try:
            run_cmd(['bcftools', 'index', '-f', vcf_path], check=True)
            return True
        except:
            return False


def create_header_file(output_dir):
    """Create header file for gnomAD annotations."""
    header_file = os.path.join(output_dir, 'gnomad_header.txt')
    with open(header_file, 'w') as f:
        f.write('##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Allele frequency in gnomAD v2.1.1 genomes">\n')
        f.write('##INFO=<ID=gnomAD_AC,Number=A,Type=Integer,Description="Allele count in gnomAD v2.1.1 genomes">\n')
        f.write('##INFO=<ID=gnomAD_AN,Number=1,Type=Integer,Description="Total alleles in gnomAD v2.1.1 genomes">\n')
        f.write('##INFO=<ID=gnomAD_nhomalt,Number=A,Type=Integer,Description="Number of homozygotes in gnomAD">\n')
    return header_file


# ============================================================================
# CHROMOSOME PROCESSING
# ============================================================================

def process_chromosome(args_tuple):
    """
    Process a single chromosome: extract, annotate, count.
    This function is designed to be called in parallel.
    """
    chrom, input_vcf, gnomad_dir, output_dir, header_file, input_has_chr = args_tuple
    
    gnomad_vcf = os.path.join(gnomad_dir, f"gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz")
    
    if not os.path.exists(gnomad_vcf):
        return {'chrom': chrom, 'status': 'no_gnomad', 'total': 0, 'annotated': 0}
    
    # Create unique temp files for this chromosome
    temp_dir = tempfile.mkdtemp(prefix=f'gnomad_chr{chrom}_')
    temp_input = os.path.join(temp_dir, f'input_chr{chrom}.vcf.gz')
    temp_gnomad = os.path.join(temp_dir, f'gnomad_chr{chrom}.vcf.gz')
    temp_output = os.path.join(temp_dir, f'output_chr{chrom}.vcf.gz')
    
    try:
        # Step 1: Extract chromosome from input VCF
        region = f"chr{chrom}" if input_has_chr else chrom
        cmd = f"bcftools view -r {region} {input_vcf} -Oz -o {temp_input}"
        result = run_cmd(cmd, shell=True)
        if result.returncode != 0:
            return {'chrom': chrom, 'status': 'extract_failed', 'error': result.stderr, 'total': 0, 'annotated': 0}
        
        run_cmd(f"tabix -f -p vcf {temp_input}", shell=True)
        
        # Step 2: If input has chr prefix, rename gnomAD chromosomes to match
        if input_has_chr:
            rename_file = os.path.join(temp_dir, 'rename.txt')
            with open(rename_file, 'w') as f:
                for c in list(range(1, 23)) + ['X', 'Y', 'M', 'MT']:
                    f.write(f"{c} chr{c}\n")
            
            cmd = f"bcftools annotate --rename-chrs {rename_file} -Oz -o {temp_gnomad} {gnomad_vcf}"
            result = run_cmd(cmd, shell=True)
            if result.returncode != 0:
                return {'chrom': chrom, 'status': 'rename_failed', 'error': result.stderr, 'total': 0, 'annotated': 0}
            
            run_cmd(f"tabix -f -p vcf {temp_gnomad}", shell=True)
            gnomad_to_use = temp_gnomad
        else:
            gnomad_to_use = gnomad_vcf
        
        # Step 3: Annotate with gnomAD
        # KEY FIX: Correct syntax is DST:=SRC
        # gnomAD_AF:=AF means "copy AF from gnomAD into gnomAD_AF in output"
        cmd = f"""bcftools annotate \
            -a {gnomad_to_use} \
            -h {header_file} \
            -c 'CHROM,POS,REF,ALT,gnomAD_AF:=AF,gnomAD_AC:=AC,gnomAD_AN:=AN,gnomAD_nhomalt:=nhomalt' \
            -Oz -o {temp_output} \
            {temp_input}"""
        
        result = run_cmd(cmd, shell=True)
        if result.returncode != 0:
            return {'chrom': chrom, 'status': 'annotate_failed', 'error': result.stderr, 'total': 0, 'annotated': 0}
        
        run_cmd(f"tabix -f -p vcf {temp_output}", shell=True)
        
        # Step 4: Count results
        result = run_cmd(f"bcftools view -H {temp_output} | wc -l", shell=True)
        total = int(result.stdout.strip()) if result.stdout.strip() else 0
        
        # Count variants with gnomAD_AF annotation
        result = run_cmd(f"bcftools view -H -i 'INFO/gnomAD_AF!=\".\"' {temp_output} | wc -l", shell=True)
        annotated = int(result.stdout.strip()) if result.stdout.strip() else 0
        
        # Step 5: Move to final output location
        final_output = os.path.join(output_dir, f"Peru.chr{chrom}.gnomad_annotated.vcf.gz")
        shutil.move(temp_output, final_output)
        if os.path.exists(f"{temp_output}.tbi"):
            shutil.move(f"{temp_output}.tbi", f"{final_output}.tbi")
        else:
            run_cmd(f"tabix -f -p vcf {final_output}", shell=True)
        
        return {
            'chrom': chrom,
            'status': 'success',
            'total': total,
            'annotated': annotated,
            'output': final_output
        }
        
    except Exception as e:
        return {'chrom': chrom, 'status': 'error', 'error': str(e), 'total': 0, 'annotated': 0}
    
    finally:
        # Cleanup temp directory
        shutil.rmtree(temp_dir, ignore_errors=True)


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Annotate VCF with gnomAD v2.1.1 allele frequencies',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00-06-00-01-gnomad-annotate.py --chr 22 --input-vcf DATA/Peru.joint.vcf.gz
  python 00-06-00-01-gnomad-annotate.py --input-vcf DATA/Peru.joint.vcf.gz --parallel 12
  python 00-06-00-01-gnomad-annotate.py --add-chr --input-vcf my_chr_prefix.vcf.gz
        """
    )
    
    parser.add_argument('--input-vcf', type=str, default=DEFAULT_INPUT_VCF,
                        help=f'Input VCF file (default: {DEFAULT_INPUT_VCF})')
    parser.add_argument('--gnomad-dir', type=str, default=DEFAULT_GNOMAD_DIR,
                        help=f'Directory containing gnomAD VCF files (default: {DEFAULT_GNOMAD_DIR})')
    parser.add_argument('--output-dir', type=str, default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--chr', type=str, default=None,
                        help='Comma-separated list of chromosomes (default: all available)')
    parser.add_argument('--parallel', type=int, default=1,
                        help=f'Number of parallel workers (default: 1, max recommended: {DEFAULT_WORKERS})')
    parser.add_argument('--add-chr', action='store_true',
                        help='Force adding "chr" prefix to gnomAD chromosomes')
    parser.add_argument('--force', action='store_true',
                        help='Force re-indexing of files')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("gnomAD ANNOTATION")
    print("=" * 70)
    print()
    
    # Check tools
    tools = check_tools()
    if not tools['bcftools']:
        print("❌ bcftools not found. Install with: brew install bcftools")
        sys.exit(1)
    if not tools['tabix']:
        print("❌ tabix not found. Install with: brew install htslib")
        sys.exit(1)
    print("✓ Required tools found (bcftools, tabix)")
    
    # Check input VCF
    if not os.path.exists(args.input_vcf):
        print(f"❌ Input VCF not found: {args.input_vcf}")
        sys.exit(1)
    print(f"✓ Input VCF: {args.input_vcf}")
    
    # Detect chromosome naming
    input_has_chr = detect_chrom_naming(args.input_vcf) or args.add_chr
    print(f"   Chromosome naming: {'chr1, chr2, ...' if input_has_chr else '1, 2, ...'}")
    
    if input_has_chr:
        print(f"   Will rename gnomAD chromosomes to match (22 → chr22)")
    
    # Check gnomAD directory
    if not os.path.exists(args.gnomad_dir):
        print(f"❌ gnomAD directory not found: {args.gnomad_dir}")
        print(f"   Run: python 00-06-00-00-gnomad-download.py --output-dir {args.gnomad_dir}")
        sys.exit(1)
    
    # Find available gnomAD files
    available_chroms = []
    for chrom in ALL_CHROMOSOMES:
        gnomad_file = os.path.join(args.gnomad_dir, f"gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz")
        if os.path.exists(gnomad_file):
            available_chroms.append(chrom)
    
    if not available_chroms:
        print(f"❌ No gnomAD VCF files found in: {args.gnomad_dir}")
        sys.exit(1)
    print(f"✓ gnomAD files found for {len(available_chroms)} chromosomes")
    
    # Parse requested chromosomes
    if args.chr:
        requested = [c.strip().upper().replace('CHR', '') for c in args.chr.split(',')]
        chromosomes = [c for c in requested if c in available_chroms]
        missing = [c for c in requested if c not in available_chroms]
        if missing:
            print(f"⚠️  gnomAD files missing for: {', '.join(missing)}")
    else:
        chromosomes = available_chroms
    
    print(f"✓ Will annotate chromosomes: {', '.join(chromosomes)}")
    print(f"✓ Parallel workers: {args.parallel}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create header file
    header_file = create_header_file(args.output_dir)
    
    print()
    print("🔄 Processing chromosomes...")
    print()
    
    start_time = datetime.now()
    results = []
    
    # Prepare arguments for each chromosome
    task_args = [
        (chrom, args.input_vcf, args.gnomad_dir, args.output_dir, header_file, input_has_chr)
        for chrom in chromosomes
    ]
    
    if args.parallel > 1 and len(chromosomes) > 1:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {executor.submit(process_chromosome, task): task[0] for task in task_args}
            
            completed = 0
            for future in as_completed(futures):
                completed += 1
                result = future.result()
                results.append(result)
                chrom = result['chrom']
                
                if result['status'] == 'success':
                    pct = 100 * result['annotated'] / result['total'] if result['total'] > 0 else 0
                    print(f"   [{completed}/{len(chromosomes)}] chr{chrom}: ✓ {result['annotated']:,}/{result['total']:,} annotated ({pct:.1f}%)")
                else:
                    print(f"   [{completed}/{len(chromosomes)}] chr{chrom}: ✗ {result['status']} - {result.get('error', '')[:60]}")
    else:
        # Sequential processing
        for i, task in enumerate(task_args, 1):
            chrom = task[0]
            print(f"   [{i}/{len(chromosomes)}] chr{chrom}...", end=' ', flush=True)
            
            result = process_chromosome(task)
            results.append(result)
            
            if result['status'] == 'success':
                pct = 100 * result['annotated'] / result['total'] if result['total'] > 0 else 0
                print(f"✓ {result['annotated']:,}/{result['total']:,} ({pct:.1f}%)")
            else:
                print(f"✗ {result['status']} - {result.get('error', '')[:40]}")
    
    elapsed = datetime.now() - start_time
    
    # Collect successful results
    successful = [r for r in results if r['status'] == 'success']
    
    if not successful:
        print()
        print("❌ No chromosomes were successfully annotated!")
        print()
        print("Debug: Run the diagnostic script:")
        print(f"   python 00-06-00-debug-gnomad.py {args.input_vcf}")
        sys.exit(1)
    
    # Merge chromosome files if multiple
    if len(successful) > 1:
        print()
        print("📦 Merging chromosome files...")
        
        # Sort by chromosome order
        chrom_order = {c: i for i, c in enumerate(ALL_CHROMOSOMES)}
        successful.sort(key=lambda x: chrom_order.get(x['chrom'], 99))
        
        # Create file list
        list_file = os.path.join(args.output_dir, 'vcf_list.txt')
        with open(list_file, 'w') as f:
            for r in successful:
                f.write(r['output'] + '\n')
        
        final_output = os.path.join(args.output_dir, "Peru.joint.gnomad_annotated.vcf.gz")
        run_cmd(f"bcftools concat -f {list_file} -Oz -o {final_output}", shell=True)
        run_cmd(f"tabix -f -p vcf {final_output}", shell=True)
        
        # Cleanup individual chromosome files
        for r in successful:
            if os.path.exists(r['output']):
                os.remove(r['output'])
            if os.path.exists(f"{r['output']}.tbi"):
                os.remove(f"{r['output']}.tbi")
        os.remove(list_file)
        
        output_file = final_output
    else:
        output_file = successful[0]['output']
    
    # Calculate totals
    total_variants = sum(r['total'] for r in successful)
    total_annotated = sum(r['annotated'] for r in successful)
    pct_annotated = 100 * total_annotated / total_variants if total_variants > 0 else 0
    
    print()
    print("=" * 70)
    print("ANNOTATION COMPLETE")
    print("=" * 70)
    print(f"   Time elapsed:        {elapsed}")
    print(f"   Chromosomes:         {len(successful)}/{len(chromosomes)}")
    print(f"   Total variants:      {total_variants:,}")
    print(f"   With gnomAD AF:      {total_annotated:,} ({pct_annotated:.1f}%)")
    print(f"   Without gnomAD AF:   {total_variants - total_annotated:,}")
    print()
    print(f"   Output: {output_file}")
    print()
    print("Next step:")
    print(f"   python 00-06-00-02-gnomad-analyze.py --input-vcf {output_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()