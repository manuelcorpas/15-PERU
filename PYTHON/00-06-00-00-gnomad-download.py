#!/usr/bin/env python3
"""
00-06-00-gnomad-download.py

Download gnomAD v2.1.1 VCF files (GRCh37/hg19) for variant annotation.

This downloads ~50GB of data. Run once, then use 00-06-00-gnomad-annotate.py

Usage:
    python 00-06-00-gnomad-download.py
    python 00-06-00-gnomad-download.py --chr 22          # Test with chr22 only (~2GB)
    python 00-06-00-gnomad-download.py --chr 1,2,22     # Specific chromosomes
    python 00-06-00-gnomad-download.py --output-dir /path/to/gnomad

Author: For Nature Health genomics paper revision
Date: December 2025
"""

import os
import sys
import argparse
import subprocess
import hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================================================================
# CONFIGURATION
# ============================================================================

# Default output directory (relative to project root)
DEFAULT_OUTPUT_DIR = "GNOMAD"

# gnomAD v2.1.1 base URL (GRCh37/hg19)
GNOMAD_BASE_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes"

# Chromosomes to download
ALL_CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X']

# Expected file sizes (approximate, in GB) for progress estimation
APPROX_SIZES = {
    '1': 5.5, '2': 5.2, '3': 4.3, '4': 4.2, '5': 3.9, '6': 3.8,
    '7': 3.5, '8': 3.3, '9': 2.7, '10': 3.0, '11': 3.0, '12': 2.9,
    '13': 2.1, '14': 2.0, '15': 1.8, '16': 2.0, '17': 1.8, '18': 1.7,
    '19': 1.4, '20': 1.4, '21': 0.9, '22': 0.9, 'X': 2.7
}


# ============================================================================
# FUNCTIONS
# ============================================================================

def check_wget():
    """Check if wget is available."""
    try:
        result = subprocess.run(['wget', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_curl():
    """Check if curl is available."""
    try:
        result = subprocess.run(['curl', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def get_file_size(filepath):
    """Get file size in GB."""
    if os.path.exists(filepath):
        return os.path.getsize(filepath) / (1024**3)
    return 0


def download_file(url, output_path, use_wget=True):
    """Download a file using wget or curl with resume support."""
    
    if use_wget:
        cmd = ['wget', '-c', '-q', '--show-progress', url, '-O', output_path]
    else:
        cmd = ['curl', '-C', '-', '-L', '-o', output_path, url]
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"      ❌ Download failed: {e}")
        return False


def download_chromosome(chrom, output_dir, use_wget=True):
    """Download VCF and index for a single chromosome."""
    
    vcf_filename = f"gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz"
    tbi_filename = f"{vcf_filename}.tbi"
    
    vcf_url = f"{GNOMAD_BASE_URL}/{vcf_filename}"
    tbi_url = f"{GNOMAD_BASE_URL}/{tbi_filename}"
    
    vcf_path = os.path.join(output_dir, vcf_filename)
    tbi_path = os.path.join(output_dir, tbi_filename)
    
    expected_size = APPROX_SIZES.get(chrom, 2.0)
    
    # Check if already downloaded
    if os.path.exists(vcf_path) and os.path.exists(tbi_path):
        current_size = get_file_size(vcf_path)
        if current_size > expected_size * 0.9:  # Allow 10% tolerance
            return {'chrom': chrom, 'status': 'exists', 'size': current_size}
    
    # Download VCF
    print(f"   Downloading chr{chrom} (~{expected_size:.1f}GB)...")
    
    success = download_file(vcf_url, vcf_path, use_wget)
    if not success:
        return {'chrom': chrom, 'status': 'failed', 'size': 0}
    
    # Download index
    print(f"   Downloading chr{chrom} index...")
    success = download_file(tbi_url, tbi_path, use_wget)
    if not success:
        return {'chrom': chrom, 'status': 'failed_index', 'size': get_file_size(vcf_path)}
    
    return {'chrom': chrom, 'status': 'downloaded', 'size': get_file_size(vcf_path)}


def verify_downloads(output_dir, chromosomes):
    """Verify downloaded files."""
    print("\n📋 Verifying downloads...")
    
    results = {'complete': [], 'missing': [], 'incomplete': []}
    
    for chrom in chromosomes:
        vcf_file = os.path.join(output_dir, f"gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz")
        tbi_file = f"{vcf_file}.tbi"
        
        if not os.path.exists(vcf_file):
            results['missing'].append(chrom)
        elif not os.path.exists(tbi_file):
            results['incomplete'].append(chrom)
        else:
            size = get_file_size(vcf_file)
            expected = APPROX_SIZES.get(chrom, 2.0)
            if size < expected * 0.8:  # Less than 80% of expected
                results['incomplete'].append(chrom)
            else:
                results['complete'].append(chrom)
    
    return results


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Download gnomAD v2.1.1 VCF files for annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00-06-00-gnomad-download.py                    # Download all chromosomes
  python 00-06-00-gnomad-download.py --chr 22           # Test with chr22 only
  python 00-06-00-gnomad-download.py --chr 1,2,22,X     # Specific chromosomes
  python 00-06-00-gnomad-download.py --verify           # Just verify existing files
        """
    )
    
    parser.add_argument('--output-dir', type=str, default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory for gnomAD files (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--chr', type=str, default=None,
                        help='Comma-separated list of chromosomes (default: all)')
    parser.add_argument('--verify', action='store_true',
                        help='Only verify existing downloads, do not download')
    parser.add_argument('--use-curl', action='store_true',
                        help='Use curl instead of wget')
    
    args = parser.parse_args()
    
    # Parse chromosomes
    if args.chr:
        chromosomes = [c.strip().upper().replace('CHR', '') for c in args.chr.split(',')]
    else:
        chromosomes = ALL_CHROMOSOMES
    
    # Calculate total size
    total_size = sum(APPROX_SIZES.get(c, 2.0) for c in chromosomes)
    
    print("=" * 70)
    print("gnomAD v2.1.1 DOWNLOAD (GRCh37/hg19)")
    print("=" * 70)
    print()
    print(f"Output directory:  {args.output_dir}")
    print(f"Chromosomes:       {', '.join(chromosomes)}")
    print(f"Estimated size:    ~{total_size:.1f} GB")
    print()
    
    # Check download tools
    use_wget = not args.use_curl and check_wget()
    use_curl = args.use_curl or (not use_wget and check_curl())
    
    if not use_wget and not use_curl:
        print("❌ Neither wget nor curl found. Please install one:")
        print("   sudo apt-get install wget")
        print("   or")
        print("   sudo apt-get install curl")
        sys.exit(1)
    
    print(f"Download tool:     {'wget' if use_wget else 'curl'}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Verify only mode
    if args.verify:
        results = verify_downloads(args.output_dir, chromosomes)
        
        print(f"\n✓ Complete:    {len(results['complete'])} chromosomes")
        print(f"⚠ Incomplete:  {len(results['incomplete'])} chromosomes")
        print(f"✗ Missing:     {len(results['missing'])} chromosomes")
        
        if results['incomplete']:
            print(f"\nIncomplete: {', '.join(results['incomplete'])}")
        if results['missing']:
            print(f"Missing: {', '.join(results['missing'])}")
        
        sys.exit(0 if not results['missing'] and not results['incomplete'] else 1)
    
    # Download
    print("\n📥 Starting downloads...")
    print("   (Use Ctrl+C to interrupt, re-run to resume)\n")
    
    downloaded = 0
    existed = 0
    failed = []
    
    for chrom in chromosomes:
        result = download_chromosome(chrom, args.output_dir, use_wget)
        
        if result['status'] == 'exists':
            print(f"   ✓ chr{chrom}: Already exists ({result['size']:.2f} GB)")
            existed += 1
        elif result['status'] == 'downloaded':
            print(f"   ✓ chr{chrom}: Downloaded ({result['size']:.2f} GB)")
            downloaded += 1
        else:
            print(f"   ✗ chr{chrom}: Failed")
            failed.append(chrom)
    
    # Summary
    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)
    print(f"   Downloaded:     {downloaded}")
    print(f"   Already existed: {existed}")
    print(f"   Failed:         {len(failed)}")
    
    if failed:
        print(f"\n⚠️  Failed chromosomes: {', '.join(failed)}")
        print("   Re-run this script to retry failed downloads.")
    
    # Verify
    results = verify_downloads(args.output_dir, chromosomes)
    
    if results['complete'] == chromosomes or len(results['complete']) == len(chromosomes):
        print(f"\n✓ All {len(chromosomes)} chromosome files ready!")
        print(f"\nNext step:")
        print(f"   python 00-06-00-gnomad-annotate.py --gnomad-dir {args.output_dir}")
    else:
        print(f"\n⚠️  {len(results['missing']) + len(results['incomplete'])} files still needed")
        print("   Re-run this script to complete downloads.")
    
    print("=" * 70)


if __name__ == "__main__":
    main()