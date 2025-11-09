#!/usr/bin/env python3
"""
Download test sequence data for DAH7PS project
"""

import os
import sys
import urllib.request
from pathlib import Path
import gzip
import shutil

# Progress bar callback
def download_progress(count, block_size, total_size):
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write(f"\r  Progress: {percent}%")
    sys.stdout.flush()

# Test organisms with DAH7PS
TEST_ORGANISMS = {
    "ecoli_k12": {
        "name": "Escherichia coli K-12",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz",
        "description": "Has aroF, aroG, aroH (all three specificities)"
    },
    "bsubtilis": {
        "name": "Bacillus subtilis 168",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000001570/UP000001570_224308.fasta.gz",
        "description": "Firmicutes representative"
    },
    "paeruginosa": {
        "name": "Pseudomonas aeruginosa PAO1",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000002438/UP000002438_208964.fasta.gz",
        "description": "Gram-negative pathogen"
    },
    "mtb": {
        "name": "Mycobacterium tuberculosis H37Rv",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000001584/UP000001584_83332.fasta.gz",
        "description": "Actinobacteria, medical importance"
    },
    "scerevisiae": {
        "name": "Saccharomyces cerevisiae S288C",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz",
        "description": "Fungi, bifunctional DAH7PS-CM (ARO3, ARO4)"
    },
    "athaliana": {
        "name": "Arabidopsis thaliana",
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006548/UP000006548_3702.fasta.gz",
        "description": "Plant, plastid-targeted enzymes"
    },
}

def main():
    print("=" * 70)
    print("DAH7PS Test Dataset Download")
    print("=" * 70)
    print()

    # Get project root
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    raw_dir = project_dir / "data" / "raw"
    proteomes_dir = raw_dir / "proteomes"

    # Create directories
    proteomes_dir.mkdir(parents=True, exist_ok=True)

    print(f"Download directory: {proteomes_dir}")
    print(f"Downloading {len(TEST_ORGANISMS)} proteomes...")
    print()

    downloaded_files = []

    for org_id, info in TEST_ORGANISMS.items():
        print(f"Downloading: {info['name']}")
        print(f"  {info['description']}")

        gz_file = proteomes_dir / f"{org_id}.faa.gz"
        faa_file = proteomes_dir / f"{org_id}.faa"

        try:
            # Download
            urllib.request.urlretrieve(info['url'], gz_file, download_progress)
            print()  # New line after progress

            # Decompress
            print(f"  Decompressing...")
            with gzip.open(gz_file, 'rb') as f_in:
                with open(faa_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove compressed file
            gz_file.unlink()

            # Count sequences
            with open(faa_file, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))

            print(f"  ✓ {seq_count} sequences downloaded")
            downloaded_files.append(faa_file)

        except Exception as e:
            print(f"  ✗ Error: {e}")
            continue

        print()

    if downloaded_files:
        # Combine all files
        print("Combining all proteomes...")
        combined_file = raw_dir / "test_proteomes.faa"

        with open(combined_file, 'w') as outfile:
            for faa_file in downloaded_files:
                with open(faa_file, 'r') as infile:
                    outfile.write(infile.read())

        # Count total sequences
        with open(combined_file, 'r') as f:
            total_seqs = sum(1 for line in f if line.startswith('>'))

        file_size = combined_file.stat().st_size / (1024 * 1024)  # MB

        print()
        print("=" * 70)
        print("✓ Download complete!")
        print("=" * 70)
        print(f"Combined file: {combined_file}")
        print(f"Total sequences: {total_seqs}")
        print(f"File size: {file_size:.1f} MB")
        print()
        print("Next steps:")
        print("  1. Update workflow/config.yaml:")
        print(f"     uniprot_ref_proteomes: \"data/raw/test_proteomes.faa\"")
        print("  2. Run: snakemake --cores 4 sequence_collection")
        print()

        return 0
    else:
        print("✗ No files downloaded successfully")
        return 1

if __name__ == "__main__":
    sys.exit(main())
