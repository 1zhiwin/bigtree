#!/usr/bin/env python3
"""
Download production-scale sequence data for DAH7PS project
Downloads reference proteomes from diverse organisms across all domains
"""

import os
import sys
import urllib.request
import urllib.error
from pathlib import Path
import gzip
import shutil
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Target organisms for comprehensive sampling
PRODUCTION_ORGANISMS = {
    # Key Bacteria - Proteobacteria
    "ecoli_k12": {
        "taxid": "83333",
        "name": "Escherichia coli K-12",
        "group": "Gammaproteobacteria"
    },
    "salmonella": {
        "taxid": "90371",
        "name": "Salmonella typhimurium LT2",
        "group": "Gammaproteobacteria"
    },
    "pseudomonas_aeruginosa": {
        "taxid": "208964", 
        "name": "Pseudomonas aeruginosa PAO1",
        "group": "Gammaproteobacteria"
    },
    "pseudomonas_putida": {
        "taxid": "160488",
        "name": "Pseudomonas putida KT2440",
        "group": "Gammaproteobacteria"
    },
    "vibrio_cholerae": {
        "taxid": "243277",
        "name": "Vibrio cholerae O1",
        "group": "Gammaproteobacteria"
    },
    "rhodobacter": {
        "taxid": "272943",
        "name": "Rhodobacter sphaeroides",
        "group": "Alphaproteobacteria"
    },
    "agrobacterium": {
        "taxid": "176299",
        "name": "Agrobacterium fabrum C58",
        "group": "Alphaproteobacteria"
    },
    "helicobacter": {
        "taxid": "85962",
        "name": "Helicobacter pylori",
        "group": "Epsilonproteobacteria"
    },
    
    # Firmicutes
    "bacillus_subtilis": {
        "taxid": "224308",
        "name": "Bacillus subtilis 168",
        "group": "Firmicutes"
    },
    "bacillus_anthracis": {
        "taxid": "198094",
        "name": "Bacillus anthracis Ames",
        "group": "Firmicutes"
    },
    "listeria": {
        "taxid": "169963",
        "name": "Listeria monocytogenes EGD-e",
        "group": "Firmicutes"
    },
    "staphylococcus": {
        "taxid": "93061",
        "name": "Staphylococcus aureus NCTC 8325",
        "group": "Firmicutes"
    },
    "streptococcus": {
        "taxid": "208435",
        "name": "Streptococcus agalactiae 2603V/R",
        "group": "Firmicutes"
    },
    "clostridium": {
        "taxid": "272563",
        "name": "Clostridium acetobutylicum ATCC 824",
        "group": "Firmicutes"
    },
    
    # Actinobacteria
    "mycobacterium_tb": {
        "taxid": "83332",
        "name": "Mycobacterium tuberculosis H37Rv",
        "group": "Actinobacteria"
    },
    "mycobacterium_smeg": {
        "taxid": "246196",
        "name": "Mycobacterium smegmatis MC2 155",
        "group": "Actinobacteria"
    },
    "corynebacterium": {
        "taxid": "196627",
        "name": "Corynebacterium glutamicum ATCC 13032",
        "group": "Actinobacteria"
    },
    "streptomyces": {
        "taxid": "100226",
        "name": "Streptomyces coelicolor A3(2)",
        "group": "Actinobacteria"
    },
    
    # Cyanobacteria
    "synechocystis": {
        "taxid": "1111708",
        "name": "Synechocystis sp. PCC 6803",
        "group": "Cyanobacteria"
    },
    "nostoc": {
        "taxid": "103690",
        "name": "Nostoc sp. PCC 7120",
        "group": "Cyanobacteria"
    },
    
    # Other bacteria
    "thermus": {
        "taxid": "300852",
        "name": "Thermus thermophilus HB8",
        "group": "Deinococcus-Thermus"
    },
    "chlamydia": {
        "taxid": "272561",
        "name": "Chlamydia trachomatis D/UW-3/CX",
        "group": "Chlamydiae"
    },
    
    # Archaea
    "methanococcus": {
        "taxid": "267377",
        "name": "Methanocaldococcus jannaschii DSM 2661",
        "group": "Archaea"
    },
    "pyrococcus": {
        "taxid": "53953",
        "name": "Pyrococcus furiosus DSM 3638",
        "group": "Archaea"
    },
    "sulfolobus": {
        "taxid": "273057",
        "name": "Sulfolobus acidocaldarius DSM 639",
        "group": "Archaea"
    },
    "halobacterium": {
        "taxid": "64091",
        "name": "Halobacterium salinarum R1",
        "group": "Archaea"
    },
    
    # Fungi
    "saccharomyces": {
        "taxid": "559292",
        "name": "Saccharomyces cerevisiae S288C",
        "group": "Fungi"
    },
    "candida": {
        "taxid": "237561",
        "name": "Candida albicans SC5314",
        "group": "Fungi"
    },
    "aspergillus_nidulans": {
        "taxid": "162425",
        "name": "Aspergillus nidulans FGSC A4",
        "group": "Fungi"
    },
    "aspergillus_fumigatus": {
        "taxid": "746128",
        "name": "Aspergillus fumigatus Af293",
        "group": "Fungi"
    },
    "neurospora": {
        "taxid": "367110",
        "name": "Neurospora crassa OR74A",
        "group": "Fungi"
    },
    "schizosaccharomyces": {
        "taxid": "284812",
        "name": "Schizosaccharomyces pombe 972h-",
        "group": "Fungi"
    },
    
    # Plants
    "arabidopsis": {
        "taxid": "3702",
        "name": "Arabidopsis thaliana",
        "group": "Plants"
    },
    "rice": {
        "taxid": "39947",
        "name": "Oryza sativa Japonica",
        "group": "Plants"
    },
    "maize": {
        "taxid": "4577",
        "name": "Zea mays",
        "group": "Plants"
    },
    "tomato": {
        "taxid": "4081",
        "name": "Solanum lycopersicum",
        "group": "Plants"
    },
    "medicago": {
        "taxid": "3880",
        "name": "Medicago truncatula",
        "group": "Plants"
    },
    "physcomitrella": {
        "taxid": "3218",
        "name": "Physcomitrella patens",
        "group": "Plants"
    },
    "chlamydomonas": {
        "taxid": "3055",
        "name": "Chlamydomonas reinhardtii",
        "group": "Green algae"
    }
}

def get_uniprot_proteome_url(taxid):
    """Get UniProt reference proteome URL for a given taxonomy ID"""
    # This would normally query UniProt API, but for simplicity using direct URLs
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    query = f"taxonomy_id:{taxid}+AND+proteome:*"
    params = "format=fasta&compressed=true"
    return f"{base_url}?query={query}&{params}"

def download_proteome(org_id, org_info, output_dir):
    """Download a single proteome"""
    try:
        url = get_uniprot_proteome_url(org_info['taxid'])
        output_file = output_dir / f"{org_id}.faa.gz"
        output_fasta = output_dir / f"{org_id}.faa"
        
        # Skip if already downloaded
        if output_fasta.exists():
            with open(output_fasta, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
            return org_id, seq_count, "exists"
        
        print(f"  Downloading {org_info['name']}...")
        
        # Download with retry logic
        max_retries = 3
        for attempt in range(max_retries):
            try:
                urllib.request.urlretrieve(url, output_file)
                break
            except urllib.error.URLError as e:
                if attempt < max_retries - 1:
                    time.sleep(5 * (attempt + 1))  # Exponential backoff
                    continue
                else:
                    raise
        
        # Decompress
        with gzip.open(output_file, 'rb') as f_in:
            with open(output_fasta, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Remove compressed file
        output_file.unlink()
        
        # Count sequences
        with open(output_fasta, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        
        return org_id, seq_count, "downloaded"
        
    except Exception as e:
        print(f"    Error downloading {org_info['name']}: {e}")
        return org_id, 0, f"error: {e}"

def download_ncbi_datasets():
    """Alternative: Download from NCBI using datasets CLI if available"""
    print("\nNote: For production use, consider using NCBI datasets CLI:")
    print("  curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'")
    print("  chmod +x datasets")
    print("  ./datasets download genome taxon bacteria --reference --include protein")
    print("  This would give you all reference bacterial proteomes at once.")

def main():
    print("=" * 70)
    print("DAH7PS Production Dataset Download")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    raw_dir = project_dir / 'data' / 'raw'
    proteomes_dir = raw_dir / 'production_proteomes'
    proteomes_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Target organisms: {len(PRODUCTION_ORGANISMS)}")
    print(f"Output directory: {proteomes_dir}")
    print()
    
    # Group organisms by domain
    groups = {}
    for org_id, info in PRODUCTION_ORGANISMS.items():
        group = info['group']
        if group not in groups:
            groups[group] = []
        groups[group].append((org_id, info))
    
    print("Organism distribution:")
    for group, orgs in sorted(groups.items()):
        print(f"  {group}: {len(orgs)} organisms")
    print()
    
    # Download proteomes with parallel processing
    print("Starting parallel downloads (max 5 concurrent)...")
    results = {}
    
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {}
        for org_id, org_info in PRODUCTION_ORGANISMS.items():
            future = executor.submit(download_proteome, org_id, org_info, proteomes_dir)
            futures[future] = org_id
        
        for future in as_completed(futures):
            org_id, seq_count, status = future.result()
            results[org_id] = {'sequences': seq_count, 'status': status}
            if status == "downloaded":
                print(f"  ✓ {org_id}: {seq_count} sequences")
            elif status == "exists":
                print(f"  • {org_id}: {seq_count} sequences (already exists)")
    
    # Combine all proteomes
    print("\nCombining all proteomes...")
    combined_file = raw_dir / 'production_proteomes.faa'
    total_sequences = 0
    
    with open(combined_file, 'w') as outfile:
        for org_id in PRODUCTION_ORGANISMS:
            fasta_file = proteomes_dir / f"{org_id}.faa"
            if fasta_file.exists():
                with open(fasta_file, 'r') as infile:
                    content = infile.read()
                    outfile.write(content)
                    total_sequences += content.count('>')
    
    # Save download summary
    summary = {
        'date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'organisms': len(PRODUCTION_ORGANISMS),
        'total_sequences': total_sequences,
        'results': results
    }
    
    summary_file = proteomes_dir / 'download_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 70)
    print("Download Summary:")
    print("=" * 70)
    
    successful = sum(1 for r in results.values() if r['sequences'] > 0)
    total_seqs = sum(r['sequences'] for r in results.values())
    
    print(f"Successfully downloaded: {successful}/{len(PRODUCTION_ORGANISMS)} organisms")
    print(f"Total sequences: {total_sequences:,}")
    print(f"Combined file: {combined_file}")
    print(f"File size: {combined_file.stat().st_size / (1024*1024):.1f} MB")
    print()
    
    # Group statistics
    print("Sequences by group:")
    for group, orgs in sorted(groups.items()):
        group_seqs = sum(results.get(org_id, {}).get('sequences', 0) 
                        for org_id, _ in orgs)
        if group_seqs > 0:
            print(f"  {group}: {group_seqs:,} sequences")
    
    print("\n✓ Production dataset ready for analysis!")
    print("\nNext step: Run production_pipeline.py")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
