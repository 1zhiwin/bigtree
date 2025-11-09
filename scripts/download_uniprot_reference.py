#!/usr/bin/env python3
"""
Download UniProt reference proteomes for production-scale DAH7PS analysis
Uses UniProt FTP server for reliable access
"""

import os
import sys
import urllib.request
from pathlib import Path
import gzip
import shutil
import time

# UniProt Reference Proteomes - Well-studied organisms with good annotation
REFERENCE_PROTEOMES = {
    # Bacteria - Proteobacteria
    "UP000000625": {"name": "Escherichia coli K-12", "taxid": "83333"},
    "UP000002438": {"name": "Pseudomonas aeruginosa PAO1", "taxid": "208964"},
    "UP000000556": {"name": "Pseudomonas putida KT2440", "taxid": "160488"},
    "UP000000535": {"name": "Salmonella typhimurium LT2", "taxid": "90371"},
    "UP000000586": {"name": "Shigella flexneri 2a", "taxid": "198214"},
    "UP000002493": {"name": "Rhodobacter sphaeroides", "taxid": "272943"},
    "UP000002671": {"name": "Caulobacter crescentus", "taxid": "190650"},
    "UP000000429": {"name": "Helicobacter pylori", "taxid": "85962"},
    
    # Bacteria - Firmicutes  
    "UP000001570": {"name": "Bacillus subtilis 168", "taxid": "224308"},
    "UP000000594": {"name": "Bacillus anthracis Ames", "taxid": "198094"},
    "UP000000817": {"name": "Listeria monocytogenes", "taxid": "169963"},
    "UP000008816": {"name": "Staphylococcus aureus", "taxid": "93061"},
    "UP000000814": {"name": "Clostridium acetobutylicum", "taxid": "272563"},
    "UP000000253": {"name": "Streptococcus pneumoniae", "taxid": "170187"},
    
    # Bacteria - Actinobacteria
    "UP000001584": {"name": "Mycobacterium tuberculosis H37Rv", "taxid": "83332"},
    "UP000000582": {"name": "Corynebacterium glutamicum", "taxid": "196627"},
    "UP000001973": {"name": "Streptomyces coelicolor", "taxid": "100226"},
    
    # Bacteria - Cyanobacteria
    "UP000001425": {"name": "Synechocystis sp. PCC 6803", "taxid": "1111708"},
    "UP000002279": {"name": "Nostoc sp. PCC 7120", "taxid": "103690"},
    
    # Bacteria - Others
    "UP000000552": {"name": "Thermus thermophilus HB8", "taxid": "300852"},
    "UP000008183": {"name": "Chlamydia trachomatis", "taxid": "272561"},
    
    # Archaea
    "UP000000805": {"name": "Methanocaldococcus jannaschii", "taxid": "267377"},
    "UP000001013": {"name": "Pyrococcus furiosus", "taxid": "53953"},
    "UP000001974": {"name": "Sulfolobus solfataricus", "taxid": "273057"},
    "UP000008546": {"name": "Halobacterium salinarum", "taxid": "64091"},
    
    # Fungi
    "UP000002311": {"name": "Saccharomyces cerevisiae S288C", "taxid": "559292"},
    "UP000000559": {"name": "Candida albicans SC5314", "taxid": "237561"},
    "UP000000560": {"name": "Aspergillus nidulans", "taxid": "162425"},
    "UP000001596": {"name": "Neurospora crassa", "taxid": "367110"},
    "UP000002485": {"name": "Schizosaccharomyces pombe", "taxid": "284812"},
    
    # Plants
    "UP000006548": {"name": "Arabidopsis thaliana", "taxid": "3702"},
    "UP000059680": {"name": "Oryza sativa Japonica", "taxid": "39947"},
    "UP000007305": {"name": "Zea mays", "taxid": "4577"},
    "UP000004994": {"name": "Solanum lycopersicum", "taxid": "4081"},
    "UP000002051": {"name": "Medicago truncatula", "taxid": "3880"},
    "UP000006727": {"name": "Physcomitrella patens", "taxid": "3218"},
    "UP000006906": {"name": "Chlamydomonas reinhardtii", "taxid": "3055"}
}

def download_reference_proteome(proteome_id, info, output_dir):
    """Download a reference proteome from UniProt FTP"""
    # Construct FTP URL based on proteome ID
    # UniProt organizes by taxonomy groups
    taxid = info['taxid']
    name = info['name']
    
    # Determine subdirectory based on taxonomy
    if taxid in ["3702", "39947", "4577", "4081", "3880", "3218", "3055"]:
        subdir = "Eukaryota"
    elif taxid in ["559292", "237561", "162425", "367110", "284812"]:
        subdir = "Eukaryota"
    elif taxid in ["267377", "53953", "273057", "64091"]:
        subdir = "Archaea"
    else:
        subdir = "Bacteria"
    
    url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{subdir}/{proteome_id}/{proteome_id}_{taxid}.fasta.gz"
    
    output_gz = output_dir / f"{proteome_id}.fasta.gz"
    output_fasta = output_dir / f"{proteome_id}.fasta"
    
    # Skip if already downloaded
    if output_fasta.exists():
        with open(output_fasta, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"  • {name}: {seq_count} sequences (exists)")
        return proteome_id, seq_count, "exists"
    
    try:
        print(f"  Downloading {name}...")
        urllib.request.urlretrieve(url, output_gz)
        
        # Decompress
        with gzip.open(output_gz, 'rb') as f_in:
            with open(output_fasta, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        output_gz.unlink()
        
        # Count sequences
        with open(output_fasta, 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        
        print(f"    ✓ {seq_count} sequences")
        return proteome_id, seq_count, "downloaded"
        
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return proteome_id, 0, f"error"

def main():
    print("=" * 70)
    print("UniProt Reference Proteomes Download")
    print("=" * 70)
    print()
    
    project_dir = Path('/home/luogy/bigtree')
    raw_dir = project_dir / 'data' / 'raw'
    proteomes_dir = raw_dir / 'reference_proteomes'
    proteomes_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Downloading {len(REFERENCE_PROTEOMES)} reference proteomes")
    print(f"Output directory: {proteomes_dir}\n")
    
    # Download each proteome
    results = []
    total_sequences = 0
    
    for proteome_id, info in REFERENCE_PROTEOMES.items():
        result = download_reference_proteome(proteome_id, info, proteomes_dir)
        results.append(result)
        total_sequences += result[1]
        time.sleep(0.5)  # Be polite to UniProt servers
    
    # Combine all proteomes
    print("\nCombining all proteomes...")
    combined_file = raw_dir / 'uniprot_reference_proteomes.faa'
    
    with open(combined_file, 'w') as outfile:
        for proteome_id, seq_count, status in results:
            if seq_count > 0:
                fasta_file = proteomes_dir / f"{proteome_id}.fasta"
                if fasta_file.exists():
                    with open(fasta_file, 'r') as infile:
                        # Add proteome ID to headers for tracking
                        for line in infile:
                            if line.startswith('>'):
                                # Add proteome ID to header
                                outfile.write(f">{proteome_id}|{line[1:]}")
                            else:
                                outfile.write(line)
    
    # Summary
    successful = sum(1 for r in results if r[1] > 0)
    
    print("\n" + "=" * 70)
    print("Download Complete!")
    print("=" * 70)
    print(f"Successfully downloaded: {successful}/{len(REFERENCE_PROTEOMES)} proteomes")
    print(f"Total sequences: {total_sequences:,}")
    print(f"Combined file: {combined_file}")
    print(f"File size: {combined_file.stat().st_size / (1024**2):.1f} MB")
    
    # Count by domain
    bacteria = sum(1 for p, i in REFERENCE_PROTEOMES.items() 
                  if not any(x in i['taxid'] for x in ["3702","39947","4577","4081","3880","3218","3055","559292","237561","162425","367110","284812","267377","53953","273057","64091"]))
    archaea = 4
    fungi = 5
    plants = 7
    
    print(f"\nOrganism distribution:")
    print(f"  Bacteria: {bacteria}")
    print(f"  Archaea: {archaea}")
    print(f"  Fungi: {fungi}")
    print(f"  Plants: {plants}")
    
    print("\n✓ Ready for production-scale HMMER search!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
