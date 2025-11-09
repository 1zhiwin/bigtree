#!/usr/bin/env python3
"""
Simplified production pipeline that uses the already filtered sequences
"""

import os
import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import subprocess
import json

def cluster_with_cdhit(input_file, output_file, identity=0.90, cpu=8):
    """Cluster sequences using CD-HIT"""
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/cd-hit',
        '-i', str(input_file),
        '-o', str(output_file),
        '-c', str(identity),
        '-n', '5',
        '-M', '8000',
        '-T', str(cpu)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"CD-HIT error: {result.stderr}")
    
    # Count results
    count = 0
    with open(output_file, 'r') as f:
        count = sum(1 for line in f if line.startswith('>'))
    
    return count

def classify_sequences(input_file, output_dir):
    """Classify sequences into Type Iα, Iβ, and II"""
    output_dir = Path(output_dir)
    
    classified = {
        'Type_Ialpha': [],
        'Type_Ibeta': [],
        'Type_II': []
    }
    
    # Read and classify
    for record in SeqIO.parse(input_file, 'fasta'):
        desc = record.description.lower()
        
        # Check for Type II markers
        if 'dahp_synth_2' in desc or 'type ii' in desc:
            classified['Type_II'].append(record)
        # Type Iα markers (Phe-sensitive)
        elif any(x in desc for x in ['arog', 'phe-sensitive', 'phenylalanine']):
            classified['Type_Ialpha'].append(record)
        # Type Iβ markers (Tyr-sensitive)
        elif any(x in desc for x in ['arof', 'aro4', 'tyr-sensitive', 'tyrosine', 'tyr-inhibited']):
            classified['Type_Ibeta'].append(record)
        # Type I Trp-sensitive
        elif any(x in desc for x in ['aroh', 'trp-sensitive', 'tryptophan']):
            classified['Type_Ibeta'].append(record)  # Group with Iβ
        else:
            # Default classification based on organism
            if 'escherichia' in desc or 'salmonella' in desc:
                classified['Type_Ialpha'].append(record)
            else:
                classified['Type_Ibeta'].append(record)
    
    # Save classified sequences
    results = {}
    for class_name, records in classified.items():
        if records:
            output_file = output_dir / f'{class_name.lower()}.faa'
            SeqIO.write(records, output_file, 'fasta')
            results[class_name] = len(records)
            print(f"  {class_name}: {len(records)} sequences")
    
    # Save combined file
    all_records = []
    for records in classified.values():
        all_records.extend(records)
    
    all_file = output_dir / 'all_classified.faa'
    SeqIO.write(all_records, all_file, 'fasta')
    
    return results

def main():
    print("=" * 70)
    print("Production Pipeline - Clustering and Classification")
    print("=" * 70)
    
    project_dir = Path('/home/luogy/bigtree')
    seqs_dir = project_dir / 'seqs_production'
    
    # Input file (already filtered)
    input_file = seqs_dir / 'dah7ps_hits_filtered.faa'
    
    # Count input
    input_count = sum(1 for line in open(input_file) if line.startswith('>'))
    print(f"\nInput: {input_count} filtered sequences")
    
    # Step 1: Cluster with CD-HIT
    print("\nClustering at 90% identity...")
    clustered_file = seqs_dir / 'dah7ps_clustered.faa'
    cluster_count = cluster_with_cdhit(input_file, clustered_file)
    print(f"  Reduced to {cluster_count} representative sequences")
    
    # Step 2: Classify sequences
    print("\nClassifying sequences...")
    classification_results = classify_sequences(clustered_file, seqs_dir)
    
    # Step 3: Create metadata
    print("\nCreating metadata...")
    metadata = []
    
    for record in SeqIO.parse(clustered_file, 'fasta'):
        # Extract organism
        desc = record.description
        organism = 'Unknown'
        if 'OS=' in desc:
            start = desc.index('OS=') + 3
            end = desc.index(' OX=') if ' OX=' in desc else len(desc)
            organism = desc[start:end]
        
        # Determine classification
        classification = 'Unknown'
        desc_lower = desc.lower()
        if 'arog' in desc_lower or 'phe-sensitive' in desc_lower:
            classification = 'Type_Ialpha'
        elif 'arof' in desc_lower or 'aro4' in desc_lower or 'tyr' in desc_lower:
            classification = 'Type_Ibeta'
        elif 'aroh' in desc_lower or 'trp-sensitive' in desc_lower:
            classification = 'Type_Ibeta'
        elif 'dahp_synth_2' in desc_lower:
            classification = 'Type_II'
        
        metadata.append({
            'seq_id': record.id,
            'description': record.description,
            'organism': organism,
            'classification': classification,
            'length': len(record.seq)
        })
    
    # Save metadata
    metadata_df = pd.DataFrame(metadata)
    metadata_file = project_dir / 'data' / 'processed' / 'production_metadata.tsv'
    metadata_df.to_csv(metadata_file, sep='\t', index=False)
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)
    print(f"Total clustered sequences: {cluster_count}")
    print("\nClassification distribution:")
    for class_name, count in classification_results.items():
        print(f"  {class_name}: {count}")
    
    print("\nOrganism distribution (top 10):")
    org_counts = metadata_df['organism'].value_counts().head(10)
    for org, count in org_counts.items():
        print(f"  {org}: {count}")
    
    print(f"\nMetadata saved to: {metadata_file}")
    print("\n✓ Production data ready for alignment and tree building!")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
