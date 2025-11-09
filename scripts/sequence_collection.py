#!/usr/bin/env python3
"""
Phase 2: Sequence Collection for DAH7PS Project
Run HMMER searches and filter sequences
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import tempfile
import shutil

# Configuration
CONFIG = {
    'pfam_hmm': '/home/luogy/bigtree/data/raw/Pfam-A.hmm',
    'proteomes': '/home/luogy/bigtree/data/raw/test_proteomes.faa',
    'type_I_profile': 'DAHP_synth_1',  # PF00793
    'type_II_profile': 'DAHP_synth_2',  # PF01474
    'evalue_threshold': 1.0e-10,
    'min_domain_coverage': 0.85,
    'max_ambiguous_fraction': 0.20,
    'min_length': 200,
    'max_length': 1500,
    'cpu': 4
}

def run_hmmer_search(proteomes_file, hmm_file, profile, output_dir, cpu=4):
    """Run hmmsearch for a specific Pfam profile"""
    print(f"Running HMMER search for {profile}...")
    
    # Extract specific HMM
    temp_hmm = output_dir / f"{profile}.hmm"
    cmd = f"hmmfetch {hmm_file} {profile} > {temp_hmm}"
    subprocess.run(cmd, shell=True, check=True)
    
    # Run hmmsearch
    out_file = output_dir / f"hmmsearch_{profile}.txt"
    domtbl_file = output_dir / f"hmmsearch_{profile}.domtbl"
    
    cmd = [
        'hmmsearch',
        '--domtblout', str(domtbl_file),
        '-o', str(out_file),
        '--cpu', str(cpu),
        '--cut_ga',  # Use gathering cutoffs
        str(temp_hmm),
        str(proteomes_file)
    ]
    
    subprocess.run(cmd, check=True)
    
    # Parse domain table
    hits = []
    with open(domtbl_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 22:
                hit = {
                    'target_name': parts[0],
                    'target_acc': parts[1],
                    'query_name': parts[3],
                    'query_acc': parts[4],
                    'evalue': float(parts[6]),
                    'domain_score': float(parts[7]),
                    'ali_from': int(parts[17]),
                    'ali_to': int(parts[18]),
                    'hmm_from': int(parts[15]),
                    'hmm_to': int(parts[16]),
                    'env_from': int(parts[19]),
                    'env_to': int(parts[20])
                }
                hits.append(hit)
    
    print(f"  Found {len(hits)} hits for {profile}")
    return pd.DataFrame(hits)

def filter_sequences(hits_df, sequences_dict, profile_type):
    """Apply quality filters to HMMER hits"""
    filtered = []
    
    for _, hit in hits_df.iterrows():
        seq_id = hit['target_name']
        if seq_id not in sequences_dict:
            continue
            
        seq_record = sequences_dict[seq_id]
        seq_len = len(seq_record.seq)
        
        # Check length filters
        if seq_len < CONFIG['min_length'] or seq_len > CONFIG['max_length']:
            continue
        
        # Check ambiguous residues
        seq_str = str(seq_record.seq).upper()
        x_count = seq_str.count('X')
        if x_count / seq_len > CONFIG['max_ambiguous_fraction']:
            continue
        
        # Check domain coverage
        domain_len = hit['ali_to'] - hit['ali_from'] + 1
        hmm_len = hit['hmm_to'] - hit['hmm_from'] + 1
        # Get expected domain length (varies by type)
        expected_len = 250 if profile_type == 'Type_I' else 200
        if domain_len / expected_len < CONFIG['min_domain_coverage']:
            continue
        
        # Add to filtered list
        filtered.append({
            'seq_id': seq_id,
            'description': seq_record.description,
            'organism': extract_organism(seq_record.description),
            'length': seq_len,
            'type': profile_type,
            'evalue': hit['evalue'],
            'domain_score': hit['domain_score'],
            'domain_start': hit['ali_from'],
            'domain_end': hit['ali_to'],
            'sequence': str(seq_record.seq)
        })
    
    print(f"  Filtered to {len(filtered)} sequences for {profile_type}")
    return filtered

def extract_organism(description):
    """Extract organism name from sequence description"""
    # Common patterns in UniProt headers
    if 'OS=' in description:
        start = description.index('OS=') + 3
        end = description.index(' OX=') if ' OX=' in description else len(description)
        return description[start:end].strip()
    elif '|' in description:
        parts = description.split('|')
        if len(parts) > 2:
            # Extract from typical UniProt format
            return parts[2].split()[1:3]  # Get genus and species
    return 'Unknown'

def cluster_sequences(sequences, identity=0.90):
    """Cluster sequences using CD-HIT or MMseqs2"""
    print(f"Clustering sequences at {identity*100}% identity...")
    
    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as f:
        for seq_info in sequences:
            f.write(f">{seq_info['seq_id']}\n{seq_info['sequence']}\n")
        temp_input = f.name
    
    # Output files
    temp_output = temp_input.replace('.faa', '_nr.faa')
    cluster_file = temp_output + '.clstr'  # CD-HIT adds .clstr to output name
    
    # Run CD-HIT
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/cd-hit',
        '-i', temp_input,
        '-o', temp_output,
        '-c', str(identity),
        '-n', '5',  # Word size
        '-M', '4000',  # Memory limit (MB)
        '-T', str(CONFIG['cpu'])
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"CD-HIT error: {e.stderr}")
        raise
    
    # Parse results
    representatives = set()
    with open(temp_output, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            representatives.add(record.id)
    
    # Parse cluster mapping
    cluster_map = {}
    current_cluster = None
    with open(cluster_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
            else:
                seq_id = line.split('>')[1].split('...')[0]
                cluster_map[seq_id] = current_cluster
    
    # Clean up temp files
    os.remove(temp_input)
    os.remove(temp_output)
    os.remove(cluster_file)
    
    # Filter to representatives
    nr_sequences = [seq for seq in sequences if seq['seq_id'] in representatives]
    print(f"  Reduced to {len(nr_sequences)} representative sequences")
    
    return nr_sequences, cluster_map

def classify_sequences(sequences):
    """Classify sequences into Type Iα, Iβ, and II"""
    classified = {
        'Type_Ialpha': [],
        'Type_Ibeta': [],
        'Type_II': []
    }
    
    for seq in sequences:
        if seq['type'] == 'Type_II':
            classified['Type_II'].append(seq)
        else:
            # Type I - further classify based on characteristics
            # This is simplified - real classification would use more features
            seq_str = seq['sequence'].upper()
            
            # Look for key motifs (simplified)
            if 'KPR' in seq_str[100:200]:  # Example motif
                classified['Type_Ialpha'].append(seq)
            else:
                classified['Type_Ibeta'].append(seq)
    
    for class_name, seqs in classified.items():
        print(f"  {class_name}: {len(seqs)} sequences")
    
    return classified

def main():
    print("=" * 70)
    print("DAH7PS Sequence Collection (Phase 2)")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    seqs_dir = project_dir / 'seqs'
    data_dir = project_dir / 'data' / 'processed'
    
    # Load proteomes
    print("Loading proteome sequences...")
    sequences_dict = {}
    with open(CONFIG['proteomes'], 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            sequences_dict[record.id] = record
    print(f"  Loaded {len(sequences_dict)} sequences")
    print()
    
    # Run HMMER searches
    all_sequences = []
    
    # Search for Type I (PF00793)
    hits_type_I = run_hmmer_search(
        CONFIG['proteomes'],
        CONFIG['pfam_hmm'],
        CONFIG['type_I_profile'],
        seqs_dir,
        CONFIG['cpu']
    )
    
    if not hits_type_I.empty:
        filtered_I = filter_sequences(hits_type_I, sequences_dict, 'Type_I')
        all_sequences.extend(filtered_I)
    
    print()
    
    # Search for Type II (PF01474)
    hits_type_II = run_hmmer_search(
        CONFIG['proteomes'],
        CONFIG['pfam_hmm'],
        CONFIG['type_II_profile'],
        seqs_dir,
        CONFIG['cpu']
    )
    
    if not hits_type_II.empty:
        filtered_II = filter_sequences(hits_type_II, sequences_dict, 'Type_II')
        all_sequences.extend(filtered_II)
    
    print()
    print(f"Total sequences after filtering: {len(all_sequences)}")
    print()
    
    # Cluster sequences
    nr_sequences, cluster_map = cluster_sequences(all_sequences)
    
    # Classify sequences
    print("\nClassifying sequences...")
    classified = classify_sequences(nr_sequences)
    
    # Save results
    print("\nSaving results...")
    
    # Save all nonredundant sequences
    all_nr_file = seqs_dir / 'dah7ps_nonredundant.faa'
    with open(all_nr_file, 'w') as f:
        for seq in nr_sequences:
            f.write(f">{seq['seq_id']} {seq['description']}\n")
            f.write(f"{seq['sequence']}\n")
    
    # Save classified sequences
    for class_name, seqs in classified.items():
        class_file = seqs_dir / f'dah7ps_{class_name.lower()}.faa'
        with open(class_file, 'w') as f:
            for seq in seqs:
                f.write(f">{seq['seq_id']} {seq['description']}\n")
                f.write(f"{seq['sequence']}\n")
    
    # Create metadata table
    metadata = pd.DataFrame(nr_sequences)
    metadata['cluster_id'] = metadata['seq_id'].map(cluster_map)
    
    # Add classification
    metadata['classification'] = 'Unknown'
    for class_name, seqs in classified.items():
        seq_ids = [s['seq_id'] for s in seqs]
        metadata.loc[metadata['seq_id'].isin(seq_ids), 'classification'] = class_name
    
    # Save metadata
    metadata_file = data_dir / 'sequence_metadata.tsv'
    metadata.to_csv(metadata_file, sep='\t', index=False)
    
    print(f"  Saved {len(nr_sequences)} nonredundant sequences to {all_nr_file}")
    print(f"  Saved metadata to {metadata_file}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Sequence Collection Summary:")
    print("=" * 70)
    print(f"Total unique sequences: {len(nr_sequences)}")
    print(f"Type Iα: {len(classified['Type_Ialpha'])}")
    print(f"Type Iβ: {len(classified['Type_Ibeta'])}")
    print(f"Type II: {len(classified['Type_II'])}")
    print(f"Average length: {metadata['length'].mean():.0f} ± {metadata['length'].std():.0f}")
    print(f"Length range: {metadata['length'].min()}-{metadata['length'].max()}")
    
    # Organism distribution
    print("\nOrganism distribution:")
    org_counts = metadata['organism'].value_counts().head(10)
    for org, count in org_counts.items():
        print(f"  {org}: {count}")
    
    print("\n✓ Phase 2 complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
