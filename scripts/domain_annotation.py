#!/usr/bin/env python3
"""
Phase 3: Domain Annotation for DAH7PS Project
Annotate regulatory domains (ACT, CM) and create comprehensive metadata
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import tempfile
import json

# Configuration
CONFIG = {
    'pfam_hmm': '/home/luogy/bigtree/data/raw/Pfam-A.hmm',
    'sequences': '/home/luogy/bigtree/seqs/dah7ps_nonredundant.faa',
    'act_domain': 'ACT',  # PF01842
    'cm_domains': ['Chorismate_mut', 'CM_2'],  # PF01817, PF04715
    'evalue_threshold': 1.0e-5,
    'coverage_threshold': 0.70,
    'min_domain_separation': 20,
    'cpu': 4
}

def run_hmmscan(sequences_file, hmm_file, output_dir, cpu=4):
    """Run hmmscan to identify all domains in sequences"""
    print("Running hmmscan for domain annotation...")
    
    # Run hmmscan
    out_file = output_dir / "hmmscan_domains.txt"
    domtbl_file = output_dir / "hmmscan_domains.domtbl"
    
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/hmmscan',
        '--domtblout', str(domtbl_file),
        '-o', str(out_file),
        '--cpu', str(cpu),
        '-E', str(CONFIG['evalue_threshold']),
        str(hmm_file),
        str(sequences_file)
    ]
    
    subprocess.run(cmd, check=True)
    
    # Parse domain table
    domains = []
    with open(domtbl_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 22:
                domain = {
                    'target_name': parts[0],  # Domain name
                    'target_acc': parts[1],   # Domain accession
                    'query_name': parts[3],   # Sequence name
                    'evalue': float(parts[6]),
                    'score': float(parts[7]),
                    'bias': float(parts[8]),
                    'domain_num': int(parts[9]),
                    'total_domains': int(parts[10]),
                    'hmm_from': int(parts[15]),
                    'hmm_to': int(parts[16]),
                    'ali_from': int(parts[17]),
                    'ali_to': int(parts[18]),
                    'env_from': int(parts[19]),
                    'env_to': int(parts[20]),
                    'accuracy': float(parts[21])
                }
                domains.append(domain)
    
    print(f"  Found {len(domains)} domain hits")
    return pd.DataFrame(domains)

def classify_domain_architecture(seq_domains):
    """Classify sequence based on domain architecture"""
    has_dahp = False
    has_act = 0
    has_cm = False
    
    dahp_coords = None
    act_coords = []
    cm_coords = None
    
    for _, domain in seq_domains.iterrows():
        domain_name = domain['target_name']
        
        # Check for DAHP synthase domain
        if 'DAHP_synth' in domain_name:
            has_dahp = True
            dahp_coords = (domain['ali_from'], domain['ali_to'])
        
        # Check for ACT domain
        elif 'ACT' in domain_name:
            has_act += 1
            act_coords.append((domain['ali_from'], domain['ali_to']))
        
        # Check for chorismate mutase domain
        elif any(cm in domain_name for cm in ['Chorismate_mut', 'CM_']):
            has_cm = True
            cm_coords = (domain['ali_from'], domain['ali_to'])
    
    # Determine architecture
    if not has_dahp:
        return 'no_dahp', {}
    
    architecture = 'core_only'
    domain_positions = {'dahp': dahp_coords}
    
    if has_act > 0 and has_cm:
        architecture = 'core+ACT+CM'
        domain_positions['act'] = act_coords
        domain_positions['cm'] = cm_coords
    elif has_act == 2:
        architecture = 'core+ACT2'
        domain_positions['act'] = act_coords
    elif has_act == 1:
        architecture = 'core+ACT'
        domain_positions['act'] = act_coords
    elif has_cm:
        architecture = 'core+CM'
        domain_positions['cm'] = cm_coords
    
    return architecture, domain_positions

def annotate_sequences(sequences_file, domains_df):
    """Annotate sequences with domain architecture"""
    print("\nAnnotating sequences with domain architecture...")
    
    # Load sequences
    sequences = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        sequences[record.id] = {
            'id': record.id,
            'description': record.description,
            'length': len(record.seq),
            'sequence': str(record.seq)
        }
    
    # Group domains by sequence
    annotated = []
    for seq_id, seq_info in sequences.items():
        seq_domains = domains_df[domains_df['query_name'] == seq_id]
        
        architecture, positions = classify_domain_architecture(seq_domains)
        
        # Count domains
        n_act = len([d for _, d in seq_domains.iterrows() if 'ACT' in d['target_name']])
        n_cm = len([d for _, d in seq_domains.iterrows() if any(cm in d['target_name'] for cm in ['Chorismate_mut', 'CM_'])])
        
        # Get all domain names
        all_domains = seq_domains['target_name'].tolist()
        
        annotated.append({
            'seq_id': seq_id,
            'description': seq_info['description'],
            'length': seq_info['length'],
            'architecture': architecture,
            'n_act_domains': n_act,
            'n_cm_domains': n_cm,
            'domain_positions': json.dumps(positions),
            'all_domains': ';'.join(all_domains) if all_domains else 'none',
            'n_total_domains': len(seq_domains)
        })
    
    return pd.DataFrame(annotated)

def add_phenotype_annotations(metadata_df):
    """Add known phenotype annotations based on organism and gene"""
    print("\nAdding phenotype annotations...")
    
    # Known phenotypes for specific organisms/genes (simplified)
    phenotypes = {
        'Escherichia coli': {
            'aroF': {'effector': 'Tyr', 'inhibition_mode': 'single'},
            'aroG': {'effector': 'Phe', 'inhibition_mode': 'single'},
            'aroH': {'effector': 'Trp', 'inhibition_mode': 'single'},
        },
        'Saccharomyces cerevisiae': {
            'ARO3': {'effector': 'Phe', 'inhibition_mode': 'single'},
            'ARO4': {'effector': 'Tyr', 'inhibition_mode': 'single'},
        },
        'Pseudomonas': {
            'default': {'effector': 'Multiple', 'inhibition_mode': 'synergistic'},
        },
        'Bacillus': {
            'default': {'effector': 'Prephenate/Chorismate', 'inhibition_mode': 'synergistic'},
        },
        'Arabidopsis': {
            'default': {'effector': 'None', 'inhibition_mode': 'none'},  # Plastid-targeted, often unregulated
        }
    }
    
    # Add phenotype columns
    metadata_df['effector'] = 'Unknown'
    metadata_df['inhibition_mode'] = 'Unknown'
    metadata_df['regulatory_type'] = 'Unknown'
    
    for idx, row in metadata_df.iterrows():
        desc = row['description']
        
        # Extract organism info
        organism = None
        gene_name = None
        
        # Parse description for organism and gene
        if 'Escherichia coli' in desc:
            organism = 'Escherichia coli'
            if 'aroF' in desc or 'AROF' in desc:
                gene_name = 'aroF'
            elif 'aroG' in desc or 'AROG' in desc:
                gene_name = 'aroG'
            elif 'aroH' in desc or 'AROH' in desc:
                gene_name = 'aroH'
        elif 'Saccharomyces' in desc:
            organism = 'Saccharomyces cerevisiae'
            if 'ARO3' in desc:
                gene_name = 'ARO3'
            elif 'ARO4' in desc:
                gene_name = 'ARO4'
        elif 'Pseudomonas' in desc:
            organism = 'Pseudomonas'
        elif 'Bacillus' in desc:
            organism = 'Bacillus'
        elif 'Arabidopsis' in desc:
            organism = 'Arabidopsis'
        
        # Assign phenotype
        if organism in phenotypes:
            if gene_name and gene_name in phenotypes[organism]:
                pheno = phenotypes[organism][gene_name]
            elif 'default' in phenotypes[organism]:
                pheno = phenotypes[organism]['default']
            else:
                continue
            
            metadata_df.at[idx, 'effector'] = pheno['effector']
            metadata_df.at[idx, 'inhibition_mode'] = pheno['inhibition_mode']
            
            # Determine regulatory type based on architecture
            if row['architecture'] == 'core_only':
                metadata_df.at[idx, 'regulatory_type'] = 'Allosteric_core'
            elif 'ACT' in row['architecture']:
                metadata_df.at[idx, 'regulatory_type'] = 'ACT-mediated'
            elif 'CM' in row['architecture']:
                metadata_df.at[idx, 'regulatory_type'] = 'CM-mediated'
    
    return metadata_df

def main():
    print("=" * 70)
    print("DAH7PS Domain Annotation (Phase 3)")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    seqs_dir = project_dir / 'seqs'
    domains_dir = seqs_dir / 'domains'
    data_dir = project_dir / 'data' / 'processed'
    
    domains_dir.mkdir(exist_ok=True)
    
    # Run hmmscan
    domains_df = run_hmmscan(
        CONFIG['sequences'],
        CONFIG['pfam_hmm'],
        domains_dir,
        CONFIG['cpu']
    )
    
    # Save raw domain results
    domains_file = domains_dir / 'all_domains.tsv'
    domains_df.to_csv(domains_file, sep='\t', index=False)
    print(f"\nSaved domain scan results to {domains_file}")
    
    # Annotate sequences with architecture
    annotations_df = annotate_sequences(CONFIG['sequences'], domains_df)
    
    # Load existing metadata
    existing_metadata_file = data_dir / 'sequence_metadata.tsv'
    if existing_metadata_file.exists():
        existing_df = pd.read_csv(existing_metadata_file, sep='\t')
        # Merge with annotations
        metadata_df = existing_df.merge(
            annotations_df[['seq_id', 'architecture', 'n_act_domains', 'n_cm_domains', 
                          'domain_positions', 'all_domains', 'n_total_domains']],
            on='seq_id',
            how='left'
        )
    else:
        metadata_df = annotations_df
    
    # Add phenotype annotations
    metadata_df = add_phenotype_annotations(metadata_df)
    
    # Save enhanced metadata
    enhanced_metadata_file = data_dir / 'sequence_metadata_annotated.tsv'
    metadata_df.to_csv(enhanced_metadata_file, sep='\t', index=False)
    print(f"Saved annotated metadata to {enhanced_metadata_file}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Domain Architecture Summary:")
    print("=" * 70)
    arch_counts = metadata_df['architecture'].value_counts()
    for arch, count in arch_counts.items():
        print(f"  {arch}: {count}")
    
    print("\nRegulatory Domain Distribution:")
    print(f"  Sequences with ACT domain(s): {(metadata_df['n_act_domains'] > 0).sum()}")
    print(f"  Sequences with CM domain: {(metadata_df['n_cm_domains'] > 0).sum()}")
    print(f"  Sequences with both ACT and CM: {((metadata_df['n_act_domains'] > 0) & (metadata_df['n_cm_domains'] > 0)).sum()}")
    
    print("\nKnown Phenotypes:")
    eff_counts = metadata_df[metadata_df['effector'] != 'Unknown']['effector'].value_counts()
    for eff, count in eff_counts.items():
        print(f"  {eff}: {count}")
    
    print("\nInhibition Modes:")
    inh_counts = metadata_df[metadata_df['inhibition_mode'] != 'Unknown']['inhibition_mode'].value_counts()
    for mode, count in inh_counts.items():
        print(f"  {mode}: {count}")
    
    # Save sequences by architecture
    print("\nSaving sequences by architecture...")
    for arch in metadata_df['architecture'].unique():
        if arch == 'no_dahp':
            continue
        arch_seqs = metadata_df[metadata_df['architecture'] == arch]['seq_id'].tolist()
        if arch_seqs:
            arch_file = seqs_dir / f'dah7ps_{arch.replace("+", "_")}.faa'
            sequences = list(SeqIO.parse(CONFIG['sequences'], 'fasta'))
            filtered = [s for s in sequences if s.id in arch_seqs]
            if filtered:
                SeqIO.write(filtered, arch_file, 'fasta')
                print(f"  Saved {len(filtered)} {arch} sequences")
    
    print("\n✓ Phase 3 complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
