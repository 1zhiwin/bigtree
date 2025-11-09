#!/usr/bin/env python3
"""
Phase 6: Ancestral Sequence Reconstruction for DAH7PS Project
Reconstruct ancestral sequences using IQ-TREE ASR
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO, Phylo, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import json

# Configuration
CONFIG = {
    'iqtree_path': '/home/luogy/miniconda3/envs/dah7ps/bin/iqtree',
    'posterior_cutoff': 0.80,  # Minimum PP for consensus
    'cpu': 4
}

def run_iqtree_asr(alignment_file, tree_file, model, prefix, cpu=4):
    """Run IQ-TREE ancestral sequence reconstruction"""
    print(f"  Running IQ-TREE ASR...")
    
    cmd = [
        CONFIG['iqtree_path'],
        '-s', str(alignment_file),
        '-te', str(tree_file),  # Use fixed tree
        '-m', model,  # Use best model from tree inference
        '--ancestral',  # Perform ASR
        '-pre', str(prefix),
        '-nt', str(cpu),
        '-quiet'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    Error: {result.stderr}")
        return False
    
    return True

def parse_asr_state_file(state_file):
    """Parse IQ-TREE state file to extract ancestral sequences"""
    ancestral_seqs = {}
    posterior_probs = {}
    
    try:
        # Read the file as tab-separated
        import pandas as pd
        df = pd.read_csv(state_file, sep='\t', comment='#')
        
        # Group by node
        for node_name in df['Node'].unique():
            node_df = df[df['Node'] == node_name].sort_values('Site')
            
            # Get the state (most probable AA) for each position
            sequence = ''.join(node_df['State'].tolist())
            
            # Get the posterior probability for the chosen state
            pps = []
            for _, row in node_df.iterrows():
                state = row['State']
                # Get the probability for this state
                pp_col = f'p_{state}'
                if pp_col in row:
                    pps.append(row[pp_col])
                else:
                    # If column doesn't exist, get max probability
                    prob_cols = [col for col in row.index if col.startswith('p_')]
                    if prob_cols:
                        max_pp = max([row[col] for col in prob_cols])
                        pps.append(max_pp)
                    else:
                        pps.append(0.0)
            
            ancestral_seqs[node_name] = sequence
            posterior_probs[node_name] = pps
        
    except Exception as e:
        print(f"    Error parsing state file: {e}")
    
    return ancestral_seqs, posterior_probs

def identify_key_nodes(tree_file):
    """Identify key nodes for ASR (LCAs, duplication points, etc.)"""
    key_nodes = {}
    
    try:
        tree = Phylo.read(tree_file, 'newick')
        
        # Find internal nodes
        internal_nodes = []
        for i, clade in enumerate(tree.find_clades()):
            if not clade.is_terminal():
                # Get terminal descendants
                terminals = [t.name for t in clade.get_terminals()]
                
                # Classify based on terminal names
                type_i_alpha = any('AROG' in t or 'aroG' in t or 'Phe' in desc 
                                  for t in terminals for desc in [t])
                type_i_beta = any('AROF' in t or 'aroF' in t or 'Tyr' in desc 
                                 for t in terminals for desc in [t])
                type_ii = any('DAHP_synth_2' in desc for t in terminals for desc in [t])
                
                # Name the node
                node_id = f"Node{i}"
                if len(terminals) == len(tree.get_terminals()):
                    key_nodes[node_id] = {
                        'name': 'Root',
                        'description': 'Root of the tree (all sequences)',
                        'terminals': terminals
                    }
                elif type_i_alpha and not type_i_beta and not type_ii:
                    key_nodes[node_id] = {
                        'name': 'LCA_Type_Ialpha',
                        'description': 'Last common ancestor of Type Iα',
                        'terminals': terminals
                    }
                elif type_i_beta and not type_i_alpha and not type_ii:
                    key_nodes[node_id] = {
                        'name': 'LCA_Type_Ibeta',
                        'description': 'Last common ancestor of Type Iβ',
                        'terminals': terminals
                    }
                elif type_ii and not type_i_alpha and not type_i_beta:
                    key_nodes[node_id] = {
                        'name': 'LCA_Type_II',
                        'description': 'Last common ancestor of Type II',
                        'terminals': terminals
                    }
                
    except Exception as e:
        print(f"    Error identifying key nodes: {e}")
    
    return key_nodes

def create_consensus_sequence(sequence, posterior_probs, cutoff=0.80):
    """Create consensus sequence with ambiguous positions marked"""
    consensus = []
    ambiguous_sites = []
    
    for i, (aa, pp) in enumerate(zip(sequence, posterior_probs)):
        if pp >= cutoff:
            consensus.append(aa)
        else:
            consensus.append('X')  # Mark as ambiguous
            ambiguous_sites.append(i + 1)  # 1-based position
    
    return ''.join(consensus), ambiguous_sites

def reconstruct_ancestors(alignment_file, tree_file, model, partition_name, asr_dir):
    """Perform ancestral sequence reconstruction for a partition"""
    print(f"\nReconstructing ancestors for: {partition_name}")
    
    # Check inputs
    if not Path(alignment_file).exists() or not Path(tree_file).exists():
        print(f"  Error: Missing input files")
        return None
    
    # Check if state file already exists
    prefix = asr_dir / f"{partition_name}_asr"
    state_file = f"{prefix}.state"
    
    if not Path(state_file).exists():
        # Run IQ-TREE ASR
        if not run_iqtree_asr(alignment_file, tree_file, model, prefix, CONFIG['cpu']):
            return None
    else:
        print(f"  Using existing state file")
    
    # Parse ASR results
    if not Path(state_file).exists():
        print(f"  Error: State file not available")
        return None
    
    ancestral_seqs, posterior_probs = parse_asr_state_file(state_file)
    
    if not ancestral_seqs:
        print(f"  Error: No ancestral sequences recovered")
        return None
    
    print(f"  Reconstructed {len(ancestral_seqs)} ancestral nodes")
    
    # Identify key nodes
    key_nodes = identify_key_nodes(tree_file)
    
    # Process ancestral sequences
    results = []
    sequences_to_save = []
    
    for node_id, sequence in ancestral_seqs.items():
        if node_id in posterior_probs:
            pps = posterior_probs[node_id]
            
            # Calculate statistics
            mean_pp = sum(pps) / len(pps) if pps else 0
            high_conf_sites = len([p for p in pps if p >= CONFIG['posterior_cutoff']])
            prop_high_conf = high_conf_sites / len(pps) if pps else 0
            
            # Create consensus
            consensus_seq, ambiguous_sites = create_consensus_sequence(
                sequence, pps, CONFIG['posterior_cutoff']
            )
            
            # Get node information
            node_info = key_nodes.get(node_id, {})
            node_name = node_info.get('name', node_id)
            node_desc = node_info.get('description', f'Ancestral node {node_id}')
            
            # Create SeqRecord
            seq_record = SeqRecord(
                Seq(consensus_seq),
                id=f"{partition_name}_{node_name}",
                description=node_desc
            )
            sequences_to_save.append(seq_record)
            
            # Store results
            result = {
                'node_id': node_id,
                'node_name': node_name,
                'description': node_desc,
                'sequence_length': len(sequence),
                'mean_posterior': mean_pp,
                'high_conf_sites': high_conf_sites,
                'prop_high_conf': prop_high_conf,
                'n_ambiguous': len(ambiguous_sites),
                'ambiguous_positions': ambiguous_sites[:10] if len(ambiguous_sites) > 10 else ambiguous_sites
            }
            results.append(result)
            
            # Print key nodes
            if node_name != node_id:
                print(f"    {node_name}: PP={mean_pp:.3f}, {prop_high_conf*100:.1f}% high confidence")
    
    # Save ancestral sequences
    if sequences_to_save:
        output_file = asr_dir / f"{partition_name}_ancestors.faa"
        SeqIO.write(sequences_to_save, output_file, 'fasta')
        print(f"  Saved {len(sequences_to_save)} ancestral sequences")
    
    return results

def main():
    print("=" * 70)
    print("DAH7PS Ancestral Sequence Reconstruction (Phase 6)")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    trees_dir = project_dir / 'trees'
    msa_dir = project_dir / 'msa'
    asr_dir = project_dir / 'asr'
    data_dir = project_dir / 'data' / 'processed'
    
    asr_dir.mkdir(exist_ok=True)
    
    # Load tree summary
    tree_summary_file = data_dir / 'tree_summary.tsv'
    if not tree_summary_file.exists():
        print("Error: No tree summary found. Run Phase 5 first.")
        return 1
    
    trees_df = pd.read_csv(tree_summary_file, sep='\t')
    
    # Perform ASR for each partition
    all_results = []
    priority_partitions = ['all_sequences', 'core_Type_I_combined', 
                          'core_Type_Ialpha', 'core_Type_Ibeta', 'core_Type_II']
    
    for partition in priority_partitions:
        if partition not in trees_df['partition'].values:
            continue
        
        row = trees_df[trees_df['partition'] == partition].iloc[0]
        
        # Get files and model
        alignment_file = msa_dir / f"{partition}.trim.faa"
        tree_file = row['rooted_tree']
        model = row['best_model']
        
        if Path(alignment_file).exists() and Path(tree_file).exists():
            results = reconstruct_ancestors(
                alignment_file, tree_file, model, partition, asr_dir
            )
            
            if results:
                for r in results:
                    r['partition'] = partition
                    all_results.extend([r])
    
    # Save ASR summary
    if all_results:
        summary_df = pd.DataFrame(all_results)
        summary_file = data_dir / 'asr_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\nSaved ASR summary to {summary_file}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Ancestral Reconstruction Summary:")
    print("=" * 70)
    
    # Group by partition
    for partition in priority_partitions:
        partition_results = [r for r in all_results if r['partition'] == partition]
        if partition_results:
            print(f"\n{partition}:")
            print(f"  Reconstructed {len(partition_results)} ancestral nodes")
            
            # Key ancestors
            for r in partition_results:
                if r['node_name'] != r['node_id']:
                    print(f"    {r['node_name']}: {r['prop_high_conf']*100:.1f}% sites with PP≥{CONFIG['posterior_cutoff']}")
    
    # Quality assessment
    print("\nQuality Assessment:")
    key_ancestors = [r for r in all_results if r['node_name'] != r['node_id']]
    
    if key_ancestors:
        good_quality = [r for r in key_ancestors if r['mean_posterior'] >= 0.70]
        print(f"  Key ancestors with mean PP≥0.70: {len(good_quality)}/{len(key_ancestors)}")
        
        if len(good_quality) < len(key_ancestors):
            print("  ⚠ Some key ancestors have low posterior probability")
        else:
            print("  ✓ All key ancestors meet quality thresholds")
    else:
        print("  ⚠ No key ancestral nodes identified")
    
    print("\n✓ Phase 6 complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
