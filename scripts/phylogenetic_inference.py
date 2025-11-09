#!/usr/bin/env python3
"""
Phase 5: Phylogenetic Inference for DAH7PS Project
Build phylogenetic trees using IQ-TREE with model selection and bootstrap support
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import Phylo, SeqIO, AlignIO
import re

# Configuration
CONFIG = {
    'iqtree_path': '/home/luogy/miniconda3/envs/dah7ps/bin/iqtree',
    'model_set': 'LG,WAG,JTT',
    'bootstrap_replicates': 1000,
    'sh_alrt_replicates': 1000,
    'cpu': 4,
    'memory': '4G',
    'rooting': 'midpoint'  # or 'outgroup'
}

def run_iqtree(alignment_file, prefix, cpu=4):
    """Run IQ-TREE with model selection and bootstrap"""
    print(f"  Running IQ-TREE...")
    
    cmd = [
        CONFIG['iqtree_path'],
        '-s', str(alignment_file),
        '-pre', str(prefix),
        '-m', 'MFP',  # ModelFinder Plus
        '-mset', CONFIG['model_set'],
        '-bb', str(CONFIG['bootstrap_replicates']),  # UFBoot
        '-alrt', str(CONFIG['sh_alrt_replicates']),  # SH-aLRT
        '-nt', str(cpu),
        '-seed', '42',  # For reproducibility
        '-quiet'
    ]
    
    # Add memory limit if specified
    if CONFIG.get('memory'):
        cmd.extend(['-mem', CONFIG['memory']])
    
    # Run IQ-TREE
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    Error: {result.stderr}")
        # Try without memory limit
        if '-mem' in cmd:
            print("    Retrying without memory limit...")
            cmd.remove('-mem')
            cmd.remove(CONFIG['memory'])
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd)
    
    # Parse output for best model
    log_file = f"{prefix}.log"
    best_model = None
    if Path(log_file).exists():
        with open(log_file, 'r') as f:
            log_content = f.read()
            model_match = re.search(r'Best-fit model: ([^\s]+)', log_content)
            if model_match:
                best_model = model_match.group(1)
    
    return best_model

def parse_tree_support(tree_file):
    """Parse tree to get support value statistics"""
    try:
        tree = Phylo.read(tree_file, 'newick')
        
        ufboot_values = []
        shalrt_values = []
        
        for clade in tree.find_clades():
            if clade.confidence is not None:
                # IQ-TREE format: SH-aLRT/UFBoot
                if '/' in str(clade.name):
                    parts = clade.name.split('/')
                    if len(parts) == 2:
                        try:
                            shalrt_values.append(float(parts[0]))
                            ufboot_values.append(float(parts[1]))
                        except:
                            pass
        
        stats = {
            'n_internal_nodes': len([c for c in tree.find_clades() if not c.is_terminal()]) - 1,
            'n_supported_ufboot95': len([v for v in ufboot_values if v >= 95]),
            'n_supported_ufboot80': len([v for v in ufboot_values if v >= 80]),
            'n_supported_shalrt80': len([v for v in shalrt_values if v >= 80]),
            'n_supported_shalrt70': len([v for v in shalrt_values if v >= 70]),
            'mean_ufboot': sum(ufboot_values) / len(ufboot_values) if ufboot_values else 0,
            'mean_shalrt': sum(shalrt_values) / len(shalrt_values) if shalrt_values else 0
        }
        
        return stats
    except Exception as e:
        print(f"    Error parsing tree: {e}")
        return None

def root_tree(tree_file, output_file, method='midpoint', outgroup=None):
    """Root tree using specified method"""
    try:
        tree = Phylo.read(tree_file, 'newick')
        
        if method == 'midpoint':
            # Midpoint rooting
            tree.root_at_midpoint()
        elif method == 'outgroup' and outgroup:
            # Outgroup rooting
            outgroup_clade = tree.find_any(name=outgroup)
            if outgroup_clade:
                tree.root_with_outgroup(outgroup_clade)
        
        # Write rooted tree
        Phylo.write(tree, output_file, 'newick')
        return True
    except Exception as e:
        print(f"    Error rooting tree: {e}")
        # Copy unrooted tree if rooting fails
        import shutil
        shutil.copy(tree_file, output_file)
        return False

def build_tree(alignment_file, partition_name, trees_dir):
    """Build phylogenetic tree for a partition"""
    print(f"\nBuilding tree for: {partition_name}")
    
    # Check alignment
    try:
        alignment = AlignIO.read(alignment_file, 'fasta')
        n_seqs = len(alignment)
        aln_len = alignment.get_alignment_length()
        print(f"  Input alignment: {n_seqs} sequences x {aln_len} positions")
        
        if n_seqs < 4:
            print(f"  Warning: Too few sequences for reliable phylogeny ({n_seqs})")
    except Exception as e:
        print(f"  Error reading alignment: {e}")
        return None
    
    # Run IQ-TREE
    prefix = trees_dir / partition_name
    best_model = run_iqtree(alignment_file, prefix, CONFIG['cpu'])
    
    if best_model:
        print(f"  Best-fit model: {best_model}")
    
    # Check output files
    tree_file = f"{prefix}.treefile"
    if not Path(tree_file).exists():
        print(f"  Error: Tree file not created")
        return None
    
    # Parse support values
    support_stats = parse_tree_support(tree_file)
    if support_stats:
        print(f"  Support values:")
        print(f"    UFBoot ≥95: {support_stats['n_supported_ufboot95']}/{support_stats['n_internal_nodes']} nodes")
        print(f"    SH-aLRT ≥80: {support_stats['n_supported_shalrt80']}/{support_stats['n_internal_nodes']} nodes")
        print(f"    Mean UFBoot: {support_stats['mean_ufboot']:.1f}")
        print(f"    Mean SH-aLRT: {support_stats['mean_shalrt']:.1f}")
    
    # Root tree
    rooted_file = f"{prefix}.rooted.treefile"
    if root_tree(tree_file, rooted_file, CONFIG['rooting']):
        print(f"  Tree rooted using {CONFIG['rooting']} method")
    
    return {
        'partition': partition_name,
        'alignment_file': str(alignment_file),
        'tree_file': tree_file,
        'rooted_tree': rooted_file,
        'best_model': best_model,
        'n_sequences': n_seqs,
        'alignment_length': aln_len,
        'support_stats': support_stats
    }

def main():
    print("=" * 70)
    print("DAH7PS Phylogenetic Inference (Phase 5)")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    msa_dir = project_dir / 'msa'
    trees_dir = project_dir / 'trees'
    data_dir = project_dir / 'data' / 'processed'
    
    trees_dir.mkdir(exist_ok=True)
    
    # Load alignment summary
    alignment_summary_file = data_dir / 'alignment_summary.tsv'
    if not alignment_summary_file.exists():
        print("Error: No alignment summary found. Run Phase 4 first.")
        return 1
    
    alignments_df = pd.read_csv(alignment_summary_file, sep='\t')
    
    # Build trees for each partition
    results = []
    priority_partitions = ['all_sequences', 'core_Type_I_combined', 'core_Type_Ialpha', 
                          'core_Type_Ibeta', 'core_Type_II']
    
    # Sort partitions by priority
    partitions = []
    for p in priority_partitions:
        if p in alignments_df['partition'].values:
            partitions.append(p)
    
    for partition in partitions:
        row = alignments_df[alignments_df['partition'] == partition].iloc[0]
        alignment_file = row['trimmed_file']
        
        if Path(alignment_file).exists():
            result = build_tree(alignment_file, partition, trees_dir)
            if result:
                results.append(result)
    
    # Save tree building summary
    if results:
        summary_data = []
        for r in results:
            summary_data.append({
                'partition': r['partition'],
                'n_sequences': r['n_sequences'],
                'alignment_length': r['alignment_length'],
                'best_model': r['best_model'],
                'n_internal_nodes': r['support_stats']['n_internal_nodes'] if r['support_stats'] else 0,
                'ufboot95_support': r['support_stats']['n_supported_ufboot95'] if r['support_stats'] else 0,
                'shalrt80_support': r['support_stats']['n_supported_shalrt80'] if r['support_stats'] else 0,
                'mean_ufboot': r['support_stats']['mean_ufboot'] if r['support_stats'] else 0,
                'mean_shalrt': r['support_stats']['mean_shalrt'] if r['support_stats'] else 0,
                'tree_file': r['tree_file'],
                'rooted_tree': r['rooted_tree']
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = data_dir / 'tree_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\nSaved tree summary to {summary_file}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Phylogenetic Inference Summary:")
    print("=" * 70)
    
    for r in results:
        print(f"\n{r['partition']}:")
        print(f"  Model: {r['best_model']}")
        if r['support_stats']:
            prop_ufboot = r['support_stats']['n_supported_ufboot95'] / r['support_stats']['n_internal_nodes'] * 100
            prop_shalrt = r['support_stats']['n_supported_shalrt80'] / r['support_stats']['n_internal_nodes'] * 100
            print(f"  Strong support: {prop_ufboot:.1f}% UFBoot≥95, {prop_shalrt:.1f}% SH-aLRT≥80")
    
    # Quality assessment
    print("\nQuality Assessment:")
    passed = True
    for r in results:
        if r['support_stats']:
            prop_ufboot = r['support_stats']['n_supported_ufboot95'] / r['support_stats']['n_internal_nodes']
            if prop_ufboot < 0.5:  # Less than 50% strong support
                print(f"  ⚠ {r['partition']}: Low bootstrap support (only {prop_ufboot*100:.1f}% UFBoot≥95)")
                passed = False
    
    if passed:
        print("  ✓ All trees meet quality thresholds")
    
    print("\n✓ Phase 5 complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
