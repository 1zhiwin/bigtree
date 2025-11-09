#!/usr/bin/env python3
"""
Phase 4: Multiple Sequence Alignment for DAH7PS Project
Generate and trim alignments for each sequence partition
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
import shutil

# Configuration
CONFIG = {
    'aligner': 'mafft',  # or 'muscle5'
    'mafft_algorithm': 'linsi',  # L-INS-i for accuracy
    'trimal_mode': 'automated1',  # automated1, gappyout, or strict
    'min_sequences': 3,  # Minimum sequences for alignment
    'cpu': 4
}

def run_mafft(input_file, output_file, algorithm='linsi', cpu=4):
    """Run MAFFT alignment"""
    print(f"  Running MAFFT ({algorithm})...")
    
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/mafft',
        '--maxiterate', '1000',
        '--thread', str(cpu),
        '--reorder',
        '--quiet'  # Suppress stderr output
    ]
    
    # Add algorithm-specific options
    if algorithm == 'linsi':
        cmd.extend(['--localpair'])
    elif algorithm == 'einsi':
        cmd.extend(['--genafpair'])
    elif algorithm == 'ginsi':
        cmd.extend(['--globalpair'])
    
    cmd.extend([str(input_file)])
    
    # Run MAFFT and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    Error: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    
    # Write output
    with open(output_file, 'w') as f:
        f.write(result.stdout)
    
    return output_file

def run_muscle5(input_file, output_file, cpu=4):
    """Run MUSCLE v5 alignment"""
    print(f"  Running MUSCLE v5...")
    
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/muscle',
        '-align', str(input_file),
        '-output', str(output_file),
        '-threads', str(cpu)
    ]
    
    subprocess.run(cmd, check=True, capture_output=True)
    
    return output_file

def run_trimal(input_file, output_file, mode='automated1'):
    """Run trimAl to trim alignment"""
    print(f"  Running trimAl ({mode})...")
    
    cmd = [
        '/home/luogy/miniconda3/envs/dah7ps/bin/trimal',
        '-in', str(input_file),
        '-out', str(output_file),
        f'-{mode}'
    ]
    
    # Add HTML output for visualization
    html_file = str(output_file).replace('.faa', '.html')
    cmd.extend(['-htmlout', html_file])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    Warning: trimAl failed: {result.stderr}")
        # Copy original if trimming fails
        shutil.copy(input_file, output_file)
    
    return output_file

def calculate_alignment_stats(alignment_file):
    """Calculate statistics for an alignment"""
    try:
        alignment = AlignIO.read(alignment_file, 'fasta')
        
        n_seqs = len(alignment)
        aln_len = alignment.get_alignment_length()
        
        # Calculate conservation
        summary = AlignInfo.SummaryInfo(alignment)
        consensus = summary.dumb_consensus(threshold=0.7)
        
        # Count gaps
        gap_counts = []
        for record in alignment:
            gaps = str(record.seq).count('-')
            gap_counts.append(gaps)
        
        avg_gaps = sum(gap_counts) / len(gap_counts) if gap_counts else 0
        gap_percent = (avg_gaps / aln_len * 100) if aln_len > 0 else 0
        
        return {
            'n_sequences': n_seqs,
            'alignment_length': aln_len,
            'avg_gaps': avg_gaps,
            'gap_percentage': gap_percent
        }
    except Exception as e:
        print(f"    Error calculating stats: {e}")
        return None

def align_partition(sequences_file, partition_name, msa_dir, method='mafft'):
    """Align a partition of sequences"""
    print(f"\nAligning partition: {partition_name}")
    
    # Check if file exists and has enough sequences
    sequences = list(SeqIO.parse(sequences_file, 'fasta'))
    if len(sequences) < CONFIG['min_sequences']:
        print(f"  Skipping - only {len(sequences)} sequences (minimum: {CONFIG['min_sequences']})")
        return None
    
    print(f"  Input: {len(sequences)} sequences")
    
    # Run alignment
    raw_alignment = msa_dir / f'{partition_name}.aln.faa'
    
    if method == 'mafft':
        run_mafft(sequences_file, raw_alignment, CONFIG['mafft_algorithm'], CONFIG['cpu'])
    elif method == 'muscle5':
        run_muscle5(sequences_file, raw_alignment, CONFIG['cpu'])
    else:
        raise ValueError(f"Unknown alignment method: {method}")
    
    # Calculate stats for raw alignment
    raw_stats = calculate_alignment_stats(raw_alignment)
    if raw_stats:
        print(f"  Raw alignment: {raw_stats['n_sequences']} seqs x {raw_stats['alignment_length']} positions")
        print(f"  Average gaps: {raw_stats['gap_percentage']:.1f}%")
    
    # Trim alignment
    trimmed_alignment = msa_dir / f'{partition_name}.trim.faa'
    run_trimal(raw_alignment, trimmed_alignment, CONFIG['trimal_mode'])
    
    # Calculate stats for trimmed alignment
    trim_stats = calculate_alignment_stats(trimmed_alignment)
    if trim_stats:
        print(f"  Trimmed: {trim_stats['alignment_length']} positions retained")
    
    return {
        'partition': partition_name,
        'sequences_file': str(sequences_file),
        'raw_alignment': str(raw_alignment),
        'trimmed_alignment': str(trimmed_alignment),
        'raw_stats': raw_stats,
        'trim_stats': trim_stats
    }

def main():
    print("=" * 70)
    print("DAH7PS Multiple Sequence Alignment (Phase 4)")
    print("=" * 70)
    print()
    
    # Setup directories
    project_dir = Path('/home/luogy/bigtree')
    seqs_dir = project_dir / 'seqs'
    msa_dir = project_dir / 'msa'
    data_dir = project_dir / 'data' / 'processed'
    
    msa_dir.mkdir(exist_ok=True)
    
    # Load metadata to get sequence classifications
    metadata_file = data_dir / 'sequence_metadata_annotated.tsv'
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    
    # Define partitions based on classification
    partitions = []
    
    # Type-based partitions
    for class_type in ['Type_Ialpha', 'Type_Ibeta', 'Type_II']:
        class_file = seqs_dir / f'dah7ps_{class_type.lower()}.faa'
        if class_file.exists():
            partitions.append({
                'name': f'core_{class_type}',
                'file': class_file,
                'type': 'core'
            })
    
    # Architecture-based partitions (for regulatory domains)
    for arch in ['core_CM', 'core_ACT', 'core_ACT2', 'core_ACT_CM']:
        arch_file = seqs_dir / f'dah7ps_{arch.lower()}.faa'
        if arch_file.exists():
            partitions.append({
                'name': arch,
                'file': arch_file,
                'type': 'regulatory'
            })
    
    # Combined partitions
    # All Type I sequences
    type_i_seqs = metadata_df[metadata_df['classification'].isin(['Type_Ialpha', 'Type_Ibeta'])]['seq_id'].tolist()
    if len(type_i_seqs) >= CONFIG['min_sequences']:
        all_seqs = list(SeqIO.parse(seqs_dir / 'dah7ps_nonredundant.faa', 'fasta'))
        type_i_records = [s for s in all_seqs if s.id in type_i_seqs]
        if type_i_records:
            type_i_file = seqs_dir / 'dah7ps_type_i_combined.faa'
            SeqIO.write(type_i_records, type_i_file, 'fasta')
            partitions.append({
                'name': 'core_Type_I_combined',
                'file': type_i_file,
                'type': 'combined'
            })
    
    # All sequences (for global tree)
    all_file = seqs_dir / 'dah7ps_nonredundant.faa'
    if all_file.exists():
        partitions.append({
            'name': 'all_sequences',
            'file': all_file,
            'type': 'global'
        })
    
    # Run alignments
    results = []
    for partition in partitions:
        result = align_partition(
            partition['file'],
            partition['name'],
            msa_dir,
            CONFIG['aligner']
        )
        if result:
            result['type'] = partition['type']
            results.append(result)
    
    # Save alignment summary
    if results:
        summary_data = []
        for r in results:
            summary_data.append({
                'partition': r['partition'],
                'type': r['type'],
                'n_sequences': r['raw_stats']['n_sequences'] if r['raw_stats'] else 0,
                'raw_length': r['raw_stats']['alignment_length'] if r['raw_stats'] else 0,
                'trimmed_length': r['trim_stats']['alignment_length'] if r['trim_stats'] else 0,
                'gap_percentage': r['raw_stats']['gap_percentage'] if r['raw_stats'] else 0,
                'raw_file': r['raw_alignment'],
                'trimmed_file': r['trimmed_alignment']
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = data_dir / 'alignment_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"\nSaved alignment summary to {summary_file}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Alignment Summary:")
    print("=" * 70)
    for r in results:
        if r['raw_stats'] and r['trim_stats']:
            print(f"{r['partition']}:")
            print(f"  Sequences: {r['raw_stats']['n_sequences']}")
            print(f"  Length: {r['raw_stats']['alignment_length']} → {r['trim_stats']['alignment_length']} (after trimming)")
            reduction = (1 - r['trim_stats']['alignment_length'] / r['raw_stats']['alignment_length']) * 100
            print(f"  Trimmed: {reduction:.1f}% of positions removed")
    
    # Check for potential issues
    print("\nQuality checks:")
    issues = []
    for r in results:
        if r['trim_stats']:
            if r['trim_stats']['n_sequences'] < 4:
                issues.append(f"  ⚠ {r['partition']}: Very few sequences ({r['trim_stats']['n_sequences']})")
            if r['trim_stats']['alignment_length'] < 100:
                issues.append(f"  ⚠ {r['partition']}: Very short alignment ({r['trim_stats']['alignment_length']} positions)")
    
    if issues:
        for issue in issues:
            print(issue)
    else:
        print("  ✓ All alignments meet quality thresholds")
    
    print("\n✓ Phase 4 complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
