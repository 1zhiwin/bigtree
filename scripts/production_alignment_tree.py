#!/usr/bin/env python3
"""
Production-scale alignment and tree building
Handles 58 representative sequences from diverse organisms
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO, AlignIO, Phylo
import time

class AlignmentTreePipeline:
    def __init__(self):
        self.project_dir = Path('/home/luogy/bigtree')
        self.seqs_dir = self.project_dir / 'seqs_production'
        self.msa_dir = self.project_dir / 'msa_production'
        self.trees_dir = self.project_dir / 'trees_production'
        self.asr_dir = self.project_dir / 'asr_production'
        self.data_dir = self.project_dir / 'data' / 'processed'
        
        # Ensure directories exist
        for dir_path in [self.msa_dir, self.trees_dir, self.asr_dir]:
            dir_path.mkdir(exist_ok=True)
        
        self.cpu = 8
    
    def run_mafft(self, input_file, output_file, algorithm='einsi'):
        """Run MAFFT alignment with E-INS-i for better accuracy"""
        print(f"  Running MAFFT ({algorithm})...")
        
        cmd = [
            '/home/luogy/miniconda3/envs/dah7ps/bin/mafft',
            '--maxiterate', '1000',
            '--thread', str(self.cpu),
            '--reorder',
            '--quiet'
        ]
        
        # E-INS-i is best for sequences with large insertions
        if algorithm == 'einsi':
            cmd.extend(['--genafpair'])
        elif algorithm == 'linsi':
            cmd.extend(['--localpair'])
        
        cmd.append(str(input_file))
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            
            # Count sequences
            aln = AlignIO.read(output_file, 'fasta')
            print(f"    Aligned {len(aln)} sequences, {aln.get_alignment_length()} positions")
        else:
            print(f"    Error: {result.stderr}")
            return False
        
        return True
    
    def run_trimal(self, input_file, output_file, mode='automated1'):
        """Trim alignment with trimAl"""
        print(f"  Running trimAl ({mode})...")
        
        cmd = [
            '/home/luogy/miniconda3/envs/dah7ps/bin/trimal',
            '-in', str(input_file),
            '-out', str(output_file),
            f'-{mode}'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            # Count positions retained
            aln = AlignIO.read(output_file, 'fasta')
            print(f"    Trimmed to {aln.get_alignment_length()} positions")
        else:
            print(f"    Warning: {result.stderr}")
            # Copy original if trimming fails
            import shutil
            shutil.copy(input_file, output_file)
        
        return True
    
    def run_iqtree(self, alignment_file, prefix, model='TEST', bootstrap=1000):
        """Run IQ-TREE with model testing and bootstrap"""
        print(f"  Running IQ-TREE...")
        
        cmd = [
            '/home/luogy/miniconda3/envs/dah7ps/bin/iqtree',
            '-s', str(alignment_file),
            '-pre', str(prefix),
            '-nt', str(self.cpu),
            '-seed', '42',
            '-quiet'
        ]
        
        if model == 'TEST':
            # Use ModelFinder Plus with mixture models
            cmd.extend(['-m', 'MFP'])
            cmd.extend(['-mset', 'LG,WAG,JTT'])
            cmd.extend(['-madd', 'LG+C20,LG+C60'])  # Add mixture models
        else:
            cmd.extend(['-m', model])
        
        # Add bootstrap
        if bootstrap > 0:
            cmd.extend(['-bb', str(bootstrap)])  # UFBoot
            cmd.extend(['-alrt', str(bootstrap)])  # SH-aLRT
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            # Parse log for best model
            log_file = f"{prefix}.log"
            if Path(log_file).exists():
                with open(log_file, 'r') as f:
                    log = f.read()
                    if 'Best-fit model:' in log:
                        model_line = [l for l in log.split('\n') if 'Best-fit model:' in l][0]
                        best_model = model_line.split('Best-fit model:')[1].strip().split()[0]
                        print(f"    Best model: {best_model}")
            
            tree_file = f"{prefix}.treefile"
            if Path(tree_file).exists():
                tree = Phylo.read(tree_file, 'newick')
                print(f"    Tree with {len(tree.get_terminals())} taxa")
        else:
            print(f"    Error: {result.stderr}")
            return False
        
        return True
    
    def run_ancestral_reconstruction(self, alignment_file, tree_file, model, prefix):
        """Run ancestral sequence reconstruction"""
        print(f"  Running ancestral reconstruction...")
        
        cmd = [
            '/home/luogy/miniconda3/envs/dah7ps/bin/iqtree',
            '-s', str(alignment_file),
            '-te', str(tree_file),
            '-m', model,
            '--ancestral',
            '-pre', str(prefix),
            '-nt', str(self.cpu),
            '-quiet'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            state_file = f"{prefix}.state"
            if Path(state_file).exists():
                # Count ancestral nodes
                node_count = 0
                with open(state_file, 'r') as f:
                    for line in f:
                        if line.startswith('Node') and '\t' not in line:
                            node_count += 1
                print(f"    Reconstructed {node_count} ancestral nodes")
        else:
            print(f"    Error: {result.stderr}")
            return False
        
        return True
    
    def process_partition(self, input_file, partition_name):
        """Process a single partition through alignment, tree, and ASR"""
        print(f"\nProcessing partition: {partition_name}")
        print("-" * 50)
        
        # Step 1: Alignment
        raw_aln = self.msa_dir / f"{partition_name}.aln.faa"
        if not self.run_mafft(input_file, raw_aln):
            return False
        
        # Step 2: Trimming
        trimmed_aln = self.msa_dir / f"{partition_name}.trim.faa"
        if not self.run_trimal(raw_aln, trimmed_aln):
            return False
        
        # Step 3: Tree building
        tree_prefix = self.trees_dir / partition_name
        if not self.run_iqtree(trimmed_aln, tree_prefix):
            return False
        
        # Step 4: Ancestral reconstruction
        tree_file = f"{tree_prefix}.treefile"
        
        # Get best model from log
        best_model = 'LG+G4'  # Default
        log_file = f"{tree_prefix}.log"
        if Path(log_file).exists():
            with open(log_file, 'r') as f:
                log = f.read()
                if 'Best-fit model:' in log:
                    model_line = [l for l in log.split('\n') if 'Best-fit model:' in l][0]
                    best_model = model_line.split('Best-fit model:')[1].strip().split()[0]
        
        asr_prefix = self.asr_dir / f"{partition_name}_asr"
        if Path(tree_file).exists():
            self.run_ancestral_reconstruction(trimmed_aln, tree_file, best_model, asr_prefix)
        
        return True
    
    def run_pipeline(self):
        """Run complete pipeline for production data"""
        print("=" * 70)
        print("Production Alignment, Tree, and ASR Pipeline")
        print("=" * 70)
        
        start_time = time.time()
        
        # Define partitions to process
        partitions = [
            ('all_classified.faa', 'all_sequences'),
            ('type_ialpha.faa', 'Type_Ialpha'),
            ('type_ibeta.faa', 'Type_Ibeta')
        ]
        
        results = []
        
        for filename, partition_name in partitions:
            input_file = self.seqs_dir / filename
            if input_file.exists():
                success = self.process_partition(input_file, partition_name)
                results.append((partition_name, success))
            else:
                print(f"\nSkipping {partition_name}: file not found")
        
        # Summary
        elapsed = time.time() - start_time
        
        print("\n" + "=" * 70)
        print("Pipeline Summary")
        print("=" * 70)
        
        for partition, success in results:
            status = "✓" if success else "✗"
            print(f"  {status} {partition}")
        
        print(f"\nTime elapsed: {elapsed/60:.1f} minutes")
        
        # List output files
        print("\nOutput files:")
        print(f"  Alignments: {self.msa_dir}")
        print(f"  Trees: {self.trees_dir}")
        print(f"  ASR: {self.asr_dir}")
        
        return True

def main():
    pipeline = AlignmentTreePipeline()
    pipeline.run_pipeline()
    return 0

if __name__ == "__main__":
    sys.exit(main())
