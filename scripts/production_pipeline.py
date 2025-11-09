#!/usr/bin/env python3
"""
Production-scale pipeline for DAH7PS analysis
Processes large proteome datasets efficiently
"""

import os
import sys
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import tempfile
import json
import time
from collections import defaultdict
import multiprocessing as mp

class ProductionPipeline:
    def __init__(self):
        self.project_dir = Path('/home/luogy/bigtree')
        self.data_dir = self.project_dir / 'data'
        self.seqs_dir = self.project_dir / 'seqs_production'
        self.msa_dir = self.project_dir / 'msa_production'
        self.trees_dir = self.project_dir / 'trees_production'
        self.asr_dir = self.project_dir / 'asr_production'
        
        # Create directories
        for dir_path in [self.seqs_dir, self.msa_dir, self.trees_dir, self.asr_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Configuration
        self.config = {
            'pfam_hmm': self.data_dir / 'raw' / 'Pfam-A.hmm',
            'proteomes': self.data_dir / 'raw' / 'uniprot_reference_proteomes.faa',
            'type_I_profile': 'DAHP_synth_1',
            'type_II_profile': 'DAHP_synth_2',
            'evalue_threshold': 1.0e-10,
            'min_domain_coverage': 0.85,
            'max_ambiguous_fraction': 0.20,
            'min_length': 250,
            'max_length': 1500,
            'clustering_identity': 0.90,
            'clustering_coverage': 0.85,
            'cpu': min(mp.cpu_count() - 1, 8)
        }
    
    def run_hmmer_search_production(self):
        """Run HMMER search on production dataset"""
        print("=" * 70)
        print("Phase 1: HMMER Search on Production Dataset")
        print("=" * 70)
        
        # Count input sequences
        print("\nCounting input sequences...")
        seq_count = 0
        with open(self.config['proteomes'], 'r') as f:
            seq_count = sum(1 for line in f if line.startswith('>'))
        print(f"  Input: {seq_count:,} sequences from reference proteomes")
        
        all_hits = []
        
        # Search for Type I
        print(f"\nSearching for Type I ({self.config['type_I_profile']})...")
        hits_I = self._run_hmmer_profile(self.config['type_I_profile'], 'Type_I')
        all_hits.extend(hits_I)
        
        # Search for Type II
        print(f"\nSearching for Type II ({self.config['type_II_profile']})...")
        hits_II = self._run_hmmer_profile(self.config['type_II_profile'], 'Type_II')
        all_hits.extend(hits_II)
        
        print(f"\nTotal hits before filtering: {len(all_hits)}")
        
        # Save raw hits
        hits_df = pd.DataFrame(all_hits)
        hits_file = self.seqs_dir / 'hmmer_hits_raw.tsv'
        hits_df.to_csv(hits_file, sep='\t', index=False)
        
        return all_hits
    
    def _run_hmmer_profile(self, profile, profile_type):
        """Run HMMER for a specific profile"""
        # Extract HMM
        temp_hmm = self.seqs_dir / f"{profile}.hmm"
        cmd = f"hmmfetch {self.config['pfam_hmm']} {profile} > {temp_hmm}"
        subprocess.run(cmd, shell=True, check=True)
        
        # Run hmmsearch
        out_file = self.seqs_dir / f"hmmsearch_{profile}.txt"
        domtbl_file = self.seqs_dir / f"hmmsearch_{profile}.domtbl"
        
        cmd = [
            '/home/luogy/miniconda3/envs/dah7ps/bin/hmmsearch',
            '--domtblout', str(domtbl_file),
            '-o', str(out_file),
            '--cpu', str(self.config['cpu']),
            '--cut_ga',  # Use gathering cutoffs
            str(temp_hmm),
            str(self.config['proteomes'])
        ]
        
        print(f"  Running hmmsearch with {self.config['cpu']} CPUs...")
        subprocess.run(cmd, check=True)
        
        # Parse results
        hits = []
        with open(domtbl_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 22:
                    hit = {
                        'seq_id': parts[0],
                        'profile': profile,
                        'type': profile_type,
                        'evalue': float(parts[6]),
                        'score': float(parts[7]),
                        'ali_from': int(parts[17]),
                        'ali_to': int(parts[18]),
                        'hmm_from': int(parts[15]),
                        'hmm_to': int(parts[16])
                    }
                    hits.append(hit)
        
        print(f"  Found {len(hits)} hits")
        return hits
    
    def filter_and_extract_sequences(self, hits):
        """Filter hits and extract sequences"""
        print("\n" + "=" * 70)
        print("Phase 2: Sequence Filtering and Extraction")
        print("=" * 70)
        
        # Load all sequences into memory (chunked for large files)
        print("\nLoading sequences...")
        sequences = {}
        
        # Process in chunks to handle large files
        chunk_size = 10000
        chunk = []
        
        with open(self.config['proteomes'], 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                chunk.append(record)
                if len(chunk) >= chunk_size:
                    for rec in chunk:
                        sequences[rec.id] = rec
                    chunk = []
                    print(f"  Loaded {len(sequences):,} sequences...", end='\r')
        
        # Process remaining chunk
        for rec in chunk:
            sequences[rec.id] = rec
        
        print(f"  Loaded {len(sequences):,} sequences total")
        
        # Extract hit sequences
        print("\nExtracting and filtering hit sequences...")
        filtered_sequences = []
        
        hits_df = pd.DataFrame(hits)
        unique_seq_ids = hits_df['seq_id'].unique()
        
        for seq_id in unique_seq_ids:
            if seq_id not in sequences:
                continue
            
            seq_record = sequences[seq_id]
            seq_len = len(seq_record.seq)
            
            # Apply filters
            if seq_len < self.config['min_length'] or seq_len > self.config['max_length']:
                continue
            
            seq_str = str(seq_record.seq).upper()
            x_fraction = seq_str.count('X') / seq_len
            if x_fraction > self.config['max_ambiguous_fraction']:
                continue
            
            # Get hit info
            seq_hits = hits_df[hits_df['seq_id'] == seq_id]
            best_hit = seq_hits.loc[seq_hits['score'].idxmax()]
            
            # Check domain coverage
            domain_len = best_hit['ali_to'] - best_hit['ali_from'] + 1
            expected_len = 250 if 'Type_I' in best_hit['type'] else 200
            if domain_len / expected_len < self.config['min_domain_coverage']:
                continue
            
            # Extract organism info from header
            organism = self._extract_organism(seq_record.description)
            proteome = seq_id.split('|')[0] if '|' in seq_id else 'unknown'
            
            filtered_sequences.append({
                'seq_id': seq_id,
                'proteome': proteome,
                'organism': organism,
                'description': seq_record.description,
                'length': seq_len,
                'type': best_hit['type'],
                'score': best_hit['score'],
                'evalue': best_hit['evalue'],
                'sequence': str(seq_record.seq),
                'seq_record': seq_record
            })
        
        print(f"  Filtered to {len(filtered_sequences)} high-quality sequences")
        
        # Save filtered sequences
        fasta_file = self.seqs_dir / 'dah7ps_hits_filtered.faa'
        with open(fasta_file, 'w') as f:
            for seq_info in filtered_sequences:
                f.write(f">{seq_info['seq_id']} {seq_info['description']}\n")
                f.write(f"{seq_info['sequence']}\n")
        
        return filtered_sequences
    
    def _extract_organism(self, description):
        """Extract organism name from sequence description"""
        if 'OS=' in description:
            start = description.index('OS=') + 3
            end = description.index(' OX=') if ' OX=' in description else len(description)
            return description[start:end].strip()
        return 'Unknown'
    
    def cluster_sequences_production(self, sequences):
        """Cluster sequences using MMseqs2 for speed"""
        print("\n" + "=" * 70)
        print("Phase 3: Sequence Clustering")
        print("=" * 70)
        
        # Separate by type
        type_I_seqs = [s for s in sequences if s['type'] == 'Type_I']
        type_II_seqs = [s for s in sequences if s['type'] == 'Type_II']
        
        print(f"\nType I sequences: {len(type_I_seqs)}")
        print(f"Type II sequences: {len(type_II_seqs)}")
        
        clustered_all = []
        
        # Cluster Type I
        if type_I_seqs:
            print(f"\nClustering Type I at {self.config['clustering_identity']*100}% identity...")
            clustered_I = self._run_mmseqs2_clustering(type_I_seqs, 'type_I')
            clustered_all.extend(clustered_I)
        
        # Cluster Type II
        if type_II_seqs:
            print(f"\nClustering Type II at {self.config['clustering_identity']*100}% identity...")
            clustered_II = self._run_mmseqs2_clustering(type_II_seqs, 'type_II')
            clustered_all.extend(clustered_II)
        
        print(f"\nTotal representative sequences: {len(clustered_all)}")
        
        # Further classify Type I into Iα and Iβ
        classified = self._classify_sequences(clustered_all)
        
        return classified
    
    def _run_mmseqs2_clustering(self, sequences, prefix):
        """Use MMseqs2 for fast clustering"""
        # Check if mmseqs2 is available, otherwise use CD-HIT
        try:
            subprocess.run(['mmseqs', '-h'], capture_output=True, check=True)
            use_mmseqs = True
        except:
            use_mmseqs = False
            print("  MMseqs2 not found, using CD-HIT instead")
        
        # Write sequences to temp file
        temp_fasta = self.seqs_dir / f'{prefix}_temp.faa'
        with open(temp_fasta, 'w') as f:
            for seq in sequences:
                f.write(f">{seq['seq_id']}\n{seq['sequence']}\n")
        
        if use_mmseqs:
            # MMseqs2 clustering
            db = self.seqs_dir / f'{prefix}_db'
            clu = self.seqs_dir / f'{prefix}_clu'
            rep = self.seqs_dir / f'{prefix}_rep.faa'
            
            cmd_create = f"mmseqs createdb {temp_fasta} {db}"
            cmd_cluster = f"mmseqs cluster {db} {clu} tmp --min-seq-id {self.config['clustering_identity']} -c {self.config['clustering_coverage']} --threads {self.config['cpu']}"
            cmd_result = f"mmseqs result2repseq {db} {clu} {rep}"
            
            subprocess.run(cmd_create, shell=True, check=True)
            subprocess.run(cmd_cluster, shell=True, check=True)
            subprocess.run(cmd_result, shell=True, check=True)
            
            # Clean up
            subprocess.run(f"rm -rf tmp {db}* {clu}*", shell=True)
            
        else:
            # CD-HIT clustering
            rep = self.seqs_dir / f'{prefix}_rep.faa'
            cmd = [
                '/home/luogy/miniconda3/envs/dah7ps/bin/cd-hit',
                '-i', str(temp_fasta),
                '-o', str(rep),
                '-c', str(self.config['clustering_identity']),
                '-n', '5',
                '-M', '8000',
                '-T', str(self.config['cpu'])
            ]
            subprocess.run(cmd, check=True)
        
        # Read representatives
        representatives = []
        seq_dict = {s['seq_id']: s for s in sequences}
        
        for record in SeqIO.parse(rep, 'fasta'):
            if record.id in seq_dict:
                representatives.append(seq_dict[record.id])
        
        print(f"  Reduced to {len(representatives)} representatives")
        
        # Clean up
        temp_fasta.unlink()
        
        return representatives
    
    def _classify_sequences(self, sequences):
        """Classify sequences into Type Iα, Iβ, and II"""
        print("\nClassifying sequences...")
        
        classified = {
            'Type_Ialpha': [],
            'Type_Ibeta': [],
            'Type_II': []
        }
        
        for seq in sequences:
            if seq['type'] == 'Type_II':
                classified['Type_II'].append(seq)
            else:
                # Type I - classify based on organism and known patterns
                # This is a simplified classification
                organism = seq['organism'].lower()
                desc = seq['description'].lower()
                
                # Known Type Iα patterns (Phe-sensitive)
                if any(x in desc for x in ['arog', 'phe-sensitive', 'phenylalanine']):
                    classified['Type_Ialpha'].append(seq)
                # Known Type Iβ patterns (Tyr-sensitive)
                elif any(x in desc for x in ['arof', 'tyr-sensitive', 'tyrosine']):
                    classified['Type_Ibeta'].append(seq)
                # Default classification based on organism
                elif 'escherichia' in organism or 'salmonella' in organism:
                    # E. coli has both, split randomly for now
                    if len(classified['Type_Ialpha']) <= len(classified['Type_Ibeta']):
                        classified['Type_Ialpha'].append(seq)
                    else:
                        classified['Type_Ibeta'].append(seq)
                else:
                    # Default to Type Iβ
                    classified['Type_Ibeta'].append(seq)
        
        for class_name, seqs in classified.items():
            print(f"  {class_name}: {len(seqs)} sequences")
        
        # Save classified sequences
        for class_name, seqs in classified.items():
            if seqs:
                fasta_file = self.seqs_dir / f'dah7ps_{class_name.lower()}_production.faa'
                with open(fasta_file, 'w') as f:
                    for seq in seqs:
                        f.write(f">{seq['seq_id']} {seq['description']}\n")
                        f.write(f"{seq['sequence']}\n")
        
        # Combine all for global analysis
        all_file = self.seqs_dir / 'dah7ps_all_production.faa'
        with open(all_file, 'w') as f:
            for class_name, seqs in classified.items():
                for seq in seqs:
                    f.write(f">{seq['seq_id']} {seq['description']}\n")
                    f.write(f"{seq['sequence']}\n")
        
        return classified
    
    def create_metadata(self, classified):
        """Create comprehensive metadata table"""
        print("\n" + "=" * 70)
        print("Phase 4: Metadata Generation")
        print("=" * 70)
        
        metadata = []
        
        for class_name, seqs in classified.items():
            for seq in seqs:
                metadata.append({
                    'seq_id': seq['seq_id'],
                    'proteome': seq['proteome'],
                    'organism': seq['organism'],
                    'classification': class_name,
                    'length': seq['length'],
                    'hmmer_score': seq['score'],
                    'evalue': seq['evalue'],
                    'description': seq['description']
                })
        
        metadata_df = pd.DataFrame(metadata)
        
        # Add organism groups
        metadata_df['domain'] = metadata_df['organism'].apply(self._assign_domain)
        
        # Save metadata
        metadata_file = self.data_dir / 'processed' / 'production_metadata.tsv'
        metadata_df.to_csv(metadata_file, sep='\t', index=False)
        
        print(f"\nMetadata saved: {len(metadata_df)} sequences")
        
        # Print statistics
        print("\nOrganism distribution:")
        domain_counts = metadata_df['domain'].value_counts()
        for domain, count in domain_counts.items():
            print(f"  {domain}: {count}")
        
        return metadata_df
    
    def _assign_domain(self, organism):
        """Assign taxonomic domain based on organism name"""
        organism_lower = organism.lower()
        
        if any(x in organism_lower for x in ['escherichia', 'pseudomonas', 'bacillus', 'salmonella',
                                              'mycobacterium', 'streptomyces', 'synechocystis',
                                              'listeria', 'staphylococcus', 'corynebacterium',
                                              'helicobacter', 'clostridium', 'streptococcus']):
            return 'Bacteria'
        elif any(x in organism_lower for x in ['saccharomyces', 'candida', 'aspergillus',
                                                'neurospora', 'schizosaccharomyces']):
            return 'Fungi'
        elif any(x in organism_lower for x in ['arabidopsis', 'oryza', 'zea', 'solanum',
                                                'medicago', 'physcomitrella', 'chlamydomonas']):
            return 'Plants'
        elif any(x in organism_lower for x in ['sulfolobus', 'pyrococcus', 'methanocaldococcus',
                                                'halobacterium']):
            return 'Archaea'
        else:
            return 'Unknown'
    
    def run_production_pipeline(self):
        """Main pipeline execution"""
        print("\n" + "=" * 80)
        print("DAH7PS PRODUCTION-SCALE ANALYSIS PIPELINE")
        print("=" * 80)
        
        start_time = time.time()
        
        # Phase 1: HMMER search
        hits = self.run_hmmer_search_production()
        
        # Phase 2: Filter and extract
        filtered_seqs = self.filter_and_extract_sequences(hits)
        
        # Phase 3: Clustering
        classified = self.cluster_sequences_production(filtered_seqs)
        
        # Phase 4: Metadata
        metadata = self.create_metadata(classified)
        
        # Summary
        elapsed = time.time() - start_time
        
        print("\n" + "=" * 80)
        print("PRODUCTION PIPELINE COMPLETE")
        print("=" * 80)
        print(f"\nTime elapsed: {elapsed/60:.1f} minutes")
        print(f"Final sequences: {len(metadata)}")
        print(f"Output directory: {self.seqs_dir}")
        
        print("\nNext steps:")
        print("1. Run domain annotation: python scripts/production_domains.py")
        print("2. Generate alignments: python scripts/production_alignments.py")
        print("3. Build trees: python scripts/production_trees.py")
        print("4. Reconstruct ancestors: python scripts/production_asr.py")
        
        return metadata

def main():
    pipeline = ProductionPipeline()
    pipeline.run_production_pipeline()
    return 0

if __name__ == "__main__":
    sys.exit(main())
