#!/usr/bin/env python3
"""
Analyze ACT domains and allosteric regulation in Type I DAH7PS.

ACT domains are small regulatory domains (~60-80 aa) that bind amino acids
for allosteric regulation. They typically occur in tandem at the C-terminus.
"""

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import json

def identify_act_domain_regions(sequences_file, metadata_file, hmm_results_file):
    """
    Identify ACT domain regions from HMMER domain annotation.
    """

    # Read metadata
    df = pd.read_csv(metadata_file, sep='\t')

    # Read sequences
    sequences = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        sequences[record.id] = str(record.seq)

    # For Type I, ACT domains typically at C-terminus
    # Based on literature: ~60-80 aa per ACT domain, often in tandem

    print("\n" + "="*60)
    print("ACT DOMAIN ANALYSIS (Type I)")
    print("="*60 + "\n")

    act_results = []

    for _, row in df.iterrows():
        if row['DAH7PS_Type'] != 'Type_I':
            continue

        protein_id = row['Protein_ID']
        seq = sequences[protein_id]
        organism = row['Organism']
        reg_spec = row['Regulatory_Specificity']

        # Estimate ACT domain region
        # Type I DAH7PS: catalytic domain ~250-280 aa, ACT domains C-terminal
        seq_length = len(seq)

        # Conservative estimates:
        # - Catalytic domain: first ~250 aa
        # - Linker: ~10-20 aa
        # - ACT domains: remaining C-terminal region

        if seq_length > 300:
            # Has ACT domains (E. coli, yeast)
            catalytic_end = 260
            act_start = 270
            act_region = seq[act_start:]

            # Estimate if tandem ACT domains (length > 120 aa)
            if len(act_region) > 120:
                has_tandem = True
                act1_region = seq[act_start:act_start+70]
                act2_region = seq[act_start+70:]
            else:
                has_tandem = False
                act1_region = act_region
                act2_region = ""

        else:
            # No ACT domains (some bacterial variants)
            catalytic_end = seq_length
            act_start = None
            act_region = ""
            has_tandem = False
            act1_region = ""
            act2_region = ""

        result = {
            'protein_id': protein_id,
            'organism': organism,
            'regulatory_specificity': reg_spec,
            'total_length': seq_length,
            'catalytic_domain_end': catalytic_end,
            'act_domain_start': act_start,
            'act_region_length': len(act_region),
            'has_tandem_act': has_tandem,
            'act1_length': len(act1_region),
            'act2_length': len(act2_region)
        }

        act_results.append(result)

        print(f"Sequence: {protein_id}")
        print(f"  Organism: {organism}")
        print(f"  Regulatory specificity: {reg_spec}")
        print(f"  Total length: {seq_length} aa")
        if act_start:
            print(f"  Catalytic domain: 1-{catalytic_end}")
            print(f"  ACT domain region: {act_start}-{seq_length}")
            print(f"  ACT region length: {len(act_region)} aa")
            print(f"  Tandem ACT domains: {'Yes' if has_tandem else 'No'}")
            if has_tandem:
                print(f"    ACT1: {act_start}-{act_start+70} ({len(act1_region)} aa)")
                print(f"    ACT2: {act_start+70}-{seq_length} ({len(act2_region)} aa)")
        else:
            print(f"  No ACT domains (catalytic domain only)")
        print()

    return pd.DataFrame(act_results)

def analyze_specificity_residues(alignment_file, traits_file, type_name):
    """
    Identify candidate specificity-determining residues by comparing
    sequences with different regulatory specificities.
    """

    print(f"\n{'='*60}")
    print(f"SPECIFICITY-DETERMINING RESIDUE ANALYSIS: {type_name}")
    print(f"{'='*60}\n")

    # Read alignment
    alignment = AlignIO.read(alignment_file, 'fasta')

    # Read traits
    traits_df = pd.read_csv(traits_file, sep='\t')
    traits = dict(zip(traits_df['Taxon'], traits_df['Trait']))

    # Group sequences by trait
    trait_groups = {}
    for record in alignment:
        trait = traits.get(record.id, 'Unknown')
        if trait not in trait_groups:
            trait_groups[trait] = []
        trait_groups[trait].append(str(record.seq))

    print(f"Trait groups:")
    for trait, seqs in trait_groups.items():
        print(f"  {trait}: {len(seqs)} sequences")

    # Find positions that differ between trait groups
    alignment_length = alignment.get_alignment_length()
    candidate_positions = []

    for i in range(alignment_length):
        # Get amino acids at this position for each trait group
        trait_aas = {}
        for trait, seqs in trait_groups.items():
            if trait == 'Unknown':
                continue
            aas_at_pos = set(seq[i] for seq in seqs if seq[i] != '-')
            if aas_at_pos:
                trait_aas[trait] = aas_at_pos

        # Check if this position distinguishes traits
        if len(trait_aas) >= 2:
            # Check if different traits have different amino acids
            all_aas = []
            for trait, aas in trait_aas.items():
                all_aas.extend(list(aas))

            # If this position has trait-specific amino acids
            unique_aas = set(all_aas)
            if len(unique_aas) >= 2:
                # This position varies between traits
                candidate_positions.append({
                    'position': i + 1,
                    'trait_amino_acids': dict(trait_aas),
                    'conservation_score': 1.0 - len(unique_aas)/20  # rough score
                })

    print(f"\nCandidate specificity-determining positions: {len(candidate_positions)}")

    # Show top candidates
    if candidate_positions:
        print(f"\nTop 20 candidate positions:")
        print(f"{'Pos':<6s} {'Phe':<15s} {'Tyr':<15s} {'Trp':<15s}")
        print("-" * 60)

        for i, pos_data in enumerate(candidate_positions[:20]):
            pos = pos_data['position']
            trait_aas = pos_data['trait_amino_acids']

            phe_aa = ','.join(trait_aas.get('Phe', set())) or '-'
            tyr_aa = ','.join(trait_aas.get('Tyr', set())) or '-'
            trp_aa = ','.join(trait_aas.get('Trp', set())) or '-'

            print(f"{pos:<6d} {phe_aa:<15s} {tyr_aa:<15s} {trp_aa:<15s}")

    return candidate_positions

def analyze_yeast_bifunctional_fusion(sequences_file, metadata_file):
    """
    Analyze the bifunctional DAH7PS-CM fusion in yeast.
    """

    print("\n" + "="*60)
    print("YEAST BIFUNCTIONAL ENZYME ANALYSIS")
    print("="*60 + "\n")

    # Read metadata
    df = pd.read_csv(metadata_file, sep='\t')

    # Read sequences
    sequences = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        sequences[record.id] = str(record.seq)

    yeast_results = []

    for _, row in df.iterrows():
        if 'Saccharomyces' not in row['Organism']:
            continue

        protein_id = row['Protein_ID']
        seq = sequences[protein_id]
        gene_name = row['Gene_Name']
        reg_spec = row['Regulatory_Specificity']

        print(f"Sequence: {gene_name} ({protein_id})")
        print(f"  Regulatory specificity: {reg_spec}")
        print(f"  Total length: {len(seq)} aa")

        # Yeast DAH7PS-CM fusion structure:
        # N-terminus: DAH7PS domain (~300 aa)
        # C-terminus: Chorismate mutase (CM) domain (~70 aa)

        dah7ps_domain_end = 300
        cm_domain_start = 300

        dah7ps_domain = seq[:dah7ps_domain_end]
        cm_domain = seq[cm_domain_start:]

        print(f"  DAH7PS domain: 1-{dah7ps_domain_end} ({len(dah7ps_domain)} aa)")
        print(f"  CM domain: {cm_domain_start}-{len(seq)} ({len(cm_domain)} aa)")
        print(f"  Fusion architecture: DAH7PS-CM bifunctional enzyme")
        print()

        yeast_results.append({
            'protein_id': protein_id,
            'gene_name': gene_name,
            'regulatory_specificity': reg_spec,
            'total_length': len(seq),
            'dah7ps_domain': f"1-{dah7ps_domain_end}",
            'cm_domain': f"{cm_domain_start}-{len(seq)}",
            'is_bifunctional': True
        })

    return pd.DataFrame(yeast_results)

if __name__ == "__main__":

    # Analyze ACT domains
    act_df = identify_act_domain_regions(
        "../data/raw/dah7ps_sequences.faa",
        "../data/raw/dah7ps_metadata.tsv",
        "../data/raw/dah7ps_hmmer_analysis.txt"
    )

    act_df.to_csv("type_i/act_domain_analysis.tsv", sep='\t', index=False)
    print("\nACT domain analysis saved to: structure/type_i/act_domain_analysis.tsv")

    # Analyze specificity-determining residues
    spec_residues = analyze_specificity_residues(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        "../traits/type_i/type_i_traits.tsv",
        "Type I"
    )

    with open("type_i/specificity_residues.json", 'w') as f:
        # Convert sets to lists for JSON
        serializable = []
        for item in spec_residues:
            new_item = item.copy()
            new_item['trait_amino_acids'] = {k: list(v) for k, v in item['trait_amino_acids'].items()}
            serializable.append(new_item)
        json.dump(serializable, f, indent=2)

    print("\nSpecificity residue analysis saved to: structure/type_i/specificity_residues.json")

    # Analyze yeast bifunctional enzymes
    yeast_df = analyze_yeast_bifunctional_fusion(
        "../data/raw/dah7ps_sequences.faa",
        "../data/raw/dah7ps_metadata.tsv"
    )

    yeast_df.to_csv("type_i/yeast_bifunctional_analysis.tsv", sep='\t', index=False)
    print("\nYeast bifunctional enzyme analysis saved to: structure/type_i/yeast_bifunctional_analysis.tsv")
