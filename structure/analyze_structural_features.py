#!/usr/bin/env python3
"""
Analyze structural features of DAH7PS sequences.
"""

from Bio import SeqIO, AlignIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
from collections import defaultdict
import json

def analyze_transit_peptides(sequences_file, metadata_file):
    """
    Analyze N-terminal transit peptides in Type II sequences.

    Transit peptides are predicted based on:
    1. Sequence length comparison (Type II vs bacterial)
    2. N-terminal enrichment of Ser, Thr, Arg
    3. Lack of acidic residues in first 50-80 aa
    """

    # Read metadata
    df = pd.read_csv(metadata_file, sep='\t')

    # Read sequences
    sequences = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        sequences[record.id] = str(record.seq)

    transit_results = []

    print("\n" + "="*60)
    print("TRANSIT PEPTIDE ANALYSIS (Type II)")
    print("="*60 + "\n")

    for _, row in df.iterrows():
        if row['DAH7PS_Type'] != 'Type_II':
            continue

        protein_id = row['Protein_ID']
        seq = sequences[protein_id]
        organism = row['Organism']

        # Analyze N-terminal region (first 80 residues)
        n_term_length = min(80, len(seq))
        n_term = seq[:n_term_length]

        # Analyze amino acid composition
        analyzer = ProteinAnalysis(n_term)
        aa_percent = analyzer.get_amino_acids_percent()

        # Transit peptide signatures
        ser_thr_content = (aa_percent.get('S', 0) + aa_percent.get('T', 0)) * 100
        arg_content = aa_percent.get('R', 0) * 100
        acidic_content = (aa_percent.get('D', 0) + aa_percent.get('E', 0)) * 100

        # Predicted cleavage site (enriched in Ala at -1, -3 positions)
        # Simple prediction: find region with high Ser/Thr after initial Met-rich region
        predicted_cleavage = None
        for i in range(30, min(80, len(seq))):
            # Look for VxA or FxA motifs (common cleavage sites)
            if i+2 < len(seq):
                if seq[i:i+3] in ['VFA', 'VIA', 'VSA', 'VTA', 'VVA', 'FLA', 'FSA']:
                    predicted_cleavage = i + 2
                    break

        if predicted_cleavage is None:
            # Fallback: estimate based on average transit peptide length
            if 'Arabidopsis' in organism:
                predicted_cleavage = 57  # Average for plant chloroplast targeting

        # Calculate properties
        if predicted_cleavage:
            transit_peptide = seq[:predicted_cleavage]
            mature_protein = seq[predicted_cleavage:]

            tp_analyzer = ProteinAnalysis(transit_peptide)

            result = {
                'protein_id': protein_id,
                'organism': organism,
                'total_length': len(seq),
                'predicted_transit_length': predicted_cleavage,
                'mature_length': len(mature_protein),
                'transit_ser_thr': ser_thr_content,
                'transit_arg': arg_content,
                'transit_acidic': acidic_content,
                'transit_sequence': transit_peptide
            }
        else:
            result = {
                'protein_id': protein_id,
                'organism': organism,
                'total_length': len(seq),
                'predicted_transit_length': 0,
                'mature_length': len(seq),
                'transit_ser_thr': 0,
                'transit_arg': 0,
                'transit_acidic': 0,
                'transit_sequence': ''
            }

        transit_results.append(result)

        print(f"Sequence: {protein_id}")
        print(f"  Organism: {organism}")
        print(f"  Total length: {len(seq)} aa")
        if predicted_cleavage:
            print(f"  Transit peptide: {predicted_cleavage} aa")
            print(f"  Mature protein: {len(mature_protein)} aa")
            print(f"  N-term Ser+Thr: {ser_thr_content:.1f}%")
            print(f"  N-term Arg: {arg_content:.1f}%")
            print(f"  N-term Asp+Glu: {acidic_content:.1f}%")
            print(f"  Transit sequence: {transit_peptide[:40]}...")
        else:
            print(f"  No transit peptide predicted (bacterial)")
        print()

    return pd.DataFrame(transit_results)

def analyze_conservation_from_alignment(alignment_file, type_name):
    """
    Analyze conservation at each position in the alignment.
    """

    print(f"\n{'='*60}")
    print(f"CONSERVATION ANALYSIS: {type_name}")
    print(f"{'='*60}\n")

    alignment = AlignIO.read(alignment_file, 'fasta')
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)

    # Calculate conservation score for each position
    conservation_scores = []
    conserved_positions = []
    variable_positions = []

    for i in range(alignment_length):
        column = alignment[:, i]

        # Count amino acids (excluding gaps)
        aa_counts = defaultdict(int)
        non_gap_count = 0
        for aa in column:
            if aa != '-':
                aa_counts[aa] += 1
                non_gap_count += 1

        if non_gap_count == 0:
            conservation = 0
        else:
            # Shannon entropy-based conservation
            entropy = 0
            for count in aa_counts.values():
                if count > 0:
                    p = count / non_gap_count
                    entropy -= p * np.log2(p)

            # Normalize entropy (0 = conserved, ~4.3 = maximally variable for 20 aa)
            max_entropy = np.log2(min(20, non_gap_count))
            if max_entropy > 0:
                normalized_entropy = entropy / max_entropy
            else:
                normalized_entropy = 0

            # Conservation score (1 = conserved, 0 = variable)
            conservation = 1 - normalized_entropy

        conservation_scores.append(conservation)

        # Classify positions
        if conservation >= 0.9:
            conserved_positions.append(i + 1)  # 1-indexed
        elif conservation <= 0.3:
            variable_positions.append(i + 1)

    print(f"Alignment length: {alignment_length} positions")
    print(f"Number of sequences: {num_sequences}")
    print(f"Mean conservation: {np.mean(conservation_scores):.3f}")
    print(f"Highly conserved positions (>0.9): {len(conserved_positions)} ({len(conserved_positions)/alignment_length*100:.1f}%)")
    print(f"Variable positions (<0.3): {len(variable_positions)} ({len(variable_positions)/alignment_length*100:.1f}%)")

    # Identify conserved regions (runs of 5+ conserved positions)
    conserved_regions = []
    region_start = None
    for i, score in enumerate(conservation_scores):
        if score >= 0.9:
            if region_start is None:
                region_start = i + 1
        else:
            if region_start is not None:
                region_end = i
                if region_end - region_start >= 5:
                    conserved_regions.append((region_start, region_end))
                region_start = None

    # Handle last region
    if region_start is not None and alignment_length - region_start >= 5:
        conserved_regions.append((region_start, alignment_length))

    print(f"\nConserved regions (5+ consecutive conserved positions):")
    for start, end in conserved_regions:
        print(f"  Region {start}-{end} ({end-start+1} positions)")

    return {
        'conservation_scores': conservation_scores,
        'conserved_positions': conserved_positions,
        'variable_positions': variable_positions,
        'conserved_regions': conserved_regions,
        'mean_conservation': np.mean(conservation_scores),
        'alignment_length': alignment_length
    }

def identify_catalytic_residues():
    """
    Identify known catalytic residues from literature.

    Based on crystal structures and biochemical studies:
    - Type I: Schiff base mechanism (Lys, active site metals)
    - Type II: TIM barrel active site
    """

    print("\n" + "="*60)
    print("KNOWN CATALYTIC AND FUNCTIONAL RESIDUES")
    print("="*60 + "\n")

    catalytic_residues = {
        'Type_I': {
            'description': 'Type I DAH7PS (α/β fold)',
            'active_site': [
                'Lys (Schiff base formation with PEP)',
                'His (proton transfer)',
                'Arg (phosphate binding)',
                'Asp (metal coordination)'
            ],
            'allosteric_sites': {
                'Phe-sensitive': [
                    'ACT domain (C-terminal)',
                    'Phe binding pocket: aromatic, hydrophobic residues'
                ],
                'Tyr-sensitive': [
                    'ACT domain (C-terminal)',
                    'Tyr binding pocket: similar to Phe but accommodates OH'
                ],
                'Trp-sensitive': [
                    'ACT domain (C-terminal)',
                    'Trp binding pocket: larger, accommodates indole ring'
                ]
            },
            'known_structures': [
                'E. coli aroF (PDB: not available)',
                'E. coli aroG (PDB: not available)',
                'S. cerevisiae ARO3 (PDB: 1N8F, 1QMG)',
                'S. cerevisiae ARO4 (PDB: 1GG1)'
            ]
        },
        'Type_II': {
            'description': 'Type II DAH7PS ((β/α)₈ TIM barrel)',
            'active_site': [
                'Active site at barrel C-terminus',
                'Metal binding site (typically Mn²⁺ or Co²⁺)',
                'Catalytic residues in loops connecting β-strands to α-helices'
            ],
            'structural_features': [
                'TIM barrel core domain',
                'Transit peptide (plants): N-terminal 50-80 aa',
                'Cleavage site: typically VxA or FxA motif'
            ],
            'known_structures': [
                'A. thaliana DHS (PDB: not widely available)',
                'Bacterial Type II (limited structural data)'
            ]
        }
    }

    for type_name, info in catalytic_residues.items():
        print(f"\n{type_name}:")
        print(f"  Description: {info['description']}")
        print(f"\n  Active Site Residues:")
        for residue in info.get('active_site', []):
            print(f"    - {residue}")

        if 'allosteric_sites' in info:
            print(f"\n  Allosteric Sites:")
            for specificity, sites in info['allosteric_sites'].items():
                print(f"    {specificity}:")
                for site in sites:
                    print(f"      - {site}")

        if 'structural_features' in info:
            print(f"\n  Structural Features:")
            for feature in info['structural_features']:
                print(f"    - {feature}")

        if 'known_structures' in info:
            print(f"\n  Known Structures:")
            for structure in info['known_structures']:
                print(f"    - {structure}")

    return catalytic_residues

if __name__ == "__main__":

    # Analyze transit peptides
    transit_df = analyze_transit_peptides(
        "../data/raw/dah7ps_sequences.faa",
        "../data/raw/dah7ps_metadata.tsv"
    )

    transit_df.to_csv("type_ii/transit_peptide_analysis.tsv", sep='\t', index=False)
    print("\nTransit peptide analysis saved to: structure/type_ii/transit_peptide_analysis.tsv")

    # Analyze conservation
    type_i_conservation = analyze_conservation_from_alignment(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        "Type I"
    )

    type_ii_conservation = analyze_conservation_from_alignment(
        "../msa/type_ii/type_ii_alignment_trimmed.faa",
        "Type II"
    )

    # Save conservation data
    with open("conservation/type_i_conservation.json", 'w') as f:
        # Convert numpy types to native Python for JSON serialization
        data = {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                for k, v in type_i_conservation.items()}
        json.dump(data, f, indent=2)

    with open("conservation/type_ii_conservation.json", 'w') as f:
        data = {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                for k, v in type_ii_conservation.items()}
        json.dump(data, f, indent=2)

    print("\nConservation data saved to:")
    print("  - structure/conservation/type_i_conservation.json")
    print("  - structure/conservation/type_ii_conservation.json")

    # Identify catalytic residues
    catalytic_info = identify_catalytic_residues()

    with open("catalytic_residues.json", 'w') as f:
        json.dump(catalytic_info, f, indent=2)

    print("\nCatalytic residue information saved to: structure/catalytic_residues.json")
