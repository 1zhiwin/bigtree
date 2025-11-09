#!/usr/bin/env python3
"""
Analyze domain interactions and allosteric coupling in Type I DAH7PS.

Focus on:
1. Linker region properties (connects catalytic and ACT domains)
2. Co-evolution between domains
3. Trait-specific residues at domain interfaces
4. Allosteric coupling mechanisms
"""

from Bio import AlignIO, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import json
from scipy.stats import spearmanr
from itertools import combinations

def analyze_linker_region(alignment_file):
    """
    Analyze the linker region connecting catalytic and ACT domains.

    Linker: positions 261-269 (9 residues)
    Key properties: flexibility, conservation, secondary structure propensity
    """

    print(f"\n{'='*60}")
    print(f"LINKER REGION ANALYSIS")
    print(f"{'='*60}\n")

    alignment = AlignIO.read(alignment_file, 'fasta')

    # Extract linker region (positions 261-269, 0-indexed: 260-268)
    linker_start = 260
    linker_end = 269

    print(f"Linker region: positions {linker_start+1}-{linker_end}")
    print(f"Length: {linker_end - linker_start} residues\n")

    linker_sequences = []

    for record in alignment:
        linker_seq = str(record.seq)[linker_start:linker_end]
        linker_sequences.append({
            'id': record.id,
            'sequence': linker_seq
        })
        print(f"{record.id}: {linker_seq}")

    # Calculate consensus
    linker_columns = []
    for i in range(linker_end - linker_start):
        column = [seq['sequence'][i] for seq in linker_sequences if i < len(seq['sequence'])]
        column_clean = [aa for aa in column if aa != '-']
        if column_clean:
            most_common = Counter(column_clean).most_common(1)[0][0]
            linker_columns.append(most_common)
        else:
            linker_columns.append('-')

    consensus = ''.join(linker_columns)
    print(f"\nConsensus: {consensus}")

    # Analyze properties
    # Remove gaps for analysis
    consensus_clean = consensus.replace('-', '')
    if consensus_clean:
        analyzer = ProteinAnalysis(consensus_clean)

        # Flexibility markers
        gly_pro = consensus_clean.count('G') + consensus_clean.count('P')
        flexibility_score = gly_pro / len(consensus_clean) * 100

        print(f"\nLinker properties:")
        print(f"  Glycine + Proline: {gly_pro}/{len(consensus_clean)} ({flexibility_score:.1f}%)")
        print(f"  GRAVY: {analyzer.gravy():.3f}")

        # Amino acid composition
        aa_counts = Counter(consensus_clean)
        print(f"  Composition: {dict(aa_counts)}")

    # Conservation analysis
    print(f"\nConservation per position:")
    for i in range(linker_end - linker_start):
        column = [seq['sequence'][i] for seq in linker_sequences if i < len(seq['sequence'])]
        column_clean = [aa for aa in column if aa != '-']
        if column_clean:
            most_common_aa, count = Counter(column_clean).most_common(1)[0]
            conservation = count / len(column_clean) * 100
            print(f"  Position {linker_start + i + 1}: {most_common_aa} ({conservation:.1f}% conserved)")

    return linker_sequences

def identify_potential_interface_residues(alignment_file, trait_sites_file):
    """
    Identify residues that may be at the catalytic-ACT domain interface.

    Strategy:
    1. Residues near domain boundary (positions 250-270 and 270-280)
    2. Trait-specific residues in these regions
    3. Conserved residues that may mediate coupling
    """

    print(f"\n{'='*60}")
    print(f"POTENTIAL INTERFACE RESIDUE IDENTIFICATION")
    print(f"{'='*60}\n")

    # Read alignment
    alignment = AlignIO.read(alignment_file, 'fasta')

    # Read trait-specific sites
    with open(trait_sites_file) as f:
        trait_sites = json.load(f)

    # Define interface regions
    # Near C-terminus of catalytic domain
    catalytic_interface = range(250, 271)  # positions 250-270
    # Near N-terminus of ACT domain
    act_interface = range(270, 291)  # positions 270-290

    interface_regions = {
        'Catalytic_C_term': catalytic_interface,
        'ACT_N_term': act_interface
    }

    interface_residues = defaultdict(list)

    for region_name, positions in interface_regions.items():
        print(f"\n{region_name} (positions {min(positions)}-{max(positions)}):")

        for pos in positions:
            # Get column at this position (0-indexed)
            column = alignment[:, pos-1]

            # Count amino acids
            aa_counts = Counter([aa for aa in column if aa != '-'])

            if not aa_counts:
                continue

            # Check if this is a trait-specific site
            is_trait_specific = any(site['position'] == pos for site in trait_sites)

            # Calculate conservation
            total = sum(aa_counts.values())
            most_common_aa, max_count = aa_counts.most_common(1)[0]
            conservation = max_count / total

            # Identify potential interface residues
            if is_trait_specific or conservation > 0.7:
                interface_residues[region_name].append({
                    'position': pos,
                    'consensus': most_common_aa,
                    'conservation': conservation,
                    'trait_specific': is_trait_specific,
                    'amino_acids': dict(aa_counts)
                })

                marker = "[TRAIT-SPECIFIC]" if is_trait_specific else "[CONSERVED]"
                print(f"  Position {pos}: {most_common_aa} (cons={conservation:.2f}) {marker}")

    return interface_residues

def analyze_coevolution_between_domains(alignment_file):
    """
    Analyze coevolution between catalytic and ACT domains.

    Uses simple mutual information / correlation approach.
    Identifies pairs of positions that coevolve.
    """

    print(f"\n{'='*60}")
    print(f"CO-EVOLUTION ANALYSIS BETWEEN DOMAINS")
    print(f"{'='*60}\n")

    alignment = AlignIO.read(alignment_file, 'fasta')
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)

    # Define domains
    catalytic_domain = range(1, 261)  # 1-260
    act_domain = range(270, 347)  # 270-346

    print(f"Analyzing coevolution between:")
    print(f"  Catalytic domain: positions 1-260")
    print(f"  ACT domain: positions 270-346")
    print(f"  Total pairwise comparisons: {len(catalytic_domain) * len(act_domain)}")
    print(f"\nThis may take a while...\n")

    # Encode amino acids as numbers for correlation
    aa_to_num = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY-')}

    # Convert alignment to numerical matrix
    alignment_matrix = []
    for record in alignment:
        seq_encoded = [aa_to_num.get(aa, 20) for aa in str(record.seq)]
        alignment_matrix.append(seq_encoded)

    alignment_array = np.array(alignment_matrix)

    # Calculate correlations between domains
    coevolving_pairs = []

    # Sample subset for performance (every 10th position)
    cat_sample = list(catalytic_domain)[::10]
    act_sample = list(act_domain)[::5]

    print(f"Sampling positions (for performance):")
    print(f"  Catalytic: {len(cat_sample)} positions")
    print(f"  ACT: {len(act_sample)} positions")
    print(f"  Total pairs: {len(cat_sample) * len(act_sample)}\n")

    for cat_pos in cat_sample:
        for act_pos in act_sample:
            # Get columns
            cat_col = alignment_array[:, cat_pos-1]
            act_col = alignment_array[:, act_pos-1]

            # Remove positions with gaps in either column
            mask = (cat_col != 20) & (act_col != 20)
            if mask.sum() < 4:  # Need at least 4 sequences
                continue

            cat_col_clean = cat_col[mask]
            act_col_clean = act_col[mask]

            # Calculate Spearman correlation
            if len(set(cat_col_clean)) > 1 and len(set(act_col_clean)) > 1:
                corr, pval = spearmanr(cat_col_clean, act_col_clean)

                # High absolute correlation suggests coevolution
                if abs(corr) > 0.6 and pval < 0.05:
                    coevolving_pairs.append({
                        'cat_position': cat_pos,
                        'act_position': act_pos,
                        'correlation': corr,
                        'p_value': pval
                    })

    # Sort by absolute correlation
    coevolving_pairs.sort(key=lambda x: abs(x['correlation']), reverse=True)

    print(f"Coevolving pairs identified (|r| > 0.6, p < 0.05): {len(coevolving_pairs)}\n")

    if coevolving_pairs:
        print("Top 10 coevolving position pairs:")
        print(f"{'Catalytic':<12s} {'ACT':<8s} {'Correlation':<12s} {'P-value':<10s}")
        print("-" * 50)
        for pair in coevolving_pairs[:10]:
            print(f"{pair['cat_position']:<12d} {pair['act_position']:<8d} {pair['correlation']:< 12.3f} {pair['p_value']:<10.4f}")
    else:
        print("No significant coevolving pairs found.")
        print("Note: Small sample size (8 sequences) limits statistical power.")

    return coevolving_pairs

def correlate_interface_with_trait_evolution(interface_residues, trait_sites):
    """
    Correlate interface residues with trait-specific evolution.
    """

    print(f"\n{'='*60}")
    print(f"INTERFACE-TRAIT CORRELATION ANALYSIS")
    print(f"{'='*60}\n")

    # Read trait-specific sites
    with open(trait_sites) as f:
        trait_sites_data = json.load(f)

    trait_positions = set(site['position'] for site in trait_sites_data)

    # Count interface residues that are trait-specific
    interface_counts = {}

    for region_name, residues in interface_residues.items():
        trait_specific_count = sum(1 for r in residues if r['trait_specific'])
        total_count = len(residues)

        interface_counts[region_name] = {
            'total': total_count,
            'trait_specific': trait_specific_count,
            'percent': trait_specific_count / total_count * 100 if total_count > 0 else 0
        }

        print(f"{region_name}:")
        print(f"  Total interface residues: {total_count}")
        print(f"  Trait-specific: {trait_specific_count} ({interface_counts[region_name]['percent']:.1f}%)")
        print()

    return interface_counts

if __name__ == "__main__":

    print("\n" + "="*60)
    print("DOMAIN INTERACTION ANALYSIS (TYPE I)")
    print("="*60)

    # Analyze linker region
    linker_results = analyze_linker_region("../msa/type_i/type_i_alignment_trimmed.faa")

    # Identify interface residues
    interface_residues = identify_potential_interface_residues(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        "../selection/type_i/trait_specific_sites.json"
    )

    # Save interface residues
    interface_dict = {k: v for k, v in interface_residues.items()}
    with open("type_i/interface_residues.json", 'w') as f:
        json.dump(interface_dict, f, indent=2)

    # Coevolution analysis
    coevolving_pairs = analyze_coevolution_between_domains("../msa/type_i/type_i_alignment_trimmed.faa")

    if coevolving_pairs:
        pd.DataFrame(coevolving_pairs).to_csv("coevolution/coevolving_pairs.tsv", sep='\t', index=False)

    # Correlate interface with trait evolution
    interface_trait_corr = correlate_interface_with_trait_evolution(
        interface_residues,
        "../selection/type_i/trait_specific_sites.json"
    )

    with open("type_i/interface_trait_correlation.json", 'w') as f:
        json.dump(interface_trait_corr, f, indent=2)

    print("\n" + "="*60)
    print("DOMAIN INTERACTION ANALYSIS COMPLETE")
    print("="*60)
    print("\nResults saved to:")
    print("  - interactions/type_i/interface_residues.json")
    print("  - interactions/type_i/interface_trait_correlation.json")
    print("  - interactions/coevolution/coevolving_pairs.tsv (if pairs found)")
