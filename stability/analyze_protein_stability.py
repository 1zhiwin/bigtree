#!/usr/bin/env python3
"""
Analyze protein stability features and predict effects of mutations.

Without 3D structures for all sequences, we use sequence-based stability predictors:
1. Amino acid composition (hydrophobic, charged, aromatic)
2. Instability index (Guruprasad method)
3. Aliphatic index (thermostability)
4. GRAVY (grand average of hydropathicity)
5. Charge distribution
6. Secondary structure propensity
"""

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import json
from collections import defaultdict

def calculate_stability_indices(sequence):
    """
    Calculate sequence-based stability indices.

    Returns dict with various stability metrics.
    """

    analyzer = ProteinAnalysis(sequence)

    # Basic properties
    length = len(sequence)
    molecular_weight = analyzer.molecular_weight()

    # Amino acid composition (get once to avoid multiple calls)
    aa_percent = analyzer.get_amino_acids_percent()

    # Instability index (< 40 = stable, > 40 = unstable)
    instability_index = analyzer.instability_index()

    # Aliphatic index (thermostability, higher = more stable)
    # Manual calculation: relative volume of aliphatic side chains (A, V, I, L)
    aliphatic_index = (aa_percent.get('A', 0) * 100 * 1.0 +
                      aa_percent.get('V', 0) * 100 * 2.9 +
                      aa_percent.get('I', 0) * 100 * 3.9 +
                      aa_percent.get('L', 0) * 100 * 3.9)

    # GRAVY (grand average of hydropathicity)
    # Negative = hydrophilic, Positive = hydrophobic
    gravy = analyzer.gravy()

    # Stability-relevant composition
    hydrophobic = sum(aa_percent.get(aa, 0) for aa in ['A', 'I', 'L', 'M', 'F', 'W', 'V']) * 100
    charged = sum(aa_percent.get(aa, 0) for aa in ['D', 'E', 'K', 'R']) * 100
    aromatic = sum(aa_percent.get(aa, 0) for aa in ['F', 'W', 'Y']) * 100
    polar = sum(aa_percent.get(aa, 0) for aa in ['S', 'T', 'N', 'Q']) * 100
    small = sum(aa_percent.get(aa, 0) for aa in ['A', 'G', 'S']) * 100

    # Proline content (can disrupt secondary structure)
    proline = aa_percent.get('P', 0) * 100

    # Cysteine content (disulfide bonds)
    cysteine = aa_percent.get('C', 0) * 100

    # Charge properties
    positive = sum(aa_percent.get(aa, 0) for aa in ['K', 'R']) * 100
    negative = sum(aa_percent.get(aa, 0) for aa in ['D', 'E']) * 100
    net_charge = positive - negative

    # Flexibility (Gly + Pro content)
    flexibility = (aa_percent.get('G', 0) + aa_percent.get('P', 0)) * 100

    return {
        'length': length,
        'molecular_weight': molecular_weight,
        'instability_index': instability_index,
        'aliphatic_index': aliphatic_index,
        'gravy': gravy,
        'hydrophobic_percent': hydrophobic,
        'charged_percent': charged,
        'aromatic_percent': aromatic,
        'polar_percent': polar,
        'small_percent': small,
        'proline_percent': proline,
        'cysteine_percent': cysteine,
        'positive_charge_percent': positive,
        'negative_charge_percent': negative,
        'net_charge': net_charge,
        'flexibility': flexibility,
        'predicted_stability': 'Stable' if instability_index < 40 else 'Unstable'
    }

def analyze_all_sequences(sequences_file, metadata_file, type_filter=None):
    """Analyze stability features for all sequences."""

    print(f"\n{'='*60}")
    print(f"SEQUENCE-BASED STABILITY ANALYSIS")
    if type_filter:
        print(f"Type: {type_filter}")
    print(f"{'='*60}\n")

    # Read metadata
    df_meta = pd.read_csv(metadata_file, sep='\t')

    # Read sequences
    sequences = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        sequences[record.id] = str(record.seq)

    results = []

    for _, row in df_meta.iterrows():
        if type_filter and row['DAH7PS_Type'] != type_filter:
            continue

        protein_id = row['Protein_ID']
        seq = sequences.get(protein_id)

        if not seq:
            continue

        # Calculate stability indices
        stability = calculate_stability_indices(seq)

        result = {
            'protein_id': protein_id,
            'organism': row['Organism'],
            'gene_name': row['Gene_Name'],
            'type': row['DAH7PS_Type'],
            'regulatory_specificity': row['Regulatory_Specificity'],
            **stability
        }

        results.append(result)

        # Print summary
        print(f"{protein_id}")
        print(f"  Organism: {row['Organism']}")
        print(f"  Specificity: {row['Regulatory_Specificity']}")
        print(f"  Length: {stability['length']} aa")
        print(f"  Instability Index: {stability['instability_index']:.2f} ({'Stable' if stability['instability_index'] < 40 else 'Unstable'})")
        print(f"  Aliphatic Index: {stability['aliphatic_index']:.2f} (thermostability)")
        print(f"  GRAVY: {stability['gravy']:.3f} ({'hydrophobic' if stability['gravy'] > 0 else 'hydrophilic'})")
        print(f"  Hydrophobic: {stability['hydrophobic_percent']:.1f}%")
        print(f"  Charged: {stability['charged_percent']:.1f}%")
        print(f"  Net charge: {stability['net_charge']:.2f}%")
        print()

    return pd.DataFrame(results)

def compare_domain_stability(sequences_file, domain_boundaries, type_name):
    """Compare stability features across domains."""

    print(f"\n{'='*60}")
    print(f"DOMAIN-SPECIFIC STABILITY: {type_name}")
    print(f"{'='*60}\n")

    # For Type I, we can analyze catalytic vs ACT domain
    # For Type II, analyze full TIM barrel

    # Read first sequence as example (they should be similar)
    sequences = list(SeqIO.parse(sequences_file, 'fasta'))

    domain_stats = {}

    for seq_record in sequences:
        seq = str(seq_record.seq)

        print(f"Sequence: {seq_record.id}")

        for domain_name, (start, end) in domain_boundaries.items():
            # Extract domain sequence (1-indexed to 0-indexed)
            domain_seq = seq[start-1:end]

            if len(domain_seq) < 10:  # Skip very short domains
                continue

            # Calculate stability for this domain
            stability = calculate_stability_indices(domain_seq)

            if domain_name not in domain_stats:
                domain_stats[domain_name] = []

            domain_stats[domain_name].append(stability)

            print(f"  {domain_name} ({start}-{end}, {len(domain_seq)} aa):")
            print(f"    Instability: {stability['instability_index']:.2f}")
            print(f"    Aliphatic: {stability['aliphatic_index']:.2f}")
            print(f"    GRAVY: {stability['gravy']:.3f}")
            print(f"    Hydrophobic: {stability['hydrophobic_percent']:.1f}%")

        print()
        break  # Just analyze first sequence for domain comparison

    return domain_stats

def predict_mutation_stability_effects(trait_specific_sites, alignment_file):
    """
    Predict stability effects of trait-specific mutations.

    Uses empirical rules:
    - Hydrophobic → Charged: Likely destabilizing (buried residue exposure)
    - Charged → Hydrophobic: Context-dependent (surface vs buried)
    - Proline introduction: Can disrupt α-helix
    - Proline removal: Can increase flexibility
    - Glycine → larger: Can create steric clash
    - Cysteine changes: Affects disulfide bonds
    """

    print(f"\n{'='*60}")
    print(f"MUTATION STABILITY EFFECT PREDICTIONS")
    print(f"{'='*60}\n")

    # Amino acid properties
    hydrophobic = set(['A', 'I', 'L', 'M', 'F', 'W', 'V'])
    charged = set(['D', 'E', 'K', 'R'])
    polar = set(['S', 'T', 'N', 'Q'])
    special = set(['G', 'P', 'C'])

    mutation_effects = []

    # Load trait-specific sites
    with open(trait_specific_sites) as f:
        sites_data = json.load(f)

    # Analyze top sites
    top_sites = sorted(sites_data, key=lambda x: x.get('specificity_score', 0), reverse=True)[:20]

    print(f"Analyzing top 20 trait-specific sites:\n")

    for site in top_sites:
        position = site['position']
        trait_aas = site['trait_amino_acids']

        print(f"Position {position}:")

        # Compare amino acids between traits
        for trait1, aas1 in trait_aas.items():
            for trait2, aas2 in trait_aas.items():
                if trait1 >= trait2:  # Avoid duplicates
                    continue

                # Find differences
                aa1_set = set(aas1)
                aa2_set = set(aas2)

                diff1 = aa1_set - aa2_set
                diff2 = aa2_set - aa1_set

                for aa1 in diff1:
                    for aa2 in diff2:
                        # Predict effect
                        effect = predict_single_mutation_effect(aa1, aa2)

                        mutation_effects.append({
                            'position': position,
                            'from_trait': trait1,
                            'to_trait': trait2,
                            'from_aa': aa1,
                            'to_aa': aa2,
                            'effect': effect['effect'],
                            'confidence': effect['confidence'],
                            'rationale': effect['rationale']
                        })

                        print(f"  {trait1} ({aa1}) → {trait2} ({aa2}): {effect['effect']} ({effect['confidence']}) - {effect['rationale']}")

        print()

    return mutation_effects

def predict_single_mutation_effect(from_aa, to_aa):
    """
    Predict stability effect of single amino acid substitution.

    Returns dict with effect, confidence, and rationale.
    """

    # Amino acid groups
    hydrophobic = set(['A', 'I', 'L', 'M', 'F', 'W', 'V'])
    charged_pos = set(['K', 'R'])
    charged_neg = set(['D', 'E'])
    charged = charged_pos | charged_neg
    polar = set(['S', 'T', 'N', 'Q'])
    aromatic = set(['F', 'W', 'Y'])
    small = set(['A', 'G', 'S'])

    # Size change (rough)
    size = {
        'G': 1, 'A': 2, 'S': 2, 'C': 2, 'P': 2, 'T': 3, 'V': 3,
        'D': 3, 'N': 3, 'E': 4, 'Q': 4, 'I': 4, 'L': 4, 'M': 4,
        'K': 5, 'R': 5, 'H': 4, 'F': 5, 'Y': 5, 'W': 6
    }

    size_change = size.get(to_aa, 3) - size.get(from_aa, 3)

    # Apply empirical rules

    # Rule 1: Hydrophobic to charged (likely destabilizing if buried)
    if from_aa in hydrophobic and to_aa in charged:
        return {
            'effect': 'Destabilizing',
            'confidence': 'Medium',
            'rationale': 'Hydrophobic→Charged (buried hydrophobic may become exposed charged)'
        }

    # Rule 2: Charged to hydrophobic (destabilizing if on surface)
    if from_aa in charged and to_aa in hydrophobic:
        return {
            'effect': 'Destabilizing',
            'confidence': 'Low',
            'rationale': 'Charged→Hydrophobic (surface charge loss, or buried hydrophobic gain)'
        }

    # Rule 3: Proline introduction (can disrupt secondary structure)
    if from_aa != 'P' and to_aa == 'P':
        return {
            'effect': 'Destabilizing',
            'confidence': 'Medium',
            'rationale': 'Proline introduction (disrupts α-helix, restricts backbone)'
        }

    # Rule 4: Proline removal (can increase flexibility)
    if from_aa == 'P' and to_aa != 'P':
        return {
            'effect': 'Stabilizing',
            'confidence': 'Low',
            'rationale': 'Proline removal (increases backbone flexibility)'
        }

    # Rule 5: Glycine to larger (can cause steric clash)
    if from_aa == 'G' and size_change > 2:
        return {
            'effect': 'Destabilizing',
            'confidence': 'Medium',
            'rationale': f'Glycine→{to_aa} (loss of flexibility, potential steric clash)'
        }

    # Rule 6: Cysteine changes (affects disulfide bonds)
    if from_aa == 'C' or to_aa == 'C':
        return {
            'effect': 'Destabilizing',
            'confidence': 'Medium',
            'rationale': 'Cysteine change (potential disulfide bond disruption)'
        }

    # Rule 7: Charge reversal (strong destabilization)
    if (from_aa in charged_pos and to_aa in charged_neg) or (from_aa in charged_neg and to_aa in charged_pos):
        return {
            'effect': 'Strongly Destabilizing',
            'confidence': 'High',
            'rationale': 'Charge reversal (electrostatic repulsion)'
        }

    # Rule 8: Conservative substitutions (same group)
    if (from_aa in hydrophobic and to_aa in hydrophobic) or \
       (from_aa in charged_pos and to_aa in charged_pos) or \
       (from_aa in charged_neg and to_aa in charged_neg) or \
       (from_aa in polar and to_aa in polar):
        return {
            'effect': 'Neutral',
            'confidence': 'High',
            'rationale': 'Conservative substitution (same physicochemical group)'
        }

    # Rule 9: Size change >3 (likely destabilizing)
    if abs(size_change) >= 3:
        return {
            'effect': 'Destabilizing',
            'confidence': 'Low',
            'rationale': f'Large size change ({from_aa}→{to_aa}, Δsize={size_change})'
        }

    # Default: Unknown
    return {
        'effect': 'Unknown',
        'confidence': 'Low',
        'rationale': f'{from_aa}→{to_aa} (insufficient information)'
    }

if __name__ == "__main__":

    print("\n" + "="*60)
    print("PROTEIN STABILITY ANALYSIS")
    print("="*60)

    # Analyze all Type I sequences
    print("\n### TYPE I SEQUENCES ###")
    type_i_stability = analyze_all_sequences(
        "../data/raw/dah7ps_sequences.faa",
        "../data/raw/dah7ps_metadata.tsv",
        type_filter="Type_I"
    )

    type_i_stability.to_csv("type_i/sequence_stability.tsv", sep='\t', index=False)

    # Domain-specific stability (Type I)
    type_i_domains = {
        'Catalytic_Domain': (1, 260),
        'ACT_Domain': (270, 346)
    }

    type_i_domain_stability = compare_domain_stability(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        type_i_domains,
        "Type I"
    )

    # Analyze all Type II sequences
    print("\n### TYPE II SEQUENCES ###")
    type_ii_stability = analyze_all_sequences(
        "../data/raw/dah7ps_sequences.faa",
        "../data/raw/dah7ps_metadata.tsv",
        type_filter="Type_II"
    )

    type_ii_stability.to_csv("type_ii/sequence_stability.tsv", sep='\t', index=False)

    # Predict mutation effects for trait-specific sites
    type_i_mutations = predict_mutation_stability_effects(
        "../selection/type_i/trait_specific_sites.json",
        "../msa/type_i/type_i_alignment_trimmed.faa"
    )

    # Save mutation predictions
    pd.DataFrame(type_i_mutations).to_csv("type_i/mutation_stability_predictions.tsv", sep='\t', index=False)

    # Summary statistics
    print(f"\n{'='*60}")
    print("STABILITY SUMMARY")
    print(f"{'='*60}\n")

    print("Type I:")
    print(f"  Mean Instability Index: {type_i_stability['instability_index'].mean():.2f}")
    print(f"  Stable sequences: {sum(type_i_stability['instability_index'] < 40)}/{len(type_i_stability)}")
    print(f"  Mean Aliphatic Index: {type_i_stability['aliphatic_index'].mean():.2f}")
    print(f"  Mean GRAVY: {type_i_stability['gravy'].mean():.3f}")

    print("\nType II:")
    print(f"  Mean Instability Index: {type_ii_stability['instability_index'].mean():.2f}")
    print(f"  Stable sequences: {sum(type_ii_stability['instability_index'] < 40)}/{len(type_ii_stability)}")
    print(f"  Mean Aliphatic Index: {type_ii_stability['aliphatic_index'].mean():.2f}")
    print(f"  Mean GRAVY: {type_ii_stability['gravy'].mean():.3f}")

    print("\nMutation Effects (Type I, top 20 sites):")
    mutation_df = pd.DataFrame(type_i_mutations)
    if not mutation_df.empty:
        effect_counts = mutation_df['effect'].value_counts()
        for effect, count in effect_counts.items():
            pct = count / len(mutation_df) * 100
            print(f"  {effect}: {count} ({pct:.1f}%)")

    print("\n" + "="*60)
    print("STABILITY ANALYSIS COMPLETE")
    print("="*60)
    print("\nResults saved to:")
    print("  - stability/type_i/sequence_stability.tsv")
    print("  - stability/type_i/mutation_stability_predictions.tsv")
    print("  - stability/type_ii/sequence_stability.tsv")
