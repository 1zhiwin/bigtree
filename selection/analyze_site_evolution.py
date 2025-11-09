#!/usr/bin/env python3
"""
Analyze site-specific evolutionary patterns in DAH7PS protein sequences.

Since we only have protein sequences (not DNA), we cannot use codon-based
dN/dS methods. Instead, we analyze:
1. Site-specific evolutionary rates based on amino acid variability
2. Conservation patterns across lineages
3. Lineage-specific substitution patterns
4. Correlation with functional domains
"""

import numpy as np
from Bio import AlignIO
from collections import defaultdict, Counter
import pandas as pd
import json

def calculate_site_specific_rates(alignment_file, type_name):
    """
    Calculate evolutionary rate for each site based on amino acid variability.

    Rate metrics:
    - Number of different amino acids (diversity)
    - Shannon entropy (information content)
    - Normalized variability score
    """

    print(f"\n{'='*60}")
    print(f"SITE-SPECIFIC EVOLUTIONARY RATE ANALYSIS: {type_name}")
    print(f"{'='*60}\n")

    alignment = AlignIO.read(alignment_file, 'fasta')
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)

    site_rates = []

    for i in range(alignment_length):
        column = alignment[:, i]

        # Count amino acids (exclude gaps)
        aa_counts = Counter([aa for aa in column if aa != '-'])
        total_non_gap = sum(aa_counts.values())

        if total_non_gap == 0:
            # All gaps
            site_rates.append({
                'position': i + 1,
                'num_amino_acids': 0,
                'entropy': 0,
                'variability': 0,
                'gap_frequency': 1.0,
                'rate_category': 'gap'
            })
            continue

        # Calculate metrics
        num_aa_types = len(aa_counts)
        gap_frequency = 1 - (total_non_gap / num_sequences)

        # Shannon entropy
        entropy = 0
        for count in aa_counts.values():
            p = count / total_non_gap
            if p > 0:
                entropy -= p * np.log2(p)

        # Normalize entropy by maximum possible (log2 of 20 amino acids)
        max_entropy = np.log2(min(20, total_non_gap))
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

        # Variability score (0-1): combines diversity and entropy
        diversity_score = (num_aa_types - 1) / 19  # Normalize to 0-1 (max 20 types)
        variability = (diversity_score + normalized_entropy) / 2

        # Categorize rate
        if variability < 0.2:
            rate_category = 'conserved'
        elif variability < 0.5:
            rate_category = 'moderate'
        elif variability < 0.7:
            rate_category = 'variable'
        else:
            rate_category = 'highly_variable'

        site_rates.append({
            'position': i + 1,
            'num_amino_acids': num_aa_types,
            'entropy': entropy,
            'normalized_entropy': normalized_entropy,
            'variability': variability,
            'gap_frequency': gap_frequency,
            'rate_category': rate_category
        })

    # Summary statistics
    variability_scores = [s['variability'] for s in site_rates if s['rate_category'] != 'gap']

    print(f"Alignment length: {alignment_length} positions")
    print(f"Number of sequences: {num_sequences}")
    print(f"\nVariability distribution:")
    print(f"  Mean variability: {np.mean(variability_scores):.3f}")
    print(f"  Median variability: {np.median(variability_scores):.3f}")
    print(f"  Std deviation: {np.std(variability_scores):.3f}")

    # Count by category
    category_counts = Counter([s['rate_category'] for s in site_rates])
    print(f"\nSite categories:")
    for category, count in sorted(category_counts.items()):
        pct = count / alignment_length * 100
        print(f"  {category:20s}: {count:4d} ({pct:5.1f}%)")

    return site_rates

def identify_rapidly_evolving_sites(site_rates, threshold=0.7):
    """Identify sites evolving rapidly (potential positive selection sites)."""

    rapidly_evolving = [s for s in site_rates
                       if s.get('variability', 0) >= threshold]

    return rapidly_evolving

def identify_conserved_sites(site_rates, threshold=0.2):
    """Identify highly conserved sites (potential functional constraints)."""

    conserved = [s for s in site_rates
                if s.get('variability', 0) <= threshold and s['rate_category'] != 'gap']

    return conserved

def compare_domain_evolution(site_rates, domain_boundaries, type_name):
    """Compare evolutionary rates between functional domains."""

    print(f"\n{'='*60}")
    print(f"DOMAIN-SPECIFIC EVOLUTIONARY RATES: {type_name}")
    print(f"{'='*60}\n")

    domain_stats = {}

    for domain_name, (start, end) in domain_boundaries.items():
        domain_sites = [s for s in site_rates
                       if start <= s['position'] <= end and s['rate_category'] != 'gap']

        if not domain_sites:
            continue

        variabilities = [s['variability'] for s in domain_sites]
        entropies = [s['normalized_entropy'] for s in domain_sites]

        domain_stats[domain_name] = {
            'num_sites': len(domain_sites),
            'mean_variability': np.mean(variabilities),
            'median_variability': np.median(variabilities),
            'std_variability': np.std(variabilities),
            'mean_entropy': np.mean(entropies),
            'conserved_sites': len([s for s in domain_sites if s['rate_category'] == 'conserved']),
            'variable_sites': len([s for s in domain_sites if s['rate_category'] in ['variable', 'highly_variable']])
        }

        print(f"{domain_name}:")
        print(f"  Positions: {start}-{end}")
        print(f"  Sites analyzed: {len(domain_sites)}")
        print(f"  Mean variability: {domain_stats[domain_name]['mean_variability']:.3f}")
        print(f"  Conserved sites: {domain_stats[domain_name]['conserved_sites']} ({domain_stats[domain_name]['conserved_sites']/len(domain_sites)*100:.1f}%)")
        print(f"  Variable sites: {domain_stats[domain_name]['variable_sites']} ({domain_stats[domain_name]['variable_sites']/len(domain_sites)*100:.1f}%)")
        print()

    return domain_stats

def analyze_lineage_specific_substitutions(alignment_file, trait_file, type_name):
    """
    Analyze substitutions specific to lineages with different regulatory traits.
    """

    print(f"\n{'='*60}")
    print(f"LINEAGE-SPECIFIC SUBSTITUTION ANALYSIS: {type_name}")
    print(f"{'='*60}\n")

    alignment = AlignIO.read(alignment_file, 'fasta')

    # Read traits
    traits_df = pd.read_csv(trait_file, sep='\t')
    traits = dict(zip(traits_df['Taxon'], traits_df['Trait']))

    # Group sequences by trait
    trait_groups = defaultdict(list)
    for record in alignment:
        trait = traits.get(record.id, 'Unknown')
        trait_groups[trait].append(str(record.seq))

    # Find trait-specific substitutions
    alignment_length = alignment.get_alignment_length()
    trait_specific_sites = []

    for i in range(alignment_length):
        # Get amino acids for each trait group at this position
        trait_aas = {}
        for trait, seqs in trait_groups.items():
            if trait == 'Unknown':
                continue
            aas = set(seq[i] for seq in seqs if seq[i] != '-')
            if aas:
                trait_aas[trait] = aas

        # Check if this site has trait-specific amino acids
        if len(trait_aas) >= 2:
            # Check if traits have distinct amino acids
            all_aas = [aa for aas in trait_aas.values() for aa in aas]
            if len(set(all_aas)) >= 2:
                # Calculate specificity score: how trait-specific is this position?
                # Perfect specificity: each trait has unique amino acid(s)
                unique_count = sum(1 for trait, aas in trait_aas.items()
                                  if not any(aa in aas for other_trait, other_aas in trait_aas.items()
                                            if other_trait != trait for aa in other_aas))

                trait_specific_sites.append({
                    'position': i + 1,
                    'trait_amino_acids': {k: list(v) for k, v in trait_aas.items()},
                    'num_traits': len(trait_aas),
                    'total_amino_acids': len(set(all_aas)),
                    'specificity_score': unique_count / len(trait_aas) if trait_aas else 0
                })

    print(f"Trait-specific sites identified: {len(trait_specific_sites)}")

    if trait_specific_sites:
        print(f"\nTop 10 trait-specific sites (by specificity score):")
        sorted_sites = sorted(trait_specific_sites, key=lambda x: x['specificity_score'], reverse=True)
        print(f"{'Position':<10s} {'Specificity':<12s} {'Traits':<10s} {'Amino Acids':<30s}")
        print("-" * 80)
        for site in sorted_sites[:10]:
            pos = site['position']
            spec = site['specificity_score']
            num_traits = site['num_traits']
            aa_summary = '; '.join([f"{trait}:{','.join(aas)}"
                                   for trait, aas in site['trait_amino_acids'].items()])
            print(f"{pos:<10d} {spec:<12.3f} {num_traits:<10d} {aa_summary[:50]:<30s}")

    return trait_specific_sites

def correlate_with_conservation(site_rates, conservation_file):
    """
    Correlate site-specific rates with conservation scores from Phase 8.
    """

    with open(conservation_file) as f:
        conservation_data = json.load(f)

    conservation_scores = conservation_data['conservation_scores']

    # Match positions
    correlations = []
    for i, site_rate in enumerate(site_rates):
        if i < len(conservation_scores):
            correlations.append({
                'position': site_rate['position'],
                'variability': site_rate.get('variability', 0),
                'conservation': conservation_scores[i],
                'rate_category': site_rate['rate_category']
            })

    # Calculate correlation coefficient
    variabilities = [c['variability'] for c in correlations if c['rate_category'] != 'gap']
    conservations = [c['conservation'] for c in correlations if c['rate_category'] != 'gap']

    if variabilities and conservations:
        correlation = np.corrcoef(variabilities, conservations)[0, 1]
        print(f"\nCorrelation between variability and conservation: {correlation:.3f}")
        print("(Expected negative correlation: high conservation = low variability)")

    return correlations

if __name__ == "__main__":

    print("\n" + "="*60)
    print("SITE-SPECIFIC EVOLUTIONARY ANALYSIS")
    print("="*60)

    # Type I analysis
    print("\n### TYPE I DAH7PS ###")

    type_i_rates = calculate_site_specific_rates(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        "Type I"
    )

    # Save results
    pd.DataFrame(type_i_rates).to_csv("type_i/site_rates.tsv", sep='\t', index=False)

    # Identify rapidly evolving sites
    rapidly_evolving_i = identify_rapidly_evolving_sites(type_i_rates, threshold=0.7)
    conserved_i = identify_conserved_sites(type_i_rates, threshold=0.2)

    print(f"\nRapidly evolving sites (variability >= 0.7): {len(rapidly_evolving_i)}")
    print(f"Highly conserved sites (variability <= 0.2): {len(conserved_i)}")

    # Domain-specific analysis
    type_i_domains = {
        'Catalytic_Domain': (1, 260),
        'Linker': (261, 269),
        'ACT_Domain': (270, 346)
    }

    type_i_domain_stats = compare_domain_evolution(type_i_rates, type_i_domains, "Type I")

    # Save domain stats
    with open("type_i/domain_evolution_stats.json", 'w') as f:
        json.dump(type_i_domain_stats, f, indent=2)

    # Lineage-specific substitutions
    type_i_trait_sites = analyze_lineage_specific_substitutions(
        "../msa/type_i/type_i_alignment_trimmed.faa",
        "../traits/type_i/type_i_traits.tsv",
        "Type I"
    )

    with open("type_i/trait_specific_sites.json", 'w') as f:
        json.dump(type_i_trait_sites, f, indent=2)

    # Correlation with conservation
    print(f"\n{'='*60}")
    print("TYPE I: Correlation with Phase 8 Conservation")
    print(f"{'='*60}")
    type_i_correlation = correlate_with_conservation(
        type_i_rates,
        "../structure/conservation/type_i_conservation.json"
    )

    # Type II analysis
    print("\n\n### TYPE II DAH7PS ###")

    type_ii_rates = calculate_site_specific_rates(
        "../msa/type_ii/type_ii_alignment_trimmed.faa",
        "Type II"
    )

    pd.DataFrame(type_ii_rates).to_csv("type_ii/site_rates.tsv", sep='\t', index=False)

    rapidly_evolving_ii = identify_rapidly_evolving_sites(type_ii_rates, threshold=0.7)
    conserved_ii = identify_conserved_sites(type_ii_rates, threshold=0.2)

    print(f"\nRapidly evolving sites (variability >= 0.7): {len(rapidly_evolving_ii)}")
    print(f"Highly conserved sites (variability <= 0.2): {len(conserved_ii)}")

    # Domain-specific analysis (Type II doesn't have ACT domains)
    type_ii_domains = {
        'TIM_Barrel_Core': (1, 473)  # Entire catalytic domain
    }

    type_ii_domain_stats = compare_domain_evolution(type_ii_rates, type_ii_domains, "Type II")

    with open("type_ii/domain_evolution_stats.json", 'w') as f:
        json.dump(type_ii_domain_stats, f, indent=2)

    # Lineage-specific (limited - only plastid vs non-plastid)
    type_ii_trait_sites = analyze_lineage_specific_substitutions(
        "../msa/type_ii/type_ii_alignment_trimmed.faa",
        "../traits/type_ii/type_ii_traits.tsv",
        "Type II"
    )

    with open("type_ii/trait_specific_sites.json", 'w') as f:
        json.dump(type_ii_trait_sites, f, indent=2)

    # Correlation with conservation
    print(f"\n{'='*60}")
    print("TYPE II: Correlation with Phase 8 Conservation")
    print(f"{'='*60}")
    type_ii_correlation = correlate_with_conservation(
        type_ii_rates,
        "../structure/conservation/type_ii_conservation.json"
    )

    print("\n" + "="*60)
    print("SITE-SPECIFIC ANALYSIS COMPLETE")
    print("="*60)
    print("\nResults saved to:")
    print("  - selection/type_i/site_rates.tsv")
    print("  - selection/type_i/domain_evolution_stats.json")
    print("  - selection/type_i/trait_specific_sites.json")
    print("  - selection/type_ii/site_rates.tsv")
    print("  - selection/type_ii/domain_evolution_stats.json")
    print("  - selection/type_ii/trait_specific_sites.json")
