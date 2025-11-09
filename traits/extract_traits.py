#!/usr/bin/env python3
"""
Extract regulatory specificity traits from metadata for ancestral trait reconstruction.
"""

import pandas as pd
from collections import defaultdict

def extract_regulatory_traits(metadata_file):
    """Extract regulatory specificity for each sequence."""

    # Read metadata
    df = pd.read_csv(metadata_file, sep='\t')

    # Create trait dictionaries
    type_i_traits = {}
    type_ii_traits = {}

    for _, row in df.iterrows():
        protein_id = row['Protein_ID']
        dah7ps_type = row['DAH7PS_Type']
        reg_spec = row['Regulatory_Specificity']

        # Simplify trait categories
        if 'Phe-sensitive' in reg_spec:
            trait = 'Phe'
        elif 'Tyr-sensitive' in reg_spec:
            trait = 'Tyr'
        elif 'Trp-sensitive' in reg_spec:
            trait = 'Trp'
        elif 'Plastid-targeted' in reg_spec:
            trait = 'Plastid'
        elif 'bifunctional' in reg_spec.lower():
            if 'Phe' in reg_spec:
                trait = 'Phe_bifunctional'
            elif 'Tyr' in reg_spec:
                trait = 'Tyr_bifunctional'
            else:
                trait = 'Unknown'
        else:
            trait = 'Unknown'

        if dah7ps_type == 'Type_I':
            type_i_traits[protein_id] = trait
        else:
            type_ii_traits[protein_id] = trait

    return type_i_traits, type_ii_traits

def write_trait_file(traits_dict, output_file):
    """Write traits in format for ancestral state reconstruction."""

    with open(output_file, 'w') as f:
        f.write("Taxon\tTrait\n")
        for taxon, trait in sorted(traits_dict.items()):
            f.write(f"{taxon}\t{trait}\n")

def summarize_traits(traits_dict, type_name):
    """Summarize trait distribution."""

    trait_counts = defaultdict(int)
    for trait in traits_dict.values():
        trait_counts[trait] += 1

    print(f"\n{type_name} Trait Distribution:")
    print(f"{'='*50}")
    for trait, count in sorted(trait_counts.items()):
        pct = count / len(traits_dict) * 100
        print(f"{trait:20s}: {count:2d} ({pct:5.1f}%)")
    print(f"{'='*50}")
    print(f"Total sequences: {len(traits_dict)}")

if __name__ == "__main__":

    metadata_file = "../data/raw/dah7ps_metadata.tsv"

    # Extract traits
    print("Extracting regulatory specificity traits...")
    type_i_traits, type_ii_traits = extract_regulatory_traits(metadata_file)

    # Write trait files
    write_trait_file(type_i_traits, "type_i/type_i_traits.tsv")
    write_trait_file(type_ii_traits, "type_ii/type_ii_traits.tsv")

    # Summarize
    summarize_traits(type_i_traits, "Type I")
    summarize_traits(type_ii_traits, "Type II")

    print("\nTrait files written:")
    print("  - traits/type_i/type_i_traits.tsv")
    print("  - traits/type_ii/type_ii_traits.tsv")
