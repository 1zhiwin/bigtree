#!/usr/bin/env python3
"""
Perform ancestral trait reconstruction for regulatory specificity using parsimony.
"""

import pandas as pd
from Bio import Phylo
from collections import defaultdict, Counter
import json

def read_traits(trait_file):
    """Read trait assignments for tip sequences."""
    df = pd.read_csv(trait_file, sep='\t')
    traits = {}
    for _, row in df.iterrows():
        traits[row['Taxon']] = row['Trait']
    return traits

def fitch_parsimony(tree, traits):
    """
    Perform Fitch parsimony to reconstruct ancestral states.

    This implements the Fitch algorithm:
    1. Post-order traversal (tips to root): compute possible states
    2. Pre-order traversal (root to tips): resolve ambiguities
    """

    # Initialize node states
    node_states = {}
    node_final_states = {}

    # Post-order traversal: compute possible states for each node
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            # Tip: assign observed trait
            trait = traits.get(clade.name, 'Unknown')
            node_states[clade] = {trait}
        else:
            # Internal node: intersect children's states
            children = list(clade.clades)
            if len(children) >= 2:
                child_states = [node_states.get(child, set()) for child in children]

                # Fitch intersection
                intersection = child_states[0].copy()
                for states in child_states[1:]:
                    intersection &= states

                if intersection:
                    # If intersection not empty, use it
                    node_states[clade] = intersection
                else:
                    # If empty, use union
                    node_states[clade] = set().union(*child_states)

    # Pre-order traversal: resolve ambiguities
    root = tree.root
    # For root, choose most parsimonious state (most common in descendants)
    if len(node_states[root]) > 1:
        # Count occurrences in tips
        tip_traits = list(traits.values())
        trait_counts = Counter(tip_traits)
        # Choose most common trait from root's possible states
        best_state = max(node_states[root], key=lambda x: trait_counts.get(x, 0))
        node_final_states[root] = best_state
    else:
        node_final_states[root] = list(node_states[root])[0]

    # Pre-order: propagate down
    for clade in tree.find_clades(order='preorder'):
        if clade == root:
            continue

        parent = tree.get_path(clade)[-2] if len(tree.get_path(clade)) > 1 else root
        parent_state = node_final_states.get(parent)

        if len(node_states[clade]) == 1:
            node_final_states[clade] = list(node_states[clade])[0]
        else:
            # If parent's state is in possible states, use it (minimize changes)
            if parent_state in node_states[clade]:
                node_final_states[clade] = parent_state
            else:
                # Otherwise, pick most common
                tip_traits = list(traits.values())
                trait_counts = Counter(tip_traits)
                best_state = max(node_states[clade], key=lambda x: trait_counts.get(x, 0))
                node_final_states[clade] = best_state

    return node_states, node_final_states

def count_trait_changes(tree, node_final_states):
    """Count the number of trait changes along the tree."""
    changes = 0
    change_events = []

    for clade in tree.find_clades():
        if clade.is_terminal() or clade == tree.root:
            continue

        parent = None
        for c in tree.find_clades():
            if clade in c.clades:
                parent = c
                break

        if parent:
            parent_state = node_final_states[parent]
            clade_state = node_final_states[clade]
            if parent_state != clade_state:
                changes += 1
                change_events.append({
                    'from': parent_state,
                    'to': clade_state,
                    'node': clade.name if clade.name else 'internal'
                })

    return changes, change_events

def analyze_trait_evolution(tree, traits, type_name):
    """Perform full ancestral trait reconstruction and analysis."""

    print(f"\n{'='*60}")
    print(f"Ancestral Trait Reconstruction: {type_name}")
    print(f"{'='*60}\n")

    # Perform Fitch parsimony
    print("Running Fitch parsimony algorithm...")
    node_possible_states, node_final_states = fitch_parsimony(tree, traits)

    # Assign names to internal nodes
    node_counter = 1
    node_names = {}
    for clade in tree.find_clades():
        if not clade.is_terminal() and not clade.name:
            clade.name = f"Node{node_counter}"
            node_counter += 1
        node_names[clade] = clade.name if clade.name else "Root"

    # Print results
    print("\nAncestral States (Fitch Parsimony):")
    print(f"{'Node':<15s} {'Possible States':<30s} {'Final State':<15s}")
    print("-" * 60)

    ancestral_results = []
    for clade in tree.find_clades():
        node_name = node_names[clade]
        possible = node_possible_states.get(clade, set())
        final = node_final_states.get(clade, 'Unknown')

        possible_str = ', '.join(sorted(possible)) if possible else 'N/A'

        print(f"{node_name:<15s} {possible_str:<30s} {final:<15s}")

        if not clade.is_terminal():
            ancestral_results.append({
                'node': node_name,
                'possible_states': list(possible),
                'final_state': final,
                'ambiguous': len(possible) > 1
            })

    # Count trait changes
    changes, change_events = count_trait_changes(tree, node_final_states)

    print(f"\n\nTrait Evolution Summary:")
    print(f"  Total trait changes: {changes}")
    print(f"  Parsimony score: {changes}")

    if change_events:
        print(f"\nTrait Change Events:")
        for i, event in enumerate(change_events, 1):
            print(f"  {i}. {event['from']} → {event['to']} (at {event['node']})")

    # Analyze trait distribution
    print(f"\n\nTrait Distribution in Tree:")
    trait_counts = Counter(node_final_states.values())
    for trait, count in sorted(trait_counts.items()):
        print(f"  {trait:<20s}: {count:2d}")

    return {
        'possible_states': {node_names[k]: list(v) for k, v in node_possible_states.items()},
        'final_states': {node_names[k]: v for k, v in node_final_states.items()},
        'ancestral_nodes': ancestral_results,
        'num_changes': changes,
        'change_events': change_events,
        'trait_distribution': dict(trait_counts)
    }

def write_annotated_tree(tree, node_final_states, output_file):
    """Write tree with ancestral state annotations."""

    # Annotate tree
    for clade in tree.find_clades():
        state = node_final_states.get(clade, 'Unknown')
        clade.comment = f"[&trait={state}]"

    # Write to file
    Phylo.write(tree, output_file, 'newick')
    print(f"\nAnnotated tree written to: {output_file}")

if __name__ == "__main__":

    # Type I analysis
    print("\n" + "="*60)
    print("TYPE I DAH7PS TRAIT EVOLUTION ANALYSIS")
    print("="*60)

    type_i_tree = Phylo.read("../trees/type_i/type_i_tree.treefile", "newick")
    type_i_traits = read_traits("type_i/type_i_traits.tsv")

    type_i_results = analyze_trait_evolution(type_i_tree, type_i_traits, "Type I")

    write_annotated_tree(type_i_tree,
                        {k: v for clade in type_i_tree.find_clades()
                         for k, v in [(clade, type_i_results['final_states'].get(clade.name if clade.name else 'Root'))]},
                        "type_i/type_i_tree_annotated.nwk")

    # Save results
    with open("type_i/type_i_trait_reconstruction.json", 'w') as f:
        json.dump(type_i_results, f, indent=2)

    print("\n" + "="*60)
    print("TYPE II DAH7PS TRAIT EVOLUTION ANALYSIS")
    print("="*60)

    type_ii_tree = Phylo.read("../trees/type_ii/type_ii_tree.treefile", "newick")
    type_ii_traits = read_traits("type_ii/type_ii_traits.tsv")

    type_ii_results = analyze_trait_evolution(type_ii_tree, type_ii_traits, "Type II")

    write_annotated_tree(type_ii_tree,
                         {k: v for clade in type_ii_tree.find_clades()
                          for k, v in [(clade, type_ii_results['final_states'].get(clade.name if clade.name else 'Root'))]},
                         "type_ii/type_ii_tree_annotated.nwk")

    # Save results
    with open("type_ii/type_ii_trait_reconstruction.json", 'w') as f:
        json.dump(type_ii_results, f, indent=2)

    print("\n" + "="*60)
    print("TRAIT RECONSTRUCTION COMPLETE")
    print("="*60)
    print("\nOutput files:")
    print("  - type_i/type_i_trait_reconstruction.json")
    print("  - type_i/type_i_tree_annotated.nwk")
    print("  - type_ii/type_ii_trait_reconstruction.json")
    print("  - type_ii/type_ii_tree_annotated.nwk")
