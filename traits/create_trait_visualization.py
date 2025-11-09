#!/usr/bin/env python3
"""
Create ASCII tree visualizations with trait annotations.
"""

from Bio import Phylo
import json
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import io

def create_ascii_tree_with_traits(tree_file, trait_file, output_file):
    """Create ASCII tree representation with trait annotations."""

    # Load tree and traits
    tree = Phylo.read(tree_file, "newick")
    with open(trait_file) as f:
        trait_data = json.load(f)

    # Create output
    output = io.StringIO()

    # Write ASCII tree
    Phylo.draw_ascii(tree, file=output)

    # Get the ASCII representation
    ascii_tree = output.getvalue()

    # Write to file with trait annotations
    with open(output_file, 'w') as f:
        f.write("Tree with Trait Annotations\n")
        f.write("=" * 80 + "\n\n")
        f.write(ascii_tree)
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("Trait Legend:\n")
        f.write("=" * 80 + "\n\n")

        # Write final states
        f.write("Final Ancestral States (Fitch Parsimony):\n")
        f.write("-" * 80 + "\n")
        for node, state in sorted(trait_data['final_states'].items()):
            f.write(f"  {node:<40s} : {state}\n")

        f.write("\n")
        f.write("Trait Changes Detected:\n")
        f.write("-" * 80 + "\n")
        if trait_data['change_events']:
            for i, event in enumerate(trait_data['change_events'], 1):
                f.write(f"  {i}. {event['from']} → {event['to']} (at {event['node']})\n")
        else:
            f.write("  No trait changes detected\n")

        f.write("\n")
        f.write(f"Parsimony Score: {trait_data['num_changes']} change(s)\n")

    print(f"ASCII tree with traits written to: {output_file}")

def create_trait_distribution_plot(trait_data, output_file, title):
    """Create bar plot of trait distribution."""

    traits = list(trait_data['trait_distribution'].keys())
    counts = list(trait_data['trait_distribution'].values())

    # Define colors for each trait
    color_map = {
        'Phe': '#FF6B6B',
        'Tyr': '#4ECDC4',
        'Trp': '#95E1D3',
        'Unknown': '#CCCCCC',
        'Plastid': '#45B7D1'
    }
    colors = [color_map.get(t, '#888888') for t in traits]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(traits, counts, color=colors, edgecolor='black', linewidth=1.5)

    plt.xlabel('Regulatory Specificity', fontsize=12, fontweight='bold')
    plt.ylabel('Number of Sequences/Nodes', fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3, linestyle='--')

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Trait distribution plot saved to: {output_file}")

def create_parsimony_analysis_plot(type_i_data, type_ii_data, output_file):
    """Create comparison plot of parsimony scores and ambiguity."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Parsimony scores
    ax1 = axes[0]
    types = ['Type I', 'Type II']
    scores = [type_i_data['num_changes'], type_ii_data['num_changes']]
    colors = ['#FF6B6B', '#4ECDC4']

    bars1 = ax1.bar(types, scores, color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Number of Trait Changes', fontsize=12, fontweight='bold')
    ax1.set_title('Parsimony Score Comparison', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')

    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontweight='bold', fontsize=12)

    # Plot 2: Ambiguous nodes
    ax2 = axes[1]
    type_i_ambiguous = sum(1 for node in type_i_data['ancestral_nodes'] if node['ambiguous'])
    type_ii_ambiguous = sum(1 for node in type_ii_data['ancestral_nodes'] if node['ambiguous'])

    type_i_total = len(type_i_data['ancestral_nodes'])
    type_ii_total = len(type_ii_data['ancestral_nodes'])

    ambiguous_counts = [type_i_ambiguous, type_ii_ambiguous]
    total_counts = [type_i_total, type_ii_total]

    x = range(len(types))
    width = 0.35

    bars2a = ax2.bar([i - width/2 for i in x], ambiguous_counts, width,
                     label='Ambiguous', color='#FFA07A', edgecolor='black', linewidth=1.5)
    bars2b = ax2.bar([i + width/2 for i in x],
                     [total_counts[i] - ambiguous_counts[i] for i in range(len(types))],
                     width, label='Unambiguous', color='#90EE90', edgecolor='black', linewidth=1.5)

    ax2.set_ylabel('Number of Ancestral Nodes', fontsize=12, fontweight='bold')
    ax2.set_title('Reconstruction Ambiguity', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(types)
    ax2.legend(fontsize=10)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Parsimony analysis plot saved to: {output_file}")

if __name__ == "__main__":

    print("\n" + "="*60)
    print("CREATING TRAIT EVOLUTION VISUALIZATIONS")
    print("="*60 + "\n")

    # Type I visualizations
    print("Type I visualizations:")
    create_ascii_tree_with_traits(
        "../trees/type_i/type_i_tree.treefile",
        "type_i/type_i_trait_reconstruction.json",
        "type_i/type_i_tree_traits.txt"
    )

    with open("type_i/type_i_trait_reconstruction.json") as f:
        type_i_data = json.load(f)

    create_trait_distribution_plot(
        type_i_data,
        "type_i/type_i_trait_distribution.png",
        "Type I DAH7PS Regulatory Specificity Distribution"
    )

    # Type II visualizations
    print("\nType II visualizations:")
    create_ascii_tree_with_traits(
        "../trees/type_ii/type_ii_tree.treefile",
        "type_ii/type_ii_trait_reconstruction.json",
        "type_ii/type_ii_tree_traits.txt"
    )

    with open("type_ii/type_ii_trait_reconstruction.json") as f:
        type_ii_data = json.load(f)

    create_trait_distribution_plot(
        type_ii_data,
        "type_ii/type_ii_trait_distribution.png",
        "Type II DAH7PS Regulatory Specificity Distribution"
    )

    # Comparative visualization
    print("\nComparative analysis:")
    create_parsimony_analysis_plot(
        type_i_data,
        type_ii_data,
        "trait_evolution_comparison.png"
    )

    print("\n" + "="*60)
    print("VISUALIZATION COMPLETE")
    print("="*60)
    print("\nGenerated files:")
    print("  - type_i/type_i_tree_traits.txt")
    print("  - type_i/type_i_trait_distribution.png")
    print("  - type_ii/type_ii_tree_traits.txt")
    print("  - type_ii/type_ii_trait_distribution.png")
    print("  - trait_evolution_comparison.png")
