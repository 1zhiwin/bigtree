#!/usr/bin/env python3
"""
Create visualizations for structure-function analysis.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import json

def create_domain_architecture_plot():
    """
    Create schematic domain architecture plots for Type I and Type II.
    """

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    # Type I architecture
    ax1.set_xlim(0, 400)
    ax1.set_ylim(0, 10)
    ax1.axis('off')
    ax1.set_title('Type I DAH7PS Domain Architecture', fontsize=16, fontweight='bold', pad=20)

    # E. coli aroF (Phe-sensitive)
    y_pos = 8
    ax1.text(-20, y_pos, 'E. coli aroF\n(Phe)', ha='right', va='center', fontsize=10)
    ax1.add_patch(patches.Rectangle((0, y_pos-0.4), 260, 0.8, facecolor='#4ECDC4', edgecolor='black', linewidth=1.5))
    ax1.add_patch(patches.Rectangle((270, y_pos-0.4), 86, 0.8, facecolor='#FF6B6B', edgecolor='black', linewidth=1.5))
    ax1.text(130, y_pos, 'Catalytic Domain', ha='center', va='center', fontsize=9, fontweight='bold')
    ax1.text(313, y_pos, 'ACT', ha='center', va='center', fontsize=9, fontweight='bold')

    # E. coli aroG (Tyr-sensitive)
    y_pos = 6
    ax1.text(-20, y_pos, 'E. coli aroG\n(Tyr)', ha='right', va='center', fontsize=10)
    ax1.add_patch(patches.Rectangle((0, y_pos-0.4), 260, 0.8, facecolor='#4ECDC4', edgecolor='black', linewidth=1.5))
    ax1.add_patch(patches.Rectangle((270, y_pos-0.4), 80, 0.8, facecolor='#FFA07A', edgecolor='black', linewidth=1.5))
    ax1.text(313, y_pos, 'ACT', ha='center', va='center', fontsize=9, fontweight='bold')

    # E. coli aroH (Trp-sensitive)
    y_pos = 4
    ax1.text(-20, y_pos, 'E. coli aroH\n(Trp)', ha='right', va='center', fontsize=10)
    ax1.add_patch(patches.Rectangle((0, y_pos-0.4), 260, 0.8, facecolor='#4ECDC4', edgecolor='black', linewidth=1.5))
    ax1.add_patch(patches.Rectangle((270, y_pos-0.4), 78, 0.8, facecolor='#95E1D3', edgecolor='black', linewidth=1.5))
    ax1.text(309, y_pos, 'ACT', ha='center', va='center', fontsize=9, fontweight='bold')

    # Yeast ARO3 (bifunctional)
    y_pos = 2
    ax1.text(-20, y_pos, 'Yeast ARO3\n(Phe, bifunc)', ha='right', va='center', fontsize=10)
    ax1.add_patch(patches.Rectangle((0, y_pos-0.4), 260, 0.8, facecolor='#4ECDC4', edgecolor='black', linewidth=1.5))
    ax1.add_patch(patches.Rectangle((270, y_pos-0.4), 30, 0.8, facecolor='#FF6B6B', edgecolor='black', linewidth=1.5))
    ax1.add_patch(patches.Rectangle((300, y_pos-0.4), 70, 0.8, facecolor='#FFD93D', edgecolor='black', linewidth=1.5))
    ax1.text(285, y_pos, 'ACT', ha='center', va='center', fontsize=8, fontweight='bold')
    ax1.text(335, y_pos, 'CM Domain', ha='center', va='center', fontsize=9, fontweight='bold')

    # Legend for Type I
    ax1.add_patch(patches.Rectangle((0, 0.2), 30, 0.4, facecolor='#4ECDC4', edgecolor='black'))
    ax1.text(35, 0.4, 'Catalytic Domain (α/β fold)', va='center', fontsize=9)

    ax1.add_patch(patches.Rectangle((200, 0.2), 15, 0.4, facecolor='#FF6B6B', edgecolor='black'))
    ax1.text(220, 0.4, 'ACT Domain', va='center', fontsize=9)

    ax1.add_patch(patches.Rectangle((300, 0.2), 15, 0.4, facecolor='#FFD93D', edgecolor='black'))
    ax1.text(320, 0.4, 'CM Domain', va='center', fontsize=9)

    # Scale bar
    ax1.plot([0, 100], [0.8, 0.8], 'k-', linewidth=2)
    ax1.text(50, 0.5, '100 aa', ha='center', va='top', fontsize=9)

    # Type II architecture
    ax2.set_xlim(0, 600)
    ax2.set_ylim(0, 8)
    ax2.axis('off')
    ax2.set_title('Type II DAH7PS Domain Architecture', fontsize=16, fontweight='bold', pad=20)

    # M. tuberculosis (bacterial, no transit peptide)
    y_pos = 6.5
    ax2.text(-20, y_pos, 'M. tuberculosis\n(bacterial)', ha='right', va='center', fontsize=10)
    ax2.add_patch(patches.Rectangle((0, y_pos-0.4), 462, 0.8, facecolor='#45B7D1', edgecolor='black', linewidth=1.5))
    ax2.text(231, y_pos, 'TIM Barrel Domain', ha='center', va='center', fontsize=10, fontweight='bold')

    # Arabidopsis DHS1 (plastid-targeted)
    y_pos = 5
    ax2.text(-20, y_pos, 'A. thaliana DHS1\n(plastid)', ha='right', va='center', fontsize=10)
    ax2.add_patch(patches.Rectangle((0, y_pos-0.4), 47, 0.8, facecolor='#90EE90', edgecolor='black', linewidth=1.5))
    ax2.add_patch(patches.Rectangle((47, y_pos-0.4), 478, 0.8, facecolor='#45B7D1', edgecolor='black', linewidth=1.5))
    ax2.text(23, y_pos, 'TP', ha='center', va='center', fontsize=8, fontweight='bold')
    ax2.text(286, y_pos, 'TIM Barrel Domain', ha='center', va='center', fontsize=10, fontweight='bold')

    # Arabidopsis DHS2
    y_pos = 3.5
    ax2.text(-20, y_pos, 'A. thaliana DHS2\n(plastid)', ha='right', va='center', fontsize=10)
    ax2.add_patch(patches.Rectangle((0, y_pos-0.4), 32, 0.8, facecolor='#90EE90', edgecolor='black', linewidth=1.5))
    ax2.add_patch(patches.Rectangle((32, y_pos-0.4), 475, 0.8, facecolor='#45B7D1', edgecolor='black', linewidth=1.5))
    ax2.text(16, y_pos, 'TP', ha='center', va='center', fontsize=8, fontweight='bold')
    ax2.text(269, y_pos, 'TIM Barrel Domain', ha='center', va='center', fontsize=10, fontweight='bold')

    # Arabidopsis DHS3
    y_pos = 2
    ax2.text(-20, y_pos, 'A. thaliana DHS3\n(plastid)', ha='right', va='center', fontsize=10)
    ax2.add_patch(patches.Rectangle((0, y_pos-0.4), 57, 0.8, facecolor='#90EE90', edgecolor='black', linewidth=1.5))
    ax2.add_patch(patches.Rectangle((57, y_pos-0.4), 470, 0.8, facecolor='#45B7D1', edgecolor='black', linewidth=1.5))
    ax2.text(28, y_pos, 'TP', ha='center', va='center', fontsize=8, fontweight='bold')
    ax2.text(292, y_pos, 'TIM Barrel Domain', ha='center', va='center', fontsize=10, fontweight='bold')

    # Legend for Type II
    ax2.add_patch(patches.Rectangle((0, 0.2), 30, 0.4, facecolor='#90EE90', edgecolor='black'))
    ax2.text(35, 0.4, 'Transit Peptide (TP)', va='center', fontsize=9)

    ax2.add_patch(patches.Rectangle((170, 0.2), 30, 0.4, facecolor='#45B7D1', edgecolor='black'))
    ax2.text(205, 0.4, 'TIM Barrel (β/α)₈', va='center', fontsize=9)

    # Scale bar
    ax2.plot([400, 500], [0.8, 0.8], 'k-', linewidth=2)
    ax2.text(450, 0.5, '100 aa', ha='center', va='top', fontsize=9)

    plt.tight_layout()
    plt.savefig('domain_architecture.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Domain architecture plot saved to: structure/domain_architecture.png")

def create_conservation_plots():
    """
    Create conservation profile plots for Type I and Type II.
    """

    # Load conservation data
    with open('conservation/type_i_conservation.json') as f:
        type_i_data = json.load(f)

    with open('conservation/type_ii_conservation.json') as f:
        type_ii_data = json.load(f)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

    # Type I conservation
    scores_i = type_i_data['conservation_scores']
    positions_i = list(range(1, len(scores_i) + 1))

    ax1.plot(positions_i, scores_i, 'b-', alpha=0.5, linewidth=0.5)
    ax1.fill_between(positions_i, scores_i, alpha=0.3)
    ax1.axhline(y=0.9, color='r', linestyle='--', linewidth=1, label='Highly conserved (>0.9)')
    ax1.axhline(y=0.3, color='orange', linestyle='--', linewidth=1, label='Variable (<0.3)')
    ax1.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Conservation Score', fontsize=12, fontweight='bold')
    ax1.set_title('Type I DAH7PS Conservation Profile', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 1)
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.legend(fontsize=10)

    # Mark ACT domain region (approximately positions 270-350 in alignment)
    ax1.axvspan(270, 346, alpha=0.1, color='red', label='ACT domain region')
    ax1.text(308, 0.95, 'ACT Domain', ha='center', fontsize=10, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Type II conservation
    scores_ii = type_ii_data['conservation_scores']
    positions_ii = list(range(1, len(scores_ii) + 1))

    ax2.plot(positions_ii, scores_ii, 'g-', alpha=0.5, linewidth=0.5)
    ax2.fill_between(positions_ii, scores_ii, alpha=0.3, color='green')
    ax2.axhline(y=0.9, color='r', linestyle='--', linewidth=1, label='Highly conserved (>0.9)')
    ax2.axhline(y=0.3, color='orange', linestyle='--', linewidth=1, label='Variable (<0.3)')
    ax2.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Conservation Score', fontsize=12, fontweight='bold')
    ax2.set_title('Type II DAH7PS Conservation Profile', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 1)
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.legend(fontsize=10)

    # Mark TIM barrel region (core catalytic domain)
    ax2.axvspan(50, 450, alpha=0.1, color='blue')
    ax2.text(250, 0.95, 'TIM Barrel Core', ha='center', fontsize=10, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.tight_layout()
    plt.savefig('conservation_profiles.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Conservation profile plots saved to: structure/conservation_profiles.png")

def create_transit_peptide_comparison():
    """
    Create comparison plot of transit peptide properties.
    """

    # Load transit peptide data
    df = pd.read_csv('type_ii/transit_peptide_analysis.tsv', sep='\t')

    # Filter only sequences with transit peptides
    df_with_tp = df[df['predicted_transit_length'] > 0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: Transit peptide lengths
    ax1 = axes[0]
    organisms = [org.split()[-1] if 'Arabidopsis' in org else org.split()[-2]
                 for org in df_with_tp['organism']]
    lengths = df_with_tp['predicted_transit_length']
    colors = ['#90EE90' if 'Arabidopsis' in org else '#FFA07A'
              for org in df_with_tp['organism']]

    bars1 = ax1.bar(range(len(organisms)), lengths, color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Transit Peptide Length (aa)', fontsize=11, fontweight='bold')
    ax1.set_title('Transit Peptide Lengths', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(organisms)))
    ax1.set_xticklabels(organisms, rotation=45, ha='right')
    ax1.grid(axis='y', alpha=0.3)

    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')

    # Plot 2: Ser+Thr content
    ax2 = axes[1]
    ser_thr = df_with_tp['transit_ser_thr']
    bars2 = ax2.bar(range(len(organisms)), ser_thr, color=colors, edgecolor='black', linewidth=1.5)
    ax2.set_ylabel('Ser + Thr Content (%)', fontsize=11, fontweight='bold')
    ax2.set_title('N-terminal Ser+Thr Enrichment', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(organisms)))
    ax2.set_xticklabels(organisms, rotation=45, ha='right')
    ax2.grid(axis='y', alpha=0.3)

    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)

    # Plot 3: Arg content
    ax3 = axes[2]
    arg_content = df_with_tp['transit_arg']
    bars3 = ax3.bar(range(len(organisms)), arg_content, color=colors, edgecolor='black', linewidth=1.5)
    ax3.set_ylabel('Arg Content (%)', fontsize=11, fontweight='bold')
    ax3.set_title('N-terminal Arg Content', fontsize=12, fontweight='bold')
    ax3.set_xticks(range(len(organisms)))
    ax3.set_xticklabels(organisms, rotation=45, ha='right')
    ax3.grid(axis='y', alpha=0.3)

    for i, bar in enumerate(bars3):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)

    plt.tight_layout()
    plt.savefig('transit_peptide_properties.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Transit peptide comparison plot saved to: structure/transit_peptide_properties.png")

if __name__ == "__main__":

    print("\n" + "="*60)
    print("CREATING STRUCTURE-FUNCTION VISUALIZATIONS")
    print("="*60 + "\n")

    create_domain_architecture_plot()
    create_conservation_plots()
    create_transit_peptide_comparison()

    print("\n" + "="*60)
    print("VISUALIZATION COMPLETE")
    print("="*60)
    print("\nGenerated files:")
    print("  - structure/domain_architecture.png")
    print("  - structure/conservation_profiles.png")
    print("  - structure/transit_peptide_properties.png")
