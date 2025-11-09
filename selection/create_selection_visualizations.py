#!/usr/bin/env python3
"""
Create visualizations for selection analysis results.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import json

def create_evolutionary_rate_profile_plot():
    """Create evolutionary rate profiles along sequences."""

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

    # Type I
    type_i_rates = pd.read_csv('type_i/site_rates.tsv', sep='\t')
    positions_i = type_i_rates['position']
    variability_i = type_i_rates['variability']

    ax1.plot(positions_i, variability_i, 'b-', alpha=0.5, linewidth=0.5)
    ax1.fill_between(positions_i, variability_i, alpha=0.3)
    ax1.axhline(y=0.7, color='r', linestyle='--', linewidth=1, label='Rapidly evolving (>0.7)')
    ax1.axhline(y=0.2, color='g', linestyle='--', linewidth=1, label='Conserved (<0.2)')

    # Mark domain boundaries
    ax1.axvspan(1, 260, alpha=0.1, color='blue')
    ax1.axvspan(270, 346, alpha=0.1, color='red')

    ax1.text(130, 0.95, 'Catalytic Domain', ha='center', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    ax1.text(308, 0.95, 'ACT Domain', ha='center', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

    ax1.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Evolutionary Variability', fontsize=12, fontweight='bold')
    ax1.set_title('Type I DAH7PS Evolutionary Rate Profile', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 1)
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.legend(fontsize=10)

    # Type II
    type_ii_rates = pd.read_csv('type_ii/site_rates.tsv', sep='\t')
    positions_ii = type_ii_rates['position']
    variability_ii = type_ii_rates['variability']

    ax2.plot(positions_ii, variability_ii, 'g-', alpha=0.5, linewidth=0.5)
    ax2.fill_between(positions_ii, variability_ii, alpha=0.3, color='green')
    ax2.axhline(y=0.7, color='r', linestyle='--', linewidth=1, label='Rapidly evolving (>0.7)')
    ax2.axhline(y=0.2, color='g', linestyle='--', linewidth=1, label='Conserved (<0.2)')

    # Mark TIM barrel
    ax2.axvspan(1, 473, alpha=0.1, color='blue')
    ax2.text(237, 0.95, 'TIM Barrel Domain', ha='center', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

    ax2.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Evolutionary Variability', fontsize=12, fontweight='bold')
    ax2.set_title('Type II DAH7PS Evolutionary Rate Profile', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 1)
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig('visualization/evolutionary_rate_profiles.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Evolutionary rate profile plot saved to: selection/visualization/evolutionary_rate_profiles.png")

def create_domain_comparison_plot():
    """Compare evolutionary rates across domains."""

    # Load domain stats
    with open('type_i/domain_evolution_stats.json') as f:
        type_i_stats = json.load(f)

    with open('type_ii/domain_evolution_stats.json') as f:
        type_ii_stats = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Mean variability by domain (Type I)
    ax1 = axes[0]
    domains_i = list(type_i_stats.keys())
    variability_i = [type_i_stats[d]['mean_variability'] for d in domains_i]
    colors_i = ['#4ECDC4', '#FFD93D', '#FF6B6B']

    bars1 = ax1.bar(range(len(domains_i)), variability_i, color=colors_i,
                    edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Mean Variability Score', fontsize=11, fontweight='bold')
    ax1.set_title('Type I: Domain-Specific Evolutionary Rates', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(domains_i)))
    ax1.set_xticklabels(['Catalytic\nDomain', 'Linker', 'ACT\nDomain'])
    ax1.set_ylim(0, 0.5)
    ax1.grid(axis='y', alpha=0.3)

    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.3f}', ha='center', va='bottom', fontweight='bold')

    # Plot 2: Conserved vs Variable sites (Type I)
    ax2 = axes[1]
    domains = list(type_i_stats.keys())
    conserved_counts = [type_i_stats[d]['conserved_sites'] for d in domains]
    variable_counts = [type_i_stats[d]['variable_sites'] for d in domains]

    x = np.arange(len(domains))
    width = 0.35

    bars2a = ax2.bar(x - width/2, conserved_counts, width, label='Conserved (<0.2)',
                     color='#90EE90', edgecolor='black', linewidth=1.5)
    bars2b = ax2.bar(x + width/2, variable_counts, width, label='Variable (>0.5)',
                     color='#FFA07A', edgecolor='black', linewidth=1.5)

    ax2.set_ylabel('Number of Sites', fontsize=11, fontweight='bold')
    ax2.set_title('Type I: Conserved vs Variable Sites by Domain', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Catalytic\nDomain', 'Linker', 'ACT\nDomain'])
    ax2.legend(fontsize=10)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('visualization/domain_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Domain comparison plot saved to: selection/visualization/domain_comparison.png")

def create_conservation_vs_variability_scatter():
    """Create scatter plot of conservation vs. variability."""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Type I
    type_i_rates = pd.read_csv('type_i/site_rates.tsv', sep='\t')
    with open('../structure/conservation/type_i_conservation.json') as f:
        type_i_conservation = json.load(f)

    conserv_scores_i = type_i_conservation['conservation_scores'][:len(type_i_rates)]
    variab_scores_i = type_i_rates['variability'].values

    # Filter out gaps
    mask_i = type_i_rates['rate_category'] != 'gap'
    conserv_i = [conserv_scores_i[i] for i in range(len(conserv_scores_i)) if mask_i.iloc[i]]
    variab_i = [variab_scores_i[i] for i in range(len(variab_scores_i)) if mask_i.iloc[i]]

    ax1.scatter(conserv_i, variab_i, alpha=0.5, s=20, color='blue')
    ax1.set_xlabel('Conservation Score (Phase 8)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Evolutionary Variability (Phase 9)', fontsize=11, fontweight='bold')
    ax1.set_title('Type I: Conservation vs. Variability', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)

    # Add trend line
    z_i = np.polyfit(conserv_i, variab_i, 1)
    p_i = np.poly1d(z_i)
    x_line_i = np.linspace(min(conserv_i), max(conserv_i), 100)
    ax1.plot(x_line_i, p_i(x_line_i), "r--", linewidth=2, label=f'Linear fit (r={np.corrcoef(conserv_i, variab_i)[0,1]:.3f})')
    ax1.legend()

    # Type II
    type_ii_rates = pd.read_csv('type_ii/site_rates.tsv', sep='\t')
    with open('../structure/conservation/type_ii_conservation.json') as f:
        type_ii_conservation = json.load(f)

    conserv_scores_ii = type_ii_conservation['conservation_scores'][:len(type_ii_rates)]
    variab_scores_ii = type_ii_rates['variability'].values

    mask_ii = type_ii_rates['rate_category'] != 'gap'
    conserv_ii = [conserv_scores_ii[i] for i in range(len(conserv_scores_ii)) if mask_ii.iloc[i]]
    variab_ii = [variab_scores_ii[i] for i in range(len(variab_scores_ii)) if mask_ii.iloc[i]]

    ax2.scatter(conserv_ii, variab_ii, alpha=0.5, s=20, color='green')
    ax2.set_xlabel('Conservation Score (Phase 8)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Evolutionary Variability (Phase 9)', fontsize=11, fontweight='bold')
    ax2.set_title('Type II: Conservation vs. Variability', fontsize=12, fontweight='bold')
    ax2.grid(alpha=0.3)

    # Add trend line
    z_ii = np.polyfit(conserv_ii, variab_ii, 1)
    p_ii = np.poly1d(z_ii)
    x_line_ii = np.linspace(min(conserv_ii), max(conserv_ii), 100)
    ax2.plot(x_line_ii, p_ii(x_line_ii), "r--", linewidth=2, label=f'Linear fit (r={np.corrcoef(conserv_ii, variab_ii)[0,1]:.3f})')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('visualization/conservation_variability_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Conservation-variability scatter plot saved to: selection/visualization/conservation_variability_correlation.png")

def create_rate_category_distribution():
    """Create pie charts showing distribution of rate categories."""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Type I
    type_i_rates = pd.read_csv('type_i/site_rates.tsv', sep='\t')
    type_i_categories = type_i_rates['rate_category'].value_counts()

    colors = ['#90EE90', '#FFD93D', '#FFA07A', '#CCCCCC']
    # Create explode tuple matching number of categories
    explode = tuple([0.05 if i == 0 else 0 for i in range(len(type_i_categories))])

    ax1.pie(type_i_categories.values, labels=type_i_categories.index, autopct='%1.1f%%',
            colors=colors[:len(type_i_categories)], explode=explode, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax1.set_title('Type I: Evolutionary Rate Categories', fontsize=13, fontweight='bold')

    # Type II
    type_ii_rates = pd.read_csv('type_ii/site_rates.tsv', sep='\t')
    type_ii_categories = type_ii_rates['rate_category'].value_counts()
    explode_ii = tuple([0.05 if i == 0 else 0 for i in range(len(type_ii_categories))])

    ax2.pie(type_ii_categories.values, labels=type_ii_categories.index, autopct='%1.1f%%',
            colors=colors[:len(type_ii_categories)], explode=explode_ii, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax2.set_title('Type II: Evolutionary Rate Categories', fontsize=13, fontweight='bold')

    plt.tight_layout()
    plt.savefig('visualization/rate_category_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Rate category distribution plot saved to: selection/visualization/rate_category_distribution.png")

if __name__ == "__main__":

    print("\n" + "="*60)
    print("CREATING SELECTION ANALYSIS VISUALIZATIONS")
    print("="*60 + "\n")

    create_evolutionary_rate_profile_plot()
    create_domain_comparison_plot()
    create_conservation_vs_variability_scatter()
    create_rate_category_distribution()

    print("\n" + "="*60)
    print("VISUALIZATION COMPLETE")
    print("="*60)
    print("\nGenerated files:")
    print("  - selection/visualization/evolutionary_rate_profiles.png")
    print("  - selection/visualization/domain_comparison.png")
    print("  - selection/visualization/conservation_variability_correlation.png")
    print("  - selection/visualization/rate_category_distribution.png")
