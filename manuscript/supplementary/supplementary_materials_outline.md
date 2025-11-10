# Supplementary Materials Outline

## Supplementary Tables

### Table S1: Sequence Information and Metadata
**Columns:**
- Protein_ID (UniProt accession)
- Organism
- Gene_name
- Enzyme_type (Type I / Type II)
- Length (amino acids)
- Regulatory_specificity (Phe/Tyr/Trp/Plastid/Unknown)
- Structure_available (PDB ID if applicable)
- Reference (literature citation)

**Rows:** 14 sequences (8 Type I, 6 Type II)

---

### Table S2: Multiple Sequence Alignment Statistics
**Columns:**
- Enzyme_type
- Number_of_sequences
- Original_alignment_length
- Trimmed_alignment_length
- Retention_percentage
- Gap_free_columns (%)
- Mean_conservation (%)
- Mean_pairwise_identity (%)
- Highly_conserved_positions (>90%, count)
- Variable_positions (<30%, count)

**Rows:** 2 (Type I, Type II)

---

### Table S3: Phylogenetic Model Selection and Parameters
**Columns:**
- Enzyme_type
- Best_model (BIC)
- LogLikelihood
- AIC
- BIC
- Tree_length (substitutions/site)
- Gamma_shape_alpha
- Invariable_sites_proportion
- Number_of_iterations
- CPU_time (seconds)
- Wall_clock_time (seconds)

**Rows:** 2 (Type I, Type II)

---

### Table S4: Bootstrap Support Statistics
**Columns:**
- Enzyme_type
- Total_internal_branches
- Mean_UFBoot (%)
- Mean_SH_aLRT (%)
- Branches_UFBoot_≥95 (count, %)
- Branches_UFBoot_≥80 (count, %)
- Branches_SH_aLRT_≥80 (count, %)
- Min_support
- Max_support

**Rows:** 2 (Type I, Type II)

---

### Table S5: Ancestral Sequence Reconstruction Quality
**Columns:**
- Enzyme_type
- Node_ID
- Mean_PP
- Median_PP
- High_confidence_sites (≥0.95, count, %)
- Medium_confidence_sites (0.80-0.95, count, %)
- Low_confidence_sites (<0.80, count, %)
- Ambiguous_positions (count)

**Rows:** 10 (6 Type I nodes, 4 Type II nodes)

---

### Table S6: Trait Evolution Reconstruction Results
**Columns:**
- Enzyme_type
- Node_ID
- Possible_states (Fitch parsimony)
- Final_state
- Ambiguous (Yes/No)
- UFBoot_support (if internal branch)
- SH_aLRT_support

**Rows:** 10 ancestral nodes

---

### Table S7: Parsimony Score Summary
**Columns:**
- Enzyme_type
- Total_trait_states
- Parsimony_score
- Trait_changes_inferred
- Tree_length (subs/site)
- Sequence_changes_per_trait_change

**Rows:** 2 (Type I, Type II)

---

### Table S8: Domain Architecture Analysis
**Columns:**
- Protein_ID
- Total_length
- Catalytic_domain_boundaries
- Regulatory_domain_type (ACT / Transit peptide / None)
- Regulatory_domain_boundaries
- Regulatory_domain_length
- Linker_region (if applicable)
- Additional_domains (e.g., CM for yeast)

**Rows:** 14 sequences

---

### Table S9: Conservation Scores (All Positions)
**Columns:**
- Enzyme_type
- Position (alignment column)
- Consensus_amino_acid
- Conservation_score (0-1)
- Shannon_entropy
- Gap_percentage
- Trait_specific (Yes/No)
- Domain (Catalytic/ACT/Linker/Transit/TIM_barrel)

**Rows:**
- Type I: 346 positions
- Type II: 473 positions

**Format:** TSV file, also provided as Excel workbook

---

### Table S10: Specificity-Determining Residues (Type I)
**Columns:**
- Position (alignment column)
- Phe_sensitive_consensus
- Tyr_sensitive_consensus
- Trp_sensitive_consensus
- Unknown_consensus
- Chi_square_statistic
- P_value
- P_value_corrected (Bonferroni)
- Significance_level
- Domain (Catalytic/ACT)
- Interface_residue (Yes/No)

**Rows:** 231 trait-specific positions

---

### Table S11: Linker Region Detailed Analysis
**Columns:**
- Position
- Consensus_AA
- Conservation_percentage
- Conservation_category
- Hydrophobicity (Kyte-Doolittle)
- Charge
- Secondary_structure_propensity (predicted)
- Functional_prediction

**Rows:** 9 positions (261-269)

---

### Table S12: Protein Stability Predictions
**Columns:**
- Protein_ID
- Full_length_instability_index
- Full_length_stability (Stable/Unstable)
- Catalytic_domain_instability (if applicable)
- Regulatory_domain_instability (if applicable)
- Aliphatic_index
- GRAVY
- Predicted_thermostability

**Rows:** 14 sequences

---

### Table S13: Mutation Stability Effect Predictions
**Columns:**
- Position
- From_amino_acid
- To_amino_acid (trait-specific substitution)
- Physicochemical_change
- Predicted_effect (Neutral/Destabilizing/Strongly destabilizing/Unknown)
- Confidence
- Rationale
- Domain

**Rows:** 133 trait-specific mutations

---

### Table S14: Evolutionary Variability by Position
**Columns:**
- Enzyme_type
- Position
- Diversity_score
- Shannon_entropy
- Normalized_entropy
- Variability_score (0-1)
- Conservation_score
- Correlation_coefficient (variability vs conservation)
- Domain

**Rows:**
- Type I: 346 positions
- Type II: 473 positions

**Format:** TSV file, also Excel

---

### Table S15: Domain-Specific Variability Comparison
**Columns:**
- Enzyme_type
- Domain
- Number_of_positions
- Mean_variability
- SD_variability
- Min_variability
- Max_variability
- Rapidly_evolving_sites (>0.7, count)

**Rows:** 4 (Type I catalytic, Type I ACT, Type II full, Type II barrel)

---

### Table S16: Interface Residue Analysis
**Columns:**
- Position
- Domain_location (Catalytic_C_term / ACT_N_term)
- Consensus_AA
- Conservation_score
- Trait_specific (Yes/No)
- Coevolving_partner (if any)
- Predicted_structural_role

**Rows:** 41 interface residues (20 catalytic C-term, 21 ACT N-term)

---

### Table S17: Coevolving Position Pairs
**Columns:**
- Catalytic_position
- ACT_position
- Spearman_correlation
- P_value
- P_value_corrected
- Significance_level
- Correlation_category (Perfect/Strong positive/Strong negative)
- Interface_pair (Yes/No)
- Predicted_interaction

**Rows:** 17 coevolving pairs

---

### Table S18: Transit Peptide Properties (Type II)
**Columns:**
- Protein_ID
- Transit_peptide_sequence
- Length (aa)
- Ser_Thr_percentage
- Asp_Glu_percentage
- Arg_percentage
- Lys_percentage
- Net_charge
- Predicted_cleavage_site
- Cleavage_motif
- TargetP_score
- ChloroP_score

**Rows:** 3 (*Arabidopsis* DHS1, DHS2, DHS3)

---

## Supplementary Figures

### Figure S1: Complete Multiple Sequence Alignments
**Panel A:** Type I alignment (346 positions, 8 sequences)
- Color-coded by amino acid properties
- Conservation histogram above alignment
- Domain boundaries indicated

**Panel B:** Type II alignment (473 positions, 6 sequences)
- Color-coded by amino acid properties
- Conservation histogram above alignment
- Transit peptide and TIM barrel regions indicated

---

### Figure S2: Phylogenetic Tree with Full Details
**Panel A:** Type I tree with:
- All branch lengths labeled
- UFBoot and SH-aLRT values at all branches
- Outgroup (if used) indicated
- Scale bar with precise units

**Panel B:** Type II tree with:
- All branch lengths labeled
- UFBoot and SH-aLRT values at all branches
- Outgroup (if used) indicated
- Scale bar with precise units

---

### Figure S3: Posterior Probability Heatmaps
**Panel A:** Type I ancestral nodes
- Heatmap: Rows = nodes, Columns = positions
- Color scale: PP from 0 (red) to 1 (green)
- Domain boundaries indicated

**Panel B:** Type II ancestral nodes
- Heatmap: Rows = nodes, Columns = positions
- Color scale: PP from 0 (red) to 1 (green)
- Domain boundaries indicated

---

### Figure S4: Domain Boundary Predictions
**Panel A:** InterProScan domain predictions
- All 14 sequences shown
- Color-coded domains (catalytic, ACT, CM, transit peptide)
- E-value thresholds indicated

**Panel B:** Multiple methods comparison
- InterProScan, Pfam, SMART predictions overlaid
- Consensus boundaries highlighted

---

### Figure S5: Extended Conservation Profiles
**Panel A:** Type I conservation per position
- Line plot with conservation score (0-1)
- Domain regions color-coded
- Highly conserved positions (>0.9) labeled
- Trait-specific positions indicated

**Panel B:** Type II conservation per position
- Line plot with conservation score (0-1)
- Domain regions color-coded
- Active site regions highlighted
- Transit peptide variability shown

---

### Figure S6: Complete Coevolution Matrix
**Panel A:** Full correlation matrix
- Heatmap: Rows = catalytic positions, Columns = ACT positions
- Color scale: Correlation from -1 (blue) to +1 (red)
- Significant pairs (p < 0.05) outlined

**Panel B:** P-value matrix
- Heatmap: Rows = catalytic positions, Columns = ACT positions
- Color scale: P-value from 0 (dark) to 1 (light)
- Bonferroni-corrected threshold indicated

---

### Figure S7: Stability Prediction Distributions
**Panel A:** Instability index distributions
- Histogram for Type I (full-length, catalytic, ACT separately)
- Histogram for Type II (full-length)
- Vertical line at II = 40 (threshold)

**Panel B:** Aliphatic index distributions
- Histogram for both types
- Higher values = more thermostable

**Panel C:** GRAVY distributions
- Histogram for both types
- Negative values expected (hydrophilic proteins)

---

### Figure S8: Trait-Specific Mutation Effects
**Panel A:** Mutation type distribution
- Pie chart of physicochemical changes
- Hydrophobic→charged, polar→hydrophobic, charge reversals, etc.

**Panel B:** Predicted stability effects by domain
- Stacked bar chart: Catalytic vs ACT domains
- Categories: Neutral, Destabilizing, Strongly destabilizing, Unknown

---

### Figure S9: Model Selection Details
**Panel A:** BIC scores for all models tested (Type I)
- Bar chart of top 10 models
- LG+G4 highlighted as best

**Panel B:** BIC scores for all models tested (Type II)
- Bar chart of top 10 models
- Q.PFAM+G4 highlighted as best

---

### Figure S10: Yeast Bifunctional Enzyme Analysis
**Panel A:** ARO3 and ARO4 domain architecture
- Schematic showing DAH7PS + ACT + CM fusion
- Domain boundaries indicated

**Panel B:** Conservation within yeast sequences
- Pairwise comparison of ARO3 vs ARO4
- Identity plot along sequence length
- Variable regions (likely specificity determinants) highlighted

---

## Supplementary Data Files

### Data File S1: All Input Sequences
**Format:** FASTA
**Contents:**
- Type I sequences (8, with headers including metadata)
- Type II sequences (6, with headers including metadata)

---

### Data File S2: Multiple Sequence Alignments
**Format:** FASTA (trimmed alignments)
**Contents:**
- Type I trimmed alignment (346 positions)
- Type II trimmed alignment (473 positions)

---

### Data File S3: Phylogenetic Trees
**Format:** Newick
**Contents:**
- Type I ML tree (.treefile)
- Type I consensus tree with bootstrap (.contree)
- Type II ML tree (.treefile)
- Type II consensus tree with bootstrap (.contree)

---

### Data File S4: Ancestral Sequences
**Format:** FASTA
**Contents:**
- Type I ancestral sequences (6 nodes, 346 aa each)
- Type II ancestral sequences (4 nodes, 473 aa each)
- Headers include node IDs and mean posterior probabilities

---

### Data File S5: Posterior Probability Matrices
**Format:** TSV (tab-separated values)
**Contents:**
- Type I: 6 files (one per node), 346 rows (positions) × 20 columns (amino acids)
- Type II: 4 files (one per node), 473 rows × 20 columns
- Values are posterior probabilities (0-1)

---

### Data File S6: Trait Reconstruction Results
**Format:** JSON
**Contents:**
- Tree topology with node IDs
- Trait assignments for all tips
- Ancestral state reconstructions (all possible states per node)
- Final parsimonious states
- Parsimony score

---

### Data File S7: Conservation and Variability Scores
**Format:** TSV
**Contents:**
- Position-wise conservation scores (Shannon entropy, normalized)
- Position-wise variability scores (diversity, entropy-based)
- Domain annotations
- Trait-specificity annotations

---

### Data File S8: Coevolution Results
**Format:** TSV
**Contents:**
- All tested position pairs (catalytic × ACT)
- Spearman correlation coefficients
- P-values (raw and corrected)
- Significant pairs highlighted

---

### Data File S9: Analysis Scripts
**Format:** Python scripts (.py), documented
**Contents:**
- 01_sequence_alignment.py (MAFFT wrapper)
- 02_alignment_trimming.py (trimAl wrapper)
- 03_phylogenetic_inference.py (IQ-TREE wrapper)
- 04_ancestral_reconstruction.py (IQ-TREE ASR wrapper)
- 05_trait_evolution.py (Fitch parsimony implementation)
- 06_domain_analysis.py (structure-function analysis)
- 07_conservation_analysis.py (conservation and variability)
- 08_selection_analysis.py (site-specific rates)
- 09_stability_prediction.py (BioPython-based)
- 10_coevolution_analysis.py (correlation detection)
- 11_interface_analysis.py (interface enrichment)
- utils.py (shared functions)
- README_scripts.md (documentation)

---

### Data File S10: Complete IQ-TREE Outputs
**Format:** Various (native IQ-TREE formats)
**Contents:**
- .iqtree (full reports)
- .log (analysis logs)
- .mldist (ML pairwise distances)
- .model.gz (model parameters)
- .splits.nex (bootstrap split frequencies)
- For both Type I and Type II analyses

---

## Supplementary Methods

### Extended Methods S1: Detailed Bioinformatics Workflow
**Contents:**
- Step-by-step commands with parameters
- Software installation instructions
- Computational environment specifications
- Quality control checkpoints
- Troubleshooting guide

---

### Extended Methods S2: Statistical Analysis Details
**Contents:**
- Detailed description of all statistical tests
- Assumptions verification
- Multiple testing correction procedures
- Power analysis
- Effect size calculations

---

### Extended Methods S3: Visualization Scripts
**Contents:**
- Python matplotlib code for all figures
- Figure styling and color schemes
- Data preprocessing for visualization
- Publication-quality export settings

---

## Supplementary Notes

### Note S1: Alternative Phylogenetic Models
**Contents:**
- Comparison of top 5 models by BIC
- Sensitivity analysis: Does model choice affect conclusions?
- Tree topology comparison across models

---

### Note S2: Alternative Trait Reconstruction Methods
**Contents:**
- Comparison of Fitch parsimony vs. ML-based methods
- Sensitivity to trait encoding (discrete vs. continuous)
- Effect of unknown trait assignments

---

### Note S3: Coevolution Method Validation
**Contents:**
- Comparison of Spearman vs. Pearson correlation
- Robustness to alignment gaps
- Comparison with mutual information methods
- False positive rate estimation

---

### Note S4: AlphaFold Structure Predictions
**Contents:**
- AlphaFold2 predictions for all 14 sequences (if generated)
- Confidence scores (pLDDT)
- Comparison with experimental structures (yeast ARO3/ARO4)
- Mapping of trait-specific and coevolving residues onto structures

---

## File Organization

```
supplementary_materials/
├── tables/
│   ├── TableS01_sequence_information.xlsx
│   ├── TableS02_alignment_statistics.xlsx
│   ├── TableS03_phylogenetic_parameters.xlsx
│   ├── TableS04_bootstrap_support.xlsx
│   ├── TableS05_asr_quality.xlsx
│   ├── TableS06_trait_reconstruction.xlsx
│   ├── TableS07_parsimony_summary.xlsx
│   ├── TableS08_domain_architecture.xlsx
│   ├── TableS09_conservation_scores.tsv
│   ├── TableS10_specificity_residues.xlsx
│   ├── TableS11_linker_analysis.xlsx
│   ├── TableS12_stability_predictions.xlsx
│   ├── TableS13_mutation_effects.xlsx
│   ├── TableS14_variability_scores.tsv
│   ├── TableS15_domain_variability.xlsx
│   ├── TableS16_interface_residues.xlsx
│   ├── TableS17_coevolution_pairs.xlsx
│   └── TableS18_transit_peptides.xlsx
├── figures/
│   ├── FigureS01_alignments.pdf
│   ├── FigureS02_phylogeny_detailed.pdf
│   ├── FigureS03_pp_heatmaps.pdf
│   ├── FigureS04_domain_predictions.pdf
│   ├── FigureS05_conservation_profiles.pdf
│   ├── FigureS06_coevolution_matrix.pdf
│   ├── FigureS07_stability_distributions.pdf
│   ├── FigureS08_mutation_effects.pdf
│   ├── FigureS09_model_selection.pdf
│   └── FigureS10_yeast_analysis.pdf
├── data/
│   ├── DataS01_sequences.fasta
│   ├── DataS02_alignments.fasta
│   ├── DataS03_trees.nwk
│   ├── DataS04_ancestral_sequences.fasta
│   ├── DataS05_posterior_probabilities/
│   ├── DataS06_trait_reconstruction.json
│   ├── DataS07_scores.tsv
│   ├── DataS08_coevolution.tsv
│   ├── DataS09_scripts/
│   └── DataS10_iqtree_outputs/
├── methods/
│   ├── ExtendedMethodsS1_workflow.pdf
│   ├── ExtendedMethodsS2_statistics.pdf
│   └── ExtendedMethodsS3_visualization.pdf
├── notes/
│   ├── NoteS1_alternative_models.pdf
│   ├── NoteS2_trait_methods.pdf
│   ├── NoteS3_coevolution_validation.pdf
│   └── NoteS4_alphafold_predictions.pdf
└── README_supplementary.md
```

---

**Total Supplementary Tables:** 18
**Total Supplementary Figures:** 10
**Total Supplementary Data Files:** 10
**Total Supplementary Methods:** 3
**Total Supplementary Notes:** 4

**Estimated Total Supplementary Size:** ~500 MB (including scripts, data, and figures)
