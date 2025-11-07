# DAH7PS Allostery Evolution & Ancestral Sequence Reconstruction

## Project Overview

**Goal:** Conduct comprehensive phylogenetic analysis and ancestral sequence reconstruction across the DAH7PS (3-deoxy-D-arabinoheptulosonate 7-phosphate synthase) family to understand the evolution of allosteric regulation mechanisms.

**Principal Investigator:** [To be filled]
**Start Date:** 2025-11-07
**Project Repository:** bigtree

---

## 1. Scientific Background

### 1.1 DAH7PS Enzyme System

DAH7PS catalyzes the first committed step of the shikimate pathway, condensing phosphoenolpyruvate (PEP) and erythrose 4-phosphate (E4P) to form DAH7P. This pathway is essential in bacteria, archaea, plants, fungi, and apicomplexan parasites for biosynthesis of aromatic amino acids.

### 1.2 Enzyme Classes

- **Type Iα**: Classical bacterial DAH7PS with α-barrel core
- **Type Iβ**: Related to Iα but with distinct regulatory features
- **Type II**: Alternative β-barrel architecture, primarily archaeal/some bacteria

### 1.3 Regulatory Mechanisms

**Allosteric Regulation:**
- **ACT domain-mediated**: Small molecule binding domains (often tandem repeats)
  - Effectors: Tyr, Phe, Trp
  - Can show single or synergistic inhibition
- **Chorismate mutase (CM) fusion**: Bifunctional enzyme with metabolic channeling
  - Effectors: Chorismate, prephenate, aromatic amino acids
- **No regulation**: Some organisms have unregulated DAH7PS

**Effector Specificity:**
- Tyr-sensitive (aroG-like)
- Phe-sensitive (aroF-like)
- Trp-sensitive (aroH-like)
- Synergistic (multiple effectors required)

---

## 2. Project Scope

### 2.1 Taxonomic Sampling

**Primary Focus:**
- **Bacteria**: Comprehensive sampling across major phyla (Proteobacteria, Firmicutes, Actinobacteria, etc.)
- **Archaea**: Representative sampling (Euryarchaeota, Crenarchaeota, etc.)
- **Plantae**: Emphasis on plastid-targeted enzymes
- **Fungi**: Ascomycota and Basidiomycota representatives

**Exclusions:**
- Eukaryotic cytoplasmic DAH7PS (unless specifically noted)
- Sequences <85% domain coverage
- Known KDO8PS contamination

### 2.2 Molecular Features of Interest

**Domain Architecture:**
- Catalytic core presence and boundaries
- ACT domain count (0, 1, 2+)
- CM domain fusion (N- or C-terminal)
- Transit peptides and localization signals

**Regulatory Traits:**
- Effector identity (Tyr/Phe/Trp/Chorismate/Prephenate)
- Inhibition mode (single vs synergistic)
- Kinetic parameters (when available)
- Literature-documented phenotypes

**Structural Features:**
- Oligomerization state
- Interface residues (core-ACT, core-CM)
- Catalytic site conservation
- Allosteric pocket composition

### 2.3 Evolutionary Questions

1. **Domain Recruitment:** When and how many times were ACT and CM domains recruited to DAH7PS cores?
2. **Effector Switching:** What substitutions underlie changes in effector specificity?
3. **Synergy Evolution:** How did synergistic inhibition evolve from single-effector ancestors?
4. **Functional Convergence:** Do independent lineages converge on similar regulatory solutions?
5. **Ancestral States:** What were the regulatory properties of ancient DAH7PS enzymes?

---

## 3. Methodology

### 3.1 Sequence Collection (Phase 1)

**Data Sources:**
- UniProt Reference Proteomes + UniRef90
- NCBI RefSeq (to fill taxonomic gaps)
- EnsemblPlants/Phytozome (plant plastid DHS)

**Search Strategy:**
- HMMER3 searches using Pfam profiles:
  - PF00793 (DAHP_synth_1; Type I)
  - PF01474 (DAHP_synth_2; Type II)
- Trusted cutoff (`--cut_ga`)
- E-value threshold: 1e-10 (backup)

**Quality Control:**
- Domain coverage ≥85-90%
- No premature stop codons
- ≤20% ambiguous residues ('X')
- Length filters based on expected domain architecture

**Redundancy Reduction:**
- Cluster at 90% identity (CD-HIT or MMseqs2)
- Keep longest representative per cluster
- Maintain cluster membership mapping

**Contamination Removal:**
- Filter KDO8PS (kdsA) via KEGG KO annotation
- Manual inspection of borderline cases
- QC gate: ≤1% suspected contamination

**Expected Output:**
- 50-200 representatives per major class
- Broad taxonomic distribution
- `seqs/dah7ps_candidates.faa`
- `seqs/dah7ps_nonredundant.faa`
- `data/processed/cluster_map.tsv`

### 3.2 Domain Annotation (Phase 2)

**Domain Detection:**
- HMMER hmmscan against Pfam-A
- Regulatory domains:
  - PF01842 (ACT domain) - record all copies and coordinates
  - AroQ/CM_2 (Chorismate mutase) - record coordinates
- Parse domain architecture per sequence

**Localization Prediction:**
- SignalP 6.0 for signal peptides
- TargetP 2.0 for chloroplast/mitochondrial transit peptides
- **Mask** targeting sequences from alignments

**Architecture Classification:**
- `core_only`
- `core+ACT` (specify count)
- `core+CM`
- `core+ACT+CM`
- `other/unknown`

**Phenotype Curation:**
- Literature mining for characterized enzymes
- Extract effector specificity, Ki values, synergy
- UniProt functional annotation review
- Create `traits/phenotypes.tsv` with citations

**Metadata Compilation:**
- Taxonomy (NCBI taxid, lineage)
- Sequence length
- DAH7PS class (Iα/Iβ/II)
- Domain architecture
- KEGG KO assignment
- Predicted localization
- Curated phenotype
- Data source
- Output: `data/processed/metadata.tsv`

### 3.3 Multiple Sequence Alignment (Phase 3)

**Alignment Strategy:**

*Partition-based approach to handle domain modularity:*

**Core Domain Alignments (separate by class):**
- Type Iα catalytic core
- Type Iβ catalytic core
- Type II catalytic core
- Masked: transit peptides, low-complexity regions

**Regulatory Domain Alignments:**
- ACT domain (all sequences possessing ACT)
- CM domain (all CM-fusion sequences)

**Alignment Tools:**
- Primary: MAFFT L-INS-i or E-INS-i (iterative refinement)
- Secondary: MUSCLE5 (for comparison)
- Ensemble alignments (MUSCLE5 `-perturb`) for uncertainty

**Quality Assessment:**
- GUIDANCE2 (optional): column and sequence confidence scores
- Visual inspection of catalytic motifs

**Trimming:**
- trimAl: automated1 or gappyout mode
- BMGE: entropy and gap-based trimming
- **Constraint**: Retain known catalytic residues and regulatory interface sites
- Document trimming parameters

**Output:**
- `msa/core_Ialpha.trim.faa`
- `msa/core_Ibeta.trim.faa`
- `msa/core_II.trim.faa`
- `msa/ACT.trim.faa`
- `msa/CM.trim.faa`
- Optional: `msa/concat.faa` + `msa/partitions.nex` (if concatenation used)

### 3.4 Phylogenetic Inference (Phase 4)

**Tree Building:**

*Per-partition trees:*
- IQ-TREE 2 for each core alignment (Iα, Iβ, II)
- Optional: concatenated core+ACT+CM tree

**Model Selection:**
- Automatic with ModelFinder (`-m MFP`)
- Allow complex models: LG+R, LG+C60+F+R, LG+C20+F+G (mixture models)
- Profile mixture models (PMSF) for large datasets

**Branch Support:**
- SH-aLRT: 1000 replicates (fast, good for long branches)
- UFBoot2: 1000 replicates (ultrafast bootstrap)
- **Acceptance threshold**: UFBoot ≥95 AND SH-aLRT ≥80 for key splits

**Rooting Strategy:**
- Within-class trees: midpoint rooting
- Global tree: outgroup with curated KDO8PS (if applicable)
- Check for long-branch attraction artifacts with mixture models

**Output:**
- `trees/core_Ialpha.treefile` (Newick)
- `trees/core_Ibeta.treefile`
- `trees/core_II.treefile`
- `trees/concat.treefile` (if used)
- `trees/support_summaries.tsv`
- iTOL annotation files for visualization

**QC Gate:**
- Majority of key splits with UFBoot ≥95, SH-aLRT ≥80

### 3.5 Ancestral Sequence Reconstruction (Phase 5)

**Node Selection:**
- Last common ancestors (LCAs) of each major class (Iα, Iβ, II)
- Duplication nodes separating aroF/aroG/aroH-like clades
- Nodes at domain gain/loss events:
  - Earliest ACT-bearing ancestor
  - Earliest CM-fusion ancestor
  - Nodes subtending clades with regulatory transitions

**Reconstruction Method:**

*Per-partition ASR:*
- IQ-TREE `--asr` or `--ancestral` (marginal reconstruction under best model)
- Core, ACT, CM domains reconstructed separately
- Optional: PAML codeml for joint reconstruction comparison
- Optional: FastML for indel reconstruction

**Uncertainty Quantification:**
- Site-wise posterior probabilities from IQ-TREE
- Consensus sequences at PP ≥0.8 (adjustable)
- Identify ambiguous positions (PP <0.8) for experimental soft randomization
- Tree uncertainty: ASR over bootstrap tree samples
- Model uncertainty: ASR under alternative models

**Output:**
- `asr/<node_id>_core.faa`
- `asr/<node_id>_ACT.faa`
- `asr/<node_id>_CM.faa`
- `asr/<node_id>_posteriors.tsv` (per-site PP)
- `asr/ambiguous_sites_list.tsv`
- `asr/asr_summary.md`

**QC Gate:**
- Median site PP ≥0.8 in catalytic core for focal ancestors

### 3.6 Trait Evolution Analysis (Phase 6)

**Trait Matrix Construction:**

*Discrete traits:*
- ACT domain presence (binary)
- CM domain presence (binary)
- Effector class (categorical: Tyr, Phe, Trp, Chorismate, Prephenate, Multiple, None)
- Inhibition mode (categorical: Single, Synergistic, None)

**Ancestral Trait Reconstruction:**
- Maximum likelihood: Mk model (PastML or phytools in R)
- Bayesian (optional): BayesTraits or RevBayes for uncertainty
- Stochastic mapping for transition rate estimation

**Analyses:**
- Count domain recruitment and loss events
- Timing of effector switching
- Correlation: domain architecture ↔ effector specificity
- Lineage-specific patterns

**Output:**
- `traits/trait_matrix.tsv`
- `traits/recon.tsv` (ancestral states)
- `figs/traits_evolution.pdf`
- `figs/domain_timeline.pdf`

### 3.7 Structural Mapping (Phase 7)

**Structure Collection:**
- Representative PDBs for Type Iα, Iβ, II
- ACT domain exemplars bound to effectors
- CM domain structures

**Substitution Mapping:**
- Map ancestral and lineage-defining substitutions onto structures
- Focus regions:
  - Core catalytic site
  - Core ↔ ACT interface
  - Core ↔ CM interface
  - Oligomerization interfaces
  - Allosteric pockets

**Analyses:**
- Highlight positions with ΔPP (ASR uncertainty hotspots)
- Identify convergent substitutions
- Rate4Site: site-wise evolutionary rates
- Interface vs non-interface substitution enrichment

**Visualization:**
- PyMOL/ChimeraX sessions
- Structure figures for manuscript

**Output:**
- `structures/mapping/substitution_maps.tsv`
- `structures/mapping/*.pse` (PyMOL sessions)
- `figs/structure_*.png`

### 3.8 Robustness & Sensitivity (Phase 8)

**Alignment Robustness:**
- Compare MAFFT vs MUSCLE5 alignments
- Compare trimAl vs BMGE trimming
- Re-run phylogeny and ASR on alternate alignments
- Quantify topology and ASR differences

**Sampling Robustness:**
- "Full" dataset (all 90% NR representatives)
- "Lite" dataset (further clustered to 70% or 50% ID)
- Outlier removal (long branches, suspected HGT)
- Re-infer trees and compare support values

**Model Robustness:**
- Site-homogeneous (e.g., LG+G) vs mixture models (LG+C60+F+R)
- Compare ASR posterior probability distributions
- Quantify model adequacy (posterior predictive checks if Bayesian)

**Documentation:**
- `docs/robustness.md`
- Summary tables of topology concordance
- ASR stability metrics

**QC Gate:**
- Key nodes retain high support (≥80 UFBoot) across settings
- Core ancestral residues stable across ≥2 tree/model settings

---

## 4. Computational Infrastructure

### 4.1 Environment Management

**Conda Environment:** `env/dah7ps.yaml`

*Core dependencies:*
- HMMER 3.3+
- CD-HIT 4.8+ OR MMseqs2
- MAFFT 7.5+
- MUSCLE 5.1+
- trimAl 1.4+
- BMGE 1.12+
- IQ-TREE 2.2+
- PAML 4.10+
- SignalP 6.0
- TargetP 2.0
- Python 3.9+ (BioPython, pandas, seaborn, ete3)
- R 4.2+ (phytools, ggtree, ape, PastML)
- PyMOL 2.5+

**Lockfiles:** `env/lockfiles/` for exact version reproducibility

### 4.2 Workflow Management

**Tool:** Snakemake or Nextflow

**Configuration:** `workflow/config.yaml`

*Key parameters:*
- Paths (input DBs, HMMs, reference proteomes)
- Thresholds (E-value, coverage, identity)
- Toggles (enable/disable optional steps)

**Features:**
- Modular rules per analysis phase
- Automatic checkpointing
- Logging to `logs/`
- Parameterized for easy re-runs

### 4.3 Data Management

**Version Control:**
- Git for code, scripts, docs, small config files
- Git LFS or external storage (Zenodo/OSF) for large data

**Reproducibility:**
- Lock environments
- Version tag releases
- DOI for data/results deposits

---

## 5. Deliverables

### 5.1 Data Products

1. **Curated Sequence Dataset**
   - Nonredundant FASTA
   - Metadata table with annotations
   - Cluster membership mapping

2. **Domain Architecture Annotations**
   - Domain coordinates
   - Masked sequences (transit peptides removed)

3. **Multiple Sequence Alignments**
   - Trimmed, partitioned alignments
   - Quality scores (GUIDANCE if used)

4. **Phylogenetic Trees**
   - Per-class high-support trees
   - Optional concatenated tree
   - Support values (UFBoot, SH-aLRT)

5. **Ancestral Sequences**
   - Reconstructed sequences for focal nodes
   - Posterior probability tables
   - Ambiguity lists

6. **Trait Reconstructions**
   - Trait evolution summaries
   - Effector switching events
   - Domain gain/loss timeline

7. **Structural Analyses**
   - Substitution maps
   - Visualization sessions (PyMOL)

### 5.2 Documentation

- `README.md`: Quickstart guide
- `docs/plan.md`: This document (project plan)
- `docs/lit_notes.md`: Literature review and references
- `docs/methods.md`: Detailed methods and parameters
- `docs/robustness.md`: Sensitivity analyses
- `docs/report.md`: Final results summary

### 5.3 Reproducible Workflow

- Complete Snakemake/Nextflow pipeline
- Environment files with locked versions
- Scripts for all custom analyses
- Example minimal dataset

### 5.4 Publication-Ready Figures

- Phylogenetic trees with trait annotations
- Domain architecture evolution timeline
- Structural maps of key substitutions
- Trait transition diagrams

---

## 6. Quality Control Gates

### Phase-Specific QC

1. **Sequence Collection**
   - Post-filter N per class: 50-200 representatives
   - Taxonomic breadth: ≥5 phyla for bacteria
   - KDO8PS contamination ≤1%

2. **Domain Annotation**
   - Manual spot-check of 10% for architecture accuracy
   - Targeting peptide prediction concordance (SignalP + TargetP)

3. **MSA Quality**
   - Catalytic residue conservation ≥95%
   - GUIDANCE column scores (if used) ≥0.5 for retained columns

4. **Phylogenetic Support**
   - ≥60% of splits with UFBoot ≥95
   - Key regulatory clades: UFBoot ≥95 AND SH-aLRT ≥80

5. **ASR Quality**
   - Median site PP ≥0.8 for catalytic core
   - Focal ancestors: ≥70% of sites with PP ≥0.9

6. **Trait Reconstruction**
   - Consistency across ≥2 tree/model settings
   - Biological plausibility checks (e.g., no ACT domain loss followed by identical ACT re-gain)

---

## 7. Timeline Estimate

*Assuming computational resources available and manual curation effort:*

- **Phase 0-1** (Setup & Environment): 1 week
- **Phase 2** (Sequence Collection): 2 weeks
- **Phase 3** (Domain Annotation): 2 weeks
- **Phase 4** (MSA): 1 week
- **Phase 5** (Phylogenetics): 1-2 weeks
- **Phase 6** (ASR): 1 week
- **Phase 7** (Trait Evolution): 1 week
- **Phase 8** (Structural Mapping): 2 weeks
- **Phase 9** (Robustness): 1 week
- **Phase 10** (Packaging & Documentation): 1 week

**Total:** ~12-14 weeks (3-3.5 months)

*Note: Timelines assume iterative refinement and manual quality checks.*

---

## 8. Key References

### Primary Literature

1. **PMID: 36041629**
   [To be summarized in `docs/lit_notes.md`]

2. **PMID: 34813062**
   [To be summarized in `docs/lit_notes.md`]

### Methods Resources

- **Pfam:** Mistry et al. (2021) Nucleic Acids Res.
- **HMMER:** Eddy (2011) PLoS Comput Biol.
- **IQ-TREE 2:** Minh et al. (2020) Mol Biol Evol.
- **MAFFT:** Katoh & Standley (2013) Mol Biol Evol.
- **trimAl:** Capella-Gutiérrez et al. (2009) Bioinformatics.
- **SignalP 6:** Teufel et al. (2022) Nat Biotechnol.
- **PastML:** Ishikawa et al. (2019) Mol Biol Evol.

---

## 9. Risk Mitigation

### Potential Challenges

1. **Taxonomic Sampling Bias**
   - *Risk:* Over-representation of model organisms
   - *Mitigation:* Actively sample underrepresented phyla; use UniRef for diversity

2. **Alignment Ambiguity in Regulatory Domains**
   - *Risk:* Low sequence identity in ACT/CM domains
   - *Mitigation:* Structure-guided alignment; partition-based analysis

3. **Long-Branch Attraction**
   - *Risk:* Fast-evolving lineages incorrectly grouped
   - *Mitigation:* Use mixture models; test alternate rootings; site-heterogeneous models

4. **ASR Uncertainty at Ancient Nodes**
   - *Risk:* Low posterior probabilities for deep ancestors
   - *Mitigation:* Report uncertainty; use ensemble reconstructions; validate with structural constraints

5. **Incomplete Phenotype Data**
   - *Risk:* Literature covers <10% of sequences
   - *Mitigation:* Focus trait analysis on well-characterized clades; predictive modeling for unknowns

6. **Computational Resources**
   - *Risk:* IQ-TREE mixture models computationally intensive
   - *Mitigation:* Use PMSF approximations; HPC access; partition data if needed

---

## 10. Success Criteria

**Project considered successful if:**

1. High-quality, nonredundant sequence dataset (50-200 per class, broad taxonomic spread)
2. Well-supported phylogenies (majority key nodes UFBoot ≥95)
3. High-confidence ancestral reconstructions (median PP ≥0.8 in core)
4. Clear evolutionary narrative for domain recruitment and effector switching
5. Reproducible workflow with full documentation
6. Publication-ready figures and data deposit

**Stretch Goals:**

- Experimental validation of 1-2 ancestral proteins
- Predictive model for effector specificity from sequence
- 3D structural modeling of ancestral enzymes (AlphaFold2)

---

## Document Revision History

| Version | Date       | Changes                              | Author |
|---------|------------|--------------------------------------|--------|
| 1.0     | 2025-11-07 | Initial project plan created         | Claude |

---

**END OF PROJECT PLAN**
