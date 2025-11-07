# DAH7PS Allostery Evolution & Ancestral Sequence Reconstruction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-TBD-blue)](https://doi.org/TBD)

## Overview

This repository contains a comprehensive computational pipeline for studying the evolutionary history of **DAH7PS** (3-deoxy-D-arabinoheptulosonate 7-phosphate synthase) enzymes, with a focus on the evolution of **allosteric regulation** mechanisms.

### Scientific Goals

- **Phylogenetic Analysis**: Reconstruct evolutionary relationships across DAH7PS Types Iα, Iβ, and II
- **Ancestral Sequence Reconstruction (ASR)**: Resurrect ancient enzyme sequences at key evolutionary transitions
- **Trait Evolution**: Map the evolution of regulatory mechanisms (ACT domains, chorismate mutase fusions)
- **Mechanistic Understanding**: Identify structural determinants of effector specificity and synergistic inhibition

### Key Features

- 🧬 **Multi-kingdom sampling**: Bacteria, Archaea, Plants, Fungi
- 🔍 **Domain architecture analysis**: ACT domains, chorismate mutase fusions
- 🌳 **Robust phylogenetics**: IQ-TREE 2 with mixture models and comprehensive support metrics
- 🧪 **High-quality ASR**: Marginal reconstructions with uncertainty quantification
- 📊 **Trait mapping**: Evolution of effector specificity (Tyr/Phe/Trp) and inhibition modes
- 🔬 **Structure mapping**: Ancestral substitutions on solved crystal structures
- ♻️ **Fully reproducible**: Conda environments, Snakemake workflow, version-controlled

---

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- Git and Git LFS
- 16+ GB RAM (32+ GB recommended)
- 50+ GB disk space

### Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/bigtree.git
cd bigtree

# Create and activate the conda environment
conda env create -f env/dah7ps.yaml
conda activate dah7ps

# Download Pfam database (required for HMMER)
cd data/raw
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
cd ../..

# Verify installation
python scripts/check_installation.py  # To be created
```

### Running the Pipeline

```bash
# Configure the workflow
# Edit workflow/config.yaml to set paths and parameters

# Run full pipeline with Snakemake
snakemake --cores 8 --use-conda all

# Or run specific phases:
snakemake --cores 8 sequence_collection
snakemake --cores 8 domain_annotation
snakemake --cores 8 alignment
snakemake --cores 8 phylogenetics
snakemake --cores 8 asr
```

### Quick Test Run

```bash
# Run with test dataset (smaller, faster)
snakemake --cores 4 --use-conda test_run
```

---

## Project Structure

```
bigtree/
├── data/
│   ├── raw/              # Input databases, proteomes
│   └── processed/        # Curated metadata, filtered sequences
├── seqs/
│   ├── *.faa             # Sequence collections
│   └── domains/          # Domain-specific sequences
├── msa/                  # Multiple sequence alignments
├── trees/                # Phylogenetic trees
├── asr/                  # Ancestral sequence reconstructions
├── traits/               # Trait evolution analyses
├── structures/
│   ├── pdbs/             # Crystal structures
│   └── mapping/          # Substitution mappings
├── figs/                 # Publication figures
├── workflow/
│   ├── config.yaml       # Main configuration
│   ├── Snakefile         # Workflow rules (to be created)
│   └── rules/            # Modular Snakemake rules
├── scripts/              # Analysis scripts
├── notebooks/            # Jupyter notebooks for exploration
├── env/
│   ├── dah7ps.yaml       # Conda environment
│   └── lockfiles/        # Exact version snapshots
├── docs/
│   ├── plan.md           # Detailed project plan
│   ├── lit_notes.md      # Literature review
│   ├── methods.md        # Methods documentation
│   └── report.md         # Final results summary
└── logs/                 # Execution logs
```

---

## Scientific Background

### DAH7PS Enzyme

**Function**: Catalyzes the first committed step of the shikimate pathway
**Reaction**: Phosphoenolpyruvate (PEP) + Erythrose 4-phosphate (E4P) → DAH7P + Pi
**EC Number**: 2.5.1.54
**Importance**: Essential for aromatic amino acid biosynthesis; absent in animals (drug target)

### Enzyme Classes

| Type | Fold | Distribution | Regulation |
|------|------|--------------|------------|
| **Iα** | (β/α)₈ barrel | Bacteria | ACT domains (Tyr/Phe/Trp) |
| **Iβ** | (β/α)₈ barrel | Bacteria | ACT domains (divergent) |
| **II** | β-barrel | Archaea, some Bacteria | Less characterized |

### Regulatory Mechanisms

1. **ACT Domain Allosteric Inhibition**
   - Small molecule binding domains
   - Effectors: Aromatic amino acids (Tyr, Phe, Trp)
   - Single vs synergistic inhibition

2. **Chorismate Mutase (CM) Fusion**
   - Bifunctional enzyme
   - Feedback from pathway intermediates
   - Metabolic channeling

### Evolutionary Questions

1. How many times were ACT domains recruited to DAH7PS?
2. What drives effector specificity (Tyr vs Phe vs Trp)?
3. How did synergistic inhibition evolve?
4. Are there convergent regulatory solutions in different lineages?
5. What were the regulatory properties of ancestral enzymes?

---

## Methodology Overview

### Phase 1: Sequence Collection
- HMMER searches against UniProt, RefSeq
- Quality filtering (coverage, ambiguity)
- Redundancy reduction (90% identity clustering)
- Taxonomic sampling across Bacteria, Archaea, Plants, Fungi

### Phase 2: Domain Annotation
- HMMER domain scans (ACT, CM)
- Localization prediction (SignalP, TargetP)
- Architecture classification
- Phenotype curation from literature

### Phase 3: Multiple Sequence Alignment
- Partition-based strategy (core, ACT, CM separate)
- MAFFT L-INS-i (primary), MUSCLE5 (comparison)
- trimAl/BMGE for quality trimming
- Conservation of catalytic residues

### Phase 4: Phylogenetic Inference
- IQ-TREE 2 with ModelFinder
- Mixture models (LG+C60+F+R) for heterogeneity
- UFBoot2 + SH-aLRT support (1000 replicates each)
- Per-class trees (Iα, Iβ, II)

### Phase 5: Ancestral Sequence Reconstruction
- Marginal reconstruction (IQ-TREE --asr)
- Per-partition ASR (core, ACT, CM)
- Uncertainty quantification (posterior probabilities)
- Focus on LCAs and domain gain/loss nodes

### Phase 6: Trait Evolution
- Discrete trait reconstruction (PastML)
- Domain presence/absence
- Effector specificity evolution
- Inhibition mode transitions

### Phase 7: Structural Mapping
- Map ancestral substitutions to crystal structures
- Interface residue analysis
- Convergent evolution detection
- PyMOL visualization

### Phase 8: Robustness Testing
- Alignment sensitivity (MAFFT vs MUSCLE, trimming methods)
- Sampling robustness (full vs lite datasets)
- Model sensitivity (site-homogeneous vs mixture)

---

## Key Dependencies

**Bioinformatics Tools:**
- HMMER 3.4 (sequence search)
- CD-HIT / MMseqs2 (clustering)
- MAFFT 7.5+ / MUSCLE 5 (alignment)
- trimAl / BMGE (alignment trimming)
- IQ-TREE 2.2+ (phylogenetics, ASR)
- PAML 4.10+ (optional ASR)
- SignalP 6 / TargetP 2 (localization)

**Data Analysis:**
- Python 3.10+ (BioPython, pandas, ete3)
- R 4.3+ (ape, phytools, ggtree)
- PyMOL (structure visualization)

**Workflow:**
- Snakemake 7.3+ (pipeline management)

See `env/dah7ps.yaml` for complete dependency list.

---

## Quality Control Gates

| Phase | Metric | Threshold |
|-------|--------|-----------|
| Sequence Collection | Representatives per class | 50–200 |
| Sequence Collection | Bacterial phyla diversity | ≥5 |
| Sequence Collection | KDO8PS contamination | ≤1% |
| Phylogenetics | UFBoot support (key nodes) | ≥95 |
| Phylogenetics | SH-aLRT support (key nodes) | ≥80 |
| ASR | Median posterior probability (core) | ≥0.80 |
| ASR | Focal ancestor site coverage (PP≥0.9) | ≥70% |

---

## Expected Outputs

### Data Products
- Curated nonredundant sequence set with metadata
- Domain architecture annotations
- High-quality trimmed alignments
- Well-supported phylogenetic trees
- Ancestral sequences with uncertainty metrics
- Trait evolution reconstructions

### Figures
- Phylogenetic trees with domain architecture annotations
- Effector specificity evolution timeline
- Structural maps of key substitutions
- Trait transition diagrams

### Reports
- Comprehensive methods documentation
- Results summary with biological interpretation
- Robustness and sensitivity analyses

---

## Documentation

- **[Project Plan](docs/plan.md)**: Detailed methodology and rationale
- **[Literature Notes](docs/lit_notes.md)**: Key references and background
- **[Methods](docs/methods.md)**: Step-by-step protocols (to be created)
- **[Report](docs/report.md)**: Final results (to be created)

---

## Citation

If you use this pipeline or data, please cite:

```
[Citation to be added upon publication]
```

### Key References

1. **PMID: 36041629** - [To be filled]
2. **PMID: 34813062** - [To be filled]

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit your changes (`git commit -am 'Add new analysis'`)
4. Push to the branch (`git push origin feature/new-analysis`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- **Pfam**: Domain annotations
- **UniProt**: Curated protein sequences
- **NCBI**: RefSeq genomes
- **IQ-TREE Team**: Phylogenetic software
- **All contributors** to open-source bioinformatics tools

---

## Contact

**Project Lead**: [Your Name]
**Email**: [your.email@institution.edu]
**Institution**: [Your Institution]

**Issues & Questions**: Please use the [GitHub Issues](https://github.com/YOUR_USERNAME/bigtree/issues) page

---

## Roadmap

### Phase 0-1: Setup ✅
- [x] Directory structure
- [x] Environment configuration
- [x] Documentation framework

### Phase 2-3: Data Collection & Annotation
- [ ] Sequence collection
- [ ] Domain annotation
- [ ] Phenotype curation

### Phase 4: Alignment & Phylogenetics
- [ ] Multiple sequence alignments
- [ ] Phylogenetic trees

### Phase 5-6: ASR & Trait Evolution
- [ ] Ancestral reconstructions
- [ ] Trait mapping

### Phase 7-8: Analysis & Validation
- [ ] Structural mapping
- [ ] Robustness testing

### Phase 9-10: Publication
- [ ] Manuscript preparation
- [ ] Data deposition
- [ ] Code release

---

## Version History

| Version | Date | Description |
|---------|------|-------------|
| 0.1.0 | 2025-11-07 | Initial project structure |

---

**Project initiated:** 2025-11-07
**Last updated:** 2025-11-07
