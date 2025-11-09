# DAH7PS Project - Installation Status

**Date:** 2025-11-09
**Environment:** dah7ps
**Status:** ✅ Core installation complete and functional

---

## Installation Summary

### ✅ Successfully Installed (Core Tools)

**Sequence Analysis:**
- ✅ HMMER 3.4 (hmmsearch, hmmscan) - Sequence search
- ✅ CD-HIT 4.8.1 - Sequence clustering
- ✅ MMseqs2 18.8cc5c - Fast sequence clustering
- ✅ BLAST 2.17.0 - Sequence alignment
- ✅ SeqKit, BioAwk, SeqTK - Sequence manipulation

**Multiple Sequence Alignment:**
- ✅ MAFFT 7.5+ - High-quality alignments
- ✅ MUSCLE 5.3 - Alternative aligner
- ✅ Clustal Omega - Additional aligner
- ✅ trimAl 1.5 - Alignment trimming
- ✅ BMGE - Entropy-based trimming

**Phylogenetic Inference:**
- ✅ IQ-TREE 3.0.1 - Main phylogenetic tool (with ASR)
- ✅ RAxML-NG 1.2.2 - Alternative phylogenetic inference
- ✅ FastTree - Fast approximation

**Ancestral Reconstruction:**
- ✅ IQ-TREE (built-in --asr)
- ✅ PAML (codeml available)

**Python Environment (3.10.19):**
- ✅ BioPython 1.86
- ✅ pandas 2.3.3
- ✅ numpy 2.2.6
- ✅ scipy 1.15.2
- ✅ matplotlib 3.10.7
- ✅ seaborn 0.13.2
- ✅ ete3 3.1.3 - Tree manipulation
- ✅ dendropy 5.0.8 - Phylogenetics
- ✅ snakemake 7.32.4 - Workflow management
- ✅ pyhmmer 0.11.2 - Python HMMER bindings
- ✅ biotite 1.2.0 - Structural bioinformatics
- ✅ scikit-bio 0.7.1 - Biology algorithms
- ✅ logomaker 0.8.7 - Sequence logos

**R Environment (4.5.2):**
- ✅ ape 5.8.1 - Phylogenetic analyses
- ✅ phytools 2.5.2 - Phylogenetic tools
- ✅ ggplot2 4.0.0 - Visualization
- ✅ phangorn - Phylogenetic reconstruction
- ✅ treeio 1.34.0 - Tree I/O (newly installed)
- ✅ aplot 0.2.9 - Plot utilities (newly installed)

**Workflow & Utilities:**
- ✅ Snakemake 7.32.4
- ✅ Jupyter Lab & Notebooks
- ✅ Git, wget, curl
- ✅ csvtk, jq - Data processing

---

## ⚠️ Optional Tools (Not Installed)

These tools are optional for the core workflow but may be useful for specific analyses:

**Protein Localization Prediction:**
- ⚠️ SignalP 6.0 - Requires manual download from DTU (academic license)
- ⚠️ TargetP 2.0 - Requires manual download from DTU (academic license)
- Alternative: DeepSig is available in conda if needed

**Visualization:**
- ⚠️ ggtree (R) - Failed due to system dependency issues with gdtools
  - Not critical: tree visualization works with ape, phytools, and ete3
  - Can be installed later if needed for publication-quality figures

**Domain Analysis:**
- ⚠️ InterProScan - Large download, optional for domain annotation
  - HMMER with Pfam is sufficient for this project

---

## System Requirements Met

- ✅ Python 3.10
- ✅ R 4.5.2
- ✅ All bioinformatics command-line tools
- ✅ Complete scientific Python stack
- ✅ Essential R phylogenetic packages
- ✅ Workflow management system

---

## Next Steps

### 1. Download Pfam Database (Required for sequence collection)

```bash
# Activate environment
conda activate dah7ps

# Download Pfam-A HMM library (~2 GB)
cd data/raw
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
cd ../..
```

This creates four index files:
- Pfam-A.hmm.h3f
- Pfam-A.hmm.h3i
- Pfam-A.hmm.h3m
- Pfam-A.hmm.h3p

### 2. Obtain Sequence Data

**UniProt Reference Proteomes:**
- Download from: https://www.uniprot.org/proteomes
- Place in `data/raw/uniprot_reference_proteomes/`

**Or use NCBI RefSeq:**
- Download specific organisms from: https://ftp.ncbi.nlm.nih.gov/refseq/

### 3. Test the Workflow

```bash
# Activate environment
conda activate dah7ps

# Run workflow dry-run to check syntax
snakemake -n

# When ready, run with test data
snakemake --cores 4 test_run
```

### 4. Optional: Install Missing Visualization Tools

If you need ggtree for publication figures:

```bash
# Install system dependencies first (Ubuntu/Debian)
sudo apt-get install libgdal-dev libcairo2-dev

# Then install ggtree
conda activate dah7ps
Rscript scripts/install_ggtree.R
```

---

## Quick Start Guide

```bash
# Activate environment (do this every time)
conda activate dah7ps

# Verify tools are working
python scripts/check_installation.py

# View project structure
ls -R

# Read documentation
cat docs/plan.md
cat workflow/README.md

# Edit configuration
nano workflow/config.yaml

# Run workflow (when ready)
snakemake --cores 8 all
```

---

## Configuration

**Main config file:** `workflow/config.yaml`

Key parameters to configure:
- Input data paths
- HMMER E-value thresholds
- Clustering identity (currently 90%)
- Alignment parameters
- Phylogenetic model settings
- ASR parameters

---

## Project Status

**Phase 0-1:** ✅ Setup complete
- [x] Environment created
- [x] Tools installed and verified
- [x] Directory structure in place
- [x] Documentation complete

**Phase 2:** ⏳ Ready to start
- [ ] Download Pfam database
- [ ] Collect sequence data
- [ ] Run HMMER searches
- [ ] Quality filtering
- [ ] Redundancy reduction

**Phase 3-8:** ⏸️ Pending
- [ ] Domain annotation
- [ ] Multiple sequence alignments
- [ ] Phylogenetic inference
- [ ] Ancestral sequence reconstruction
- [ ] Trait evolution analysis
- [ ] Structural mapping
- [ ] Robustness testing

---

## Troubleshooting

### Issue: Command not found

```bash
# Make sure environment is activated
conda activate dah7ps

# Verify tool is installed
which iqtree
which hmmsearch
```

### Issue: Out of memory

Edit `workflow/config.yaml` to:
- Reduce clustering threshold
- Use smaller subsets for testing
- Limit parallel jobs

### Issue: Snakemake errors

```bash
# Remove cached metadata
rm -rf .snakemake

# Dry run to check syntax
snakemake -n all

# Verbose output
snakemake -p --cores 1 all
```

---

## Support

- **Documentation:** `docs/` directory
- **Workflow help:** `workflow/README.md`
- **Environment help:** `env/README.md`
- **Project plan:** `docs/plan.md`

---

## Environment Details

**Location:** `/home/luogy/miniconda3/envs/dah7ps`

**To activate:**
```bash
conda activate dah7ps
```

**To deactivate:**
```bash
conda deactivate
```

**To update:**
```bash
conda env update -f env/dah7ps_minimal.yaml --prune
```

**To export (for reproducibility):**
```bash
conda list --explicit > env/lockfiles/dah7ps_$(date +%Y%m%d).txt
```

---

## Resources Required

**Disk Space:**
- Environment: ~5 GB (installed)
- Pfam database: ~2 GB
- Sequence data: 10-50 GB (depends on sampling)
- Working data: 10-100 GB (alignments, trees, etc.)
- **Total estimate:** 30-160 GB

**Memory:**
- Minimum: 16 GB RAM
- Recommended: 32 GB RAM
- IQ-TREE with mixture models may need 64+ GB for large datasets

**CPU:**
- All tools support multi-threading
- Recommended: 8+ cores for efficient analysis

---

## Installation Log

Installation logs available in:
- `logs/conda_install_minimal.log` - Full conda installation log
- Check for any warnings or errors if issues arise

---

**Status:** Ready to begin Phase 2 (Sequence Collection)!
