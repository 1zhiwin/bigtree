# Environment Setup

## Quick Start

```bash
# Create the conda environment
conda env create -f dah7ps.yaml

# Activate the environment
conda activate dah7ps

# Verify installation
python scripts/check_installation.py  # To be created
```

## Environment Management

### Creating a Lockfile

For reproducibility, create explicit lockfiles:

```bash
# Create platform-specific lockfile
conda list --explicit > lockfiles/dah7ps_$(date +%Y%m%d)_$(uname -s).txt

# Create cross-platform environment export
conda env export > lockfiles/dah7ps_$(date +%Y%m%d).yaml
```

### Updating the Environment

```bash
# Update from YAML file
conda env update -f dah7ps.yaml --prune

# Update specific package
conda update <package_name>
```

### Recreating from Lockfile

```bash
# From explicit lockfile (platform-specific)
conda create --name dah7ps --file lockfiles/dah7ps_YYYYMMDD_Linux.txt

# From exported YAML
conda env create -f lockfiles/dah7ps_YYYYMMDD.yaml
```

## Manual Installations Required

Some tools require manual installation:

### 1. SignalP 6.0

```bash
# Download from DTU Health Tech (requires academic license)
# https://services.healthtech.dtu.dk/service.php?SignalP-6.0

# Follow installation instructions
# Add to PATH or specify in workflow/config.yaml
```

### 2. TargetP 2.0

```bash
# Download from DTU Health Tech
# https://services.healthtech.dtu.dk/service.php?TargetP-2.0

# Follow installation instructions
# Add to PATH or specify in workflow/config.yaml
```

### 3. Pfam Database

```bash
# Download latest Pfam-A HMM library
cd data/raw/
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# Press the HMM database for hmmsearch/hmmscan
hmmpress Pfam-A.hmm

# This creates .h3f, .h3i, .h3m, .h3p files
```

### 4. InterProScan (Optional)

```bash
# Large download; optional for this project
# https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html
```

### 5. ChimeraX (Optional structural visualization)

Download from: https://www.cgl.ucsf.edu/chimerax/download.html

## Testing the Installation

Run the installation check script:

```bash
python scripts/check_installation.py
```

This will verify:
- Python packages
- Bioinformatics tools (HMMER, MAFFT, IQ-TREE, etc.)
- R packages
- Optional tools

## Troubleshooting

### Common Issues

**1. Conda solver takes too long**
```bash
# Use mamba as a faster solver
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

**2. R packages fail to install**
```bash
# Activate environment and install manually
conda activate dah7ps
R
> install.packages("phytools")
> BiocManager::install("ggtree")
```

**3. PyMOL not working**
```bash
# Try open-source version
conda install -c conda-forge pymol-open-source

# Or install proprietary version separately
```

**4. SignalP/TargetP license issues**
```bash
# Use alternative tools already in environment:
# - DeepSig for signal peptides
# - DeepLoc for subcellular localization
```

## Environment Validation

After installation, verify key tools:

```bash
# HMMER
hmmsearch -h
hmmscan -h

# Sequence alignment
mafft --version
muscle -version

# Phylogenetics
iqtree2 --version
raxmlHPC-NG --version

# Python environment
python -c "import Bio; print(Bio.__version__)"
python -c "import ete3; print(ete3.__version__)"

# R environment
Rscript -e "library(ape); library(phytools); cat('R packages OK\n')"
```

## Resource Requirements

**Disk Space:**
- Conda environment: ~5-10 GB
- Pfam database: ~2 GB
- Working data (depends on project): 10-100 GB

**Memory:**
- Minimum: 16 GB RAM
- Recommended: 32+ GB RAM for large phylogenetic analyses
- IQ-TREE mixture models: may require 64+ GB

**CPU:**
- All major tools support multi-threading
- Recommended: 8+ cores for parallel processing

## Reproducibility Notes

- Environment created: 2025-11-07
- Conda version: 23.x+
- Platform: Linux (primary), macOS (compatible)
- Python: 3.10
- R: 4.3.1

For maximum reproducibility, use the lockfiles in `lockfiles/` directory.
