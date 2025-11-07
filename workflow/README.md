# Workflow Documentation

## Overview

This directory contains the Snakemake workflow for the DAH7PS allostery evolution and ASR project. The workflow is organized into modular phases that can be run independently or as a complete pipeline.

## Quick Start

```bash
# Activate the conda environment
conda activate dah7ps

# Configure the workflow
# Edit config.yaml to set your paths and parameters

# Dry run to see what will be executed
snakemake -n

# Run the full pipeline
snakemake --cores 8 --use-conda all

# Run specific phase
snakemake --cores 8 sequence_collection
```

## Workflow Structure

### Configuration

**Main config file:** `config.yaml`

This file contains all parameters for the analysis:
- Input/output paths
- Tool parameters
- Quality thresholds
- Taxonomic sampling criteria

### Pipeline Phases

The workflow is organized into sequential phases:

#### Phase 1: Sequence Collection
```bash
snakemake --cores 8 seqs/dah7ps_nonredundant.faa
```

**Steps:**
1. `download_pfam` - Download Pfam-A HMM database
2. `sequence_search` - HMMER search for DAH7PS candidates
3. `filter_sequences` - Quality filtering
4. `cluster_sequences` - Redundancy reduction at 90% identity

**Outputs:**
- `seqs/dah7ps_candidates.faa`
- `seqs/dah7ps_nonredundant.faa`
- `data/processed/cluster_map.tsv`

#### Phase 2: Domain Annotation
```bash
snakemake --cores 8 data/processed/metadata.tsv
```

**Steps:**
1. `annotate_domains` - HMMER scan for ACT and CM domains
2. `predict_localization` - Signal peptide and organellar targeting
3. `create_metadata` - Compile comprehensive metadata table

**Outputs:**
- `data/processed/domain_architecture.tsv`
- `data/processed/localization_predictions.tsv`
- `data/processed/metadata.tsv`

#### Phase 3: Multiple Sequence Alignment
```bash
snakemake --cores 8 msa/
```

**Steps:**
1. `extract_partitions` - Separate core, ACT, CM domains
2. `align_mafft` - MAFFT alignment (L-INS-i)
3. `trim_alignment` - trimAl quality trimming

**Outputs:**
- `msa/core_Ialpha.trim.faa`
- `msa/core_Ibeta.trim.faa`
- `msa/core_II.trim.faa`
- `msa/ACT.trim.faa`
- `msa/CM.trim.faa`

#### Phase 4: Phylogenetic Inference
```bash
snakemake --cores 8 trees/
```

**Steps:**
1. `infer_tree` - IQ-TREE 2 with ModelFinder, UFBoot, SH-aLRT

**Outputs:**
- `trees/*.treefile` - Best ML trees
- `trees/*.log` - IQ-TREE logs
- `trees/*.iqtree` - Full IQ-TREE output

#### Phase 5: Ancestral Sequence Reconstruction
```bash
snakemake --cores 8 asr/
```

**Steps:**
1. `ancestral_reconstruction` - IQ-TREE marginal ASR
2. `extract_ancestral_sequences` - Parse and format results

**Outputs:**
- `asr/*.state` - IQ-TREE ancestral state files
- `asr/*_ancestral.faa` - Reconstructed sequences
- `asr/asr_summary.md`

#### Phase 6: Trait Evolution
```bash
snakemake --cores 8 traits/
```

**Steps:**
1. `prepare_trait_matrix` - Create discrete trait table
2. `reconstruct_traits` - PastML or phytools reconstruction

**Outputs:**
- `traits/trait_matrix.tsv`
- `traits/trait_reconstruction.tsv`

#### Phase 7: Reporting
```bash
snakemake --cores 8 docs/report.md
```

**Steps:**
1. `generate_report` - Compile results summary

**Outputs:**
- `docs/report.md`

## Running the Workflow

### Basic Execution

```bash
# Run all phases
snakemake --cores 8 all

# Run with specific number of cores
snakemake --cores 16 all

# Dry run (show what would be executed)
snakemake -n all

# Show execution plan as DAG
snakemake --dag all | dot -Tpdf > workflow_dag.pdf
```

### Targeted Execution

```bash
# Run only sequence collection
snakemake --cores 8 seqs/dah7ps_nonredundant.faa

# Run only phylogenetics for one partition
snakemake --cores 8 trees/core_Ialpha.treefile

# Force re-run of specific rule
snakemake --cores 8 --forcerun infer_tree trees/core_Ialpha.treefile
```

### Cluster Execution

For HPC environments:

```bash
# SLURM cluster
snakemake --cluster "sbatch -p normal -t 24:00:00 -c {threads} --mem=32G" \
          --jobs 10 \
          --cores 80 \
          all

# PBS cluster
snakemake --cluster "qsub -l nodes=1:ppn={threads},mem=32gb,walltime=24:00:00" \
          --jobs 10 \
          all
```

Or use cluster profiles (recommended):

```bash
# Create SLURM profile
mkdir -p ~/.config/snakemake/slurm

# Run with profile
snakemake --profile slurm all
```

## Quality Control

The workflow includes QC checkpoints:

```bash
# Check sequence quality
snakemake --cores 1 logs/qc_sequences.txt

# Check tree support
snakemake --cores 1 logs/qc_tree_support.txt
```

## Cleaning Up

```bash
# Remove intermediate files only
snakemake clean

# Remove ALL generated files (careful!)
snakemake clean_all
```

## Troubleshooting

### Common Issues

**1. Rule fails due to missing input**
```bash
# Check what inputs are missing
snakemake -n --reason all
```

**2. Out of memory errors**
```bash
# Increase memory allocation in cluster submission
# Or edit config.yaml to reduce dataset size
```

**3. Tool not found**
```bash
# Verify environment
conda activate dah7ps
python scripts/check_installation.py
```

**4. Checkpoint corruption**
```bash
# Remove Snakemake metadata and restart
rm -rf .snakemake
snakemake --cores 8 all
```

### Debugging

```bash
# Verbose output
snakemake --cores 8 -p all

# Print shell commands
snakemake --cores 8 -p all

# Detailed debugging
snakemake --cores 8 --debug-dag all
```

## Customization

### Adding New Rules

Create a new rule in `Snakefile` or in a separate file in `rules/`:

```python
rule my_new_rule:
    input:
        "path/to/input"
    output:
        "path/to/output"
    log:
        "logs/my_new_rule.log"
    script:
        "scripts/my_script.py"
```

### Modifying Parameters

Edit `config.yaml`:

```yaml
sequence_collection:
  clustering:
    identity_threshold: 0.95  # Change from 0.90 to 0.95
```

### Adding New Scripts

1. Create script in `scripts/` directory
2. Reference in Snakefile rule:
   ```python
   script:
       "scripts/my_script.py"
   ```

## Performance Optimization

### Parallel Execution

```bash
# Run up to 4 jobs in parallel
snakemake --cores 32 --jobs 4 all
```

### Resource Specification

In Snakefile:
```python
rule resource_intensive:
    input: ...
    output: ...
    threads: 16
    resources:
        mem_mb=64000,
        runtime=1440  # 24 hours in minutes
```

### Conda Integration

Use per-rule environments:
```python
rule my_rule:
    input: ...
    output: ...
    conda:
        "envs/special_env.yaml"
```

## Best Practices

1. **Always dry run first**: `snakemake -n all`
2. **Use version control**: Commit `config.yaml` changes
3. **Document parameters**: Add comments in `config.yaml`
4. **Check logs**: Review `logs/` directory after runs
5. **Validate outputs**: Use QC rules before proceeding
6. **Back up results**: Copy key outputs to permanent storage

## Workflow Reproducibility

### Exact Reproduction

```bash
# Use locked environment
conda create --name dah7ps_locked --file ../env/lockfiles/dah7ps_YYYYMMDD.txt

# Run with exact parameters
snakemake --cores 8 --use-conda all

# Generate reproducibility report
snakemake --report report.html
```

### Archiving

```bash
# Archive entire workflow state
snakemake --archive workflow_archive.tar.gz all
```

## Support

For issues or questions:
1. Check logs in `logs/` directory
2. Review Snakemake documentation: https://snakemake.readthedocs.io/
3. Open an issue on the project repository

## References

- **Snakemake**: Köster & Rahmann (2012) Bioinformatics
- **IQ-TREE**: Minh et al. (2020) Mol Biol Evol
- **HMMER**: Eddy (2011) PLoS Comput Biol
