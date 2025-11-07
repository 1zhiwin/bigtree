# Analysis Scripts

This directory contains Python, R, and shell scripts used by the Snakemake workflow.

## Structure

- `check_installation.py` - Verify all dependencies are installed
- `filter_sequences.py` - Quality filtering of HMMER hits (to be created)
- `cluster_sequences.py` - Redundancy reduction with CD-HIT/MMseqs2 (to be created)
- `predict_localization.py` - Signal peptide and organellar targeting prediction (to be created)
- `create_metadata.py` - Compile comprehensive metadata table (to be created)
- `extract_partitions.py` - Extract domain-specific sequences (to be created)
- `extract_ancestral_sequences.py` - Parse IQ-TREE ASR output (to be created)
- `prepare_trait_matrix.py` - Create trait matrix for evolution analysis (to be created)
- `reconstruct_traits.R` - Ancestral trait reconstruction in R (to be created)
- `generate_report.py` - Compile final results summary (to be created)
- `qc_sequences.py` - Quality control checks (to be created)
- `qc_tree_support.py` - Tree support validation (to be created)

## Usage

Scripts are called by Snakemake rules. They can also be run standalone:

```bash
# Example: Check installation
python scripts/check_installation.py

# Example: Run with arguments
python scripts/filter_sequences.py \
    --input data/processed/hmmer_hits.tblout \
    --proteome data/raw/uniprot.faa \
    --output seqs/candidates.faa \
    --min-coverage 0.85
```

## Development

When creating new scripts:

1. Add argparse for command-line arguments
2. Include logging to `logs/` directory
3. Add docstrings and comments
4. Test with small datasets first
5. Update this README with script description

## Testing

```bash
# Run installation check
python scripts/check_installation.py

# Add unit tests here as they're developed
pytest tests/
```
