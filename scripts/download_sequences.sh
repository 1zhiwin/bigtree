#!/bin/bash
#
# Download sequence data for DAH7PS project
# This script provides multiple options for obtaining protein sequences
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
RAW_DIR="$PROJECT_DIR/data/raw"

echo "============================================================"
echo "DAH7PS Project - Sequence Data Download"
echo "============================================================"
echo ""

# Create directory if needed
mkdir -p "$RAW_DIR/proteomes"
cd "$RAW_DIR"

echo "Available options:"
echo ""
echo "1. Test dataset (small, quick start) - 100 MB"
echo "2. Representative bacteria (medium) - ~5 GB"
echo "3. Comprehensive multi-kingdom (large) - ~20 GB"
echo "4. Custom download instructions"
echo ""

# Function to download test dataset
download_test() {
    echo "Downloading test dataset (E. coli + selected model organisms)..."
    echo ""

    # E. coli K12 (has aroF, aroG, aroH)
    echo "Downloading E. coli K-12..."
    wget -O proteomes/ecoli_k12.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz"

    # Bacillus subtilis
    echo "Downloading B. subtilis..."
    wget -O proteomes/bsubtilis.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000001570/UP000001570_224308.fasta.gz"

    # Pseudomonas aeruginosa
    echo "Downloading P. aeruginosa..."
    wget -O proteomes/paeruginosa.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000002438/UP000002438_208964.fasta.gz"

    # Mycobacterium tuberculosis
    echo "Downloading M. tuberculosis..."
    wget -O proteomes/mtb.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000001584/UP000001584_83332.fasta.gz"

    # Yeast (S. cerevisiae - has bifunctional enzymes)
    echo "Downloading S. cerevisiae..."
    wget -O proteomes/scerevisiae.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz"

    # Arabidopsis thaliana (plant)
    echo "Downloading A. thaliana..."
    wget -O proteomes/athaliana.faa.gz \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006548/UP000006548_3702.fasta.gz"

    echo ""
    echo "Decompressing files..."
    gunzip proteomes/*.gz

    echo ""
    echo "Combining into single file..."
    cat proteomes/*.faa > test_proteomes.faa

    echo ""
    echo "✓ Test dataset ready: $RAW_DIR/test_proteomes.faa"
    echo "  Size: $(du -h test_proteomes.faa | cut -f1)"
    echo "  Sequences: $(grep -c '^>' test_proteomes.faa)"
}

# Function to download representative bacteria
download_bacteria() {
    echo "This will download ~200 bacterial reference proteomes (~5 GB)"
    echo "Estimated time: 30-60 minutes depending on connection"
    echo ""
    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled."
        exit 1
    fi

    echo ""
    echo "Downloading UniProt bacterial reference proteomes list..."
    wget -O bacteria_proteomes_list.txt \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"

    echo ""
    echo "For a comprehensive download, use the UniProt FTP:"
    echo "  rsync -avz --progress \\"
    echo "    rsync://ftp.uniprot.org/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/ \\"
    echo "    proteomes/bacteria/"
    echo ""
    echo "Or visit: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/"
}

# Show custom instructions
show_custom() {
    cat << 'EOF'

============================================================
Custom Download Instructions
============================================================

Option 1: UniProt Reference Proteomes (Recommended)
----------------------------------------------------

Browse and download from:
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/

Directory structure:
- Bacteria/     - Bacterial proteomes
- Archaea/      - Archaeal proteomes
- Eukaryota/    - Eukaryotic proteomes (includes plants, fungi)

Example using rsync for specific phyla:
```bash
# Download all Proteobacteria
rsync -avz --progress \
  rsync://ftp.uniprot.org/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/ \
  --include='*Proteobacteria*' --include='*/' --exclude='*' \
  data/raw/proteomes/

# Download all Actinobacteria
rsync -avz --progress \
  rsync://ftp.uniprot.org/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/ \
  --include='*Actinobacteria*' --include='*/' --exclude='*' \
  data/raw/proteomes/
```

Option 2: NCBI RefSeq (Alternative)
------------------------------------

For specific organisms:
1. Go to: https://www.ncbi.nlm.nih.gov/genome/
2. Search for organism
3. Download protein FASTA files

For bulk download:
```bash
# Install NCBI datasets tool
conda install -c conda-forge ncbi-datasets-cli

# Download bacterial genomes by taxonomy
datasets download genome taxon "Bacteria" \
  --reference \
  --include protein \
  --filename bacteria_genomes.zip
```

Option 3: EnsemblBacteria/Plants
---------------------------------

Plants (for plastid-targeted DAH7PS):
- EnsemblPlants: http://plants.ensembl.org/
- Phytozome: https://phytozome-next.jgi.doe.gov/

Bacteria:
- EnsemblBacteria: http://bacteria.ensembl.org/

Option 4: Use existing data
----------------------------

If you have your own FASTA files:
```bash
# Copy to raw data directory
cp /path/to/your/sequences.faa data/raw/my_sequences.faa

# Or create symbolic link
ln -s /path/to/your/sequences.faa data/raw/my_sequences.faa
```

Then update workflow/config.yaml:
```yaml
databases:
  uniprot_ref_proteomes: "data/raw/my_sequences.faa"
```

============================================================
Taxonomic Sampling Strategy
============================================================

For comprehensive DAH7PS evolution analysis:

Bacteria (Priority: HIGH)
- Proteobacteria (many paralogs: aroF, aroG, aroH)
- Firmicutes (diverse regulation)
- Actinobacteria (industrial importance)
- Cyanobacteria (photosynthetic)
- Bacteroidetes
- Spirochaetes

Archaea (Priority: MEDIUM)
- Euryarchaeota (Type II DAH7PS)
- Crenarchaeota
- Thaumarchaeota

Eukaryota - Plants (Priority: HIGH)
- Arabidopsis thaliana (model)
- Oryza sativa (rice)
- Zea mays (maize)
- Other angiosperms

Eukaryota - Fungi (Priority: MEDIUM)
- Saccharomyces cerevisiae (model, bifunctional)
- Aspergillus species
- Other ascomycetes and basidiomycetes

Target: 50-200 representative sequences per DAH7PS class
Total sequences after clustering: ~500-1000

============================================================

EOF
}

# Main menu
case "${1:-}" in
    test|1)
        download_test
        ;;
    bacteria|2)
        download_bacteria
        ;;
    custom|4)
        show_custom
        ;;
    *)
        echo "Usage: $0 [option]"
        echo ""
        echo "Options:"
        echo "  test      - Download small test dataset (recommended to start)"
        echo "  bacteria  - Download bacterial reference proteomes"
        echo "  custom    - Show custom download instructions"
        echo ""
        echo "Example:"
        echo "  $0 test"
        echo ""
        show_custom
        ;;
esac

echo ""
echo "============================================================"
echo "Next steps:"
echo "  1. Check downloaded files in: $RAW_DIR"
echo "  2. Update workflow/config.yaml with sequence file path"
echo "  3. Run: snakemake --cores 8 sequence_collection"
echo "============================================================"
