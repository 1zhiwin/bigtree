# Sequence Data Collection Guide

## Overview

For the DAH7PS project, you need protein sequences from diverse organisms to:
1. Search for DAH7PS candidates using HMMER
2. Analyze taxonomic distribution
3. Study regulatory domain evolution
4. Reconstruct ancestral sequences

---

## Quick Start (Recommended for Testing)

### Option 1: Small Test Dataset (~100 MB, 5 minutes)

Download a curated set of model organisms:

```bash
cd ~/bigtree
bash scripts/download_sequences.sh test
```

This downloads:
- **E. coli K-12** - Has aroF, aroG, aroH (Tyr/Phe/Trp specific)
- **B. subtilis** - Firmicutes representative
- **P. aeruginosa** - Gram-negative pathogen
- **M. tuberculosis** - Actinobacteria, medical importance
- **S. cerevisiae** - Fungi, bifunctional DAH7PS-CM
- **A. thaliana** - Plant, plastid-targeted enzymes

**Result:** `data/raw/test_proteomes.faa`

**Pros:**
- Fast download
- Known DAH7PS sequences
- Good for testing workflow
- Covers main enzyme classes

**Cons:**
- Limited taxonomic diversity
- Not suitable for publication-quality phylogeny

---

## Comprehensive Datasets

### Option 2: UniProt Reference Proteomes (Recommended for Publication)

**What:** High-quality, manually curated proteomes representing species diversity

**Size:**
- Bacteria: ~5-10 GB (200+ species)
- Archaea: ~500 MB (50+ species)
- Eukaryotes: ~2-5 GB (100+ species)

**Download Methods:**

#### A. Browser Download (Selective)

1. Visit: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/

2. Navigate directories:
   - `Bacteria/` - All bacterial proteomes
   - `Archaea/` - Archaeal proteomes
   - `Eukaryota/` - Plants, fungi, etc.

3. Download specific organism folders

4. Files are in format: `UP000XXXXXX_TAXID.fasta.gz`

#### B. Command-Line Download (Bulk)

**Download all bacteria (careful: ~10 GB!):**

```bash
cd data/raw
mkdir -p proteomes/bacteria

# Download with wget (recursive)
wget -r -np -nH --cut-dirs=6 -R "index.html*" \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/

# Or use rsync (better for large downloads)
rsync -avz --progress \
  rsync://ftp.uniprot.org/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/ \
  proteomes/bacteria/
```

**Download specific phyla:**

```bash
# Proteobacteria only
rsync -avz --progress \
  --include='*/' --include='*.fasta.gz' --exclude='*' \
  rsync://ftp.uniprot.org/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/ \
  proteomes/proteobacteria/ \
  | grep -i proteobacteria
```

#### C. Selected Representative Organisms

Create a list of target organisms and download:

```bash
cd data/raw

# Create list of UniProt IDs
cat > organism_list.txt << 'EOF'
UP000000625  # E. coli K-12
UP000001570  # Bacillus subtilis
UP000002438  # Pseudomonas aeruginosa
UP000001584  # Mycobacterium tuberculosis
UP000000579  # Salmonella typhimurium
UP000000535  # Streptococcus pneumoniae
UP000002311  # Saccharomyces cerevisiae
UP000006548  # Arabidopsis thaliana
UP000059680  # Oryza sativa (rice)
EOF

# Download each
while read id rest; do
  [[ "$id" =~ ^# ]] && continue
  echo "Downloading $id..."
  wget -r -np -nH --cut-dirs=6 -A "*.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/" \
    --accept-regex ".*$id.*"
done < organism_list.txt
```

---

### Option 3: NCBI RefSeq

**What:** All annotated genomes, very comprehensive

**Pros:**
- More organisms than UniProt
- Includes draft genomes
- Good for gap-filling rare taxa

**Cons:**
- Variable quality
- More redundancy
- Larger files

**Download using NCBI Datasets:**

```bash
# Install if not already in environment
conda install -c conda-forge ncbi-datasets-cli

cd data/raw

# Download by taxonomy
datasets download genome taxon "Proteobacteria" \
  --reference \
  --include protein \
  --filename proteobacteria.zip

# Extract
unzip proteobacteria.zip

# Proteins are in: ncbi_dataset/data/GCF_*/protein.faa
```

**Download specific organisms:**

```bash
# E. coli K-12
datasets download genome accession GCF_000005845.2 \
  --include protein \
  --filename ecoli_k12.zip
```

---

### Option 4: EnsemblBacteria / EnsemblPlants

**For Plants (plastid-targeted enzymes):**

1. Visit: http://plants.ensembl.org/
2. Select organism
3. Download → Protein sequences (FASTA)

**For Bacteria:**

1. Visit: http://bacteria.ensembl.org/
2. Browse genomes
3. Download protein FASTA

---

## Recommended Sampling Strategy

### Minimal Dataset (for testing)
- 5-10 model organisms
- Covers main classes (Iα, Iβ, II)
- ~100 MB
- **Time to analyze:** 1-2 hours

### Standard Dataset (for initial analysis)
- 50-100 representative organisms
- Major phyla represented
- ~1-2 GB
- **Time to analyze:** 1-2 days

### Comprehensive Dataset (for publication)
- 200-500 organisms
- Broad taxonomic sampling
- Multiple strains of key species
- ~5-10 GB
- **Time to analyze:** 3-7 days

---

## Taxonomic Coverage Checklist

For comprehensive DAH7PS evolution study:

### Bacteria (Priority groups)

**Proteobacteria** (aroF/aroG/aroH paralogs)
- [ ] Escherichia coli (model)
- [ ] Salmonella enterica
- [ ] Pseudomonas species
- [ ] Vibrio species
- [ ] Neisseria species

**Firmicutes** (diverse regulation)
- [ ] Bacillus subtilis (model)
- [ ] Staphylococcus aureus
- [ ] Lactobacillus species
- [ ] Clostridium species

**Actinobacteria**
- [ ] Mycobacterium tuberculosis
- [ ] Streptomyces species
- [ ] Corynebacterium species

**Cyanobacteria** (photosynthetic)
- [ ] Synechocystis sp.
- [ ] Prochlorococcus marinus

**Others:**
- [ ] Bacteroidetes representatives
- [ ] Spirochaetes
- [ ] Aquificae (deep-branching)

### Archaea

**Euryarchaeota** (Type II DAH7PS)
- [ ] Methanocaldococcus jannaschii
- [ ] Pyrococcus furiosus
- [ ] Thermococcus kodakarensis

**Crenarchaeota**
- [ ] Sulfolobus solfataricus
- [ ] Pyrobaculum aerophilum

### Eukaryotes - Plants

**Model organisms:**
- [ ] Arabidopsis thaliana
- [ ] Oryza sativa (rice)
- [ ] Zea mays (maize)
- [ ] Solanum lycopersicum (tomato)

**Diversity:**
- [ ] Bryophytes (moss)
- [ ] Gymnosperms
- [ ] Basal angiosperms

### Eukaryotes - Fungi

**Ascomycota:**
- [ ] Saccharomyces cerevisiae (model, ARO3/ARO4)
- [ ] Neurospora crassa
- [ ] Aspergillus species

**Basidiomycota:**
- [ ] Ustilago maydis
- [ ] Cryptococcus neoformans

---

## Data Processing After Download

### 1. Combine Files

If you downloaded multiple files:

```bash
cd data/raw

# Decompress all
gunzip proteomes/**/*.gz

# Combine into single FASTA
cat proteomes/**/*.fasta > all_proteomes.faa

# Or keep compressed and use zcat
zcat proteomes/**/*.fasta.gz > all_proteomes.faa
```

### 2. Check File Quality

```bash
# Count sequences
grep -c '^>' all_proteomes.faa

# Check for issues
grep -c '^$' all_proteomes.faa  # Empty lines
head -100 all_proteomes.faa     # Visual inspection

# Get statistics
seqkit stats all_proteomes.faa
```

### 3. Remove Duplicates (if needed)

```bash
# Use seqkit to remove exact duplicates
seqkit rmdup -s all_proteomes.faa > all_proteomes_unique.faa

# Check reduction
echo "Before: $(grep -c '^>' all_proteomes.faa)"
echo "After: $(grep -c '^>' all_proteomes_unique.faa)"
```

### 4. Update Configuration

Edit `workflow/config.yaml`:

```yaml
databases:
  pfam_hmm: "data/raw/Pfam-A.hmm"
  uniprot_ref_proteomes: "data/raw/all_proteomes.faa"  # ← Update this
```

---

## Example: Step-by-Step Test

Let's do a complete test run:

```bash
# 1. Download test dataset
cd ~/bigtree
bash scripts/download_sequences.sh test

# 2. Verify download
ls -lh data/raw/test_proteomes.faa
grep -c '^>' data/raw/test_proteomes.faa

# 3. Test HMMER search (manual)
hmmsearch --domtblout data/raw/test_search.out \
  --cut_ga \
  data/raw/Pfam-A.hmm \
  data/raw/test_proteomes.faa

# 4. Check results
grep -c 'DAHP_synth' data/raw/test_search.out

# 5. Update config
sed -i 's|uniprot_reference_proteomes.faa|test_proteomes.faa|g' workflow/config.yaml

# 6. Run workflow
snakemake --cores 4 sequence_collection
```

---

## Troubleshooting

### Download fails
- Check internet connection
- Try rsync instead of wget
- Download smaller chunks
- Use UniProt FTP mirror: ftp.ebi.ac.uk

### Files too large
- Start with test dataset
- Download specific phyla only
- Use UniRef90 instead of full proteomes

### Disk space issues
```bash
# Check available space
df -h data/raw

# Compress old files
gzip data/raw/proteomes/*.faa

# Remove intermediate files
rm -rf data/raw/proteomes/temp/
```

---

## Data Quality Checklist

Before running analysis:
- [ ] Files are in FASTA format
- [ ] Headers are properly formatted (>ID description)
- [ ] No empty sequences
- [ ] File size is reasonable (not corrupted)
- [ ] Taxonomic diversity adequate for goals
- [ ] Known positive controls present (E. coli, yeast)

---

## Next Steps

After downloading sequences:

1. **Verify data:** Check file sizes and sequence counts
2. **Update config:** Edit `workflow/config.yaml`
3. **Test search:** Run manual HMMER search
4. **Start workflow:** `snakemake --cores 8 sequence_collection`

---

## Resources

**UniProt:**
- Reference Proteomes: https://www.uniprot.org/proteomes
- FTP: ftp://ftp.uniprot.org/pub/databases/uniprot/

**NCBI:**
- Genome Browser: https://www.ncbi.nlm.nih.gov/genome/
- Datasets: https://www.ncbi.nlm.nih.gov/datasets/

**Ensembl:**
- Plants: http://plants.ensembl.org/
- Bacteria: http://bacteria.ensembl.org/

**Documentation:**
- UniProt manual: https://www.uniprot.org/help/
- NCBI datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/
