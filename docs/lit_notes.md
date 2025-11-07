# Literature Notes: DAH7PS Allostery & Evolution

## Overview

This document contains summaries and key findings from papers relevant to DAH7PS evolution, allosteric regulation, and ancestral sequence reconstruction.

---

## Key Primary References

### PMID: 36041629

**Status:** To be reviewed and summarized

**Citation:** [To be filled after retrieval]

**Key Topics:**
- [To be filled]

**Main Findings:**
- [To be filled]

**Relevant Methods:**
- [To be filled]

**Implications for This Project:**
- [To be filled]

**Notes:**
- [To be filled]

---

### PMID: 34813062

**Status:** To be reviewed and summarized

**Citation:** [To be filled after retrieval]

**Key Topics:**
- [To be filled]

**Main Findings:**
- [To be filled]

**Relevant Methods:**
- [To be filled]

**Implications for This Project:**
- [To be filled]

**Notes:**
- [To be filled]

---

## DAH7PS Biochemistry & Structure

### General Background

**Enzyme Function:**
- DAH7PS (3-deoxy-D-arabinoheptulosonate 7-phosphate synthase)
- EC 2.5.1.54
- Catalyzes: PEP + E4P → DAH7P + Pi
- First committed step of shikimate pathway
- Essential for aromatic amino acid biosynthesis

**Phylogenetic Distribution:**
- Present in bacteria, archaea, plants, fungi, apicomplexans
- Absent in animals (pathway lost, making it a drug target)

**Structural Classes:**

1. **Type I (α-barrel)**
   - Subtype Iα: Most bacterial DAH7PS
   - Subtype Iβ: Divergent regulation
   - (β/α)₈-barrel (TIM barrel) fold
   - Often forms dimers or tetramers

2. **Type II (β-barrel)**
   - Predominantly archaeal, some bacterial
   - Different fold architecture
   - Distinct oligomerization

---

## Allosteric Regulation Mechanisms

### ACT Domain Regulation

**ACT Domain Properties:**
- Pfam: PF01842
- Small molecule-binding domain (~60-80 residues)
- Named after: Aspartate kinase, Chorismate mutase, TyrA
- Often present as tandem repeats in DAH7PS
- Binds aromatic amino acids (Tyr, Phe, Trp) in interdomain cleft

**Regulatory Modes:**

1. **Single effector inhibition**
   - E. coli aroG (Phe-sensitive)
   - E. coli aroF (Tyr-sensitive)
   - E. coli aroH (Trp-sensitive)

2. **Synergistic inhibition**
   - Requires two effectors simultaneously
   - Example: Some organisms require Tyr + Phe
   - Potentially more sophisticated metabolic control

**Structural Basis:**
- Effector binding → conformational change in ACT domain
- Transmitted to catalytic core via interdomain interface
- Can affect substrate binding or catalytic turnover

### Chorismate Mutase Fusion

**CM Domain Properties:**
- Pfam: AroQ_class_II (CM_2), or other CM families
- Bifunctional enzyme: DAH7PS-CM
- Catalyzes: Chorismate → Prephenate
- Metabolic channeling potential

**Regulatory Features:**
- Allosteric regulation by pathway intermediates (chorismate, prephenate)
- Feedback from aromatic amino acids
- Coordinated regulation of pathway flux

**Distribution:**
- Common in certain bacterial lineages
- Fungi (yeast ARO3, ARO4)
- Some plants

---

## Evolutionary Questions from Literature

### Domain Recruitment

**Questions:**
- How many independent ACT domain fusions?
- Timing of CM domain fusion events?
- Ancestral state: regulated or unregulated?

**Hypotheses:**
- Multiple independent recruitments of ACT domains
- Regulatory domains co-opted from other metabolic enzymes
- Positive selection for feedback inhibition in amino acid-rich environments

### Effector Specificity Evolution

**Questions:**
- What substitutions determine Tyr vs Phe vs Trp specificity?
- Can we predict effector from sequence?
- How many switches between effector types?

**Known Specificity Determinants:**
- Residues lining effector-binding pocket in ACT domain
- Shape and charge complementarity
- Examples from mutagenesis studies (to be compiled)

### Synergistic Inhibition

**Questions:**
- Did synergy evolve from single-effector ancestors?
- Structural basis of cooperativity?
- Adaptive advantage of synergistic control?

**Mechanistic Hypotheses:**
- Two effector-binding sites in tandem ACT domains
- Positive cooperativity between sites
- Fine-tuned metabolic regulation

---

## Ancestral Sequence Reconstruction (ASR) - Methods

### General ASR Principles

**Approaches:**
- Maximum likelihood (ML): IQ-TREE, PAML
- Bayesian: MrBayes, BEAST (with ancestral state sampling)
- Parsimony (less common, less accurate)

**Uncertainty:**
- Site-wise posterior probabilities
- Alternative ancestral states
- Topology uncertainty (bootstrap/Bayesian tree sampling)
- Model uncertainty

**Validation:**
- Structural plausibility (fold prediction, stability)
- Functional constraints (catalytic residues)
- Experimental resurrection (when feasible)

### DAH7PS-Specific Considerations

**Challenges:**
- Long evolutionary timescales (billions of years for deep nodes)
- Modular domain architecture (partition-based ASR)
- Recombination/domain shuffling (requires careful tree interpretation)

**Opportunities:**
- Well-characterized extant enzymes for validation
- Crystal structures for mapping
- Functional assays for resurrected ancestors (kinetics, effector response)

---

## Comparative Studies

### Related ASR Projects (Examples)

**Protein families studied:**
- Steroid receptors (Thornton lab)
- β-lactamase (antibiotic resistance evolution)
- Fluorescent proteins (ancestral GFP)
- Metabolic enzymes (various)

**Lessons Learned:**
- Importance of alignment quality
- Site-heterogeneous models improve accuracy
- Deep ancestors often generalists vs specialized descendants
- Allostery can be gained/lost multiple times

---

## Structural Data

### Available PDB Structures

**Type Iα Examples:**
- E. coli aroG (Phe-sensitive) with ACT domains
- E. coli aroF (Tyr-sensitive)
- T. maritima DAH7PS
- [To be compiled: PDB IDs and citations]

**Type II Examples:**
- Archaeal DAH7PS structures
- [To be compiled]

**CM Fusion Examples:**
- S. cerevisiae ARO3, ARO4
- [To be compiled]

**Key Structural Features to Map:**
- Active site residues (PEP and E4P binding)
- Metal coordination (if present)
- Oligomerization interfaces (dimer, tetramer)
- ACT domain effector-binding pockets
- Core-ACT interdomain interface
- CM domain active site (in fusions)

---

## Functional Assays & Kinetics

### Biochemical Characterization

**Standard Assays:**
- DAH7PS activity: coupled enzyme assay (lactate dehydrogenase, pyruvate kinase)
- Product detection: thiobarbituric acid assay
- Kinetic parameters: Km (PEP), Km (E4P), kcat, Ki (effectors)

**Allosteric Characterization:**
- Dose-response curves for effectors
- Synergy testing (Hill coefficients, combination inhibition)
- Structural probing (circular dichroism, thermal shift assays)

**Data Compilation Needed:**
- Literature survey for kinetic parameters
- Create database: organism, enzyme, effector, Ki, mode
- Link to sequence IDs in our dataset

---

## Taxonomic Distribution Patterns

### Bacteria

**Proteobacteria:**
- Often multiple paralogs (aroF, aroG, aroH)
- ACT domain regulation common
- Diverse effector specificities

**Firmicutes:**
- [To be characterized]

**Actinobacteria:**
- [To be characterized]

**Cyanobacteria:**
- [To be characterized]

### Archaea

**Type II Predominance:**
- Most archaea use Type II
- Regulatory mechanisms less characterized
- Potential for novel allosteric modes

### Eukaryotes

**Plants:**
- Plastid-localized DAH7PS
- Transit peptides (to be masked in alignments)
- Often CM fusions

**Fungi:**
- ARO3, ARO4 (S. cerevisiae)
- Bifunctional DAH7PS-CM
- Tyr/Phe feedback inhibition

---

## Bioinformatics Resources

### Databases

- **Pfam:** PF00793, PF01474, PF01842
- **UniProt/SwissProt:** Curated DAH7PS entries
- **KEGG:** K01626 (DAH7PS), K01735 (kdsA/KDO8PS - to exclude)
- **InterPro:** Integrated domain annotations
- **PDB:** Structural models
- **BRENDA:** Kinetic and functional data (EC 2.5.1.54)

### Tools

- **Sequence Search:** HMMER, BLAST, MMseqs2
- **MSA:** MAFFT, MUSCLE, Clustal Omega
- **Trimming:** trimAl, BMGE, GUIDANCE2
- **Phylogenetics:** IQ-TREE, RAxML-NG, PhyML
- **ASR:** IQ-TREE (--asr), PAML (codeml), FastML, GRASP
- **Trait Evolution:** PastML, phytools, BayesTraits, Mesquite
- **Localization:** SignalP, TargetP, ChloroP
- **Structure:** PyMOL, ChimeraX, AlphaFold2, ESMFold

---

## Project-Specific Notes

### Decision Log

**2025-11-07:**
- Project initiated
- Scope defined: Bacteria, Archaea, Plantae, Fungi; Types Iα, Iβ, II
- Focus on ACT and CM-mediated regulation
- Partition-based alignment strategy chosen

### Open Questions

1. Should we include Type II in the same analysis or separate?
   - **Decision:** Separate core alignments, but include in global trait analysis

2. How to handle fusion proteins (DAH7PS-CM)?
   - **Decision:** Partition-based ASR; reconstruct core and CM separately

3. Outgroup for rooting?
   - **Decision:** Consider KDO8PS; test midpoint rooting; use mixture models

4. Experimental validation targets?
   - **Decision:** Defer to Phase 10; prioritize ancestors at key regulatory transitions

### Data Availability

- **UniProt Reference Proteomes:** Freely available
- **Pfam HMMs:** Freely available
- **PDB Structures:** Freely available
- **SignalP 6 / TargetP 2:** May require academic license or DTU Health Tech webserver
- **Commercial tools:** None planned

---

## To Do: Literature Review Tasks

- [ ] Retrieve and read PMID: 36041629
- [ ] Retrieve and read PMID: 34813062
- [ ] Compile list of all available DAH7PS crystal structures (PDB IDs)
- [ ] Survey kinetic parameters from BRENDA and primary literature
- [ ] Create effector specificity table from literature
- [ ] Review ASR methodology papers (Thornton, Gaucher, etc.)
- [ ] Check for recent reviews on shikimate pathway evolution
- [ ] Identify any existing DAH7PS phylogenies for comparison

---

## References Template

### Format

```
Author(s). (Year). Title. Journal Volume(Issue): Pages. PMID: XXXXXXXX. DOI: XX.XXXX/XXXXXX
```

### Papers to Add

1. [PMID: 36041629]
2. [PMID: 34813062]
3. [Classic DAH7PS biochemistry papers]
4. [Structural biology papers]
5. [ASR methodology papers]
6. [Shikimate pathway reviews]

---

**Last Updated:** 2025-11-07
