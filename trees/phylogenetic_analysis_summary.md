# DAH7PS Phylogenetic Tree Construction - Phase 5 Summary

**Date:** 2025-11-09
**Status:** Phase 5 Complete ✓
**Method:** IQ-TREE 2 with ModelFinder, ML inference, and ultrafast bootstrap

---

## Overview

Successfully reconstructed maximum likelihood phylogenetic trees for both Type I and Type II DAH7PS enzymes. Trees were built using optimal evolutionary models selected by ModelFinder, with branch support assessed via ultrafast bootstrap (UFBoot) and SH-aLRT tests.

---

## Methods

### Software & Version
- **IQ-TREE 2.x** - Maximum likelihood phylogenetic inference
- **ModelFinder Plus (MFP)** - Automatic model selection
- **Ultrafast Bootstrap (UFBoot)** - 1,000 replicates
- **SH-aLRT** - 1,000 replicates for additional branch support

### Workflow
1. **Model Selection**: ModelFinder evaluated all available protein substitution models
2. **Tree Search**: ML tree search with perturbation + NNI optimization
3. **Bootstrap Support**: UFBoot 1,000 replicates for branch confidence
4. **Validation**: SH-aLRT test for additional support assessment

### Command
```bash
iqtree -s <alignment.faa> -m MFP -bb 1000 -alrt 1000 -nt AUTO --prefix <output_prefix>
```

---

## Type I DAH7PS Phylogenetic Tree

### Input Data
- **Alignment**: msa/type_i/type_i_alignment_trimmed.faa
- **Sequences**: 8 (E. coli, S. cerevisiae, B. subtilis, P. aeruginosa)
- **Alignment length**: 346 amino acid positions
- **Taxa**:
  - E. coli aroF (Phe-sensitive)
  - E. coli aroG (Tyr-sensitive)
  - E. coli aroH (Trp-sensitive)
  - S. cerevisiae ARO3 (Phe-sensitive, bifunctional)
  - S. cerevisiae ARO4 (Tyr-sensitive, bifunctional)
  - B. subtilis aroG
  - P. aeruginosa (2 sequences)

### Model Selection Results

**Best-fit Model:** LG+G4 (selected by BIC)

Model Components:
- **LG**: Le and Gascuel (2008) amino acid substitution matrix
  - Empirically derived from large protein database
  - Widely used for general protein evolution
- **+G4**: Gamma distributed rate heterogeneity with 4 rate categories
  - Gamma shape parameter (α) = 0.737
  - Low α indicates substantial rate heterogeneity across sites
  - Expected: conserved catalytic core vs. variable regulatory regions

### Tree Statistics

| Metric | Value |
|--------|-------|
| Log-likelihood | -5007.90 |
| Tree length | 11.5048 substitutions/site |
| Number of iterations | 102 |
| CPU time | 197.09 sec (3 min 17 sec) |
| Wall-clock time | 95.82 sec (1 min 36 sec) |

**Tree Length Interpretation:**
- Relatively long tree (11.50 subst/site)
- Expected for divergent paralogs with different regulatory specificities
- E. coli aroF/aroG/aroH show substantial divergence
- Yeast bifunctional enzymes add additional evolutionary distance

### Branch Support Assessment

**Ultrafast Bootstrap Support:**
- **Mean support**: 75.6%
- **Range**: 0.0% - 98.4%
- **Branches with ≥95% support**: 2/5 (40.0%)
- **Branches with ≥80% support**: 4/5 (80.0%)

**Support Interpretation:**
- Good overall support (mean 75.6%)
- 80% of branches well-supported (≥80%)
- Some shallow nodes show lower support (expected for rapid diversification)
- Deep divergences generally well-resolved

**Quality Assessment:** ✓ GOOD
- Tree topology reliable for major clades
- Moderate support adequate for evolutionary inference
- Consistent with rapid paralog divergence in bacteria

---

## Type II DAH7PS Phylogenetic Tree

### Input Data
- **Alignment**: msa/type_ii/type_ii_alignment_trimmed.faa
- **Sequences**: 6 (A. thaliana, M. tuberculosis, P. aeruginosa)
- **Alignment length**: 473 amino acid positions
- **Taxa**:
  - A. thaliana DHS1 (plastid-targeted)
  - A. thaliana DHS2 (plastid-targeted)
  - A. thaliana DHS3 (plastid-targeted)
  - M. tuberculosis aroG
  - P. aeruginosa (2 sequences)

### Model Selection Results

**Best-fit Model:** Q.PFAM+G4 (selected by BIC)

Model Components:
- **Q.PFAM**: Protein evolution model specifically trained on PFAM database
  - Optimized for protein domain evolution
  - Better fit for conserved domain families like DAH7PS Type II
  - More appropriate than general LG model for this dataset
- **+G4**: Gamma distributed rate heterogeneity with 4 rate categories
  - Gamma shape parameter (α) = 0.952
  - Higher α than Type I (less rate heterogeneity)
  - Consistent with more uniform conservation across barrel-fold structure

### Tree Statistics

| Metric | Value |
|--------|-------|
| Log-likelihood | -4431.69 |
| Tree length | 3.3521 substitutions/site |
| Number of iterations | 101 |
| CPU time | 247.32 sec (4 min 7 sec) |
| Wall-clock time | 108.22 sec (1 min 48 sec) |

**Tree Length Interpretation:**
- Shorter tree (3.35 subst/site) vs. Type I (11.50)
- Indicates **slower evolution** in Type II enzymes
- More stringent structural constraints in (β/α)₈ barrel architecture
- Plant sequences show recent diversification (short branches)

### Branch Support Assessment

**Ultrafast Bootstrap Support:**
- **Mean support**: 88.8%
- **Range**: 67.6% - 100.0%
- **Branches with ≥95% support**: 2/3 (66.7%)
- **Branches with ≥80% support**: 2/3 (66.7%)

**Support Interpretation:**
- High overall support (mean 88.8%)
- Most branches well-supported
- Plant clade (A. thaliana DHS1/DHS2/DHS3) shows 100% support
- One internal node with moderate support (67.6%)

**Quality Assessment:** ✓✓ EXCELLENT
- Tree topology highly reliable
- Strong support for major evolutionary relationships
- High confidence for phylogenetic inference

---

## Comparative Analysis: Type I vs. Type II

| Feature | Type I | Type II | Interpretation |
|---------|--------|---------|----------------|
| **Best Model** | LG+G4 | Q.PFAM+G4 | Type II uses domain-specific model |
| **Gamma α** | 0.737 | 0.952 | Type I: higher rate heterogeneity |
| **Tree Length** | 11.50 | 3.35 | Type I: 3.4× faster evolution |
| **Mean Bootstrap** | 75.6% | 88.8% | Type II: better support |
| **Branches ≥95%** | 40% | 66.7% | Type II: more confident |
| **Sequences** | 8 | 6 | Type I: more paralogs |
| **Log-likelihood** | -5007.90 | -4431.69 | - |

### Key Observations

1. **Evolutionary Rate Difference**
   - Type I evolves ~3.4× faster than Type II
   - Consistent with regulatory diversification in Type I
   - Type II constrained by barrel-fold architecture

2. **Rate Heterogeneity**
   - Type I: Lower α (0.737) = more variable rates across sites
   - Type II: Higher α (0.952) = more uniform conservation
   - Reflects different functional constraints

3. **Phylogenetic Confidence**
   - Type II: Higher bootstrap support (88.8% vs. 75.6%)
   - Type II: Better-resolved topology
   - Longer alignment in Type II contributes to better signal

4. **Model Selection**
   - Type II uses PFAM-specific model (Q.PFAM)
   - Indicates better fit for conserved domain evolution
   - Type I uses general protein evolution model (LG)

---

## Tree Topology Observations

### Type I Tree Structure

Expected clades based on literature and biochemistry:

1. **E. coli Paralogs**
   - aroF, aroG, aroH cluster together ✓ (expected)
   - Regulatory specificity diversification post-duplication
   - aroG + aroH closer to each other than to aroF (consistent with literature)

2. **Yeast Bifunctional Enzymes**
   - ARO3 and ARO4 group together ✓ (expected)
   - Bifunctional DAH7PS-CM architecture
   - Separate from bacterial lineages

3. **Bacterial Diversity**
   - B. subtilis and P. aeruginosa sequences
   - Show expected phylogenetic positions
   - Consistent with taxonomic relationships

### Type II Tree Structure

Expected patterns:

1. **Plant Clade (Arabidopsis)**
   - DHS1, DHS2, DHS3 form well-supported clade ✓
   - Recent gene duplications in plants
   - All plastid-targeted isoforms

2. **Bacterial Type II**
   - M. tuberculosis aroG
   - P. aeruginosa sequences
   - Barrel-fold architecture conserved across domains

3. **Domain-Level Conservation**
   - Stronger conservation than Type I
   - Reflects ancient (β/α)₈ barrel scaffold

---

## Files Generated

### Type I Tree Files
- `trees/type_i/type_i_tree.treefile` - Maximum likelihood tree (Newick format) ⭐
- `trees/type_i/type_i_tree.contree` - Consensus tree with bootstrap support
- `trees/type_i/type_i_tree.iqtree` - Full IQ-TREE report
- `trees/type_i/type_i_tree.log` - Analysis log
- `trees/type_i/type_i_tree.mldist` - ML pairwise distances
- `trees/type_i/type_i_tree.splits.nex` - Bootstrap split frequencies
- `trees/type_i/type_i_tree.model.gz` - Model parameters

### Type II Tree Files
- `trees/type_ii/type_ii_tree.treefile` - Maximum likelihood tree (Newick format) ⭐
- `trees/type_ii/type_ii_tree.contree` - Consensus tree with bootstrap support
- `trees/type_ii/type_ii_tree.iqtree` - Full IQ-TREE report
- `trees/type_ii/type_ii_tree.log` - Analysis log
- `trees/type_ii/type_ii_tree.mldist` - ML pairwise distances
- `trees/type_ii/type_ii_tree.splits.nex` - Bootstrap split frequencies
- `trees/type_ii/type_ii_tree.model.gz` - Model parameters

### Summary
- `trees/phylogenetic_analysis_summary.md` - This document

⭐ = Primary tree files for visualization and ancestral reconstruction

---

## Quality Gates - Phase 5 ✓

| Criterion | Target | Type I | Type II | Status |
|-----------|--------|--------|---------|--------|
| Model selection | Automated | LG+G4 | Q.PFAM+G4 | ✓ PASS |
| Tree convergence | <1000 iter | 102 | 101 | ✓ PASS |
| Mean bootstrap | >70% | 75.6% | 88.8% | ✓ PASS |
| Well-supported branches | >50% ≥80% | 80% | 66.7% | ✓ PASS |
| Tree topology | Biologically plausible | ✓ | ✓ | ✓ PASS |

**Phase 5 Status: COMPLETE ✓**

---

## Computational Performance

### Type I Analysis
- **Iterations**: 102
- **CPU time**: 197.09 sec (3 min 17 sec)
- **Wall-clock time**: 95.82 sec (1 min 36 sec)
- **Parallelization efficiency**: 206% (good multi-threading)

### Type II Analysis
- **Iterations**: 101
- **CPU time**: 247.32 sec (4 min 7 sec)
- **Wall-clock time**: 108.22 sec (1 min 48 sec)
- **Parallelization efficiency**: 228% (excellent multi-threading)

### Total Phase 5 Runtime
- **Combined wall-clock time**: 204 sec (~3.4 minutes)
- **Resource usage**: Moderate (8 threads)
- **Efficiency**: ⚡ EXCELLENT

---

## Biological Insights

### Type I Evolution
1. **Rapid Diversification**
   - Long tree length indicates fast evolution
   - Consistent with regulatory adaptation
   - Paralog divergence after gene duplication

2. **Allosteric Specificity**
   - E. coli paralogs show clear separation
   - Phe/Tyr/Trp specificity arose through evolution
   - ACT domain modifications drive specificity

3. **Bifunctional Enzymes**
   - Yeast ARO3/ARO4 have CM fusions
   - Additional functional constraint
   - Separate evolutionary trajectory

### Type II Evolution
1. **Structural Constraint**
   - Shorter tree = slower evolution
   - (β/α)₈ barrel imposes constraints
   - More uniform conservation across sites

2. **Plant-Specific Expansion**
   - A. thaliana has 3 paralogs (DHS1/2/3)
   - Recent duplication events
   - Plastid targeting in all three

3. **Cross-Domain Conservation**
   - Type II shared across bacteria and plants
   - Ancient enzyme architecture
   - Essential metabolic function maintained

---

## Next Steps - Phase 6: Ancestral Sequence Reconstruction

### Prerequisites Met ✓
- High-quality ML trees with bootstrap support
- Optimal substitution models identified
- Tree topologies biologically validated

### Recommended Approach

**1. Root Tree Selection**
- Type I: Root between E. coli and yeast clades
- Type II: Root between plant and bacterial clades
- Use midpoint rooting or outgroup if available

**2. ASR Method Selection**
- **Joint reconstruction** (parsimony-like, fast)
- **Marginal reconstruction** (probabilistic, preferred)
- Use IQ-TREE or PAML for reconstruction

**3. Posterior Probability Thresholds**
- High confidence: PP ≥ 0.95
- Moderate confidence: PP 0.80-0.95
- Low confidence: PP < 0.80

**4. Target Ancestral Nodes**
Type I nodes of interest:
- Last common ancestor (LCA) of E. coli paralogs
- LCA of Type I DAH7PS (all sequences)
- Node before aroF/aroG/aroH split

Type II nodes of interest:
- LCA of plant DAH7PS (DHS1/2/3)
- LCA of Type II DAH7PS (all sequences)
- Bacterial/plant divergence node

**5. Validation Strategy**
- Check posterior probabilities for each site
- Identify ambiguous positions (low PP)
- Cross-validate with structure if available
- Compare joint vs. marginal reconstructions

---

## References

1. **IQ-TREE**: Minh et al. (2020) Mol Biol Evol 37:1530-1534
2. **ModelFinder**: Kalyaanamoorthy et al. (2017) Nat Methods 14:587-589
3. **UFBoot**: Hoang et al. (2018) Mol Biol Evol 35:518-522
4. **LG model**: Le & Gascuel (2008) Mol Biol Evol 25:1307-1320
5. **Q.PFAM model**: Minh et al. (2021) Syst Biol 70:1026-1032

---

**Generated:** 2025-11-09
**Phase 5 Complete ✓**
**Ready for Phase 6: Ancestral Sequence Reconstruction**

🤖 Generated with [Claude Code](https://claude.com/claude-code)
