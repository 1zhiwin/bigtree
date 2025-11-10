# Main Manuscript Figure Descriptions

## Figure 1: Phylogenetic Analysis and Evolutionary Rates

**Panel A: Type I Phylogenetic Tree**
- Maximum likelihood tree of 8 Type I DAH7PS sequences
- Branch labels show UFBoot/SH-aLRT support values
- Tip labels colored by regulatory specificity:
  - Red: Phe-sensitive (aroF, ARO3)
  - Blue: Tyr-sensitive (aroG, ARO4)
  - Green: Trp-sensitive (aroH)
  - Gray: Unknown
- Scale bar shows substitutions per site
- Tree rooted at midpoint

**Panel B: Type II Phylogenetic Tree**
- Maximum likelihood tree of 6 Type II DAH7PS sequences
- Branch labels show UFBoot/SH-aLRT support values
- Tip labels colored by localization:
  - Green: Plastid-targeted (DHS1, DHS2, DHS3)
  - Gray: Cytoplasmic (bacterial)
- Scale bar shows substitutions per site
- Tree rooted at midpoint

**Panel C: Evolutionary Rate Comparison**
- Bar chart comparing tree lengths
- X-axis: Enzyme type (Type I, Type II)
- Y-axis: Tree length (substitutions per site)
- Type I: 11.50 ± SE
- Type II: 3.35 ± SE
- **P < 0.0001, 3.4-fold difference
- Error bars represent bootstrap variation

**Panel D: Model Parameters**
- Table showing best-fit substitution models:
  - Type I: LG+G4 (α = 0.737)
  - Type II: Q.PFAM+G4 (α = 0.952)
- Gamma shape parameter (α) comparison
- Lower α = more rate heterogeneity

**Figure Legend:**
Phylogenetic analysis reveals 3.4-fold faster evolution in Type I DAH7PS enzymes. (A) Type I tree (LG+G4 model) showing regulatory specificity diversity (Phe/Tyr/Trp sensitivities). (B) Type II tree (Q.PFAM+G4 model) showing compartmentalization in plants. (C) Tree length comparison demonstrates significantly faster Type I evolution (11.50 vs. 3.35 substitutions/site, p < 0.0001). (D) Gamma shape parameters indicate greater rate heterogeneity in Type I (α = 0.737) versus Type II (α = 0.952). UFBoot and SH-aLRT values shown at branches.

---

## Figure 2: Ancestral Sequence Reconstruction

**Panel A: Type I Ancestral Nodes**
- Phylogenetic tree with 6 ancestral nodes labeled (Node1-Node6)
- Nodes shown as filled circles, sized by mean posterior probability
- Node5 highlighted (deepest common ancestor, highest confidence)
- Color gradient: High PP (green) to Low PP (red)

**Panel B: Type II Ancestral Nodes**
- Phylogenetic tree with 4 ancestral nodes labeled (Node1-Node4)
- Nodes shown as filled circles, sized by mean posterior probability
- Node3 and Node4 highlighted (plant ancestors, exceptional confidence)
- Color gradient: High PP (green) to Low PP (red)

**Panel C: Reconstruction Confidence Distribution**
- Stacked bar charts for each ancestral node
- Segments: High confidence (PP ≥ 0.95, green)
           Medium confidence (PP 0.80-0.95, yellow)
           Low confidence (PP < 0.80, red)
- Separate panels for Type I (6 nodes) and Type II (4 nodes)
- Mean PP values labeled above each bar

**Panel D: Conservation vs. Reconstruction Confidence**
- Scatter plot: X-axis = sequence conservation (%)
                Y-axis = mean PP
- Points represent ancestral nodes
- Type I nodes: blue circles
- Type II nodes: red squares
- Regression line showing positive correlation
- R² and p-value displayed

**Figure Legend:**
Ancestral sequence reconstruction achieves high confidence, especially for Type II. (A) Type I tree showing 6 reconstructed ancestral nodes (circles), with Node5 representing the deepest common ancestor (mean PP = 0.8608). (B) Type II tree with 4 ancestral nodes, including exceptionally well-resolved plant ancestors Node3 and Node4 (PP > 0.95). (C) Reconstruction confidence distributions show Type II has more high-confidence sites (66.6%) than Type I (53.1%). (D) Higher sequence conservation predicts better reconstruction confidence (R² = 0.84, p < 0.001).

---

## Figure 3: Trait Evolution Analysis

**Panel A: Type I Trait Mapping**
- Phylogenetic tree with tip labels showing regulatory specificities
- Pie charts at tips: Phe (red), Tyr (blue), Trp (green), Unknown (gray)
- Ancestral nodes show reconstructed states (colored circles)
- Node 89.9/71 highlighted as Tyr-ancestral state
- Single trait change indicated with asterisk

**Panel B: Type II Trait Mapping**
- Phylogenetic tree with tip labels showing localization
- Squares at tips: Plastid (green), Cytoplasmic (gray)
- Ancestral nodes show reconstructed states
- Node 100/100 highlighted as plastid acquisition point
- Single trait change indicated with asterisk

**Panel C: Parsimony Score Comparison**
- Bar chart: X-axis = Enzyme type
            Y-axis = Parsimony score
- Both bars at height = 1 (minimal trait evolution)
- Inset table showing:
  - Tree length (substitutions/site)
  - Parsimony score
  - Ratio (sequence changes per trait change)
  - Type I: 11.50 subs/site ÷ 1 trait change = 11.50
  - Type II: 3.35 subs/site ÷ 1 trait change = 3.35

**Panel D: Trait State Transitions**
- Sankey diagram showing trait evolution
- Type I: Unknown → Tyr → (Phe, Trp)
- Type II: Cytoplasmic → Plastid
- Width of flows proportional to number of sequences
- Colors match trait categories

**Figure Legend:**
Conservative trait evolution despite extensive sequence divergence. (A) Type I trait reconstruction shows single transition to Tyr-sensitivity (node 89.9/71), with subsequent diversification to Phe and Trp. Parsimony score = 1. (B) Type II trait reconstruction reveals single plastid-targeting acquisition (node 100/100, maximum support). Parsimony score = 1. (C) Both enzyme types show minimal trait changes (score = 1) despite rapid sequence evolution, particularly in Type I (11.50 substitutions/site). (D) Trait transition diagrams illustrate single-step innovations in both lineages.

---

## Figure 4: Stability-Evolvability Trade-off

**Panel A: Overall Stability Comparison**
- Box plots: X-axis = Enzyme type/domain
            Y-axis = Instability index
- Categories: Type I full, Type II full,
              Type I catalytic, Type I ACT
- Horizontal line at II = 40 (stability threshold)
- Type I ACT domain significantly above threshold (II = 70.14)
- Statistical comparisons: *** p < 0.001

**Panel B: Domain-Specific Instability (E. coli aroF)**
- Split protein diagram showing catalytic (1-260) and ACT (270-356) domains
- Color-coded by instability: Blue (stable, II = 33.92) to Red (unstable, II = 70.14)
- Linker region (261-269) highlighted
- Instability index values labeled for each domain
- 2.1-fold difference indicated

**Panel C: Variability vs. Stability**
- Scatter plot: X-axis = Instability index
                Y-axis = Evolutionary variability
- Points represent sequence positions
- Color by domain: Blue (catalytic), Red (ACT)
- Regression line showing negative correlation
- R² = -0.65, p < 0.001
- Interpretation: Less stable → More variable

**Panel D: Mutation Stability Predictions**
- Pie chart of 133 trait-specific mutations
- Segments:
  - Neutral: 27.1% (gray)
  - Destabilizing: 19.5% (orange)
  - Strongly destabilizing: 3.0% (red)
  - Unknown: 50.4% (light gray)
- Total destabilizing: 22.5% (labeled)
- Examples of each category shown in table

**Figure Legend:**
ACT regulatory domains show stability-evolvability trade-off. (A) Overall stability comparison shows Type I more stable than Type II, but domain analysis reveals heterogeneity. (B) Domain-specific stability in *E. coli* aroF: catalytic domain stable (II = 33.92), ACT domain unstable (II = 70.14), representing 2.1-fold difference. (C) Correlation between instability and evolutionary variability (R² = -0.65) demonstrates that less stable domains evolve faster. (D) Analysis of 133 trait-specific mutations predicts 22.5% reduce stability, supporting regulatory evolution cost model. Instability index < 40 indicates stable protein.

---

## Figure 5: Domain Interface Analysis

**Panel A: Interface Residue Identification**
- Schematic of Type I architecture
- Catalytic domain (1-260): blue
- Linker (261-269): gray
- ACT domain (270-356): red
- Interface regions highlighted:
  - Catalytic C-terminus (250-270): yellow outline
  - ACT N-terminus (270-290): orange outline
- Overlap at 270 indicates linker is part of both interfaces

**Panel B: Interface Trait-Specificity Enrichment**
- Bar chart comparing trait-specific percentages
- Categories: Genome-wide (66.8%)
             Catalytic C-term interface (80.0%)
             ACT N-term interface (81.0%)
- Enrichment factor: 1.2× (indicated)
- Fisher's exact test: p < 0.01
- Error bars represent 95% confidence intervals

**Panel C: Linker Region Conservation**
- Sequence logo for positions 261-269
- Consensus: L-M-V-D-C-S-H-G-N
- Height represents conservation
- D264 and H267 maximally conserved (100%, indicated with asterisks)
- S266 and N269 highly conserved (87.5%)
- Conservation percentages labeled below

**Panel D: Interface Residue Heatmap**
- Heatmap showing conservation and trait-specificity
- Rows: Interface positions (250-290)
- Columns: Conservation score, Trait-specific (yes/no), Variability
- Color scale: Blue (conserved) to Red (variable)
- D264 and H267 highlighted (100% conservation, invariant)

**Figure Legend:**
Domain interfaces enriched for trait-specific residues. (A) Type I architecture showing interface regions at catalytic C-terminus (250-270) and ACT N-terminus (270-290), connected by structured linker (261-269). (B) Interface residues show 1.2-fold enrichment for trait-specific changes (80-81%) compared to genome-wide average (66.8%, p < 0.01). (C) Linker conservation analysis reveals two invariant residues, D264 (Asp) and H267 (His), likely critical for allosteric signal transmission. (D) Heatmap of interface residue properties identifies conserved scaffold residues (blue) versus variable specificity-determining residues (red).

---

## Figure 6: Coevolution and Allosteric Pathway

**Panel A: Coevolution Network**
- Network diagram: Nodes = sequence positions
                   Edges = coevolving pairs (|r| > 0.6)
- Catalytic domain positions: blue circles
- ACT domain positions: red squares
- Edge thickness proportional to |correlation|
- Perfect correlations (r = ±1.0): thick edges, labeled
- Hub residue 280 (ACT): large node, central position
- 17 total edges shown

**Panel B: Perfect Coevolution Pairs**
- Scatter plots showing amino acid encodings
- Panel B1: Position 91 (catalytic) vs. 280 (ACT)
  - Perfect positive correlation (r = +1.000, p < 0.0001)
- Panel B2: Position 171 (catalytic) vs. 280 (ACT)
  - Perfect negative correlation (r = -1.000, p < 0.0001)
- Each point represents a sequence
- Color by regulatory trait
- Regression lines shown

**Panel C: Interface Coevolution**
- Highlighting interface pair: 251 (catalytic C-term) ↔ 290 (ACT N-term)
- r = +0.812, p = 0.0143
- Both positions mapped onto domain architecture schematic
- Distance indicated (~40 Å, structural modeling)
- Interpretation: Direct interface communication

**Panel D: Allosteric Pathway Model**
- Cartoon representation of signal transmission
- Step 1: Aromatic amino acid binds ACT domain pocket
- Step 2: ACT conformational change (interface residues 270-290, 280 highlighted)
- Step 3: Signal through linker (D264, H267 indicated)
- Step 4: Catalytic C-terminus response (interface residues 250-270, 251 highlighted)
- Step 5: Active site modulation (substrate affinity reduced)
- Arrows show signal flow direction
- Coevolving positions (91, 171, 280, 251, 290) highlighted in pathway

**Figure Legend:**
Coevolution identifies interdomain communication network maintaining allosteric coupling. (A) Network of 17 coevolving position pairs (|r| > 0.6, p < 0.05) connects catalytic and ACT domains. Hub residue V280 (ACT N-terminus, 88% conserved) coevolves with multiple catalytic positions. (B) Perfect correlations (r = ±1.000) between catalytic positions 91, 171 and ACT position 280 indicate obligate coupling—these positions cannot evolve independently. (C) Interface coevolution (251↔290, r = +0.812) provides direct evidence for interdomain coupling at boundaries. (D) Integrated allosteric pathway model showing signal transmission: ligand binding (ACT) → interface (280, 290) → structured linker (D264, H267) → catalytic C-terminus (251, 91, 171) → active site modulation.

---

## Figure 7: Type II Compartmentalization Evolution

**Panel A: Transit Peptide Architecture**
- Schematic of three *Arabidopsis* DHS proteins
- DHS1: Transit peptide (47 aa, green) + TIM barrel (478 aa, gray)
- DHS2: Transit peptide (32 aa, green) + TIM barrel (475 aa, gray)
- DHS3: Transit peptide (57 aa, green) + TIM barrel (470 aa, gray)
- Predicted cleavage sites indicated (VFA, FxA motifs)
- Length scale bar

**Panel B: Transit Peptide Properties**
- Bar charts comparing compositional features:
  - Ser+Thr content (%): DHS1 (27.5%), DHS2 (25.0%), DHS3 (27.5%)
  - Asp+Glu content (%): DHS1 (5.0%), DHS2 (7.5%), DHS3 (2.5%)
  - Arg content (%): DHS1 (5.0%), DHS2 (3.8%), DHS3 (3.8%)
- Reference lines for typical chloroplast transit peptides
- All three within expected ranges

**Panel C: Conservation Profile**
- Line plot: X-axis = Position (1-500)
            Y-axis = Conservation score (0-1)
- Separate line for Type II alignment
- Low conservation in N-terminal region (transit peptides, variable)
- High conservation in TIM barrel core (positions ~100-470)
- Active site regions highlighted (368-374, 442-447)
- Mean conservation: 62.7% (labeled)

**Panel D: Trait Evolution and Compartmentalization**
- Phylogenetic tree (Type II) with cellular localization icons
- Bacterial sequences: cytoplasm icon (gray)
- Plant sequences (DHS1/2/3): chloroplast icon (green)
- Node 100/100 highlighted (plastid acquisition, maximum support)
- Parsimony score = 1 (single transition)
- Inset: Compartmentalization strategy schematic
  - Shikimate pathway sequestered in plastid
  - Physical separation from cytosolic amino acids

**Figure Legend:**
Type II evolution through subcellular compartmentalization in plants. (A) *Arabidopsis* DHS enzymes possess variable-length N-terminal transit peptides (32-57 aa) directing proteins to plastids, with predicted cleavage sites (VxA/FxA motifs). (B) Despite length variation, all three transit peptides show characteristic chloroplast targeting composition: enriched Ser+Thr (25-27.5%), depleted acidic residues (2.5-7.5%). (C) Conservation profile shows low conservation in transit peptide region but high conservation (62.7%) in TIM barrel core, with peaks at active site regions. (D) Phylogenetic mapping reveals single plastid-targeting acquisition (node 100/100, maximum support, parsimony = 1), followed by gene duplications producing three plant paralogs. Compartmentalization strategy provides spatial regulation alternative to Type I allostery.

---

**Total Main Figures: 7**

**Supporting Information Figures (to be generated):**
- S1: Full alignment visualizations
- S2: Bootstrap support distributions
- S3: Posterior probability heatmaps
- S4: Domain boundary predictions
- S5: Extended conservation profiles
- S6: Complete coevolution matrix
- S7: Stability prediction distributions

**Figure File Formats:**
- Publication quality: 300 DPI, TIFF or PDF
- Web display: 150 DPI, PNG
- Source files: Matplotlib/Python scripts, Adobe Illustrator (optional)
