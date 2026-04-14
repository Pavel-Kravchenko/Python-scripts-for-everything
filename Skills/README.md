# Bioinformatics Skills for Claude Code

208 compressed skills extracted from a 197-notebook, 6-tier bioinformatics course. Each skill is a dense reference card: code patterns, decision tables, API references, and pitfalls — no tutorial prose.

**Compressed from 70,710 → 31,560 lines (55% reduction).** Max 300 lines per skill. All frontmatter validated (name, description, tool_type, primary_tool). All code fences balanced.

---

## How to Use

### In Claude Code
Reference any skill by name — Claude activates it automatically:
- *"Help me preprocess scRNA-seq data"* → `bio-applied-scrna-preprocessing` activates
- *"Implement KMP string matching"* → `algo-kmp-algorithm` activates
- *"Build a GWAS Manhattan plot"* → `bio-applied-gwas` activates

### Installation
```bash
# Copy all skills to Claude Code
cp -r Skills/*/SKILL.md ~/.claude/skills/

# Or sync with the provided script
python3 scripts/curate_skills.py
```

---

## Skill Index (208 skills)

### Foundations (10 skills)
Linux, Git, Bash, R, statistics, and probability.

| Skill | Use When... |
|-------|-------------|
| `foundations-linux-fundamentals` | Shell commands, file permissions, SSH, piping |
| `foundations-git-version-control` | Branching, merging, conflict resolution |
| `foundations-bash-scripting` | Writing automation scripts, loops, conditionals |
| `foundations-character-encodings` | UTF-8/ASCII issues in FASTA/VCF files |
| `foundations-r-fundamentals` | R data types, vectors, data frames, ggplot2 |
| `foundations-biostatistics-fundamentals` | Descriptive stats, distributions, hypothesis testing |
| `foundations-probability` | Probability distributions, Bayes' theorem |
| `foundations-statistics-python` | scipy.stats, statsmodels, Python statistics |
| `foundations-r-hypothesis-testing-and-nonparametrics` | Sign test, Wilcoxon, Kruskal-Wallis in R |
| `foundations-r-regression-correlation-and-diagnostics` | Linear regression, correlation, residual analysis in R |

### Python for Bioinformatics (28 skills)
Complete Python from basics to advanced, with bioinformatics examples throughout.

| Skill | Use When... |
|-------|-------------|
| `python-bio-python-introduction` | Setting up Python, first scripts |
| `python-bio-variables` | Variable naming, assignment, scope |
| `python-bio-data-types` | int, float, str, bool, None, type conversions |
| `python-bio-operators` | Arithmetic, comparison, logical operators |
| `python-bio-expressions` | Complex expressions, operator precedence |
| `python-bio-strings` | DNA/RNA string manipulation, complement, GC content |
| `python-bio-sequences` | Sequence slicing, indexing, immutability |
| `python-bio-control-flow` | if/elif/else, for/while loops, break/continue |
| `python-bio-functions` | Function design, arguments, return values |
| `python-bio-file-operations` | Reading/writing FASTA, CSV, JSON files |
| `python-bio-lists` | List operations for sequence collections |
| `python-bio-tuples` | Immutable sequences, named tuples |
| `python-bio-dictionaries` | Gene→function mappings, codon tables |
| `python-bio-sets` | Unique gene sets, set operations for Venn diagrams |
| `python-bio-comprehensions` | List/dict/set comprehensions for data processing |
| `python-bio-iterators` | Custom iterators for FASTA parsing, lazy processing |
| `python-bio-generators` | Yield-based generators for large datasets |
| `python-bio-regular-expressions` | Regex for motif finding, sequence patterns |
| `python-bio-classes` | OOP for Sequence, Gene, ProteinRecord classes |
| `python-bio-oop` | Inheritance, polymorphism, abstract classes |
| `python-bio-decorators` | Function decorators for timing, caching |
| `python-bio-context-managers` | File handling with `with` statements |
| `python-bio-error-handling` | Try/except for robust bioinformatics scripts |
| `python-bio-numpy` | Arrays, vectorized operations for biological data |
| `python-bio-pandas` | DataFrames for gene expression, annotation tables |
| `python-bio-data-wrangling` | Merging, pivoting, cleaning biological datasets |
| `python-bio-data-visualization` | Matplotlib/seaborn for volcano plots, heatmaps |
| `python-bio-sql-for-bioinformatics` | SQL queries for biological databases |

### Core Bioinformatics (17 skills)
Fundamental bioinformatics: databases, alignment, phylogenetics, structure.

| Skill | Use When... |
|-------|-------------|
| `bio-core-biological-databases` | Searching NCBI, UniProt, Ensembl |
| `bio-core-biopython-essentials` | SeqIO, Seq objects, GenBank parsing |
| `bio-core-pairwise-sequence-alignment` | Needleman-Wunsch, Smith-Waterman, gap penalties |
| `bio-core-blast-searching` | BLAST setup, e-value interpretation, result parsing |
| `bio-core-multiple-sequence-alignment` | ClustalW, MUSCLE, alignment quality |
| `bio-core-phylogenetics` | Tree construction (NJ, UPGMA, ML), bootstrap |
| `bio-core-protein-structure` | PDB parsing, Ramachandran, DSSP |
| `bio-core-nucleic-acid-structure` | A/B/Z-DNA, RNA structure |
| `bio-core-chromatogram-analysis` | Sanger .ab1 files, base calling |
| `bio-core-sequence-motifs` | PWM construction, motif scanning |
| `bio-core-domains` | Pfam, InterPro, domain architecture |
| `bio-core-gene-ontology` | GO enrichment, term hierarchy |
| `bio-core-pathways` | KEGG, Reactome pathway analysis |
| `bio-core-comparative-genomics` | Synteny, ortholog detection |
| `bio-core-computational-genetics` | HWE, linkage, codon usage |
| `bio-core-hic-analysis` | Hi-C contact maps, TADs, compartments |
| `bio-core-motif-discovery` | De novo motif finding, MEME |

### Applied Bioinformatics (76 skills)
NGS pipelines, omics analysis, clinical genomics, and specialized methods.

| Category | Skills |
|----------|--------|
| **NGS & Variants** | `bio-applied-ngs-fundamentals` · `bio-applied-bio-data-formats` · `bio-applied-variant-calling-and-snp-analysis` · `bio-applied-snp-calling-pipeline` |
| **RNA-seq** | `bio-applied-rna-seq-analysis` |
| **Single Cell** | `bio-applied-scrna-preprocessing` · `bio-applied-dimensionality-reduction` · `bio-applied-cell-type-annotation` · `bio-applied-trajectory-analysis` · `bio-applied-single-cell-scanpy` · `bio-applied-scatac-chromatin` · `bio-applied-cite-seq-integration` · `bio-applied-sc-integration` |
| **Epigenomics** | `bio-applied-chipseq-pipeline` · `bio-applied-differential-binding` · `bio-applied-tf-footprinting` · `bio-applied-wgbs-bismark` · `bio-applied-dmr-analysis` · `bio-applied-epigenetic-clocks` |
| **Metagenomics** | `bio-applied-microbial-diversity` · `bio-applied-taxonomic-profiling` · `bio-applied-functional-annotation` · `bio-applied-assembly-binning` · `bio-applied-qiime2-16s` |
| **Genomics** | `bio-applied-genome-assembly` · `bio-applied-advanced-ngs` · `bio-applied-gwas` · `bio-applied-copy-number-analysis` · `bio-applied-spatial-transcriptomics` · `bio-applied-ont-processing` · `bio-applied-assembly-sv` · `bio-applied-isoform-analysis` |
| **Immunogenomics** | `bio-applied-vdj-biology` · `bio-applied-immune-repertoire` · `bio-applied-hla-typing` |
| **Cancer** | `bio-applied-cancer-transcriptomics` |
| **CRISPR** | `bio-applied-mageck-gene-essentiality` · `bio-applied-screen-qc-normalization` |
| **Small RNA** | `bio-applied-mirna-seq-pipeline` · `bio-applied-lncrna-classification` |
| **Metabolomics** | `bio-applied-lc-ms-preprocessing` · `bio-applied-metabolite-identification` · `bio-applied-metabolic-flux` |
| **Proteomics** | `bio-applied-proteomics` · `bio-applied-structural-methods` |
| **Multi-Omics** | `bio-applied-data-harmonization` · `bio-applied-mofa2` · `bio-applied-mixomics` |
| **Networks** | `bio-applied-ppi-networks` · `bio-applied-network-modules` · `bio-applied-gene-regulatory-networks` |
| **Cheminformatics** | `bio-applied-rdkit-basics` · `bio-applied-qsar-modeling` · `bio-applied-virtual-screening` · `bio-applied-molecular-gnn` |
| **Virology** | `bio-applied-viral-genome-assembly` · `bio-applied-phylodynamics` · `bio-applied-variant-surveillance` |
| **Other** | `bio-applied-statistics-for-bioinformatics` · `bio-applied-machine-learning-for-biology` · `bio-applied-deep-learning-for-biology` · `bio-applied-clinical-genomics` · `bio-applied-bayesian-statistics-python` · `bio-applied-molecular-modeling` · `bio-applied-docking` · `bio-applied-biochemistry` · `bio-applied-enzyme-kinetics` · `bio-applied-genetic-engineering-in-silico` · `bio-applied-population-genetics` · `bio-applied-molecular-evolution` · `bio-applied-numerical-methods-for-bioinformatics` · `bio-applied-promoter` · `bio-applied-regulatory-analysis` · `bio-applied-workflow-engines` · `bio-applied-testing-cicd` · `bio-applied-capstone-project` |

### Algorithms & Data Structures (29 skills)
CS fundamentals with bioinformatics applications (DP = alignment, suffix trees = genome indexing).

| Category | Skills |
|----------|--------|
| **Complexity** | `algo-complexity-analysis` · `algo-basic-algorithms` |
| **Sorting** | `algo-comparison-sorts` · `algo-linear-sorts` |
| **Searching** | `algo-linear-binary-search` |
| **Linear Structures** | `algo-linked-lists` · `algo-stacks-queues` · `algo-dynamic-arrays` |
| **Trees** | `algo-binary-search-trees` · `algo-avl-trees` · `algo-red-black-trees` |
| **Hashing** | `algo-hash-tables-bloom` |
| **String Matching** | `algo-naive-pattern-matching` · `algo-kmp-algorithm` · `algo-rabin-karp` · `algo-dfa-matching` |
| **String Structures** | `algo-tries` · `algo-aho-corasick` · `algo-suffix-arrays` · `algo-suffix-trees` |
| **Graphs** | `algo-graph-representations` · `algo-bfs-dfs` · `algo-dijkstra` · `algo-mst-kruskal-prim` · `algo-topological-sort` |
| **Dynamic Programming** | `algo-intro-memoization` · `algo-tabulation` · `algo-knapsack` · `algo-sequence-alignment` |

### Modern AI for Science (13 skills)
LLMs, diffusion models, AlphaFold, protein language models, genomic foundation models.

| Skill | Use When... |
|-------|-------------|
| `ai-science-llm-finetuning` | LoRA, quantization, SFTTrainer |
| `ai-science-llm-training-systems` | Distributed training, FSDP, DeepSpeed |
| `ai-science-vision-rag` | Vision-language models, ColPali, RAG pipelines |
| `ai-science-diffusion-generative-models` | DDIM, noise schedules, score matching |
| `ai-science-alphafold-protein-design` | AF2 architecture, pLDDT/PAE, RFdiffusion |
| `ai-science-genomic-llms` | Nucleotide Transformer, DNABERT-2, HyenaDNA |
| `ai-science-enformer-regulatory` | Enformer, regulatory element prediction |
| `ai-science-splicing-models` | SpliceAI, splice variant scoring |
| `ai-science-epigenomic-sequence-models` | Epiformer, chromatin accessibility |
| `ai-science-variant-to-structure-models` | AlphaMissense, variant → structure impact |
| `ai-science-esm2-embeddings` | ESM2 protein embeddings, ESMFold |
| `ai-science-zero-shot-mutation` | ESM-1v mutation effect prediction |
| `ai-science-geneformer-scgpt` | Geneformer, scGPT for single-cell |

### Legacy Skills (35 skills)
Consolidated skills from earlier course versions, migrated to directory format.

---

## Quality

- All 208 skills have valid YAML frontmatter (name, description, tool_type, primary_tool)
- All code fences properly balanced
- 201/208 skills have a Pitfalls section
- Max skill size: 300 lines, avg ~152 lines

---

## Course Origin

Extracted from the **Bioinformatics with Python** course — 197 notebooks across 6 tiers. Full course: [Course/README.md](../Course/README.md)
