# Bioinformatics Skills for Claude Code

38 skill files compressing a 106-notebook, 6-tier bioinformatics course into actionable Claude Code references. Each skill provides key patterns, code templates, complexity tables, and common pitfalls for a focused topic area.

**Maximum content. Maximum yield. Minimum tokens.**

---

## How to Use

### In Claude Code
Reference any skill by name in your prompt — Claude will activate the relevant skill automatically:
- *"Help me write a BLAST search"* → `sequence-alignment` skill activates
- *"Parse this VCF file"* → `ngs-variant-calling` skill activates
- *"Build a phylogenetic tree"* → `phylogenetics-evolution` skill activates

### Installation
Copy the `Skills/` directory into your project, or install as a Claude Code plugin:
```bash
# Option 1: Copy skills into your project
cp -r Skills/ /path/to/your/project/.claude/skills/

# Option 2: Reference from Claude Code settings
# Add the skills directory path to your Claude Code configuration
```

### Combining Skills
Skills are modular — combine them for complex workflows:
- `sequence-alignment` + `string-algorithms` = biological tool usage + underlying algorithmic understanding
- `ngs-variant-calling` + `linux-git-bash` = variant pipeline + shell automation
- `ml-deep-learning-bio` + `numpy-pandas-wrangling` = ML models + data preparation

---

## Skill Index

### Foundations (Skills 1–3)

| Skill | Use When... |
|-------|-------------|
| [`linux-git-bash`](linux-git-bash.md) | Writing shell pipelines, git workflows, processing FASTA/BAM/VCF on CLI |
| [`biostatistics-r`](biostatistics-r.md) | Choosing statistical tests, multiple testing correction, R code for bio data |
| [`advanced-r-statistics`](advanced-r-statistics.md) | Nonparametric tests in R, regression diagnostics, Hodges-Lehman, power analysis |

### Python for Bioinformatics (Skills 3–7)

| Skill | Use When... |
|-------|-------------|
| [`python-core-bio`](python-core-bio.md) | String manipulation for DNA/RNA, file I/O for FASTA/GenBank, function design |
| [`python-collections-regex`](python-collections-regex.md) | k-mer counting, gene set operations, generators for large files, regex for motifs |
| [`python-advanced-sql`](python-advanced-sql.md) | OOP for bio classes, decorators, error handling, SQL for biological databases |
| [`numpy-pandas-wrangling`](numpy-pandas-wrangling.md) | Expression matrices, annotation merges, tidy data, vectorized operations |
| [`data-visualization-bio`](data-visualization-bio.md) | Volcano plots, heatmaps, genome tracks, publication-quality figures |

### Core Bioinformatics (Skills 8–11)

| Skill | Use When... |
|-------|-------------|
| [`biopython-databases`](biopython-databases.md) | Fetching sequences from NCBI, parsing FASTA/GenBank, BioPython Seq/SeqIO |
| [`sequence-alignment`](sequence-alignment.md) | Pairwise/multiple alignment, BLAST searching, substitution matrices |
| [`phylogenetics-evolution`](phylogenetics-evolution.md) | Building trees (NJ, UPGMA, ML), bootstrap, comparative genomics, synteny |
| [`structural-bioinformatics`](structural-bioinformatics.md) | PDB parsing, protein structure, chromatograms, motifs/PWM, GO enrichment |

### Applied Bioinformatics (Skills 12–16)

| Skill | Use When... |
|-------|-------------|
| [`ngs-variant-calling`](ngs-variant-calling.md) | FASTQ QC, read alignment, SAM/BAM, variant calling pipelines, VCF parsing |
| [`rnaseq-metagenomics`](rnaseq-metagenomics.md) | Differential expression (DESeq2), 16S amplicon analysis, promoter scanning |
| [`ml-deep-learning-bio`](ml-deep-learning-bio.md) | scikit-learn classifiers, CNNs for genomics, protein language models |
| [`clinical-modeling-workflows`](clinical-modeling-workflows.md) | Molecular docking, ACMG classification, Scanpy, Snakemake, CI/CD |
| [`genome-assembly-proteomics`](genome-assembly-proteomics.md) | De novo assembly metrics, mass spectrometry data, protein identification |

### Algorithms & Data Structures (Skills 16–20)

| Skill | Use When... |
|-------|-------------|
| [`complexity-sorting-searching`](complexity-sorting-searching.md) | Big-O analysis, implementing sorts, binary search variants |
| [`linear-tree-hash-structures`](linear-tree-hash-structures.md) | Linked lists, BST/AVL/RB trees, hash tables, Bloom filters |
| [`string-algorithms`](string-algorithms.md) | KMP, Rabin-Karp, DFA matching for sequence pattern search |
| [`advanced-string-structures`](advanced-string-structures.md) | Tries, Aho-Corasick, suffix arrays/trees for genome indexing |
| [`graphs-dynamic-programming`](graphs-dynamic-programming.md) | Graph traversal, shortest paths, MST, DP for sequence alignment |

### Biology & Computation (Skills 22–28)

| Skill | Use When... |
|-------|-------------|
| [`probability-statistics-python`](probability-statistics-python.md) | Statistical tests in Python (scipy/statsmodels), regression, power analysis, distributions |
| [`genetics-computational`](genetics-computational.md) | Codon tables, restriction mapping, genetic crosses, codon usage bias, HWE |
| [`biochemistry-enzymology`](biochemistry-enzymology.md) | Enzyme kinetics fitting, inhibition analysis, EC classification, metabolic pathways |
| [`genetic-engineering-insilico`](genetic-engineering-insilico.md) | Primer design, CRISPR guides, cloning workflows, codon optimization |
| [`population-genetics-evolution`](population-genetics-evolution.md) | Drift simulation, selection, dN/dS, molecular clock, Fst, neutrality tests |
| [`numerical-methods-bio`](numerical-methods-bio.md) | Interpolation, curve fitting, optimization, FFT for bioinformatics data |

### Tier 2 Depth (Skills 29–30)

| Skill | Use When... |
|-------|-------------|
| [`hic-analysis`](hic-analysis.md) | Analyzing 3D genome organization: loading cooler/mcool files, computing A/B compartments, TAD boundaries, pileup plots with cooltools |
| [`motif-discovery`](motif-discovery.md) | Building PWMs from aligned binding sites, IC/KDIC scoring, Fisher enrichment testing, TomTom database matching |

### Applied Bioinformatics Depth (Skills 31–35)

| Skill | Use When... |
|-------|-------------|
| [`gwas-population-genetics`](gwas-population-genetics.md) | GWAS study design, QC pipelines, PCA for stratification, Manhattan/QQ plots, LD clumping, fine-mapping |
| [`spatial-transcriptomics`](spatial-transcriptomics.md) | AnnData with spatial coordinates, spatially variable gene detection, Squidpy neighborhood graphs, cell-type deconvolution, NMF programs |
| [`bayesian-python`](bayesian-python.md) | Bayesian linear and generalized linear models with pymc/bambi, model comparison via arviz LOO-CV/WAIC, prior specification |
| [`copy-number-analysis`](copy-number-analysis.md) | Log2 ratio computation, CBS segmentation, CN state calling, read depth normalization, gene annotation with pyranges, tumor purity |
| [`tf-footprinting-atac`](tf-footprinting-atac.md) | ATAC-seq QC (NFR fraction, nucleosomal ladder), Tn5 offset correction, footprint scoring, pybedtools operations, meta-profile plots |

### Modern AI for Science (Skills 36–38)

| Skill | Use When... |
|-------|-------------|
| [`llm-finetuning`](llm-finetuning.md) | Fine-tuning LLMs with LoRA and quantization, configuring SFTTrainer, designing chat templates, building instruction datasets |
| [`vision-rag`](vision-rag.md) | Vision-language models, ColPali document retrieval, RAG pipeline construction, Qwen2-VL inference patterns |
| [`diffusion-generative`](diffusion-generative.md) | DDIM sampling, linear/cosine noise schedules, score matching intuition, inverse imaging problems, scientific image restoration |

---

## Skill File Format

Each skill follows a consistent structure for fast scanning:

| Section | Purpose |
|---------|---------|
| **When to Use** | Activation scenarios (2–3 lines) |
| **Quick Reference** | Tables of facts, formulas, complexity — maximum density |
| **Key Patterns** | Numbered patterns with minimal code (≤10 lines each) |
| **Code Templates** | Copy-paste-ready functions with biological examples |
| **Common Pitfalls** | Bulleted mistakes and fixes |
| **Related Skills** | Cross-references to complementary skills |

---

## Course Origin

Based on the **Bioinformatics with Python** course — 106 notebooks across 6 tiers, built from materials by:

- **[Kodomo Bioinformatics Program](https://kodomo.fbb.msu.ru/wiki/2017)** — Moscow State University (10-semester curriculum)
- **[An Introduction to Applied Bioinformatics](https://readiab.org/)** — J. Gregory Caporaso, Northern Arizona University
- **Summer School of Bioinformatics** — Statistical methods, NGS analysis, promoter research

Full course: [Course/README.md](../Course/README.md) | Attribution: [Course/CREDITS.md](../Course/CREDITS.md)

---

## License

Educational use only. All intellectual property rights for original materials remain with their respective authors and institutions. See the main repository [disclaimer](../README.md#disclaimer).
