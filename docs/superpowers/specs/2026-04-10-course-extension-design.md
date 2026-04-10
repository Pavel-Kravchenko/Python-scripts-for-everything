# Course Extension Design
**Date:** 2026-04-10  
**Status:** Approved

---

## Overview

Extend the existing 5-tier bioinformatics curriculum with new modules across two parallel tracks:

1. **Bioinformatics depth track** — new modules in Tier 2 and Tier 3 drawn from external workshop and school materials, converted research code patterns, and a Bayesian statistics course
2. **Modern AI track** — new Tier 5 covering LLM fine-tuning, vision RAG, and diffusion/generative models

All new content follows strict data policy: no real research data copied, no sensitive identifiers (researcher names, protein targets, unpublished sequences, project-specific gene lists), public datasets only. Code patterns and algorithmic concepts extracted; research context discarded.

---

## Source Materials & Attribution

| Source | Attribution | Use |
|---|---|---|
| `Python_EDA_workshop` | MPIB workshop | Hi-C analysis module patterns |
| `MOTIF_DISCOVERY` package | TotipotencyLab, Max Planck Institute of Biochemistry (MIT license) | Motif discovery code patterns only |
| `NGSchool2023` | NGSchool 2023 | GWAS, spatial transcriptomics, copy number patterns |
| `MPIB/Python_scripts` | TotipotencyLab, MPIB | TF footprinting, chromatin accessibility patterns |
| `ML2023` / `ML2023_2` | — | Classical ML → deep learning examples (already covered in Tier 3; no new module added — patterns may inform minor improvements to existing Tier 3 ML notebooks) |
| `Monday–Friday` R stats | Fränzi Korner-Nievergelt (original R author) | Bayesian stats converted to Python |
| `Mardi` / `FinetuningSmallLanguageModels_VisionRag` | Inspired by Unsloth AI & Manuel Faysse | LLM fine-tuning, vision RAG |
| `Lundi` | Adapted from DDRM (Kawar et al.) | Diffusion models, inverse problems |

---

## Tier 2 — Core Bioinformatics: 2 New Modules

### Module: Hi-C Analysis
- **Placement:** After existing Tier 2 modules (exact number assigned during implementation plan after checking current Tier 2 notebook count)
- **Source:** `Python_EDA_workshop` (7-notebook workshop series)
- **Content:**
  - cooler format: loading, inspecting, slicing contact matrices
  - cooltools API: expected contact frequency, scaling curves
  - Eigenvector decomposition (A/B compartments)
  - Insulation score calculation (TAD boundary detection)
  - Saddle plots
  - Pileup (snipping) analysis
- **Data:** Public Hi-C datasets from 4DN Data Portal or ENCODE (`.cool` format)
- **Deliverable:** `XX_Hi-C_Analysis.ipynb` + `hic-analysis.md` skill

### Module: Motif Discovery
- **Placement:** After Hi-C module
- **Source:** `MOTIF_DISCOVERY` code patterns (MIT license)
- **Content:**
  - Position weight matrices (PWM/PPM): construction, normalization
  - Information content per position; total IC
  - KDIC scoring (mean IC normalized to [0,1])
  - IUPAC consensus sequence generation
  - Score distributions: exact enumeration (motif length ≤ 10), Monte Carlo for longer
  - Motif enrichment: Fisher's exact test + Benjamini-Hochberg correction
  - TomTom matching concept (JASPAR/HOCOMOCO public databases)
  - Pipeline design patterns: abstract tool interface, dataclass configuration
  - BED/FASTA I/O patterns; TSS region extraction; peak loading/filtering
- **Data:** JASPAR 2024 public motif database; ENCODE ChIP-seq peaks (public)
- **Bugs to avoid:** PPM shape ambiguity (L×4 vs 4×L) — teach explicit shape convention; Monte Carlo convergence — always report sample size and CI
- **Deliverable:** `XX_Motif_Discovery.ipynb` + `motif-discovery.md` skill

---

## Tier 3 — Applied Bioinformatics: 5 New Modules

### Module: GWAS
- **Source:** `NGSchool2023/GWAS_Michail` (patterns only, no bundled data)
- **Content:**
  - GWAS study design: cases/controls, phenotype definition, confounders
  - Quality control: SNP and sample filtering criteria, Hardy-Weinberg, MAF thresholds
  - Population stratification: PCA on genotype data
  - Association testing: logistic/linear regression per SNP
  - Multiple testing in GWAS: genome-wide significance threshold (5×10⁻⁸)
  - Manhattan and QQ plots
  - LD and clumping concepts
  - Downstream: fine-mapping, GWAS catalog lookup
- **Data:** 1000 Genomes Project public VCF subset or simulated genotype data
- **Attribution note:** NGSchool 2023
- **Deliverable:** `XX_GWAS.ipynb` + `gwas-population-genetics.md` skill (extends existing `population-genetics-evolution.md`)

### Module: Spatial Transcriptomics
- **Source:** `NGSchool2023` spatial notebooks (patterns only; no original data)
- **Content:**
  - Spatial data formats: AnnData with spatial coordinates, Visium/Xenium layout
  - Quality control for spatial data
  - Normalization and dimensionality reduction in spatial context
  - Spatial neighborhood graphs
  - Spatially variable gene detection
  - Cell-type deconvolution concept (RCTD, cell2location patterns)
  - Visualization: spatial scatter plots, expression overlays
- **Data:** Public 10x Visium datasets (mouse brain, available via squidpy/scanpy)
- **Attribution note:** NGSchool 2023
- **Deliverable:** `XX_Spatial_Transcriptomics.ipynb` + `spatial-transcriptomics.md` skill

### Module: DNA Copy Number Analysis
- **Source:** `NGSchool2023/DNA_copy_number` (patterns only)
- **Content:**
  - Copy number variation (CNV) concepts: gains, losses, LOH
  - Read depth normalization approaches
  - Segmentation algorithms (CBS concept)
  - Calling copy number states from segments
  - Visualization: genome-wide CN profiles
  - Downstream: gene-level CN annotation
- **Data:** Public cancer genome data (TCGA WGS subset or simulated)
- **Attribution note:** NGSchool 2023
- **Deliverable:** `XX_Copy_Number_Analysis.ipynb` (no separate skill — extends existing `ngs-variant-calling.md`)

### Module: TF Footprinting & Chromatin Accessibility
- **Source:** `MPIB/Python_scripts` code patterns (no original data, no project-specific context)
- **Content:**
  - ATAC-seq data recap: fragment size distribution, nucleosome-free regions
  - TF footprinting concept: insertion bias around motif sites
  - Computing expected vs. observed cut-site profiles
  - Footprint score calculation
  - Motif-set comparison across conditions (accessibility enrichment)
  - Intersection of BED files: genomic interval arithmetic with pybedtools
  - Accumulation plots: signal around genomic features
- **Data:** Public ATAC-seq from ENCODE (FASTQ or processed BAM)
- **Attribution:** TotipotencyLab, MPIB (code pattern inspiration)
- **Deliverable:** `XX_TF_Footprinting.ipynb` (no separate skill — extends existing `ngs-variant-calling.md` and new `motif-discovery.md`)

### Module: Bayesian Statistics in Python
- **Source:** Fränzi Korner-Nievergelt's Monday–Friday R statistics course — converted to Python
- **Content (7 sub-sections, one notebook):**
  1. Linear models: frequentist vs. Bayesian framing, parameter estimation, credible intervals (pymc + statsmodels)
  2. Prior specification: informative vs. weakly informative, prior predictive checks
  3. Multiple regression: collinearity diagnostics, VIF, model interpretation
  4. Model comparison: WAIC, LOO-CV, information criteria (arviz)
  5. Linear mixed-effects models: random intercepts/slopes (pymc + bambi)
  6. GLMs: Bernoulli, Binomial, Poisson, Negative Binomial (statsmodels + pymc)
  7. Advanced: GLMM, zero-inflated models, GAM/GAMM concepts, meta-analysis
- **Data:** Public ecological datasets — palmerpenguins, NEON ecological surveys, or equivalent open CSV; no data from Fränzi Korner-Nievergelt's course materials copied
- **Attribution:** Fränzi Korner-Nievergelt (original R course author)
- **Deliverable:** `XX_Bayesian_Statistics_Python.ipynb` + `bayesian-python.md` skill

---

## Tier 5 — Modern AI for Science (New Tier)

**Philosophy:** GPU-optional (Colab-first). Each notebook runs on free-tier Colab. Theory-first, then hands-on pattern. No production deployment scope.

**Tier 5 Structure:**
```
Course/Tier_5_Modern_AI_for_Science/
├── README.md
├── 01_LLM_Finetuning/
│   └── 01_LLM_Finetuning.ipynb
├── 02_Vision_RAG/
│   └── 02_Vision_RAG.ipynb
└── 03_Diffusion_Generative_Models/
    └── 03_Diffusion_Generative_Models.ipynb
```

### Module: LLM Fine-tuning
- **Source:** `Mardi` / `FinetuningSmallLanguageModels_VisionRag` patterns
- **Attribution:** Inspired by Unsloth AI & Manuel Faysse
- **Content:**
  - Base vs. instruction/chat models: what changes during fine-tuning
  - LoRA: low-rank adapter math, rank selection, target modules
  - Quantization: 4-bit, 8-bit; trade-offs with quality
  - Chat templating: system/user/assistant structure
  - SFTTrainer workflow: dataset prep, training loop, evaluation
  - Synthetic data generation for instruction tuning
  - Practical tips: gradient checkpointing, batch size, learning rate
- **Data:** Public instruction datasets (Alpaca, OpenHermes, or similar HuggingFace Hub datasets)
- **Deliverable:** `01_LLM_Finetuning.ipynb` + `llm-finetuning.md` skill

### Module: Vision RAG
- **Source:** Same as above
- **Attribution:** Inspired by Unsloth AI & Manuel Faysse
- **Content:**
  - Vision-language models: architecture overview (encoder + LLM decoder)
  - Document understanding: page-level vs. token-level retrieval
  - ColPali: late interaction for document retrieval (concept + API pattern)
  - RAG pipeline: retrieval → context injection → generation
  - Qwen2-VL inference pattern
  - Evaluation: retrieval recall, generation faithfulness
- **Data:** Public PDF documents (arXiv papers, open-access reports)
- **Deliverable:** `02_Vision_RAG.ipynb` + `vision-rag.md` skill

### Module: Diffusion & Generative Models
- **Source:** `Lundi` (diffusion.py, imaging.py, inv_problems.py patterns)
- **Attribution:** Adapted from DDRM (Kawar et al.), bahjat-kawar/ddrm
- **Content:**
  - Score-based generative models: intuition and math
  - DDIM: deterministic sampling, noise schedule, reverse process
  - Inverse problems in imaging: denoising, inpainting, colorization as special cases
  - SVD-based degradation operators
  - Linear scheduler implementation pattern
  - Score field visualization
  - Scientific applications: cryo-EM denoising, medical image restoration concepts
- **Data:** Public image datasets (MNIST, CIFAR-10, or STL-10 via torchvision)
- **Bugs to address from source:** `diffusion.py` linear scheduler uses implicit shape assumptions — teach explicit tensor shape annotation
- **Deliverable:** `03_Diffusion_Generative_Models.ipynb` + `diffusion-generative.md` skill

---

## New Skill Files (8 total)

| Skill File | Tier | Covers |
|---|---|---|
| `hic-analysis.md` | Tier 2 | cooler, cooltools, contact matrices, TADs, compartments |
| `motif-discovery.md` | Tier 2 | PWMs, IC, KDIC, Fisher enrichment, TomTom, pipeline design |
| `gwas-population-genetics.md` | Tier 3 | GWAS design, QC, PCA, Manhattan plots, fine-mapping |
| `spatial-transcriptomics.md` | Tier 3 | AnnData spatial, SVGs, deconvolution, squidpy |
| `bayesian-python.md` | Tier 3 | pymc, bambi, arviz, GLMs, LMEs, GAMs, priors |
| `llm-finetuning.md` | Tier 5 | LoRA, quantization, SFT, chat templates |
| `vision-rag.md` | Tier 5 | VLMs, ColPali, RAG pipeline, document retrieval |
| `diffusion-generative.md` | Tier 5 | DDIM, score matching, inverse problems, noise schedules |

All follow the existing 5-section format: **When to Use → Quick Reference → Key Patterns → Code Templates → Pitfalls**

---

## README Updates

| File | Change |
|---|---|
| `README.md` (root) | Update tier count (5→6 tiers including new Tier 5); update topic summary |
| `Course/README.md` | Add all 10 new modules to TOC with descriptions and prerequisites |
| `Skills/README.md` | Add 8 new skill entries with one-line descriptions |
| `Course/Tier_5_Modern_AI_for_Science/README.md` | New file: tier overview, GPU requirements, Colab links, attribution |

---

## Data Policy (Non-Negotiable)

- **No copying:** zero data files from source directories transferred
- **No real assets:** no images, model weights, gene lists, FASTA sequences, BED files from research sources
- **No identifiers:** no researcher names in notebooks, no protein targets, no project names, no unpublished sequences
- **Public datasets only:** every notebook links to a public, reproducible data source
- **Attribution in metadata cells:** trainer/author credits go in a dedicated markdown cell at notebook top, not in code

---

## Parallel Execution Plan (7 Sonnet agents)

| Agent | Track | Notebooks | Skills | Notes |
|---|---|---|---|---|
| 1 | Hi-C (Tier 2) | 1 | 1 | — |
| 2 | Motif Discovery (Tier 2) | 1 | 1 | — |
| 3 | GWAS + Spatial + Copy Number (Tier 3) | 3 | 2 | Heaviest bioinformatics agent; split further if needed |
| 4 | Bayesian Stats Python (Tier 3) | 1 | 1 | Long notebook (7 sub-sections); use clear section headers |
| 5 | TF Footprinting (Tier 3) | 1 | 0 | — |
| 6 | Tier 5 (LLM + Vision + Diffusion) | 3 | 3 | Heavy; Colab-first constraint must be enforced per notebook |
| 7 | README updates across all tiers | 0 | 0 | Runs after Agents 1–6 complete |

Total deliverables: **10 new notebooks, 8 new skill files, 4 updated READMEs**

---

## Success Criteria

- All new notebooks run top-to-bottom with public data only
- Every notebook has: motivation cell, theory cell, code pattern cells, exercise cell
- Every new skill file follows existing 5-section format and cross-references related skills
- No sensitive data, names, or research-specific content present anywhere
- READMEs accurately reflect new content with correct module numbering
- Attribution present for all external sources
