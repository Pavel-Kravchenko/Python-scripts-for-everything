# Module 06: Protein Language Models

**Tier 5 — Modern AI for Science | Module 06**

Protein sequence language models: ESM2 embeddings, ESMFold, zero-shot mutation scoring, and structure-conditioned design.

## What You'll Learn
- Protein language model pre-training: masked language modeling on UniRef90
- ESM2 model family (8M → 15B parameters): per-residue and per-sequence embeddings
- Downstream tasks: function prediction, localization, thermostability regression
- ESMFold: MSA-free protein structure prediction in seconds
- Evolutionary Scale Modeling attention maps as contact predictors
- Zero-shot mutation effect prediction using log-likelihood ratios (ESM-1v)
- ProteinGym benchmark for mutation effect models
- ESM-IF1 inverse folding for structure-conditioned sequence design
- Antibody language models: IgLM, AbDiffuser

## Prerequisites
- Tier 3 Module 07 (Protein Structure) — 3D structure, PDB format
- Tier 5 Module 04 (AlphaFold & Protein Design) — structure prediction concepts
- Tier 5 Module 01 (LLM Fine-tuning) — transformer fine-tuning

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [ESM2 Embeddings & ESMFold](01_esm2_embeddings.ipynb) | Embeddings, downstream tasks, ESMFold, contact maps |
| 2 | [Zero-Shot Mutation Prediction](02_zero_shot_mutation.ipynb) | ESM-1v scoring, ProteinGym benchmark, inverse folding |

## GPU Requirements
| Notebook | GPU | Notes |
|----------|-----|-------|
| ESM2 Embeddings | T4 required | A100 for ESM2-3B/15B |
| Zero-Shot Scoring | T4 recommended | CPU feasible for ESM2-650M |

## Key Models & Resources
- [ESM GitHub (Meta)](https://github.com/facebookresearch/esm)
- [ProteinGym benchmark](https://proteingym.org/)
- [bio_embeddings documentation](https://docs.bioembeddings.com/)
- [EVE (Fraternali et al. 2021)](https://www.nature.com/articles/s41586-021-04043-8)
- [ESM-IF1 paper (Hsu et al.)](https://www.biorxiv.org/content/10.1101/2022.04.10.487779)

---

[← Module 05: Genomic Foundation Models](../05_Genomic_Foundation_Models/) | [Course Overview](../../README.md) | [Next: 5.07 Foundation Models for Single-Cell →](../07_Foundation_Models_Single_Cell/)
