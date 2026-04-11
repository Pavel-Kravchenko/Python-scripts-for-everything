# Module 07: Foundation Models for Single-Cell Biology

**Tier 5 — Modern AI for Science | Module 07**

Pre-trained single-cell foundation models: Geneformer, scGPT, and CellxGene Census for cell type annotation and perturbation prediction.

## What You'll Learn
- Gene rank-based tokenization: representing single-cell transcriptomes as ordered gene sequences
- Geneformer pre-training on 30M cells from Genecorpus-30M
- Fine-tuning Geneformer for cell type classification and gene network inference
- scGPT: generative transformer for multi-task single-cell analysis
- Perturbation prediction: forecasting transcriptional responses to genetic/chemical perturbations
- Universal Cell Embeddings (scFoundation, UCE) for cross-study integration
- CellxGene Census: programmatic access to 50M+ human and mouse cells
- Zero-shot cell type annotation with foundation model embeddings

## Prerequisites
- Module 30 (Single-Cell RNA-seq) — scRNA-seq preprocessing, AnnData, clustering
- Tier 5 Module 01 (LLM Fine-tuning) — HuggingFace Trainer, LoRA

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Geneformer & scGPT](01_geneformer_scgpt.ipynb) | Tokenization, annotation, perturbation prediction, CellxGene Census |

## GPU Requirements
| Notebook | GPU | Notes |
|----------|-----|-------|
| Geneformer & scGPT | T4 required | A100 for full fine-tuning; LoRA feasible on T4 |

## Key Models & Resources
- [Geneformer (Theodoris et al. 2023, Nature)](https://www.nature.com/articles/s41586-023-06139-9)
- [scGPT (Cui et al. 2024, Nature Methods)](https://www.nature.com/articles/s41592-024-02201-0)
- [scFoundation](https://github.com/biomap-research/scFoundation)
- [CellxGene Census API](https://chanzuckerberg.github.io/cellxgene-census/)
- [Universal Cell Embeddings (UCE)](https://github.com/snap-stanford/UCE)

---

[← Module 06: Protein Language Models](../06_Protein_Language_Models/) | [Course Overview](../../README.md)
