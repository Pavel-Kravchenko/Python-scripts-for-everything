# Tier 5 Module 04: AlphaFold & Protein Design

**Tier 5 — Modern AI for Science | Module 04**

Practical protein structure prediction and computational protein design using AlphaFold2, ESMFold, RoseTTAFold, and generative design tools.

## What You'll Learn
- AlphaFold2 architecture: MSA, evoformer, structure module
- Running AlphaFold via Google Colab (free tier)
- Interpreting confidence metrics: pLDDT, PAE (Predicted Aligned Error)
- Multimer prediction and protein complex modeling
- ESMFold: language-model-only structure prediction (no MSA required)
- RoseTTAFold and RFdiffusion for protein backbone design
- ProteinMPNN for inverse folding (sequence design for a given backbone)
- Evaluating predicted structures: TM-score, RMSD, MolProbity
- Downstream analysis: binding site prediction, pocket detection
- AlphaFold DB: searching and interpreting precomputed structures

## Prerequisites
- Tier 2 Module 06 (Protein Structure) — PDB format, structural concepts
- Module T5-01 (LLM Fine-tuning) — transformer architectures

## Planned Notebooks

| # | Notebook | Topics | GPU |
|---|----------|--------|-----|
| 1 | [Structure Prediction with AlphaFold & ESMFold](01_alphafold_esmfold.ipynb) | Colab AF2, ESMFold API, pLDDT/PAE interpretation | Optional |
| 2 | [Protein Complex Modeling](02_complex_modeling.ipynb) | AlphaFold-Multimer, ColabFold, interface analysis | Optional |
| 3 | [Generative Protein Design](03_protein_design.ipynb) | RFdiffusion backbone design, ProteinMPNN inverse folding | Yes (T4) |

## Key Tools
| Tool | Purpose |
|------|---------|
| AlphaFold2 (Colab) | Structure prediction |
| ColabFold | Fast AF2 via MMseqs2 MSA |
| ESMFold | MSA-free structure prediction |
| ProteinMPNN | Sequence design for fixed backbone |
| RFdiffusion | Generative backbone design |
| BioPython | Structure parsing and analysis |
| py3Dmol / nglview | 3D structure visualization |
| TM-align | TM-score structural comparison |

## Resources
- [AlphaFold Colab notebook](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)
- [ColabFold (fast AlphaFold2)](https://github.com/sokrypton/ColabFold)
- [ESMFold — ESM Atlas](https://esmatlas.com/resources?action=fold)
- [RoseTTAFold tutorials](https://www.rosettacommons.org/)
- [ProteinMPNN paper and code](https://github.com/dauparas/ProteinMPNN)
- [RFdiffusion paper and code](https://github.com/RosettaCommons/RFdiffusion)
- [AlphaFold DB](https://alphafold.ebi.ac.uk/)

## Running on Colab

All notebooks are designed to run on free-tier Google Colab. Structure prediction cells require a T4 or A100 GPU; visualization and analysis cells run on CPU.

---

[← Previous: T5-03 Diffusion & Generative Models](../03_Diffusion_Generative_Models/) | [Course Overview](../../README.md)
