# Tier 5: Modern AI for Science

GPU-optional modules covering contemporary AI methods applied to scientific research. Every notebook is designed to run on free-tier Google Colab. Theory and code-pattern cells run without a GPU; hands-on training cells require a T4 or A100 runtime.

## Prerequisites

- Tier 3, Module 10 (Deep Learning for Biology) — transformer architectures and PyTorch basics
- Comfort with numpy, pandas, and matplotlib

## Modules

| Module | Topic | GPU Required |
|--------|-------|-------------|
| [01 LLM Fine-tuning](01_LLM_Finetuning/01_LLM_Finetuning.ipynb) | LoRA, quantization, SFTTrainer, instruction datasets, experiment tracking, ablations | Yes (T4 or better) |
| [02 Vision RAG](02_Vision_RAG/02_Vision_RAG.ipynb) | VLMs, ColPali, document retrieval, RAG pipeline | Optional (CPU feasible for inference) |
| [03 Diffusion & Generative Models](03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb) | Score matching, DDIM, inverse problems, noise schedules | Optional (CPU feasible for small examples) |
| [04 AlphaFold & Protein Design](04_AlphaFold_Protein_Design/04_AlphaFold_Protein_Design.ipynb) | AlphaFold2 architecture, pLDDT/PAE, ESMFold, RFdiffusion, ProteinMPNN | Optional (T4 for design) |
| [05 Genomic Foundation Models](05_Genomic_Foundation_Models/README.md) | Nucleotide Transformer, Enformer, Borzoi, SpliceAI, AlphaGenome, Epiformer | Optional (T4 recommended for long-context models) |
| [06 Protein Language Models](06_Protein_Language_Models/README.md) | ESM2 embeddings, ESMFold, mutation effect prediction, inverse folding | Yes (T4 recommended) |
| [07 Foundation Models for Single-Cell](07_Foundation_Models_Single_Cell/README.md) | Geneformer, scGPT, perturbation prediction, Census-scale analysis | Yes (T4 recommended) |

## Running on Colab

Each notebook contains a setup cell that installs all required packages. Recommended runtime: **T4 GPU** (free tier). Switch to CPU runtime for theory-only reading.

```python
# Standard Colab setup (at top of each notebook)
!pip install -q unsloth trl peft bitsandbytes transformers accelerate
```

## Attribution

| Module | Source | Attribution |
|--------|--------|-------------|
| LLM Fine-tuning | FinetuningSmallLanguageModels workshop patterns | Inspired by Unsloth AI & Manuel Faysse |
| Vision RAG | VisionRag workshop patterns | Inspired by Unsloth AI & Manuel Faysse |
| Diffusion & Generative Models | Lundi workshop patterns | Adapted from DDRM (Kawar et al., 2022), github.com/bahjat-kawar/ddrm |

All content uses public datasets (HuggingFace Hub, torchvision). No proprietary research data.
