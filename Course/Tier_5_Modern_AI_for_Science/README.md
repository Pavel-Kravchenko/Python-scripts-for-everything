# Tier 5: Modern AI for Science

GPU-optional modules covering contemporary AI methods applied to scientific research. Every notebook is designed to run on free-tier Google Colab. Theory and code-pattern cells run without a GPU; hands-on training cells require a T4 or A100 runtime.

## Prerequisites

- Tier 3, Module 10 (Deep Learning for Biology) — transformer architectures and PyTorch basics
- Comfort with numpy, pandas, and matplotlib

## Modules

| Module | Topic | GPU Required |
|--------|-------|-------------|
| [01 LLM Fine-tuning](01_LLM_Finetuning/01_LLM_Finetuning.ipynb) | LoRA, quantization, SFTTrainer, instruction datasets, experiment tracking, ablations | Yes (T4 or better) |
| [01b LLM Training Systems](01_LLM_Finetuning/02_llm_training_systems.ipynb) | Distributed training, FSDP, DeepSpeed, gradient checkpointing | Yes (multi-GPU) |
| [02 Vision RAG](02_Vision_RAG/02_Vision_RAG.ipynb) | VLMs, ColPali, document retrieval, RAG pipeline | Optional (CPU feasible for inference) |
| [03 Diffusion & Generative Models](03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb) | Score matching, DDIM, inverse problems, noise schedules | Optional (CPU feasible for small examples) |
| [04 AlphaFold & Protein Design](04_AlphaFold_Protein_Design/04_AlphaFold_Protein_Design.ipynb) | AlphaFold2 architecture, pLDDT/PAE, ESMFold, RFdiffusion, ProteinMPNN | Optional (T4 for design) |
| [05a Genomic LLMs](05_Genomic_Foundation_Models/01_genomic_llms.ipynb) | Nucleotide Transformer, DNABERT-2, HyenaDNA, k-mer tokenization | Optional (T4 recommended) |
| [05b Enformer & Regulatory](05_Genomic_Foundation_Models/02_enformer_regulatory.ipynb) | Enformer architecture, regulatory element prediction | Optional (T4 recommended) |
| [05c Splicing Models](05_Genomic_Foundation_Models/03_splicing_models.ipynb) | SpliceAI, Pangolin, splice variant effect prediction | Optional |
| [05d Epigenomic Models](05_Genomic_Foundation_Models/04_epigenomic_sequence_models.ipynb) | Epiformer, chromatin accessibility prediction | Optional (T4 recommended) |
| [05e Variant-to-Structure](05_Genomic_Foundation_Models/05_variant_to_structure_models.ipynb) | AlphaMissense, variant effect → structural impact | Optional |
| [06a ESM2 Embeddings](06_Protein_Language_Models/01_esm2_embeddings.ipynb) | ESM2 architecture, protein embeddings, ESMFold | Yes (T4 recommended) |
| [06b Zero-shot Mutation](06_Protein_Language_Models/02_zero_shot_mutation.ipynb) | ESM-1v scoring, masked log-likelihood, mutation effects | Yes (T4 recommended) |
| [07 Geneformer & scGPT](07_Foundation_Models_Single_Cell/01_geneformer_scgpt.ipynb) | Rank-value tokenization, cell-type annotation, perturbation prediction | Yes (T4 recommended) |

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
