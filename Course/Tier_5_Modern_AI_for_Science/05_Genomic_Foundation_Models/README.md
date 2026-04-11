# Module 05: Genomic Foundation Models

**Tier 5 — Modern AI for Science | Module 05**

DNA/RNA sequence language models: Nucleotide Transformer, HyenaDNA, Evo, and Enformer for regulatory genomics.

## What You'll Learn
- DNA tokenization strategies: k-mer, BPE, single-nucleotide character-level
- Pre-trained genomic LMs: Nucleotide Transformer (2.5B), DNABERT-2, HyenaDNA, Evo
- Extracting sequence embeddings for downstream classification (splice sites, promoters)
- HyenaDNA: Hyena convolution operators for 1M-token genomic context
- Evo: prokaryotic genome model for zero-shot functional prediction and sequence generation
- Enformer architecture: 196 kb receptive field, predicting 5313 regulatory tracks
- In-silico mutagenesis (ISM) with Enformer for variant effect prediction
- Regulatory sequence design with gradient-guided optimization

## Prerequisites
- Tier 3 Module 10 (Deep Learning for Biology) — transformer architectures, PyTorch
- Tier 5 Module 01 (LLM Fine-tuning) — LoRA, quantization, HuggingFace Trainer

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Genomic LLMs: NT, HyenaDNA, Evo](01_genomic_llms.ipynb) | Tokenization, embeddings, fine-tuning |
| 2 | [Enformer: Regulatory Prediction](02_enformer_regulatory.ipynb) | ISM, eQTL scoring, sequence design |

## GPU Requirements
| Notebook | GPU | Notes |
|----------|-----|-------|
| Genomic LLMs | T4 required (NT-500M) | A100 for NT-2.5B/Evo-7B |
| Enformer | Optional | CPU feasible for inference on short sequences |

## Key Models & Resources
- [Nucleotide Transformer (InstaDeep)](https://github.com/instadeep/nucleotide-transformer)
- [HyenaDNA (Hazyresearch)](https://github.com/HazyResearch/hyena-dna)
- [Evo (Arc Institute)](https://github.com/evo-design/evo)
- [Enformer (DeepMind)](https://github.com/google-deepmind/deepmind-research/tree/master/enformer)
- [DNABERT-2](https://github.com/Zhihan1996/DNABERT_2)

---

[← Module 04: AlphaFold & Protein Design](../04_AlphaFold_Protein_Design/) | [Course Overview](../../README.md) | [Next: 5.06 Protein Language Models →](../06_Protein_Language_Models/)
