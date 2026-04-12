# Module 05: Genomic Foundation Models

**Tier 5 — Modern AI for Science | Module 05**

DNA/RNA sequence language models: Nucleotide Transformer, HyenaDNA, Evo, Enformer, Borzoi, SpliceAI, AlphaGenome, and Epiformer for regulatory genomics.

## What You'll Learn
- DNA tokenization strategies: k-mer, BPE, single-nucleotide character-level
- Pre-trained genomic LMs: Nucleotide Transformer (2.5B), DNABERT-2, HyenaDNA, Evo
- Extracting sequence embeddings for downstream classification (splice sites, promoters)
- HyenaDNA: Hyena convolution operators for 1M-token genomic context
- Evo: prokaryotic genome model for zero-shot functional prediction and sequence generation
- Enformer and Borzoi for large-context regulatory and RNA-seq prediction
- SpliceAI for splice-specific variant effect scoring
- AlphaGenome as an emerging unified model for expression, splicing, chromatin, and contact prediction
- Epiformer for chromatin accessibility prediction from sequence plus conservation features
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
| 3 | [Splicing Models: SpliceAI & AlphaGenome](03_splicing_models.ipynb) | Delta scores, splice usage, variant prioritization |
| 4 | [Epigenomic Sequence Models: Borzoi & Epiformer](04_epigenomic_sequence_models.ipynb) | RNA-seq coverage, chromatin accessibility, modality-aware benchmarking |
| 5 | [From DNA Variants to Protein Structure](05_variant_to_structure_models.ipynb) | AlphaFold2/3, RoseTTAFold, coding variant follow-up |

## GPU Requirements
| Notebook | GPU | Notes |
|----------|-----|-------|
| Genomic LLMs | T4 required (NT-500M) | A100 for NT-2.5B/Evo-7B |
| Enformer | Optional | CPU feasible for inference on short sequences |
| Splicing Models | Optional | CPU feasible for scoring APIs and small batch inference |
| Epigenomic Sequence Models | T4 recommended | Large-window model inference benefits from GPU memory |
| Variant → Structure | Optional | AF3/RoseTTAFold full runs often need dedicated compute |

## New Sections in This Module

1. **Splicing-first variant interpretation** with SpliceAI and AlphaGenome splice outputs.
2. **Epigenomic signal modeling** with Borzoi and Epiformer for expression/accessibility tasks.
3. **Genomics-to-structure handoff** using AlphaFold2, AlphaFold3, and RoseTTAFold for coding variants.

## Key Models & Resources
- [Nucleotide Transformer (InstaDeep)](https://github.com/instadeep/nucleotide-transformer)
- [HyenaDNA (Hazyresearch)](https://github.com/HazyResearch/hyena-dna)
- [Evo (Arc Institute)](https://github.com/evo-design/evo)
- [Enformer (DeepMind)](https://github.com/google-deepmind/deepmind-research/tree/master/enformer)
- [Borzoi (Calico)](https://github.com/calico/borzoi)
- [SpliceAI (Illumina)](https://github.com/Illumina/SpliceAI)
- [AlphaGenome (DeepMind)](https://github.com/google-deepmind/alphagenome_research)
- [Epiformer](https://github.com/yal054/epiformer)
- [DNABERT-2](https://github.com/Zhihan1996/DNABERT_2)

## Assignments and Solutions
- [Assignment Set 2: Genomic Foundation Models](../../../Assignments/Tier_5_Modern_AI/02_genomic_foundation_models.ipynb)
- [Solutions Set 2: Genomic Foundation Models](../../../Solutions/Tier_5_Modern_AI/02_genomic_foundation_models_solutions.ipynb)

---

[← Module 04: AlphaFold & Protein Design](../04_AlphaFold_Protein_Design/) | [Course Overview](../../README.md) | [Next: 5.06 Protein Language Models →](../06_Protein_Language_Models/)
