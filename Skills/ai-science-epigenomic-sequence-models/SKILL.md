---
name: ai-science-epigenomic-sequence-models
description: "Borzoi vs Epiformer model selection for RNA-seq coverage vs chromatin accessibility prediction"
tool_type: python
primary_tool: NumPy
---

# Epigenomic Sequence Models: Borzoi and Epiformer

## Key Distinctions

- **Borzoi**: sequence-only input → RNA-seq coverage at 32 bp resolution. Captures splicing, UTR usage, isoform signals
- **Epiformer**: sequence + conservation track input → chromatin accessibility. Requires precomputed PhyloP scores
- **AlphaGenome** (2025): unified model for expression, splicing, chromatin, 3D contacts — avoids model-selection but heavier compute

## Model Choice Matrix

| Goal | Preferred model |
|---|---|
| RNA-seq coverage from sequence | Borzoi |
| ATAC/chromatin accessibility | Epiformer |
| Broad multi-modal variant effects | AlphaGenome |
| Splicing triage | SpliceAI |

## Gotchas

- Enformer predicts 128 bp bin averages of ChIP/CAGE/DNase. Borzoi predicts RNA-seq at 32 bp — different modality, not interchangeable
- Epiformer cannot be used without precomputed conservation track for the genome region
- Conservation helps accessibility prediction because conserved non-coding regions are enriched for functional regulatory elements
- Always benchmark models on the output modality you need — cross-task mismatch degrades performance

## Sources

- [Borzoi](https://github.com/calico/borzoi)
- [Epiformer](https://github.com/yal054/epiformer)
