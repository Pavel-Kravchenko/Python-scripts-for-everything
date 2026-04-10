# Module 14: Hi-C Analysis

**Tier 2 — Core Bioinformatics | Module 14**

3D genome organization analysis using Hi-C data.

## What You'll Learn

- cooler format: loading, inspecting, slicing contact matrices
- Contact decay curves: expected contacts and P(s) normalization
- A/B compartment detection via eigenvector decomposition
- TAD boundary calling with insulation score
- Pileup (aggregate) analysis around genomic features
- Saddle plot interpretation for compartment strength

## Prerequisites

- Module 01 (Biological Databases) — data access and public repositories
- Module 07 (Protein Structure) — 3D data structures and visualization
- Tier 1 numpy/matplotlib for matrix visualization
- Module 11 (ChIP-seq) — familiarity with genomic interval formats (helpful)

## Data

Public Hi-C contact matrices from the [4DN Data Portal](https://data.4dnucleome.org/) or [ENCODE](https://www.encodeproject.org/) (`.cool` / `.mcool` format). The notebook uses a synthetic demo cooler when real data is unavailable.

## Related Skill

`hic-analysis.md` — Quick reference tables, patterns, and pitfalls for cooler/cooltools

---

[<- Previous Module](../13_Computational_Genetics/) | [Back to Course Overview](../../README.md) | [Next Module ->](../15_Motif_Discovery/)
