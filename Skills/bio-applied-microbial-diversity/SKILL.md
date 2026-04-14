---
name: bio-applied-microbial-diversity
description: "Microbial diversity analysis: alpha/beta diversity metrics, OTU/ASV methods, taxonomy assignment, and community comparison. Use when analyzing 16S amplicon or microbiome data."
tool_type: python
primary_tool: Matplotlib
---

# Microbial Diversity Analysis

## 16S rRNA Primer Targets

| Region | Primers | Length | Notes |
|--------|---------|--------|-------|
| V4 | 515F/806R | ~250 bp | Earth Microbiome Project standard |
| V3-V4 | 341F/785R | ~460 bp | Good taxonomic resolution |
| V1-V3 | 27F/534R | ~500 bp | Alternative target |

Other markers: 18S (eukaryotes), ITS (fungi), shotgun metagenomics (full functional + taxonomic profiling).

## OTU vs ASV

| Feature | OTUs (97%) | ASVs |
|---------|-----------|------|
| Resolution | Species-level (approx) | Single-nucleotide |
| Reproducibility | Depends on algorithm | Consistent across studies |
| Tools | UCLUST, VSEARCH, CD-HIT | DADA2, Deblur, UNOISE3 |

**ASVs are now recommended** (Callahan et al., 2017). DADA2 learns error rates per nucleotide transition, tests whether each unique sequence could be explained as an error from a more abundant sequence.

## Taxonomy Assignment

| Method | Approach |
|--------|----------|
| Naive Bayes classifier | k-mer frequencies (QIIME2 classify-sklearn) |
| BLAST/VSEARCH | Sequence similarity top hit/consensus |
| Phylogenetic placement | pplacer, EPA-ng |

Reference databases: SILVA (~2M seqs, most comprehensive), Greengenes2 (genome phylogeny), RDP, GTDB (genome-based gold standard).

## Alpha Diversity

Measures diversity within a single sample:
- **Richness**: observed species count
- **Evenness**: how evenly distributed (Shannon, Simpson)
- **Phylogenetic**: Faith's PD incorporates evolutionary tree

Qualitative (presence/absence) vs quantitative (abundance-weighted).

## Feature Table

Central data structure: rows = OTUs/ASVs, columns = samples, values = read counts. Simulate with Dirichlet-multinomial model per body site.

## Pitfalls

- **Rarefaction depth**: Unequal library sizes bias richness; rarefy or use normalized metrics
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
