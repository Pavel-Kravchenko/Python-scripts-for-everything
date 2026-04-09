# 3.17 Genome Assembly and Advanced NGS
**Tier 3: Applied Bioinformatics**

This module builds directly on [3.01 NGS Fundamentals](../01_NGS_Fundamentals/) — which covered sequencing platforms, FASTQ format, basic QC with FastQC, adapter trimming, SAM/BAM format, and alignment with BWA/HISAT2. Here we go deeper: assembling genomes *de novo* when no reference exists, understanding the mathematical foundations of read mapping (BWT, FM-index), assessing assembly quality rigorously, and working with long-read technologies (Oxford Nanopore, PacBio HiFi) that are now the standard for high-quality genome projects.

## Topics Covered

- **Why de novo assembly**: when reference mapping is not possible or insufficient
- **Overlap-Layout-Consensus (OLC)**: the classical approach for long reads; Hamiltonian path problem
- **De Bruijn graph assembly**: k-mer decomposition, Eulerian path, k-mer size trade-offs — implemented from scratch in Python
- **Real assemblers**: SPAdes (short reads, multi-k), Flye, Canu (long reads), hifiasm (HiFi, haplotype-resolved)
- **Assembly quality metrics**: N50, L50, NG50, N90 — implemented from scratch
- **BUSCO**: gene-completeness benchmarking with universal single-copy orthologs
- **QUAST**: comprehensive assembly statistics and reference-based evaluation
- **Scaffolding**: paired-end, mate-pair, Hi-C, and long-read scaffolding; GFA assembly graphs
- **Gap filling**: TGS-GapCloser, LR_Gapcloser
- **Burrows-Wheeler Transform (BWT)**: construction, invertibility, compression property — implemented in Python
- **FM-index and backward search**: O(m) exact pattern matching using LF-mapping — implemented in Python
- **BWA-MEM internals**: SMEM seeding, seed chaining, banded Smith-Waterman extension, MAPQ
- **MultiQC**: aggregated QC reports across many samples
- **K-mer frequency histograms**: GenomeScope for genome size and heterozygosity estimation
- **Contamination detection**: Kraken2, blobtools, NCBI FCS-GX
- **Coverage depth analysis**: identifying assembly errors, collapsed repeats, and gaps
- **Oxford Nanopore Technology**: basecalling (Guppy/Dorado), error profiles, adaptive sampling, ultra-long reads
- **PacBio HiFi (CCS)**: circular consensus sequencing, Q30+ accuracy, comparison with CLR
- **Hybrid assembly strategies**: combining short reads, HiFi, ONT, and Hi-C for optimal results
- **Polishing**: Medaka, Pilon, Racon for error correction of long-read draft assemblies

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_genome_assembly_and_advanced_ngs.ipynb](01_genome_assembly_and_advanced_ngs.ipynb) | Complete coverage of de novo assembly algorithms, BWT/FM-index, assembly QC, scaffolding, long-read technologies, and advanced QC methods. Includes Python implementations of OLC assembler, de Bruijn graph assembler, N50/L50/NG50 calculator, BWT construction, and FM-index backward search. |

## Prerequisites

- **3.01 NGS Fundamentals** — sequencing platforms, FASTQ, FastQC, trimming, SAM/BAM, BWA/HISAT2
- **Tier 4 String Algorithms** (recommended for BWT/FM-index theory) — suffix arrays, formal BWT construction algorithms, FM-index data structure proofs

---
[← Previous: 3.16 Numerical Methods](../16_Numerical_Methods_for_Bioinformatics/) | [Course Overview](../../README.md) | [Next: 3.18 Proteomics →](../18_Proteomics_and_Structural_Methods/)
