# 1. NGS Fundamentals

**Tier 3: Applied Bioinformatics**

The complete next-generation sequencing pipeline from raw reads to aligned BAM files. Covers how Illumina, PacBio, and Nanopore platforms work, quality control, read trimming, and SAM/BAM format.

## Topics Covered

- Illumina, PacBio SMRT, and Oxford Nanopore sequencing principles and error profiles
- FASTQ format: structure, quality encoding, parsing
- Quality control with FastQC: per-base quality, adapter content, duplication
- Read trimming with Trimmomatic and fastp
- SAM/BAM format: header, FLAG fields, CIGAR strings
- Alignment with BWA and HISAT2

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_ngs_fundamentals.ipynb](01_ngs_fundamentals.ipynb) | NGS platforms, FASTQ/SAM/BAM formats, QC, trimming, and alignment |

## Prerequisites

- Tier 2: Core Bioinformatics (especially BioPython Essentials)
- Tier 1 Module 7: File I/O

---

[<- Previous Module](../00_Skills_Check/) | [Back to Course Overview](../../README.md) | [Next Module ->](../02_Variant_Calling_and_SNP_Analysis/)
