# Module 15: Motif Discovery

**Tier 2 — Core Bioinformatics | Module 15**

Quantitative motif analysis from position frequency matrices to enrichment testing.

## What You'll Learn

- Building PPM and PWM from aligned binding-site sequences
- Computing information content (IC) and KDIC quality metrics
- Generating IUPAC consensus sequences
- Scanning genomic sequences and selecting thresholds from score distributions
- Motif enrichment testing with Fisher's exact test and Benjamini-Hochberg FDR
- TomTom concept: matching against JASPAR 2024 and HOCOMOCO databases
- Pipeline design patterns: abstract interfaces and BED/FASTA I/O

## Prerequisites

- Module 10 (Sequence Motifs and Domains) — PWM basics and PROSITE
- Module 02 (BioPython Essentials) — sequence I/O
- Tier 1 numpy/pandas for matrix operations
- Modules 10–12 (ChIP-seq, ATAC-seq, peak calling) — helpful context

## Data

JASPAR 2024 public motif database; ENCODE ChIP-seq peaks (public). The notebook uses synthetic sequences for demonstration.

## Related Skill

`motif-discovery.md` — Quick reference tables, patterns, and pitfalls
