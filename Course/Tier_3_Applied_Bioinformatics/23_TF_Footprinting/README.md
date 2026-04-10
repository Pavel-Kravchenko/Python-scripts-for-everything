# Module 23: TF Footprinting & Chromatin Accessibility

**Tier 3 — Applied Bioinformatics | Module 23**

Transcription factor footprinting from ATAC-seq data — from fragment QC to accumulation plots.

## What You'll Learn

- ATAC-seq data quality assessment (fragment size distribution, nucleosomal ladder)
- Tn5 insertion bias and offset correction (+4/-5 bp)
- Computing insertion profiles around TF motif sites
- Footprint score calculation (flanking/central ratio)
- Genomic interval arithmetic with pybedtools (intersection, subtraction, slop)
- Accumulation (aggregate) plots: mean signal around genomic features
- Comparing footprint scores across conditions

## Prerequisites

- Module 01 (NGS Fundamentals) — read alignment, BAM format
- Module 15 (Motif Discovery) — PWM scanning to define motif sites
- Module 05 (Promoter Analysis) — regulatory elements background

## Data

Simulated ATAC-seq fragment data and insertion profiles. Real data: ENCODE ATAC-seq BAM files (public).

## Related Skills

- `motif-discovery.md` — PWM scanning for TF motif sites
- `ngs-variant-calling.md` — upstream ATAC-seq alignment pipeline
