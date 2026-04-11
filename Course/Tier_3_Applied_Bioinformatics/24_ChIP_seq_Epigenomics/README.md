# Module 24: ChIP-seq & Epigenomics

**Tier 3 — Applied Bioinformatics | Module 24**

ChIP-seq peak calling, differential binding analysis, and epigenomic data integration.

## What You'll Learn
- ChIP-seq experimental design and controls (Input, IgG)
- Read trimming, alignment (Bowtie2), and duplicate removal (Picard)
- Peak calling with MACS3 for transcription factors and histone marks
- Quality metrics: FRiP score, cross-correlation, IDR reproducibility
- Differential binding analysis with DiffBind
- Functional annotation of peaks with ChIPseeker
- Integration with RNA-seq (overlap enriched peaks with DEGs)
- Visualization: heatmaps with deepTools, genome browser tracks

## Prerequisites
- Module 01 (NGS Fundamentals) — alignment, SAM/BAM, QC
- Module 03 (RNA-seq Analysis) — differential analysis concepts
- Module 23 (TF Footprinting & ATAC-seq) — chromatin accessibility basics

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [ChIP-seq Processing Pipeline](01_chipseq_pipeline.ipynb) | Trimming → alignment → peak calling → QC |
| 2 | [Differential Binding & Annotation](02_differential_binding.ipynb) | DiffBind, ChIPseeker, GO enrichment |

## Key Tools
| Tool | Purpose |
|------|---------|
| MACS3 | Peak calling |
| deepTools | bamCompare, plotHeatmap, computeMatrix |
| DiffBind | Differential binding (R/Bioconductor) |
| ChIPseeker | Peak annotation to nearest gene |
| Bowtie2 | Short-read alignment |
| Picard | Duplicate marking |
| FastQC / MultiQC | Quality control |

## Resources
- [ENCODE ChIP-seq pipeline](https://www.encodeproject.org/chip-seq/)
- [Harvard HBC ChIP-seq training](https://hbctraining.github.io/Intro-ChIP-seq/)
- [Bioconductor ChIPseeker vignette](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)
- [Bioconductor DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf)
- [deepTools documentation](https://deeptools.readthedocs.io/)
- [MACS3 documentation](https://macs3-project.github.io/MACS/)

## Related Skill
`chipseq-epigenomics.md` *(planned)*

---

[← Previous: 3.23 TF Footprinting](../23_TF_Footprinting/) | [Course Overview](../../README.md) | [Next: 3.25 Long-Read Sequencing →](../25_Long_Read_Sequencing/)
