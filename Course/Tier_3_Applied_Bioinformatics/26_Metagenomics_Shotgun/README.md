# Module 26: Shotgun Metagenomics

**Tier 3 — Applied Bioinformatics | Module 26**

Whole-community sequencing analysis: taxonomic profiling, functional annotation, and comparative metagenomics. Extends Module 04 (Microbial Diversity / 16S) to shotgun approaches.

## What You'll Learn
- Shotgun vs 16S amplicon metagenomics: when to choose each
- Read QC and host decontamination (Bowtie2 against host genome)
- Taxonomic profiling with Kraken2 + Bracken; Bracken abundance re-estimation
- Functional annotation with HUMAnN3 (UniRef/MetaCyc pathways)
- De novo metagenomic assembly (MEGAHIT) and binning (MetaBAT2)
- Bin quality assessment (CheckM)
- Metagenome-assembled genome (MAG) annotation (Prokka)
- Alpha/beta diversity for shotgun data; statistical comparison (PERMANOVA)
- 16S-based analysis with QIIME2 (complementary approach)

## Prerequisites
- Module 01 (NGS Fundamentals) — read QC, alignment
- Module 04 (Microbial Diversity) — ecological diversity concepts
- Module 17 (Genome Assembly) — assembly graph concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Taxonomic Profiling](01_taxonomic_profiling.ipynb) | Kraken2, Bracken, Krona visualization |
| 2 | [Functional Annotation](02_functional_annotation.ipynb) | HUMAnN3 pathways, enzyme coverage |
| 3 | [Assembly & Binning](03_assembly_binning.ipynb) | MEGAHIT assembly, MetaBAT2 binning, CheckM, MAG annotation |
| 4 | [QIIME2 16S Workflow](04_qiime2_16s.ipynb) | Paired-end denoising, diversity metrics, taxonomy classification |

## Key Tools
| Tool | Purpose |
|------|---------|
| Kraken2 + Bracken | Taxonomic classification and abundance |
| HUMAnN3 | Functional pathway profiling |
| MEGAHIT | De novo metagenomic assembly |
| MetaBAT2 | Contig binning into MAGs |
| CheckM | Bin completeness / contamination |
| Prokka | MAG annotation |
| QIIME2 | 16S amplicon analysis |

## Resources
- [QIIME2 tutorials](https://docs.qiime2.org/2024.10/tutorials/) — Moving Pictures, 16S workflows
- [Galaxy Training Network — Metagenomics](https://training.galaxyproject.org/training-material/topics/metagenomics/)
- [Kraken2 documentation](https://github.com/DerrickWood/kraken2)
- [HUMAnN3 documentation](https://github.com/biobakery/humann)
- [MetaBAT2 documentation](https://bitbucket.org/berkeleylab/metabat)
- [bioBakery tutorials](https://github.com/biobakery/biobakery/wiki)

## Related Skill
`metagenomics-shotgun.md` *(planned)*

---

[← Previous: 3.25 Long-Read Sequencing](../25_Long_Read_Sequencing/) | [Course Overview](../../README.md) | [Next: 3.27 Multi-Omics Integration →](../27_Multi_Omics_Integration/)
