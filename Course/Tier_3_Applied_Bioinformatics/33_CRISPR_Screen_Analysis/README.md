# Module 33: CRISPR Screen Analysis

**Tier 3 — Applied Bioinformatics | Module 33**

Pooled CRISPR library screens for gene essentiality, drug resistance, and regulatory element discovery.

## What You'll Learn
- Pooled screen design: sgRNA library selection, lentiviral delivery, MOI control
- sgRNA counting with MAGeCK count from FASTQ
- Gene essentiality scoring with MAGeCK RRA algorithm (robust rank aggregation)
- Hit calling: log fold-change, p-value, FDR thresholds for positive and negative selection
- Screen quality metrics: Gini index, read depth, replicate correlation
- Copy-number bias correction with CRISPRcleanR
- CRISPRi / CRISPRa screens for regulatory element discovery
- Integration with DepMap, BioGRID, and drug-target databases

## Prerequisites
- Module 01 (NGS Fundamentals) — read processing and alignment concepts
- Module 14 (Genetic Engineering In Silico) — CRISPR guide RNA design

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [MAGeCK Gene Essentiality](01_mageck_gene_essentiality.ipynb) | sgRNA counting, RRA algorithm, volcano plots |
| 2 | [Screen QC & Advanced Methods](02_screen_qc_normalization.ipynb) | Quality metrics, CRISPRcleanR, CRISPRi/a screens |

## Key Tools
| Tool | Purpose |
|------|---------|
| MAGeCK | Counting + robust rank aggregation |
| MAGeCK-VISPR | Interactive result visualization |
| CRISPRcleanR | Copy-number bias correction |
| BAGEL2 | Bayesian gene essentiality scoring |

## Resources
- [MAGeCK documentation](https://sourceforge.net/p/mageck/wiki/Home/)
- [DepMap portal](https://depmap.org/portal/)
- [Addgene sgRNA library guide](https://www.addgene.org/crispr/libraries/)
- [CRISPRcleanR paper (Iorio et al. 2018)](https://www.nature.com/articles/s41598-018-29pompier)

## Related Skill
`crispr-screen-analysis.md` *(planned)*

---

[← Previous: 3.32 DNA Methylation Analysis](../32_DNA_Methylation_Analysis/) | [Course Overview](../../README.md) | [Next: 3.34 Immunogenomics →](../34_Immunogenomics/)
