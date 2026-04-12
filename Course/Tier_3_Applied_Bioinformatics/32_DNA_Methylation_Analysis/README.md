# Module 32: DNA Methylation Analysis

**Tier 3 — Applied Bioinformatics | Module 32**

Bisulfite sequencing, differentially methylated region (DMR) calling, and epigenetic clock analysis.

## What You'll Learn
- 5-methylcytosine biology: CpG islands, gene silencing, imprinting
- Bisulfite conversion chemistry; WGBS, RRBS, and EPIC array comparisons
- Bismark alignment pipeline: bisulfite index → alignment → methylation extraction
- Differential methylation: DMPs with methylKit logistic regression, DMRs with BSmooth
- Integrating methylation with gene expression (promoter hypermethylation vs silencing)
- Epigenetic clocks: Horvath 353-CpG clock, GrimAge, PhenoAge biological age
- Epigenetic age acceleration analysis in population cohorts

## Prerequisites
- Module 01 (NGS Fundamentals) — read trimming, alignment, BAM format
- Module 24 (ChIP-seq & Epigenomics) — chromatin and epigenetic concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [WGBS/RRBS Processing with Bismark](01_wgbs_bismark.ipynb) | Bisulfite index, alignment, methylation extraction |
| 2 | [Differential Methylation (DMRs)](02_dmr_analysis.ipynb) | methylKit DMPs, BSmooth DMRs, genomic annotation |
| 3 | [Epigenetic Clocks & Aging](03_epigenetic_clocks.ipynb) | Horvath clock, GrimAge, age acceleration |

## Key Tools
| Tool | Purpose |
|------|---------|
| Bismark | Bisulfite alignment and methylation extraction |
| methylKit | Differential methylation testing (R) |
| DSS | Dispersion shrinkage DMR calling (R) |
| minfi | EPIC/450K array analysis (R) |
| methylclock | Epigenetic age estimation (R) |
| genomation | Annotation and visualization |

## Resources
- [Bismark documentation](https://felixkrueger.github.io/Bismark/)
- [methylKit vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/)
- [Horvath clock paper (2013)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115)
- [ENCODE WGBS pipeline](https://www.encodeproject.org/wgbs/)

## Related Skill
`dna-methylation.md` *(planned)*

---

[← Previous: 3.31 Single-Cell Multi-Omics](../31_Single_Cell_Multi_Omics/) | [Course Overview](../../README.md) | [Next: 3.33 CRISPR Screen Analysis →](../33_CRISPR_Screen_Analysis/)
