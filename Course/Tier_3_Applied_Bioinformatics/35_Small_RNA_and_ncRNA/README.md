# Module 35: Small RNA & Non-Coding RNA

**Tier 3 — Applied Bioinformatics | Module 35**

miRNA-seq analysis, lncRNA discovery, circular RNA detection, and regulatory ncRNA networks.

## What You'll Learn
- miRNA biogenesis: Drosha/DGCR8 → pre-miRNA → Dicer → RISC complex
- Small RNA-seq library design; 3' adapter trimming, size selection
- miRNA quantification with miRBase annotation and differential expression (DESeq2)
- miRNA target prediction: TargetScan, miRTarBase, CLASH
- ncRNA classification: lncRNA, circRNA, piRNA, snoRNA, snRNA
- De novo lncRNA discovery with StringTie and coding potential filtering
- Circular RNA detection from back-splice junctions (CIRI2, find_circ)
- Co-expression networks for lncRNA function inference (WGCNA)

## Prerequisites
- Module 01 (NGS Fundamentals) — short-read sequencing, alignment
- Module 03 (RNA-seq Analysis) — differential expression concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [miRNA-seq Pipeline](01_mirna_seq_pipeline.ipynb) | Trimming, alignment, miRBase, DESeq2, target prediction |
| 2 | [lncRNA & ncRNA Classification](02_lncrna_classification.ipynb) | StringTie assembly, CPC2, circRNA (CIRI2), co-expression |

## Key Tools
| Tool | Purpose |
|------|---------|
| Cutadapt | 3' adapter trimming for small RNA |
| miRDeep2 | miRNA quantification and novel discovery |
| featureCounts | miRNA count matrix |
| StringTie | Transcript assembly for lncRNA discovery |
| CPC2 / CPAT | Coding potential assessment |
| CIRI2 | Circular RNA detection |
| WGCNA | Co-expression network analysis |

## Resources
- [miRBase database](https://www.mirbase.org/)
- [LNCipedia database](https://lncipedia.org/)
- [circBase](http://www.circbase.org/)
- [GENCODE lncRNA annotation](https://www.gencodegenes.org/)
- [miRDeep2 documentation](https://github.com/rajewsky-lab/mirdeep2)

## Related Skill
`small-rna-ncrna.md` *(planned)*

---

[← Previous: 3.34 Immunogenomics](../34_Immunogenomics/) | [Course Overview](../../README.md) | [Next: 3.36 Metabolomics →](../36_Metabolomics/)
