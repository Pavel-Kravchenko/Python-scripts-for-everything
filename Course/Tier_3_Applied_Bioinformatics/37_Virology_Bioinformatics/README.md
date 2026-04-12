# Module 37: Virology Bioinformatics

**Tier 3 — Applied Bioinformatics | Module 37**

Viral genome assembly, intra-host variation, phylodynamics, and real-time genomic surveillance.

## What You'll Learn
- Viral sequencing strategies: ARTIC amplicon, metagenomic, nanopore rapid sequencing
- Reference-guided viral genome assembly with iVar (SARS-CoV-2, influenza)
- Intra-host variant calling at < 5% allele frequency with LoFreq
- Viral genome annotation with VADR and Nextclade clade assignment
- Molecular clock models and timed phylogenetic trees with Nextstrain / Augur
- BEAST2 for Bayesian coalescent effective population size estimation
- Wastewater-based epidemiology: variant deconvolution with Freyja
- Recombination detection and lineage nomenclature

## Prerequisites
- Module 01 (NGS Fundamentals) — alignment, SAM/BAM, variant calling
- Module 06 (Phylogenetics) — phylogenetic methods, tree interpretation
- Module 17 (Genome Assembly) — assembly graph concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Viral Genome Assembly & Variant Analysis](01_viral_genome_assembly.ipynb) | iVar, LoFreq, VADR, Nextclade |
| 2 | [Phylodynamics & Molecular Epidemiology](02_phylodynamics.ipynb) | Nextstrain Augur, BEAST2, transmission clusters |
| 3 | [Variant Surveillance & Wastewater](03_variant_surveillance.ipynb) | Freyja, real-time dashboards, functional annotation |

## Key Tools
| Tool | Purpose |
|------|---------|
| iVar | Amplicon trimming, consensus, variants |
| LoFreq | Low-frequency variant calling |
| Minimap2 | Long + short read alignment |
| Nextstrain / Augur | Timed phylogenetics pipeline |
| TreeTime | Molecular clock dating |
| BEAST2 | Bayesian phylodynamics |
| Freyja | Wastewater variant deconvolution |
| Nextclade | Real-time clade assignment |
| pangolin | Pango lineage classification |

## Resources
- [Nextstrain documentation](https://docs.nextstrain.org/)
- [ARTIC Network protocols](https://artic.network/)
- [BEAST2 tutorials](https://www.beast2.org/tutorials/)
- [GISAID EpiCoV](https://www.gisaid.org/)
- [Freyja documentation](https://github.com/andersen-lab/Freyja)

## Related Skill
`virology-bioinformatics.md` *(planned)*

---

[← Previous: 3.36 Metabolomics](../36_Metabolomics/) | [Course Overview](../../README.md) | [Next: Tier 5 →](../../Tier_5_Modern_AI_for_Science/)
