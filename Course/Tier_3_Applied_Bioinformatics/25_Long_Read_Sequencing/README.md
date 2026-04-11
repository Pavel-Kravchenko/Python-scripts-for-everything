# Module 25: Long-Read Sequencing

**Tier 3 — Applied Bioinformatics | Module 25**

Oxford Nanopore and PacBio sequencing: basecalling, assembly, variant calling, and transcript isoform analysis.

## What You'll Learn
- Long-read technology comparison: ONT (R9/R10 chemistry) vs PacBio (HiFi/CLR)
- Basecalling with Dorado (ONT) and quality assessment
- Long-read alignment with Minimap2
- De novo assembly: Flye for metagenomes and genomes; Hifiasm for HiFi data
- Structural variant detection (Sniffles2, PBSV)
- Full-length transcript isoform analysis (FLAMES, IsoSeq)
- Methylation / base modification detection from ONT signal
- Hybrid assembly strategies (combining short + long reads)

## Prerequisites
- Module 01 (NGS Fundamentals) — sequencing concepts, SAM/BAM, QC
- Module 17 (Genome Assembly & Advanced NGS) — assembly graph concepts

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [ONT Data Processing](01_ont_processing.ipynb) | Basecalling, QC (NanoStat), alignment, methylation |
| 2 | [Long-Read Assembly & SV](02_assembly_sv.ipynb) | Flye/Hifiasm assembly, Sniffles2 SVs, QUAST evaluation |
| 3 | [Isoform Analysis](03_isoform_analysis.ipynb) | Full-length transcripts, novel isoforms, differential isoform usage |

## Key Tools
| Tool | Purpose |
|------|---------|
| Dorado | ONT basecalling (replaces Guppy) |
| Minimap2 | Long-read alignment |
| Flye | De novo genome/metagenome assembly |
| Hifiasm | HiFi de novo assembly |
| Sniffles2 | Structural variant calling |
| NanoStat / NanoPlot | ONT QC |
| QUAST | Assembly quality evaluation |
| FLAMES / IsoSeq | Full-length isoform analysis |

## Resources
- [Oxford Nanopore community tutorials](https://community.nanoporetech.com/)
- [Awesome Nanopore — GenomicsAotearoa training](https://github.com/GenomicsAotearoa/Genomics-Aotearoa-Nanopore-training)
- [Bioinformatics Data Skills (Vince Buffalo)](https://www.oreilly.com/library/view/bioinformatics-data-skills/9781449367480/) — O'Reilly, covers long reads
- [Minimap2 documentation](https://github.com/lh3/minimap2)
- [Flye documentation](https://github.com/fenderglass/Flye)
- [Sniffles2 documentation](https://github.com/fritzsedlazeck/Sniffles)

## Related Skill
`long-read-sequencing.md` *(planned)*

---

[← Previous: 3.24 ChIP-seq & Epigenomics](../24_ChIP_seq_Epigenomics/) | [Course Overview](../../README.md) | [Next: 3.26 Shotgun Metagenomics →](../26_Metagenomics_Shotgun/)
