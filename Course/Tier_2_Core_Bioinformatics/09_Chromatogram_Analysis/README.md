# 9. Chromatogram Analysis

**Tier 2: Core Bioinformatics**

Reading and interpreting Sanger sequencing chromatogram files with BioPython and matplotlib. Covers Phred quality scores, base calling, and building a consensus sequence from forward and reverse reads.

## Topics Covered

- How Sanger (dideoxy chain termination) sequencing works
- Reading `.ab1` chromatogram trace files with BioPython
- Plotting the four fluorescence channels with matplotlib
- Phred quality scores and base calling confidence
- Identifying mixed bases, heterozygous positions, and low-quality regions
- Quality trimming strategies
- Building a consensus from forward and reverse reads

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_chromatogram_analysis.ipynb](01_chromatogram_analysis.ipynb) | Parsing .ab1 files, plotting trace data, quality scoring, and consensus building |

## Prerequisites

- Tier 2 Module 2: BioPython Essentials
- Tier 1 Module 18: Data Visualization (matplotlib)

---

[<- Previous Module](../08_Nucleic_Acid_Structure/) | [Back to Course Overview](../../README.md) | [Next Module ->](../10_Sequence_Motifs_and_Domains/)
