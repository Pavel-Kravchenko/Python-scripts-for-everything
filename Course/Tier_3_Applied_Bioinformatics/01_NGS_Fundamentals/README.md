# 1. NGS Fundamentals

**Tier 3: Applied Bioinformatics**

The complete next-generation sequencing pipeline from raw reads to aligned BAM files, plus a
comprehensive survey of every major bioinformatics data format you will encounter in real-world work.

## Topics Covered

### Notebook 1 – NGS Fundamentals
- Illumina, PacBio SMRT, and Oxford Nanopore sequencing principles and error profiles
- FASTQ format: structure, quality encoding, parsing
- Quality control with FastQC: per-base quality, adapter content, duplication
- Read trimming with Trimmomatic and fastp
- SAM/BAM format: header, FLAG fields, CIGAR strings
- Alignment with BWA and HISAT2

### Notebook 2 – Bioinformatics Data Formats
- **FASTA**: structure, indexing (.fai), random-access fetching
- **FASTQ**: quality encoding, Phred scores, paired-end reads
- **SAM / BAM / CRAM**: FLAG decoding, CIGAR parsing, compression trade-offs
- **VCF / BCF**: variant records, genotype fields, INFO/FORMAT parsing
- **BED**: 0-based intervals, overlap detection, distance calculation
- **GFF3 / GTF**: gene annotations, attribute parsing, hierarchy
- **WIG / BedGraph / BigWig / BigBed**: coverage tracks, binary indexed formats
- **PDB**: fixed-width ATOM/HETATM parsing, backbone extraction, distance calculation
- **mmCIF / PDBx**: modern structure format, AlphaFold pLDDT confidence scores
- **FAST5 / POD5**: Oxford Nanopore raw signal, basecalling, format conversion
- **Newick / Nexus / NHX**: phylogenetic tree serialization, recursive parsing

## Notebooks

| Notebook | Description | Cells |
|----------|-------------|-------|
| [01_ngs_fundamentals.ipynb](01_ngs_fundamentals.ipynb) | NGS platforms, FASTQ/SAM/BAM formats, QC, trimming, and alignment | ~49 |
| [02_bio_data_formats.ipynb](02_bio_data_formats.ipynb) | Comprehensive guide to all NGS, structure, and phylogeny formats with hands-on parsing | ~35 |

## Prerequisites

- Tier 2: Core Bioinformatics (especially BioPython Essentials)
- Tier 1 Module 7: File I/O

---

[← Previous Module](../00_Skills_Check/) | [Back to Course Overview](../../README.md) | [Next Module →](../02_Variant_Calling_and_SNP_Analysis/)
