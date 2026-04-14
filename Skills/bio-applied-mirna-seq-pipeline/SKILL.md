---
name: bio-applied-mirna-seq-pipeline
description: "miRNA-seq pipeline: adapter trimming, alignment to miRBase, quantification, DE analysis, and target prediction. Use when processing small RNA sequencing data."
tool_type: python
primary_tool: Python
---

# miRNA-seq Processing and Analysis

- [miRBase database](https://www.mirbase.org/)
- [miRDeep2 documentation](https://github.com/rajewsky-lab/mirdeep2)
- [TargetScan](https://www.targetscan.org/)
- [mirTarBase](https://mirtarbase.cuhk.edu.cn/)

## miRNA Naming Convention

- `miR-21` / `miR-21-5p`: mature guide, 5' arm
- `miR-21-3p`: 3' arm (was miR-21*)
- `hsa-miR-21-5p`: human, 5' arm
- `mir-21` (lowercase): the gene/precursor

Seed region: positions 2-7 of mature miRNA, critical for target recognition.

## Library Design

Standard RNA-seq fails (miRNAs lack poly-A, are already short). Use 3' adapter ligation-based libraries. 50 bp single-end sequencing is sufficient.

## Processing Pipeline

### Cutadapt Trimming
```bash
cutadapt \
  -a TGGAATTCTCGGGTGCCAAGG \  # 3' adapter
  -m 16 -M 28 \                # length filter
  --discard-untrimmed \
  -j 8 -o trimmed.fastq.gz sample.fastq.gz
```

### Bowtie Alignment to miRBase
```bash
bowtie-build mature.fa mirbase_index
bowtie -x mirbase_index --norc -v 1 -m 5 -p 8 trimmed.fastq.gz -S aligned.sam
```

### Quantification
```bash
featureCounts -a hg38_mirbase_v22.gff3 -o raw_counts.txt -t miRNA -g Name -s 1 aligned.bam
```

## QC Thresholds

| Metric | Acceptable |
|--------|-----------|
| Total raw reads | >5M per sample |
| Adapter trimming rate | >60% |
| Alignment rate to miRBase | >50% |
| Top expressed miRNA | <50% total |

## Differential Expression

DESeq2 (same NB model as mRNA-seq). Normalization: TMM or DESeq2 median-of-ratios preferred over RPM/CPM.

```python
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# CPM normalize, log2 transform, t-test per miRNA, BH correction
reject, padj, _, _ = multipletests(pvalues, method='fdr_bh')
significant = (padj < 0.05) & (abs(log2fc) > 1.0)
```

## Target Prediction

### Seed Match Types
| Site type | Definition | Repression |
|-----------|-----------|-----------|
| 8mer | Seed + match pos 8 + A at pos 1 | Strongest |
| 7mer-m8 | Seed + match at position 8 | Strong |
| 7mer-A1 | Seed + A at position 1 | Moderate |
| 6mer | Seed only | Weak |

### Databases
- **TargetScan 8.0**: context++ score (seed type, accessibility, conservation)
- **miRDB**: ML-based MirTarget score
- **miRTarBase**: experimentally validated (CLASH, luciferase)

## Pitfalls

- **Adapter contamination**: Most reads must contain adapter; discard untrimmed reads
- **isomiR variation**: 5'/3' end variants are common; decide whether to collapse or keep
- **Multiple testing**: Apply FDR correction when testing hundreds of miRNAs
