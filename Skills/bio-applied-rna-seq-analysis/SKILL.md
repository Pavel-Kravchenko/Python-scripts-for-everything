---
name: bio-applied-rna-seq-analysis
description: RNA-seq workflow — experimental design, alignment vs pseudo-alignment, count normalization (TPM), differential expression setup, and key pitfalls
tool_type: python
primary_tool: Matplotlib
---

# RNA-seq Analysis: Reads to Differential Expression

## Workflow Overview

```
FASTQ → QC (FastQC/fastp) → Alignment (STAR/HISAT2) → Counting (featureCounts)
                          → Pseudo-align (Salmon/kallisto) → tximport
                          → Count matrix → Normalization → DESeq2/edgeR → GSEA
```

## Experimental Design Rules

- **Replicates**: ≥3 biological replicates per condition; more replicates > deeper sequencing
- **Depth**: 10–30M reads for gene-level DE; 50–100M for isoform discovery
- **Read type**: SE 50–75 bp for gene-level; PE 100–150 bp for splicing/isoforms
- **Strandedness**: use strand-specific protocol (dUTP) — necessary for overlapping genes
- **Batches**: balance conditions across batches; record batch for correction (limma `removeBatchEffect` or DESeq2 design formula)

## Alignment Tools

| Tool | Memory | Speed | Notes |
|------|--------|-------|-------|
| STAR | 30+ GB | Very fast | Most used for human/mouse; outputs GeneCounts |
| HISAT2 | ~8 GB | Fast | Lower memory; successor to TopHat2 |

```bash
# STAR index (once)
STAR --runMode genomeGenerate --genomeDir star_index/ \
     --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf

# STAR align
STAR --genomeDir star_index/ --readFilesIn R1.fastq R2.fastq \
     --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

# featureCounts
featureCounts -a genes.gtf -o counts.txt -T 4 -p --countReadPairs \
              sample1.bam sample2.bam sample3.bam
```

## Pseudo-alignment (Salmon/kallisto)

10–100× faster than alignment-based; comparable accuracy for gene-level DE.

```bash
# Salmon
salmon index -t transcriptome.fa -i salmon_index
salmon quant -i salmon_index -l A -1 R1.fastq -2 R2.fastq \
             -o sample_quant --validateMappings
```

Outputs TPM per transcript → use `tximport` (R) or sum per gene (Python) before DESeq2.

## Normalization Units

| Unit | Formula | Cross-sample comparable? | Use for DE? |
|------|---------|--------------------------|------------|
| Raw counts | — | No | Yes (DESeq2/edgeR) |
| RPKM/FPKM | counts / (len_kb × total_M) | No | No |
| TPM | (counts/len) / sum(counts/len) × 1e6 | Yes (sums to 1M) | No |

```python
def counts_to_tpm(counts, lengths_bp):
    """Convert raw counts array to TPM. lengths_bp in base pairs."""
    rate = counts / lengths_bp
    return rate / rate.sum() * 1e6
```

**Rule**: use raw counts for DE testing; use TPM only for visualization/reporting.

## Count Matrix

- Shape: (n_genes × n_samples), integer values
- Distribution: negative binomial (overdispersed Poisson) — DESeq2/edgeR model this directly
- Zeros are common for lowly expressed genes — do not impute

```python
import pandas as pd
import numpy as np

# Load featureCounts output
counts = pd.read_csv('counts.txt', sep='\t', comment='#', index_col=0)
# Drop annotation columns, keep sample columns
counts = counts.iloc[:, 5:]  # featureCounts has 5 meta-columns before samples
```

## DE Analysis Quick Setup (pydeseq2)

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.stats import DeseqStats

dds = DeseqDataSet(counts=counts_df, metadata=sample_info, design_factors="condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition", "Treatment", "Control"])
stat_res.summary()
res = stat_res.results_df  # log2FoldChange, pvalue, padj
```

## Pitfalls

- **Coordinate systems**: BED is 0-based; GTF/GFF are 1-based — off-by-one errors when mixing
- **Batch effects**: PCA or hierarchical clustering before DE to detect; correct with design formula or `removeBatchEffect` — never correct counts and feed to DESeq2
- **Multiple testing**: always use `padj` (Benjamini-Hochberg); never filter on raw p-value alone
- **Library size confounding**: do not use TPM/FPKM as input to DESeq2/edgeR — they have already normalized in a way that breaks count-based models
- **Poly-A vs rRNA depletion**: poly-A misses non-coding RNA and degraded RNA; use rRNA depletion for FFPE or total RNA experiments
- **Strandedness mismatch**: wrong `--strandedness` in featureCounts halves your counts; always verify with RSeQC `infer_experiment.py`
- **Outlier samples**: one bad sample can dominate DE results; always run PCA/hclust QC before DE
