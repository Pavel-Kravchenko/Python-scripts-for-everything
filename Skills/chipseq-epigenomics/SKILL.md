---
name: chipseq-epigenomics
description: ChIP-seq processing pipeline, peak calling with MACS3, differential binding with DiffBind, peak annotation with ChIPseeker, deepTools visualization
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# ChIP-seq & Epigenomics

## When to Use

Use this skill when:
- Processing ChIP-seq data (TF binding or histone modifications)
- Calling narrow peaks (TF) or broad peaks (H3K4me3, H3K27ac, H3K27me3)
- Running differential binding analysis between conditions
- Annotating peaks to genomic features (promoter, intron, intergenic)
- Visualizing signal enrichment with heatmaps and profile plots

## Quick Reference

| Task | Tool | Key Command / Function |
|------|------|----------------------|
| Read trimming | Trim Galore | `trim_galore --paired` |
| Alignment | Bowtie2 | `bowtie2 -x index -1 R1 -2 R2` |
| Duplicate removal | Picard | `MarkDuplicates REMOVE_DUPLICATES=true` |
| Peak calling (TF) | MACS3 | `macs3 callpeak -f BAMPE -q 0.05` |
| Peak calling (histone) | MACS3 | `macs3 callpeak --broad --broad-cutoff 0.1` |
| FRiP score | deepTools | `featureCounts -a peaks.narrowPeak -F SAF` |
| BigWig normalization | deepTools | `bamCoverage --normalizeUsing RPKM` |
| TSS heatmap | deepTools | `computeMatrix reference-point` + `plotHeatmap` |
| Differential binding | DiffBind (R) | `dba.analyze(dba)` |
| Peak annotation | ChIPseeker (R) | `annotatePeak(peaks, TxDb=txdb)` |

## Key Patterns

**Pattern 1: Full ChIP-seq pipeline**
```bash
# Trim adapters
trim_galore --paired --fastqc sample_R1.fq.gz sample_R2.fq.gz -o trimmed/

# Align with Bowtie2
bowtie2 -x hg38 -1 trimmed/R1 -2 trimmed/R2 | samtools sort -o sample.bam
samtools index sample.bam

# Remove duplicates
picard MarkDuplicates I=sample.bam O=sample_dedup.bam M=metrics.txt REMOVE_DUPLICATES=true

# Call narrow peaks (TF)
macs3 callpeak -t sample_dedup.bam -c input_dedup.bam \
    -f BAMPE -g hs -n sample --outdir peaks/ -q 0.05
```

**Pattern 2: deepTools visualization**
```bash
# Normalize to BigWig
bamCoverage -b sample_dedup.bam -o sample.bw --normalizeUsing RPKM --binSize 10

# Heatmap at TSS (±3 kb)
computeMatrix reference-point -S sample.bw -R genes.bed \
    --referencePoint TSS -b 3000 -a 3000 -o matrix.gz
plotHeatmap -m matrix.gz -out heatmap.png --colorMap RdYlBu_r
```

**Pattern 3: DiffBind differential binding (R)**
```r
library(DiffBind)
samples <- read.csv('samplesheet.csv')  # columns: SampleID, Condition, BAM, Peaks
dba <- dba(sampleSheet=samples)
dba <- dba.count(dba)
dba <- dba.normalize(dba)
dba <- dba.contrast(dba, categories=DBA_CONDITION, minMembers=2)
dba <- dba.analyze(dba)
db_peaks <- dba.report(dba)
dba.plotVolcano(dba)
```

**Pattern 4: ChIPseeker peak annotation (R)**
```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile('peaks/sample_peaks.narrowPeak')
anno <- annotatePeak(peaks, tssRegion=c(-2000, 200),
                     TxDb=txdb, annoDb='org.Hs.eg.db')
plotAnnoPie(anno)
plotDistToTSS(anno)

# Get promoter peaks (within 2 kb of TSS)
promo_peaks <- as.data.frame(anno)[abs(as.data.frame(anno)$distanceToTSS) < 2000, ]
```

## Quality Metrics

| Metric | Recommendation | How to Compute |
|--------|---------------|----------------|
| FRiP (TF) | ≥ 5% | reads in peaks / total mapped |
| FRiP (histone) | ≥ 1% | reads in peaks / total mapped |
| NSC (normalized strand coefficient) | > 1.05 | SPP / deepTools |
| RSC (relative strand correlation) | > 0.8 | SPP |
| % duplicates | < 30% (for TF) | Picard metrics |

## Common Pitfalls

- **Input control is mandatory** — always subtract input for peak calling; IgG is an alternative
- **Paired-end vs single-end** — use `-f BAMPE` for PE and `-f BAM` for SE in MACS3
- **Broad vs narrow peaks** — H3K4me3/H3K27ac = narrow; H3K9me3/H3K27me3 = broad
- **Blacklist regions** — filter ENCODE blacklisted regions before analysis
- **Normalization** — use spike-in normalization (Orlando method) for histone mark quantification across conditions

## ENCODE Standards

```bash
# Filter blacklisted regions
bedtools intersect -v -a peaks.narrowPeak -b encode_blacklist.bed > peaks_filtered.bed

# Minimum recommended depth
# TF ChIP-seq: 20M uniquely mapped reads
# Histone ChIP-seq: 40-80M uniquely mapped reads
```

## Code Templates

### Parse narrowPeak File
```python
import pandas as pd

def read_narrowpeak(path):
    cols = ['chrom','start','end','name','score','strand',
            'signal','pvalue','qvalue','peak']
    df = pd.read_csv(path, sep='\t', header=None, names=cols)
    return df

peaks = read_narrowpeak('sample_peaks.narrowPeak')
# Filter by q-value (column is -log10 q-value; 1.3 ≈ q<0.05)
sig_peaks = peaks[peaks['qvalue'] >= 1.3]
print(f"{len(sig_peaks)} significant peaks")
```

### Overlap Two Peak Sets
```python
import pybedtools

a = pybedtools.BedTool('condition_A.narrowPeak')
b = pybedtools.BedTool('condition_B.narrowPeak')

# Peaks shared in both conditions (reciprocal 50% overlap)
shared = a.intersect(b, f=0.5, r=True)
a_only  = a.intersect(b, f=0.5, r=True, v=True)
b_only  = b.intersect(a, f=0.5, r=True, v=True)
print(f"Shared: {len(shared)}, A-only: {len(a_only)}, B-only: {len(b_only)}")
```

### Compute FRiP Score
```python
import subprocess

def frip_score(bam, peaks_bed):
    """Fraction of Reads in Peaks."""
    total = int(subprocess.check_output(
        ['samtools', 'view', '-c', '-F', '4', bam]).decode().strip())
    in_peaks = int(subprocess.check_output(
        ['bedtools', 'intersect', '-a', bam, '-b', peaks_bed,
         '-u', '-f', '0.5']).decode().count('\n'))
    # Faster with featureCounts; this is a quick estimate
    return in_peaks / total

frip = frip_score('sample_dedup.bam', 'sample_peaks.narrowPeak')
print(f"FRiP = {frip:.3f}  ({'PASS' if frip >= 0.05 else 'FAIL'} for TF)")
```

### Load BigWig Signal in Python
```python
import pyBigWig
import numpy as np

bw = pyBigWig.open('sample.bw')
# Mean signal over a region
signal = bw.stats('chr1', 1_000_000, 1_100_000, type='mean', nBins=100)
signal = np.array([v if v is not None else 0 for v in signal])
bw.close()
```

## Related Skills
- `ngs-variant-calling` — FASTQ QC, alignment with BWA, samtools operations
- `rnaseq` — differential expression, DESeq2 normalization
- `tf-footprinting-atac` — ATAC-seq, open chromatin, TF footprinting
- `data-visualization-bio` — heatmaps, genome browser tracks, volcano plots
