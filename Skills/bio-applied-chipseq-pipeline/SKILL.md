---
name: bio-applied-chipseq-pipeline
description: "ChIP-seq pipeline — QC, alignment, deduplication, peak calling with MACS2, and signal normalization with deepTools"
tool_type: python
primary_tool: Matplotlib
---

# ChIP-seq Processing Pipeline

- [ENCODE ChIP-seq pipeline](https://www.encodeproject.org/chip-seq/)
- [MACS3 documentation](https://macs3-project.github.io/MACS/)
- [deepTools documentation](https://deeptools.readthedocs.io/)

## Experiment Design Reference

| Experiment Type | Peak Shape | Recommended Reads | Example Marks |
|---|---|---|---|
| TF | Narrow (< 500 bp) | 20-40 M | CTCF, GATA1, p53 |
| Active histone | Narrow | 30-50 M | H3K4me3, H3K9ac |
| Broad histone | Broad (> 5 kb) | 40-80 M | H3K27me3, H3K9me3 |
| Enhancer mark | Narrow/Mixed | 30-50 M | H3K27ac, H3K4me1 |

**Antibody**: must be ChIP-grade (not just WB-validated). **Replicates**: ENCODE requires >= 2 biological replicates + IDR analysis. **Sequencing**: paired-end strongly preferred (accurate fragment size, better dedup).

## Step 1: QC and Trimming

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc/raw/ -t 4
multiqc qc/raw/ -o qc/raw/multiqc/

# Trim Galore (auto-detects adapters)
trim_galore --paired --fastqc --cores 4 --quality 20 --length 20 \
    --output_dir trimmed/ sample_R1.fastq.gz sample_R2.fastq.gz

# Alternative: fastp (faster)
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o trimmed/sample_R1_trimmed.fastq.gz -O trimmed/sample_R2_trimmed.fastq.gz \
    --detect_adapter_for_pe --thread 4 -h trimmed/fastp_report.html
```

**Warning signs**: adapter > 20% = short fragments (re-evaluate sonication); Q < 20 beyond 10bp at 3' = aggressive trim; bimodal GC = contamination; > 30% pairs discarded = re-prep library.

## Step 2: Alignment (Bowtie2)

```bash
bowtie2 -x /data/indices/hg38 \
    -1 trimmed/sample_R1_val_1.fq.gz -2 trimmed/sample_R2_val_2.fq.gz \
    -p 8 --no-mixed --no-discordant \
    2> logs/sample_bowtie2.log | \
    samtools view -bS -F 4 -F 256 -q 30 -f 2 | \
    samtools sort -@ 4 -o aligned/sample.bam

samtools index aligned/sample.bam
samtools flagstat aligned/sample.bam | tee logs/sample_flagstat.txt

# Remove ENCODE blacklisted regions
bedtools intersect -a aligned/sample.bam -b hg38-blacklist.v2.bed -v > aligned/sample_filtered.bam
samtools index aligned/sample_filtered.bam
```

Target: > 80% alignment rate. Filter flags: `-F 4` (unmapped), `-F 256` (secondary), `-q 30` (MAPQ), `-f 2` (properly paired).

## Step 3: Duplicate Removal (Picard)

```bash
picard MarkDuplicates \
    I=aligned/sample_filtered.bam O=dedup/sample_dedup.bam \
    M=dedup/sample_dup_metrics.txt REMOVE_DUPLICATES=true \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \  # 100 for HiSeq, 2500 for NovaSeq
    VALIDATION_STRINGENCY=SILENT

samtools index dedup/sample_dedup.bam
```

**ENCODE QC thresholds**: PBC1 > 0.7 (acceptable), > 0.9 (ideal). NRF > 0.8 (acceptable), > 0.9 (ideal). Total dup rate > 50% = library complexity problem; TF ChIP should be < 30%.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
