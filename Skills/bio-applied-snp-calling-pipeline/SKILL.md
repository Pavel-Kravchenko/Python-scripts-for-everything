---
name: bio-applied-snp-calling-pipeline
description: "SNP calling pipeline: Trimmomatic, HISAT2/BWA-MEM2 alignment, samtools/bcftools variant calling, ANNOVAR annotation"
tool_type: python
primary_tool: Matplotlib
---

# SNP Calling Pipeline

## Pipeline Architecture

```
Raw FASTQ reads
    → Trimmomatic (TRAILING:20, MINLEN:50)
    → HISAT2/BWA-MEM2 (alignment)
    → samtools (SAM→BAM→sort→index→depth)
    → samtools mpileup → bcftools call -cv (raw VCF)
    → ANNOVAR (annotation: dbSNP138, refGene, 1000G, GWAS, ClinVar)
    → Final annotated report
```

## Tool Reference

| Tool | Version | Purpose |
|------|---------|---------|
| Trimmomatic | 0.36 | Read quality trimming |
| HISAT2 | 2.1.0 | Read alignment (RNA-seq aligner used with DNA flags) |
| SAMtools | >=1.9 | BAM manipulation & pileup |
| BCFtools | >=1.9 | Variant calling from pileup |
| ANNOVAR | latest | Variant annotation |

**For new DNA projects**: use **BWA-MEM2** instead of HISAT2.

```bash
conda create -n snp-pipeline python=3.10
conda install -c bioconda trimmomatic hisat2 samtools bcftools
```

## Step 1 — Read Trimming

```bash
java -jar Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
    input.fastq chr_outfile.fastq TRAILING:20 MINLEN:50
```

| Parameter | Meaning |
|-----------|---------|
| `TRAILING:20` | Remove trailing bases with Phred < 20 |
| `MINLEN:50` | Discard reads < 50 bp after trimming |

## Step 2 — Alignment

HISAT2 with DNA-mode flags (`--no-spliced-alignment --no-softclip`) behaves like Bowtie2 end-to-end mode.

| Flag | Purpose |
|------|---------|
| `--no-spliced-alignment` | DNA reads don't span introns |
| `--no-softclip` | Prevents masking true variants at read ends |

### Alignment QC thresholds
- Overall rate >= 90%: GOOD
- 70-90%: ACCEPTABLE (check contamination)
- < 70%: POOR (wrong reference or contamination)
- Unique alignment >= 80%: GOOD

## Pitfalls

- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive
- **HISAT2 for DNA**: disable spliced alignment and soft-clipping, otherwise false split-read alignments occur
- **BWA read groups**: always include `-R` flag when aligning — GATK requires read groups
