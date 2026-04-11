---
name: long-read-sequencing
description: Oxford Nanopore and PacBio long-read sequencing — basecalling (Dorado), QC (NanoStat), alignment (Minimap2), assembly (Flye, Hifiasm), SV calling (Sniffles2), methylation, isoform analysis
---

## When to Use

Use this skill when:
- Processing Oxford Nanopore (ONT) or PacBio HiFi reads
- Performing de novo genome or metagenome assembly
- Detecting structural variants (deletions, insertions, inversions, duplications)
- Analyzing CpG methylation from ONT signal
- Performing full-length transcript isoform analysis

## Quick Reference

| Task | Tool | Key Command |
|------|------|-------------|
| ONT basecalling | Dorado | `dorado basecaller hac pod5/` |
| Read QC | NanoStat | `NanoStat --fastq reads.fq.gz` |
| Read filtering | NanoFilt | `NanoFilt -q 10 -l 1000` |
| Alignment | Minimap2 | `minimap2 -ax map-ont ref.fa reads.fq.gz` |
| HiFi alignment | Minimap2 | `minimap2 -ax map-hifi ref.fa reads.fq.gz` |
| cDNA alignment | Minimap2 | `minimap2 -ax splice reads.fq.gz` |
| De novo assembly (ONT) | Flye | `flye --nano-hq reads.fq.gz --genome-size 5m` |
| De novo assembly (HiFi) | Hifiasm | `hifiasm -o out.asm -t 16 reads.fq.gz` |
| Assembly polishing | Medaka | `medaka_consensus -i reads.fq.gz -d assembly.fa` |
| Assembly QC | QUAST | `quast.py assembly.fa -r ref.fa` |
| SV calling | Sniffles2 | `sniffles --input aligned.bam --vcf svs.vcf` |
| Methylation | modkit | `modkit pileup aligned.bam meth.bed --cpg` |
| Isoform analysis | bambu (R) | `bambu(reads='aln.bam', annotations=gtf)` |

## Key Patterns

**Pattern 1: ONT processing pipeline**
```bash
# Basecalling with Dorado (high accuracy model)
dorado basecaller hac pod5_data/ | samtools fastq > reads.fastq.gz

# QC
NanoStat --fastq reads.fastq.gz --outdir nanostat/
NanoFilt -q 10 -l 1000 reads.fastq.gz > reads_filtered.fastq.gz

# Alignment
minimap2 -ax map-ont hg38.fa reads_filtered.fastq.gz | \
    samtools sort -o aligned.bam && samtools index aligned.bam
```

**Pattern 2: De novo assembly**
```bash
# Flye for ONT genome assembly
flye --nano-hq reads_filtered.fastq.gz --genome-size 3g \
    --out-dir flye_out/ --threads 16

# Hifiasm for PacBio HiFi
hifiasm -o sample.asm -t 16 hifi_reads.fastq.gz
awk '/^S/{print ">"$2"\n"$3}' sample.asm.bp.p_ctg.gfa > assembly.fasta

# Polish ONT assembly with Medaka
medaka_consensus -i reads_filtered.fastq.gz -d flye_out/assembly.fasta \
    -o medaka/ -t 8 -m r1041_e82_400bps_hac_v4.2.0
```

**Pattern 3: Structural variant calling**
```bash
sniffles --input aligned.bam --vcf svs.vcf \
    --reference hg38.fa --threads 8 --minsupport 5
```

**Pattern 4: Methylation detection**
```bash
# Basecall with methylation model
dorado basecaller hac,5mCG_5hmCG pod5/ > calls_mod.bam

# Extract CpG methylation
modkit pileup calls_mod.bam methylation.bed \
    --ref hg38.fa --cpg --threads 8
```

**Pattern 5: Isoform analysis (R)**
```r
library(bambu)
se <- bambu(reads='aligned.bam',
            annotations=gencode_gtf,
            genome=hg38_fa)
writeBambuOutput(se, path='bambu_output/')
# se@rowRanges contains isoform coordinates
# assay(se, 'CPM') contains isoform-level expression
```

## Technology Comparison

| Property | ONT R9/R10 | PacBio HiFi |
|----------|-----------|-------------|
| Read length | Typically 5–50 kb | Typically 15–25 kb |
| Raw accuracy | 97–99% (R10) | >99.9% |
| Throughput | High (P2 Solo: 80 Gb) | Moderate (Sequel II: 160 Gb) |
| Methylation | Direct (native DNA) | 5mC with Kinetics |
| Cost per Gb | Low | Higher |

## Common Pitfalls

- **Basecall model selection** — match model to flow cell (R9.4 vs R10.4) and kit chemistry
- **Assembly genome size** — always provide `--genome-size` to Flye for ploidy-aware assembly
- **Coverage for assembly** — aim for ≥50× for Flye; ≥30× for Hifiasm HiFi
- **SV minimum support** — default `--minsupport 5` for Sniffles2; lower for low-coverage data
- **Medaka model** — use the correct model matching your basecaller version and flow cell
