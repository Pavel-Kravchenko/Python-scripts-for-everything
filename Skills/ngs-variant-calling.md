---
name: ngs-variant-calling
description: NGS sequencing pipelines — FASTQ QC, read alignment (BWA/HISAT2), SAM/BAM processing, variant calling (GATK/bcftools), VCF parsing and annotation
---

# NGS & Variant Calling

## When to Use
- Processing raw sequencing data (FASTQ → BAM → VCF)
- Calling germline or somatic variants from WGS/WES
- Filtering, annotating, or interpreting VCF files
- Parsing SAM/BAM flags, CIGAR strings, or coverage stats
- Population genetics: allele frequencies, HWE, GWAS

## Quick Reference

### Sequencing Platform Comparison
| Feature | Illumina | PacBio HiFi | Oxford Nanopore |
|---------|----------|-------------|-----------------|
| Read length | 50–300 bp | 10–25 kb | 10 kb–1 Mb+ |
| Accuracy | ~99.9% | ~99.9% (Q30+) | ~95–99% |
| Error type | Substitutions | Insertions (random) | Indels |
| Throughput | Up to 6 Tb/run | ~30 Gb/cell | 50–200 Gb/cell |
| Best for | WGS, RNA-seq, ChIP-seq | Assembly, SVs | Field work, ultra-long |

### FASTQ Format
```
@READ_ID                  # Line 1: header
GATCGATCGATC              # Line 2: sequence
+                         # Line 3: separator
IIIIHHHGGG##              # Line 4: quality (Phred+33 ASCII)
```
**Phred+33 encoding**: `quality = ord(char) - 33`; `P_error = 10^(-Q/10)`

| Q score | Error prob | ASCII char |
|---------|-----------|------------|
| 10 | 10% | + |
| 20 | 1% | 5 |
| 30 | 0.1% | ? |
| 40 | 0.01% | I |

**Paired-end**: `sample_R1.fastq.gz` (forward) + `sample_R2.fastq.gz` (reverse); reads at matching line positions are from the same fragment.

### FastQC Key Metrics
| Module | Target | Common Issue |
|--------|--------|-------------|
| Per-base quality | Q28+ across all positions | Quality drops at 3' end (normal) |
| Per-sequence quality | Peak Q30+ | Bimodal = some reads failed |
| GC content | Normal distribution | Shifted = contamination |
| Duplication | <20% | High in PCR-heavy libraries |
| Adapter content | <5% | >10% → trim |
| Overrepresented seqs | None | Adapter dimers, rRNA |

Note: RNA-seq normally fails "Per-base sequence content" (hexamer bias); amplicon data fails duplication modules — not all warnings are errors.

### SAM/BAM Format
11 mandatory fields per alignment record:
`QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL`

**FLAG bits** (bitwise OR):
| Value | Meaning |
|-------|---------|
| 1 | Read is paired |
| 2 | Proper pair |
| 4 | Read unmapped |
| 8 | Mate unmapped |
| 16 | Reverse strand |
| 32 | Mate reverse |
| 64 | R1 (first in pair) |
| 128 | R2 (second in pair) |
| 256 | Secondary alignment |
| 1024 | PCR/optical duplicate |
| 2048 | Supplementary alignment |

Common values: 99 = R1 forward proper pair; 147 = R2 reverse proper pair; 4 = unmapped.

**CIGAR operations**:
| Op | Consumes query | Consumes ref | Meaning |
|----|---------------|-------------|---------|
| M | Yes | Yes | Match/mismatch |
| I | Yes | No | Insertion |
| D | No | Yes | Deletion |
| N | No | Yes | Intron skip (RNA-seq) |
| S | Yes | No | Soft clip |
| H | No | No | Hard clip |

Examples: `100M` → 100 bp aligned; `50M2I48M` → insertion; `50M10000N50M` → RNA-seq exon junction.

### VCF Format
```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS  ID   REF  ALT  QUAL  FILTER  INFO         FORMAT  SAMPLE1
chr1    10000 rs1  A    G    5000  PASS    DP=200;AF=0.5 GT:DP   0/1:100
```

Genotype encodings: `0/0` hom-ref; `0/1` het; `1/1` hom-alt; `1/2` two alt alleles; `./.` missing; `0|1` phased.

## Key Patterns

### Complete Pipeline (BWA → samtools → GATK)
```bash
# 1. QC
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc/ -t 4

# 2. Trim adapters
fastp -i R1.fastq.gz -I R2.fastq.gz \
      -o R1.trimmed.fastq.gz -O R2.trimmed.fastq.gz \
      --detect_adapter_for_pe --cut_tail --cut_tail_mean_quality 20 \
      --length_required 36 --html report.html

# 3. Index reference (once)
bwa index reference.fasta

# 4. Align + sort
bwa mem -t 8 reference.fasta R1.trimmed.fastq.gz R2.trimmed.fastq.gz \
  | samtools sort -o aligned.sorted.bam
samtools index aligned.sorted.bam

# 5. Mark duplicates
gatk MarkDuplicates -I aligned.sorted.bam -O dedup.bam -M metrics.txt

# 6. Base quality recalibration (BQSR)
gatk BaseRecalibrator -I dedup.bam -R ref.fasta \
     --known-sites dbsnp.vcf -O recal.table
gatk ApplyBQSR -I dedup.bam -R ref.fasta \
     --bqsr-recal-file recal.table -O recal.bam

# 7. Call variants
gatk HaplotypeCaller -I recal.bam -R ref.fasta -O raw.g.vcf -ERC GVCF

# 8. Joint genotyping (cohort)
gatk GenotypeGVCFs -R ref.fasta -V raw.g.vcf -O genotyped.vcf

# 9. Hard filter variants
gatk VariantFiltration -R ref.fasta -V genotyped.vcf \
  --filter-expression "QD < 2.0"   --filter-name "LowQD" \
  --filter-expression "FS > 60.0"  --filter-name "StrandBias" \
  --filter-expression "MQ < 40.0"  --filter-name "LowMQ" \
  -O filtered.vcf
```

### RNA-seq Alignment (HISAT2/STAR)
```bash
# HISAT2 (splice-aware, lighter)
hisat2-build reference.fasta ref_index
hisat2 -x ref_index -1 R1.fastq.gz -2 R2.fastq.gz \
       --dta -p 8 | samtools sort -o rna.sorted.bam

# STAR (faster for large datasets)
STAR --runMode genomeGenerate --genomeDir star_index \
     --genomeFastaFiles reference.fasta --sjdbGTFfile annotation.gtf
STAR --genomeDir star_index --readFilesIn R1.fastq.gz R2.fastq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix sample_
```

### samtools Common Operations
```bash
samtools flagstat sorted.bam               # alignment summary
samtools stats sorted.bam | grep ^SN       # detailed stats
samtools depth sorted.bam > depth.txt      # per-position depth
samtools view -b -q 30 sorted.bam -o hq.bam   # filter MAPQ>=30
samtools view -b -f 2 sorted.bam -o paired.bam # only proper pairs
samtools view -b -F 4 sorted.bam -o mapped.bam  # exclude unmapped
samtools view -c sorted.bam chr1:1-1000000  # count reads in region
samtools markdup sorted.bam dedup.bam       # mark duplicates
```

### bcftools Quick Variant Calling
```bash
bcftools mpileup -f ref.fasta sorted.bam \
  | bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz
bcftools stats variants.vcf.gz | grep ^SN   # summary stats
bcftools filter -e 'QUAL<30 || DP<10' variants.vcf.gz -o filtered.vcf
```

### Variant Annotation
```bash
# VEP (Ensembl)
vep -i variants.vcf -o annotated.vcf \
    --cache --assembly GRCh38 \
    --sift b --polyphen b --af --af_gnomad \
    --appris --canonical --vcf

# SnpEff
snpEff GRCh38.86 variants.vcf > annotated.vcf

# ANNOVAR
table_annovar.pl variants.vcf humandb/ -buildver hg38 \
    -protocol refGene,clinvar_20220320,gnomad30_genome \
    -operation g,f,f -vcfinput
```

## Code Templates

### Parse FASTQ
```python
def parse_fastq(filepath):
    with open(filepath) as f:
        while True:
            header = f.readline().strip()
            if not header: break
            seq  = f.readline().strip()
            f.readline()  # '+'
            qual = f.readline().strip()
            yield header[1:], seq, qual

def phred_scores(qual_str):
    return [ord(c) - 33 for c in qual_str]
```

### Decode SAM FLAG
```python
FLAG_BITS = {
    1: 'paired', 2: 'proper_pair', 4: 'unmapped', 8: 'mate_unmapped',
    16: 'reverse', 32: 'mate_reverse', 64: 'R1', 128: 'R2',
    256: 'secondary', 1024: 'duplicate', 2048: 'supplementary',
}

def decode_flag(flag):
    return [desc for bit, desc in FLAG_BITS.items() if flag & bit]
```

### Parse CIGAR
```python
import re

def parse_cigar(cigar):
    return [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def ref_length(cigar):
    """Bases consumed on reference."""
    return sum(n for n, op in parse_cigar(cigar) if op in 'MDN=X')
```

### Parse VCF
```python
def parse_vcf(path):
    with open(path) as f:
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#'):
                cols = line[1:].strip().split('\t')
                continue
            fields = line.strip().split('\t')
            rec = dict(zip(cols, fields))
            # Parse INFO
            rec['INFO_dict'] = dict(
                kv.split('=', 1) if '=' in kv else (kv, True)
                for kv in rec['INFO'].split(';')
            )
            yield rec
```

### Hardy-Weinberg Test
```python
from scipy.stats import chi2

def hwe_test(genotypes):
    """genotypes: list of (0,0)/(0,1)/(1,1) tuples"""
    n = len(genotypes)
    q = sum(a for g in genotypes for a in g) / (2 * n)
    p = 1 - q
    obs = [sum(1 for g in genotypes if sum(g) == k) for k in (0, 1, 2)]
    exp = [p**2 * n, 2*p*q * n, q**2 * n]
    chi2_stat = sum((o - e)**2 / e for o, e in zip(obs, exp) if e > 0)
    return chi2_stat, chi2.sf(chi2_stat, df=1)
```

## Common Pitfalls
- **Aligner choice**: Use BWA-MEM for DNA, HISAT2/STAR for RNA-seq (BWA cannot handle intron skips).
- **Reference genome mismatch**: chr1 vs 1 naming breaks pipelines; ensure FASTQ, BAM, VCF, and reference all use the same chromosome naming convention.
- **Skip BQSR only when justified**: Required for GATK best practices with WGS/WES; skip for amplicon or low-coverage data.
- **Duplicate marking before variant calling**: Always mark duplicates; do not remove them (downstream tools need the flags).
- **Hard filters vs VQSR**: Use VQSR for large cohorts (>30 WGS samples); use hard filters for small cohorts or targeted panels.
- **GQ vs QUAL**: QUAL is variant-level confidence; GQ (Genotype Quality) is per-sample. Filter both — a high QUAL with low GQ per sample can still yield wrong genotypes.
- **Allele balance for hets**: Expected ~0.5 (0.2–0.8 acceptable); extreme values (0.1 or 0.9) suggest false positives.
- **VCF multi-allelic sites**: ALT field can be comma-separated; genotype index 2 refers to the second ALT, not REF.
- **Soft-clipped reads**: CIGAR `S` bases are in SEQ but not aligned; include them in read length but not reference position calculations.
- **Population stratification in GWAS**: Lambda (genomic inflation factor) >1.05 suggests uncorrected stratification; use principal components as covariates.

---

## Related Skills
- `rnaseq` — differential expression after alignment (DESeq2/edgeR)
- `biostatistics-r` — statistical tests, HWE, chi-squared, logistic regression
- `sequence-alignment` — pairwise alignment, BLAST, motif finding
