---
name: bio-applied-variant-calling-and-snp-analysis
description: "This notebook covers the complete variant analysis pipeline: from aligned reads to annotated variants, including VCF format, variant filtering, annotation databases, population genetics, and GWAS conc"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/01_variant_calling_and_snp_analysis.ipynb"
---

# Variant Calling and SNP Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/01_variant_calling_and_snp_analysis.ipynb`*

# Variant Calling and SNP Analysis

## Tier 3 - Applied Bioinformatics

This notebook covers the complete variant analysis pipeline: from aligned reads to annotated variants, including VCF format, variant filtering, annotation databases, population genetics, and GWAS concepts.

### Learning Objectives
- Classify variant types: SNPs, indels, structural variants, CNVs
- Understand the variant calling pipeline
- Parse and work with VCF format in depth
- Apply variant quality filters
- Annotate variants using databases (ClinVar, dbSNP, gnomAD)
- Understand population genetics concepts (allele frequency, HWE)
- Interpret GWAS results: Manhattan plots, QQ plots

---

## 1. Types of Genetic Variants

### 1.1 Small Variants

| Type | Description | Size | Example |
|------|------------|------|---------|
| **SNP/SNV** | Single nucleotide change | 1 bp | A -> G |
| **MNP/MNV** | Multiple adjacent nucleotide changes | 2+ bp | AT -> GC |
| **Insertion** | Bases added | 1-50 bp | A -> ATCG |
| **Deletion** | Bases removed | 1-50 bp | ATCG -> A |
| **Indel** | Insertion or deletion | 1-50 bp | Complex |

### 1.2 Structural Variants (SVs)

| Type | Description | Size | Detection |
|------|------------|------|-----------|
| **Large deletion** | Loss of genomic segment | >50 bp | Split reads, depth changes |
| **Large insertion** | Gain of genomic segment | >50 bp | Split reads |
| **Duplication** | Segment copied | >50 bp | Increased depth |
| **Inversion** | Segment reversed | >50 bp | Discordant read pairs |
| **Translocation** | Segment moved between chromosomes | Variable | Discordant read pairs |

### 1.3 Copy Number Variants (CNVs)

```
Normal:       [====A====][====B====][====C====]
Deletion:     [====A====]           [====C====]   (loss of B)
Duplication:  [====A====][====B====][====B====][====C====]  (gain of B)
Amplification:[====A====][=B=][=B=][=B=][=B=] [====C====]  (multiple copies)

CNV detection: read depth analysis
  Normal region:  |||||||||||||||  (expected depth)
  Deletion:       ||| ||| ||       (reduced depth)
  Duplication:    |||||||||||||||||||||||  (increased depth)
```

```python
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import re
import math
import random

random.seed(42)
np.random.seed(42)

def classify_variant(ref, alt):
    """Classify a variant based on REF and ALT alleles."""
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    elif len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
        return 'MNV'
    elif len(ref) < len(alt):
        return 'insertion'
    elif len(ref) > len(alt):
        return 'deletion'
    else:
        return 'complex'

# Examples
variants = [
    ('A', 'G', 'rs334'),
    ('AT', 'GC', 'rs12345'),
    ('A', 'ATCG', '.'),
    ('ATCG', 'A', '.'),
    ('C', 'T', 'rs7903146'),
]

print(f"{'REF':>6} {'ALT':>6} {'Type':>12} {'ID':>12}")
print("-" * 42)
for ref, alt, vid in variants:
    vtype = classify_variant(ref, alt)
    print(f"{ref:>6} {alt:>6} {vtype:>12} {vid:>12}")
```

---

## 2. Variant Calling Pipeline

### 2.1 Pipeline Overview

```
Aligned BAM file
       |
       v
  Mark duplicates (Picard / samtools markdup)
       |
       v
  Base quality score recalibration (GATK BQSR) [optional]
       |
       v
  Variant calling (GATK HaplotypeCaller / bcftools / FreeBayes)
       |
       v
  Raw VCF
       |
       v
  Variant filtering (GATK VQSR / hard filters)
       |
       v
  Filtered VCF
       |
       v
  Annotation (VEP / ANNOVAR / SnpEff)
       |
       v
  Annotated variants for interpretation
```

### 2.2 Variant Calling Tools

| Tool | Algorithm | Best For | Notes |
|------|-----------|----------|-------|
| **GATK HaplotypeCaller** | Local de novo assembly | WGS, WES | Gold standard, part of GATK Best Practices |
| **bcftools mpileup/call** | Pileup-based | WGS, quick analysis | Fast, lightweight |
| **FreeBayes** | Bayesian haplotype-based | WGS, small cohorts | Good for low-frequency variants |
| **DeepVariant** | Deep learning | WGS, WES | Google's neural network caller |
| **Strelka2** | Statistical model | Somatic variants | Good for tumor-normal pairs |

### 2.3 How GATK HaplotypeCaller Works

1. **Identify active regions**: Find regions where reads suggest variants
2. **Local reassembly**: Build a de Bruijn graph of reads in the region
3. **Haplotype determination**: Extract candidate haplotypes from the graph
4. **Pairwise alignment**: Align reads to each candidate haplotype
5. **Genotype assignment**: Use Bayesian model to assign genotypes

```bash
# GATK Best Practices pipeline
# NOTE: GATK requires read groups in the BAM file.
# When aligning with BWA, always include the -R flag:
#   bwa mem -R '@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1' \
#       ref.fa R1.fq.gz R2.fq.gz | samtools sort -o sorted.bam
# Without read groups, MarkDuplicates and HaplotypeCaller will fail.

# 1. Mark duplicates
gatk MarkDuplicates -I sorted.bam -O dedup.bam -M metrics.txt

# 2. Base quality score recalibration
gatk BaseRecalibrator -I dedup.bam -R ref.fa --known-sites dbsnp.vcf -O recal.table
gatk ApplyBQSR -I dedup.bam -R ref.fa --bqsr-recal-file recal.table -O recal.bam

# 3. Call variants
gatk HaplotypeCaller -I recal.bam -R ref.fa -O raw.g.vcf -ERC GVCF

# 4. Joint genotyping (for cohorts)
gatk GenotypeGVCFs -R ref.fa -V raw.g.vcf -O genotyped.vcf

# 5. Filter variants
gatk VariantFiltration -R ref.fa -V genotyped.vcf \
  --filter-expression "QD < 2.0" --filter-name "LowQD" \
  --filter-expression "FS > 60.0" --filter-name "StrandBias" \
  --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
  -O filtered.vcf
```

---

## 3. VCF Format in Depth

### 3.1 VCF Structure

```
## Meta-information lines (start with ##)
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Allele Number">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
##FILTER=<ID=LowQual,Description="Low quality">

## Header line (starts with #)
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1  SAMPLE2

## Data lines
chr1  10000  rs123  A  G  5000  PASS  DP=200;AF=0.50;AC=1;AN=2  GT:DP:GQ:AD  0/1:100:99:50,50
```

### 3.2 Column Details

| Column | Description | Notes |
|--------|-------------|-------|
| CHROM | Chromosome | chr1, chr2, ... |
| POS | 1-based position | Position of REF allele |
| ID | Variant ID | rsID or . if novel |
| REF | Reference allele | Must match reference genome |
| ALT | Alternate allele(s) | Comma-separated for multi-allelic |
| QUAL | Phred-scaled quality | -10*log10(P[no variant]) |
| FILTER | Filter status | PASS or semicolon-separated filters |
| INFO | Variant-level annotations | Key=Value pairs, semicolon-separated |
| FORMAT | Sample field format | Colon-separated field names |
| SAMPLE | Sample genotype data | Colon-separated values matching FORMAT |

### 3.3 Genotype Encodings

| Genotype | Meaning | Zygosity |
|----------|---------|----------|
| 0/0 | Homozygous reference | - |
| 0/1 | Heterozygous | Het |
| 1/1 | Homozygous alternate | Hom-alt |
| 1/2 | Heterozygous for two different alts | Het (multi-allelic) |
| ./. | Missing / no call | - |
| 0|1 | Phased heterozygous | Het (phased) |

```python
def parse_vcf(vcf_text):
    """
    Parse VCF text into structured data.
    
    Returns:
        meta: list of meta-information lines
        header: list of column names
        variants: list of variant dictionaries
    """
    meta = []
    header = []
    variants = []
    
    for line in vcf_text.strip().split('\n'):
        if line.startswith('##'):
            meta.append(line)
        elif line.startswith('#'):
            header = line[1:].split('\t')
        else:
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            
            # Parse INFO field
            info = {}
            for item in fields[7].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info[key] = value
                else:
                    info[item] = True  # Flag fields
            
            variant = {
                'chrom': fields[0],
                'pos': int(fields[1]),
                'id': fields[2],
                'ref': fields[3],
                'alt': fields[4].split(','),  # Multiple alts possible
                'qual': float(fields[5]) if fields[5] != '.' else None,
                'filter': fields[6],
                'info': info,
                'type': classify_variant(fields[3], fields[4].split(',')[0]),
            }
            
            # Parse sample genotypes if present
            if len(fields) > 9:
                fmt_fields = fields[8].split(':')
                samples = {}
                for i, sample_data in enumerate(fields[9:]):
                    sample_values = sample_data.split(':')
                    sample_dict = {}
                    for j, key in enumerate(fmt_fields):
                        if j < len(sample_values):
                            sample_dict[key] = sample_values[j]
                    samples[f'SAMPLE{i+1}'] = sample_dict
                variant['samples'] = samples
            
            variants.append(variant)
    
    return meta, header, variants


# Example VCF with multiple samples and rich annotations
vcf_example = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3
chr1\t10000\trs123456\tA\tG\t5000\tPASS\tDP=300;AF=0.33;AC=2;AN=6\tGT:DP:GQ:AD\t0/1:100:99:50,50\t0/0:95:99:95,0\t0/1:105:99:55,50
chr1\t15000\trs789012\tC\tT\t3000\tPASS\tDP=250;AF=0.50;AC=3;AN=6\tGT:DP:GQ:AD\t0/1:80:99:40,40\t1/1:85:99:0,85\t0/0:85:99:85,0
chr1\t20000\t.\tATG\tA\t1500\tPASS\tDP=180;AF=0.17;AC=1;AN=6\tGT:DP:GQ:AD\t0/0:60:80:60,0\t0/0:55:75:55,0\t0/1:65:90:40,25
chr1\t25000\trs345678\tG\tGA\t2000\tPASS\tDP=220;AF=0.33;AC=2;AN=6\tGT:DP:GQ:AD\t0/1:75:99:40,35\t0/0:70:99:70,0\t0/1:75:95:45,30
chr2\t50000\trs111222\tT\tC\t800\tLowQual\tDP=50;AF=0.17;AC=1;AN=6\tGT:DP:GQ:AD\t0/0:15:30:15,0\t0/0:20:40:20,0\t0/1:15:20:10,5
chr2\t60000\trs333444\tA\tG,T\t4500\tPASS\tDP=280;AF=0.33,0.17;AC=2,1;AN=6\tGT:DP:GQ:AD\t0/1:90:99:45,45,0\t0/2:95:99:50,0,45\t0/0:95:99:95,0,0
chr3\t100000\t.\tC\tT\t200\tLowQual;LowDP\tDP=20;AF=0.50;AC=1;AN=2\tGT:DP:GQ:AD\t./.:.:.:.\t0/1:20:15:12,8\t./.:.:.:.\t"""

meta, header, variants = parse_vcf(vcf_example)

print(f"Meta-information lines: {len(meta)}")
print(f"Columns: {header}")
print(f"Variants: {len(variants)}\n")

for v in variants:
    alt_str = ','.join(v['alt'])
    print(f"{v['chrom']}:{v['pos']} {v['id']:>12}  {v['ref']}>{alt_str}  "
          f"QUAL={v['qual']}  {v['filter']}  Type={v['type']}  "
          f"AF={v['info'].get('AF', 'N/A')}  DP={v['info'].get('DP', 'N/A')}")
```

```python
def decode_genotype(gt_string, ref, alts):
    """Decode a VCF genotype string into alleles."""
    alleles = [ref] + alts
    separator = '|' if '|' in gt_string else '/'
    indices = gt_string.split(separator)
    
    decoded = []
    for idx in indices:
        if idx == '.':
            decoded.append('.')
        else:
            decoded.append(alleles[int(idx)])
    
    phased = separator == '|'
    
    # Determine zygosity
    if '.' in indices:
        zygosity = 'missing'
    elif len(set(indices)) == 1:
        zygosity = 'hom-ref' if indices[0] == '0' else 'hom-alt'
    else:
        zygosity = 'het'
    
    return '/'.join(decoded), zygosity, phased

# Show genotypes for each variant
print("Variant Genotypes Across Samples")
print("=" * 80)

for v in variants:
    print(f"\n{v['chrom']}:{v['pos']} ({v['id']}) {v['ref']}>{','.join(v['alt'])}")
    if 'samples' in v:
        for sample_name, sample_data in v['samples'].items():
            gt = sample_data.get('GT', './.')
            dp = sample_data.get('DP', '.')
            gq = sample_data.get('GQ', '.')
            ad = sample_data.get('AD', '.')
            alleles, zygosity, phased = decode_genotype(gt, v['ref'], v['alt'])
            print(f"  {sample_name}: GT={gt} ({alleles}, {zygosity})  DP={dp}  GQ={gq}  AD={ad}")
```
