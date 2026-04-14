---
name: bio-applied-variant-calling-and-snp-analysis
description: Variant calling pipeline, VCF format, genotype decoding, and SNP analysis with GATK/bcftools
tool_type: python
primary_tool: NumPy
---

# Variant Calling and SNP Analysis

## Variant Type Reference

| Type | Size | Example |
|------|------|---------|
| SNV | 1 bp | A→G |
| MNV | 2+ bp equal length | AT→GC |
| Insertion | 1–50 bp | A→ATCG |
| Deletion | 1–50 bp | ATCG→A |
| Large SV | >50 bp | detected by split reads / depth |
| CNV | variable | depth-based (deletion=low depth, dup=high depth) |

## GATK Best Practices Pipeline

```bash
# BWA alignment — always set read groups or MarkDuplicates/HaplotypeCaller will fail
bwa mem -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1' ref.fa R1.fq R2.fq | samtools sort -o sorted.bam

# Mark duplicates
gatk MarkDuplicates -I sorted.bam -O dedup.bam -M metrics.txt

# BQSR (optional but recommended)
gatk BaseRecalibrator -I dedup.bam -R ref.fa --known-sites dbsnp.vcf -O recal.table
gatk ApplyBQSR -I dedup.bam -R ref.fa --bqsr-recal-file recal.table -O recal.bam

# Call variants (GVCF mode for cohort)
gatk HaplotypeCaller -I recal.bam -R ref.fa -O raw.g.vcf -ERC GVCF

# Joint genotyping
gatk GenotypeGVCFs -R ref.fa -V raw.g.vcf -O genotyped.vcf

# Hard filters
gatk VariantFiltration -R ref.fa -V genotyped.vcf \
  --filter-expression "QD < 2.0" --filter-name "LowQD" \
  --filter-expression "FS > 60.0" --filter-name "StrandBias" \
  --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
  -O filtered.vcf
```

## Variant Calling Tool Comparison

| Tool | Algorithm | Best For |
|------|-----------|----------|
| GATK HaplotypeCaller | Local de novo assembly | WGS/WES cohorts, gold standard |
| bcftools mpileup/call | Pileup-based | Fast single-sample WGS |
| FreeBayes | Bayesian haplotype | Low-frequency variants |
| DeepVariant | Deep learning (CNN) | High accuracy WGS/WES |
| Strelka2 | Statistical model | Tumor-normal somatic |

## VCF Format

```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS    ID       REF  ALT  QUAL  FILTER  INFO          FORMAT  SAMPLE1
chr1    10000  rs123456 A    G    5000  PASS    DP=200;AF=0.5 GT:DP   0/1:100
```

### Genotype Encodings

| GT | Meaning |
|----|---------|
| 0/0 | Hom-ref |
| 0/1 | Het |
| 1/1 | Hom-alt |
| 1/2 | Het multi-allelic |
| ./. | Missing |
| 0\|1 | Phased het |

## Key Code Patterns

```python
def classify_variant(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 'SNV'
    elif len(ref) == len(alt) and len(ref) > 1:
        return 'MNV'
    elif len(ref) < len(alt):
        return 'insertion'
    elif len(ref) > len(alt):
        return 'deletion'
    return 'complex'

def parse_vcf_info(info_str):
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            info[k] = v
        else:
            info[item] = True
    return info

def decode_genotype(gt_string, ref, alts):
    """Returns (allele_string, zygosity, is_phased)."""
    alleles = [ref] + alts
    sep = '|' if '|' in gt_string else '/'
    indices = gt_string.split(sep)
    decoded = [alleles[int(i)] if i != '.' else '.' for i in indices]
    if '.' in indices:
        zyg = 'missing'
    elif len(set(indices)) == 1:
        zyg = 'hom-ref' if indices[0] == '0' else 'hom-alt'
    else:
        zyg = 'het'
    return '/'.join(decoded), zyg, sep == '|'
```

## Pitfalls

- **Coordinate systems**: BED is 0-based half-open; VCF/GFF are 1-based inclusive — off-by-one errors happen when mixing them
- **Missing read groups**: GATK MarkDuplicates and HaplotypeCaller require `@RG` tags — always pass `-R` to `bwa mem`
- **VQSR vs hard filters**: VQSR requires ≥30 WGS samples or ≥10 WES samples; use hard filters for smaller cohorts
- **Multi-allelic sites**: `ALT` may be comma-separated; split before per-allele analysis
- **Batch effects**: Systematic differences between sequencing runs confound variant calling — batch-correct or include batch as covariate
- **Multiple testing**: Apply Benjamini-Hochberg FDR when testing thousands of variants
