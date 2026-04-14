---
name: bio-applied-viral-genome-assembly
description: Viral genome assembly pipeline — ARTIC amplicon sequencing, iVar/LoFreq variant calling, quasispecies/minority variant detection, QC thresholds, and key pitfalls
tool_type: python
primary_tool: Matplotlib
---

# Viral Genome Assembly and Variant Analysis

## Sequencing Strategy Selection

| Strategy | Description | Use case |
|----------|-------------|----------|
| Amplicon (ARTIC) | Tiling PCR from known reference | SARS-CoV-2, HIV, Influenza |
| Metagenomic (shotgun) | Unbiased, needs host depletion | Novel/unknown pathogens |
| Hybrid capture | Probe enrichment | Low-titer (hepatitis, low-level viremia) |

**Key distinction from bacterial/human assembly**: viruses exist as **quasispecies** — a cloud of related variants, not a single clone. Consensus ≠ population.

## Coverage Thresholds

| Depth | Sufficient for |
|-------|----------------|
| ≥20× | Consensus calling (iVar) |
| ≥100× | Low-confidence minority variants |
| ≥200× | Reliable minority variant calling (≥1%) |
| ≥1000× | Very low-frequency variants (0.1–0.5%) |

## Reference-Guided Assembly Pipeline

```bash
# 1. Quality + adapter trim
fastp -i R1.fastq -I R2.fastq -o R1.trim.fastq -O R2.trim.fastq --thread 4

# 2. Align (Illumina short reads)
minimap2 -ax sr ref.fa R1.trim.fastq R2.trim.fastq | samtools sort -o sorted.bam
samtools index sorted.bam

# 3. Trim ARTIC primers (soft-clips primer regions)
ivar trim -i sorted.bam -b primer_scheme.bed -p trimmed -m 20 -q 20

# 4. Consensus (min freq 0.75, min depth 10)
samtools mpileup -aa -A -d 0 -Q 0 trimmed.bam | ivar consensus -p consensus -t 0.75 -m 10

# 5. Variant calling
lofreq call --call-indels -f ref.fa -o variants.vcf trimmed.bam
# OR iVar variants
samtools mpileup -aa -A -d 0 --reference ref.fa -Q 0 trimmed.bam | \
    ivar variants -p variants -q 20 -t 0.03 -m 10

# 6. Clade assignment + QC
nextclade run --input-fasta consensus.fa --input-dataset sars-cov-2 --output-tsv qc.tsv
```

## Key Tool Summary

| Step | Tool | Notes |
|------|------|-------|
| Short read alignment | `minimap2 -ax sr` | BWA-MEM2 also works |
| Nanopore alignment | `minimap2 -ax ont2d` | |
| Primer trimming | `ivar trim` | Required for ARTIC amplicons |
| Consensus | `ivar consensus -t 0.75` | 75% frequency threshold |
| Minority variants | `lofreq call` | Poisson + Bonferroni; best for Illumina |
| Minority variants | `ivar variants` | ARTIC standard |
| Annotation/QC | `nextclade run` | Clade + frameshifts + quality flags |

## Intra-Host Variant Calling Filters

Apply **all four** filters:

```python
filters = {
    'min_depth':    100,    # reads covering the site
    'min_af':       0.01,   # 1% alternate allele frequency
    'strand_bias':  0.001,  # Fisher exact p-value > 0.001
    'min_qual':     25,     # base quality at variant position
}
```

## Minority Variant Tool Comparison

| Tool | Model | Min freq | Notes |
|------|-------|----------|-------|
| LoFreq | Poisson + Bonferroni | 0.5% | Best accuracy on Illumina amplicons |
| iVar | Binomial exact | Configurable | ARTIC pipeline standard |
| DeepVariant | CNN | ~5% | High accuracy, slow |
| VarScan2 | Fisher exact | Configurable | Somatic/viral mode |

## QC Pass Criteria

```python
# Per-sample QC thresholds
QC_THRESHOLDS = {
    'mapped_pct':       80.0,   # % reads mapping to reference
    'genome_completeness': 90.0, # % genome covered at ≥20×
    'median_depth':     200,    # for variant calling
}
```

## Python Patterns

```python
import numpy as np
import pandas as pd
from scipy import stats

def apply_variant_filters(df_var, min_depth=100, min_af=0.01, sb_p_cutoff=0.001):
    df_var['pass_depth']  = df_var['Depth'] >= min_depth
    df_var['pass_af']     = df_var['Alt_freq'] >= min_af
    df_var['pass_strand'] = df_var['Strand_bias_p'] > sb_p_cutoff
    df_var['PASS'] = df_var['pass_depth'] & df_var['pass_af'] & df_var['pass_strand']
    return df_var

def strand_bias_pvalue(depth, alt_freq, fwd_alt_frac):
    """Fisher exact test for strand bias."""
    fwd = int(depth * alt_freq * fwd_alt_frac)
    rev = int(depth * alt_freq * (1 - fwd_alt_frac))
    fwd_ref = int(depth * (1 - alt_freq) * fwd_alt_frac)
    rev_ref = int(depth * (1 - alt_freq) * (1 - fwd_alt_frac))
    _, p = stats.fisher_exact([[fwd, rev], [fwd_ref, rev_ref]])
    return p
```

## Pitfalls

- **Primer dimers inflate low-quality variants**: always trim ARTIC primers with `ivar trim` before variant calling — un-trimmed primer sequences appear as high-frequency variants at amplicon boundaries
- **Amplicon dropout from primer mismatches**: divergent lineages can fail to amplify → dropout regions show false N-masking in consensus; update primer scheme each variant wave
- **Coverage uniformity**: ARTIC alternating pool design means some regions get systematically lower depth — check per-amplicon coverage, not just genome average
- **Quasispecies ≠ co-infection**: intermediate frequencies (30–70%) can be drift, technical noise, or genuine mixed infection — require independent validation before claiming co-infection
- **LoFreq requires realigned indels**: run `lofreq indelqual` and `samtools indelrealign` before `lofreq call --call-indels` for accurate indel variant calls
- **iVar consensus Ns**: positions below `min_depth` threshold are written as N — a 90% complete genome still has 3kb masked; downstream analysis (Pangolin, Nextclade) tolerates ≤30% Ns
- **Coordinate systems**: BED primer files are 0-based; iVar uses 0-based internally; VCF output is 1-based — off-by-one errors when cross-referencing
