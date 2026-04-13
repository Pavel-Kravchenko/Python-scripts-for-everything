---
name: bio-applied-snp-calling-pipeline
description: "This notebook walks through a complete, real-world SNP calling pipeline: from raw FASTQ reads to a final annotated variant report. We cover each tool, why it is used, what it produces, and how to inte"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/02_snp_calling_pipeline.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# SNP Calling Pipeline: A Practical Hands-On Workflow

*Source: Course notebook `Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/02_snp_calling_pipeline.ipynb`*

# SNP Calling Pipeline: A Practical Hands-On Workflow

## Tier 3 - Applied Bioinformatics

This notebook walks through a complete, real-world SNP calling pipeline: from raw FASTQ reads to a final annotated variant report. We cover each tool, why it is used, what it produces, and how to interpret the results.

The pipeline is based on: https://github.com/Pavel-Kravchenko/snp-calling-pipeline

### Learning Objectives
- Understand the end-to-end SNP calling workflow
- Apply Trimmomatic for quality-based read trimming
- Align reads with HISAT2 using genomic DNA settings
- Process BAM files and call variants with samtools / bcftools
- Annotate variants with ANNOVAR against RefGene, dbSNP, 1000 Genomes, GWAS Catalog, and ClinVar
- Parse and summarise ANNOVAR output in Python

---

## 1. Pipeline Architecture

The pipeline accepts single-end FASTQ reads and produces a variant report annotated against multiple databases.

> **Historical context:** This pipeline was built around HISAT2, a splice-aware RNA-seq aligner, which is used here with flags to disable spliced alignment (making it suitable for DNA). Modern pipelines use **BWA-MEM2** for DNA alignment. The biological results from either approach are comparable for SNP calling; see Section 4 for details.

```
Raw FASTQ reads
       |
       v
  Trimmomatic  (quality trimming: TRAILING:20, MINLEN:50)
       |
       v
  HISAT2       (alignment to reference genome)
               --no-spliced-alignment --no-softclip
       |
       v
  samtools     SAM → BAM → sort → index → idxstats → depth
       |
       v
  samtools mpileup
       |
       v
  bcftools call -cv
       |        raw VCF
       v
  ANNOVAR convert2annovar.pl   (VCF → avinput)
       |
       v
  ANNOVAR annotate_variation.pl
       ├── dbSNP138    (flag known SNPs)
       ├── refGene     (gene/exon annotation)
       ├── 1000G       (population allele frequency)
       ├── GWAS Catalog (trait associations)
       └── ClinVar     (clinical significance)
       |
       v
  Final annotated report
```

### Pipeline File Structure

```
snp-calling-pipeline/
├── reads/                 ← put your .fastq files here
├── Human/                 ← HISAT2 index + reference FASTA per chromosome
│   ├── chr1.fasta
│   ├── chr1.1.ht2 ... chr1.8.ht2   ← HISAT2 index files
│   └── builder.sh         ← builds the HISAT2 index
├── annovar/
│   ├── humandb/           ← ANNOVAR databases (downloaded separately)
│   ├── convert2annovar.pl
│   └── annotate_variation.pl
├── Trimmomatic-0.36/
│   └── trimmomatic-0.36.jar
├── hisat2-2.1.0/          ← HISAT2 binaries
├── output/                ← created at runtime, one sub-dir per sample
└── script.sh              ← main pipeline driver
```

---

## 2. Prerequisites: Tool Installation

Before running the pipeline install the required tools.

| Tool | Version used | Purpose | Install |
|------|-------------|---------|--------|
| **Trimmomatic** | 0.36 | Read quality trimming | http://www.usadellab.org/cms/?page=trimmomatic |
| **HISAT2** | 2.1.0 | Read alignment | https://ccb.jhu.edu/software/hisat2/ |
| **SAMtools** | ≥ 1.9 | BAM manipulation & pileup | http://samtools.sourceforge.net/ |
| **BCFtools** | ≥ 1.9 | Variant calling from pileup | bundled with SAMtools |
| **ANNOVAR** | latest | Variant annotation | http://annovar.openbioinformatics.org/ |

Quick conda install:

```bash
conda create -n snp-pipeline python=3.10
conda activate snp-pipeline
conda install -c bioconda trimmomatic hisat2 samtools bcftools
# ANNOVAR requires manual registration at http://annovar.openbioinformatics.org/
```

---

## 3. Step 1 — Read Trimming with Trimmomatic

Raw sequencing reads often have:
- **Low-quality bases** at the 3' end (quality decreases along the read)
- **Adapter contamination** from library preparation
- **Very short reads** that align unreliably

Trimmomatic removes these issues before alignment.

### 3.1 Pipeline Command

```bash
java -jar $home/Trimmomatic-0.36/trimmomatic-0.36.jar \
    SE \
    -phred33 \
    input.fastq \
    chr_outfile.fastq \
    TRAILING:20 \
    MINLEN:50
```

### 3.2 Parameter Breakdown

| Parameter | Value | Meaning |
|-----------|-------|--------|
| `SE` | — | Single-End mode (one file per sample) |
| `-phred33` | — | Illumina 1.8+ quality encoding (ASCII offset 33) |
| `TRAILING:20` | 20 | Remove trailing bases with Phred quality < 20 |
| `MINLEN:50` | 50 bp | Discard reads shorter than 50 bp after trimming |

**Why TRAILING rather than SLIDINGWINDOW?**  
Quality degradation in Illumina reads is worst at the 3' end. `TRAILING` is fast, simple, and effective for removing these low-quality tails. `SLIDINGWINDOW` is more thorough but slower.

**Phred quality score recap:**

| Phred Q | Error probability | Accuracy |
|---------|-----------------|----------|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1,000 | 99.9% |
| 40 | 1 in 10,000 | 99.99% |

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import random
import re

random.seed(42)
np.random.seed(42)


def phred_to_prob(q):
    """Convert Phred score to error probability."""
    return 10 ** (-q / 10)


def simulate_read_qualities(length=150, seed=42):
    """Simulate a typical Illumina read quality profile."""
    rng = np.random.default_rng(seed)
    # Good quality at start, gradual decline at 3' end
    base_quality = np.linspace(38, 15, length)
    noise = rng.normal(0, 3, length)
    qualities = np.clip(base_quality + noise, 2, 40).astype(int)
    return qualities


def trimmomatic_trailing(qualities, min_q=20):
    """Simulate Trimmomatic TRAILING: remove 3' bases below min_q."""
    trimmed_len = len(qualities)
    for i in range(len(qualities) - 1, -1, -1):
        if qualities[i] >= min_q:
            trimmed_len = i + 1
            break
    return trimmed_len


# Simulate 500 reads and show effect of trimming
read_lengths_before = []
read_lengths_after = []
kept = 0
discarded = 0
MINLEN = 50

for i in range(500):
    q = simulate_read_qualities(150, seed=i)
    read_lengths_before.append(150)
    trimmed = trimmomatic_trailing(q, min_q=20)
    if trimmed >= MINLEN:
        read_lengths_after.append(trimmed)
        kept += 1
    else:
        discarded += 1

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Quality profile of one read
q_example = simulate_read_qualities(150, seed=0)
trim_point = trimmomatic_trailing(q_example, min_q=20)
ax = axes[0]
ax.plot(range(1, 151), q_example, color='steelblue', lw=1.2, label='Quality score')
ax.axhline(y=20, color='red', linestyle='--', lw=1.5, label='TRAILING threshold (Q20)')
ax.axvline(x=trim_point, color='orange', linestyle='-', lw=2, label=f'Trim point ({trim_point} bp)')
ax.fill_betweenx([0, 40], trim_point, 150, alpha=0.15, color='red', label='Removed')
ax.set_xlabel('Position in read (bp)')
ax.set_ylabel('Phred quality score')
ax.set_title('Quality Profile: TRAILING:20 Trimming')
ax.legend(fontsize=8)
ax.set_ylim(0, 42)

# Distribution of read lengths after trimming
ax2 = axes[1]
ax2.hist(read_lengths_after, bins=30, color='steelblue', edgecolor='white', alpha=0.8)
ax2.axvline(x=MINLEN, color='red', linestyle='--', lw=2, label=f'MINLEN={MINLEN}')
ax2.set_xlabel('Read length after trimming (bp)')
ax2.set_ylabel('Number of reads')
ax2.set_title(f'Read Length Distribution After Trimming\n(kept={kept}, discarded={discarded})')
ax2.legend()

plt.tight_layout()
plt.show()

print(f"Reads kept: {kept} ({100*kept/500:.1f}%)")
print(f"Reads discarded (< {MINLEN} bp): {discarded} ({100*discarded/500:.1f}%)")
print(f"Mean length after trimming: {np.mean(read_lengths_after):.1f} bp")
```

---

## 4. Step 2 — Alignment with HISAT2

> **Note on tool choice:** HISAT2 is primarily a splice-aware aligner designed for RNA-seq. For new WGS/WES SNP-calling projects, **BWA-MEM2** (successor to BWA-MEM) or **Bowtie2** in end-to-end mode are the standard choices because they are optimised for DNA alignment and are fully integrated with GATK Best Practices. This pipeline uses HISAT2 with RNA-seq features disabled as a practical exercise; the HISAT2 flags  make it behave similarly to Bowtie2 in  mode, and the results are valid for learning purposes.

### 4.1 Building the Reference Index (builder.sh)

Before the first run, build a HISAT2 index for each reference chromosome:



This produces  index files used by the main alignment step.

For **new projects** using BWA-MEM2, the equivalent index build is:



### 4.2 Alignment Command



### 4.3 Key Flags Explained

| Flag | Purpose | Why it matters for SNP calling |
|------|---------|-------------------------------|
|  | Index prefix | Points to the per-chromosome HISAT2 index |
|  | Unpaired (single-end) reads | Matches our SE Trimmomatic output |
|  | Disables exon-junction spanning | DNA reads do not span introns; disabling this prevents false split-read alignments |
|  | Disables soft-clipping | Forces hard alignment; prevents unaligned bases near read ends from masking true variants |

### 4.4 Why  Matters

Without this flag, HISAT2 would allow reads to span introns (e.g., a read spanning a 10 kb intron boundary), which is biologically correct for RNA-seq but produces spurious gap-containing alignments in genomic DNA. Such alignments would appear as false insertions or deletions in the variant calling step.

| Aligner | Designed for | For DNA SNP calling |
|---------|-------------|---------------------|
| HISAT2 (default) | RNA-seq (splice-aware) | Use  |
| HISAT2 (with flags) | DNA-like | Acceptable; functionally similar to Bowtie2 |
| Bowtie2 --end-to-end | DNA | Preferred for learning; standard for ChIP-seq |
| BWA-MEM / BWA-MEM2 | DNA | Industry standard for WGS/WES |
| STAR | RNA-seq (splice-aware) | Not suitable for DNA SNP calling |

```python
def parse_hisat2_summary(summary_text):
    """Parse HISAT2 alignment summary statistics from stderr output."""
    stats = {}
    lines = summary_text.strip().split('\n')
    for line in lines:
        line = line.strip()
        m = re.match(r'(\d+) reads; of these:', line)
        if m:
            stats['total_reads'] = int(m.group(1))
        m = re.match(r'(\d+) \((.*?)%\) aligned exactly 1 time', line)
        if m:
            stats['uniquely_aligned'] = int(m.group(1))
            stats['unique_pct'] = float(m.group(2))
        m = re.match(r'(\d+) \((.*?)%\) aligned >1 times', line)
        if m:
            stats['multi_aligned'] = int(m.group(1))
            stats['multi_pct'] = float(m.group(2))
        m = re.match(r'(.*?)% overall alignment rate', line)
        if m:
            stats['overall_rate'] = float(m.group(1))
    return stats


# Typical HISAT2 summary from a WGS sample
hisat2_output = """
2500000 reads; of these:
  2500000 (100.00%) were unpaired; of these:
    62500 (2.50%) aligned 0 times
    2225000 (89.00%) aligned exactly 1 time
    212500 (8.50%) aligned >1 times
97.50% overall alignment rate
"""

stats = parse_hisat2_summary(hisat2_output)
total = stats.get('total_reads', 0)
unique = stats.get('uniquely_aligned', 0)
multi = stats.get('multi_aligned', 0)
unaligned = total - unique - multi

print("HISAT2 Alignment Summary")
print("=" * 40)
print(f"Total reads:           {total:>10,}")
print(f"Uniquely aligned:      {unique:>10,}  ({100*unique/total:.1f}%)")
print(f"Multi-mapped:          {multi:>10,}  ({100*multi/total:.1f}%)")
print(f"Unaligned:             {unaligned:>10,}  ({100*unaligned/total:.1f}%)")
print(f"Overall rate:          {stats.get('overall_rate', 0):>9.1f}%")
print()

# Visualise
fig, ax = plt.subplots(figsize=(7, 4))
categories = ['Uniquely\naligned', 'Multi-\nmapped', 'Unaligned']
counts = [unique, multi, unaligned]
colors = ['#2196F3', '#FF9800', '#F44336']
bars = ax.bar(categories, [c/total*100 for c in counts], color=colors, edgecolor='white', width=0.5)
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f'{count/total*100:.1f}%\n({count:,})', ha='center', va='bottom', fontsize=9)
ax.set_ylabel('Percentage of reads (%)')
ax.set_title('HISAT2 Alignment Rate')
ax.set_ylim(0, 100)
plt.tight_layout()
plt.show()

# Quality thresholds
print("Alignment quality assessment:")
if stats.get('overall_rate', 0) >= 90:
    print("  ✓ Overall rate >= 90%: GOOD")
elif stats.get('overall_rate', 0) >= 70:
    print("  ⚠ Overall rate 70-90%: ACCEPTABLE (check for contamination)")
else:
    print("  ✗ Overall rate < 70%: POOR (wrong reference or severe contamination)")
if unique / total >= 0.80:
    print("  ✓ Unique alignment rate >= 80%: GOOD")
else:
    print("  ⚠ High multi-mapping rate: consider duplicate-region filtering")
```
