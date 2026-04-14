---
name: bio-applied-wgbs-bismark
description: WGBS/RRBS Processing with Bismark with Bismark
tool_type: python
primary_tool: Bismark
---

# WGBS/RRBS Processing with Bismark

- [Bismark documentation](https://felixkrueger.github.io/Bismark/)
- [methylKit vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/)
- [Encode WGBS pipeline](https://www.encodeproject.org/wgbs/)

## DNA Methylation Biology

**5-methylcytosine (5mC)** is the most studied epigenetic mark in mammals. A methyl group is added to the 5-carbon position of cytosine by DNA methyltransferases (DNMTs): DNMT3A and DNMT3B establish de novo methylation, while DNMT1 maintains methylation through replication. The reverse reaction is catalysed by TET enzymes (TET1/2/3), which oxidize 5mC to 5-hydroxymethylcytosine (5hmC) and beyond, ultimately leading to passive or active demethylation.

### Genomic context of cytosine methylation

In mammals, methylation occurs almost exclusively at **CpG dinucleotides**. Roughly 70–80% of all CpGs in the human genome are methylated. The exceptions are:

| Feature | Methylation state | Function |
|---|---|---|
| CpG islands (CGIs) | Mostly **unmethylated** | Protect gene promoters from silencing |
| CGI shores (±2 kb) | Variable, tissue-specific | Major site of differential methylation |
| Gene bodies | Moderately methylated | Associated with active transcription |
| Repeats / transposons | Heavily methylated | Silences parasitic elements |
| Imprinted regions | Hemi-methylated (parent-of-origin) | Monoallelic gene expression |

**CpG islands** are defined as regions ≥200 bp with CpG observed/expected ratio > 0.6 and GC content > 50%. Despite representing only ~1% of the genome, they overlap ~70% of gene promoters. Their unmethylated state allows transcription factor binding.

### Methylation in development and disease

- During early embryogenesis, methylation is globally erased and re-established in a lineage-specific manner — explaining why different tissues have distinct methylation profiles.
- In cancer, CpG island promoter **hypermethylation** silences tumor suppressor genes (e.g., *BRCA1*, *MLH1*, *CDKN2A*), while global hypomethylation destabilizes the genome.
- Aging is associated with a gradual drift in methylation — the basis for epigenetic clocks (Notebook 3).

### Why measure methylation at single-CpG resolution?

Array-based methods (Illumina 450K/EPIC) measure 450,000–850,000 CpGs but miss most of the genome. **Whole-Genome Bisulfite Sequencing (WGBS)** provides base-resolution methylation across all ~28 million CpGs in the human genome, at the cost of sequencing depth (~30× per strand).

## Bisulfite Sequencing Chemistry

The key insight behind bisulfite sequencing is a simple chemical reaction:

> **Unmethylated cytosine + sodium bisulfite → uracil → reads as T**
> **Methylated cytosine + sodium bisulfite → protected → reads as C**

This allows methylation status to be inferred at single-base resolution by comparing the bisulfite-treated sequence to the reference genome: a C in the read means the site was methylated; a T means it was unmethylated.

### Step-by-step chemistry

1. **Denaturation**: DNA is denatured to single strands (required for bisulfite access).
2. **Sulfonation**: Bisulfite adds to the C5–C6 double bond of cytosine, forming cytosine-sulfonate.
3. **Deamination**: The amino group is removed, yielding uracil-sulfonate.
4. **Desulfonation**: Under alkaline conditions, the sulfonate is removed, leaving uracil.
5. **PCR amplification**: Uracil is copied as thymine.

Methylated cytosines resist sulfonation because the methyl group at C5 sterically blocks the reaction.

### Bisulfite conversion efficiency

Incomplete conversion (unconverted unmethylated C) is a critical quality metric. It is estimated using:
- **Lambda spike-in**: Bacteriophage lambda DNA (known to be unmethylated) is spiked in and used to compute non-conversion rate.
- **CHH/CHG sites**: In mammalian cells, non-CpG methylation is rare, so CHH/CHG C→T conversion rate reflects bisulfite efficiency.
- **Target**: > 99.5% conversion efficiency.

### WGBS vs RRBS vs EPIC array

| Method | Coverage | Cost | Best for |
|---|---|---|---|
| WGBS | All ~28M CpGs | High (>$500/sample) | Comprehensive studies |
| RRBS | ~5M CpGs (CpG-enriched) | Medium | Cost-effective targeted |
| EPIC array | 850K CpGs | Low (~$200) | Large cohorts |

**RRBS** uses MspI restriction enzyme (cuts at CCGG) to enrich CpG-dense regions before bisulfite treatment, reducing sequencing cost ~6× while covering most CpG islands.

### Effect on alignment

Bisulfite treatment creates four distinct strands (OT, CTOT, OB, CTOB), each with a different C→T conversion pattern. Standard aligners cannot handle this complexity — specialized aligners like **Bismark** are required.

## Bismark Alignment Pipeline

Bismark is the standard tool for aligning bisulfite-converted reads. It works by:

1. **Index preparation**: Converting all cytosines in the reference genome to thymine (CT conversion) and also converting all guanines on the reverse complement to adenine (GA conversion), then building Bowtie2 indices for both strands.
2. **Read alignment**: For each read, Bismark performs a CT conversion and aligns against the CT-converted genome; it also performs a GA conversion and aligns against the GA-converted genome. The best alignment wins.
3. **Methylation extraction**: Comparing aligned reads to the original (unconverted) reference, calling each CpG site as methylated (C in read) or unmethylated (T in read).

### Step 1: Build the bisulfite genome index

```bash
# Prepare genome index (only needs to be done once per genome)
bismark_genome_preparation /path/to/genome/hg38/

# This creates CT_conversion/ and GA_conversion/ subdirectories
# with Bowtie2 indices
```python

### Step 2: Quality trimming

Bisulfite-treated libraries require aggressive adapter trimming (Trim Galore wraps cutadapt with bisulfite-specific settings):

```bash
trim_galore --paired --fastqc \
    sample_R1.fastq.gz sample_R2.fastq.gz
```python

### Step 3: Alignment

```bash
bismark --genome /path/to/genome/hg38/ \
    -1 sample_R1_val_1.fq.gz \
    -2 sample_R2_val_2.fq.gz \
    --output_dir bismark_output/
```python

Key output: `sample_R1_val_1_bismark_bt2_pe.bam`

### Step 4: Deduplication

PCR duplicates are problematic in bisulfite sequencing because two reads mapping to the same position could represent the same molecule (since bisulfite treatment reduces sequence complexity). Deduplication is essential for WGBS (less critical for RRBS, where many reads genuinely overlap).

```bash
deduplicate_bismark -p sample_R1_val_1_bismark_bt2_pe.bam
```python

### Step 5: Methylation extraction

```bash
bismark_methylation_extractor \
    --paired-end \
    --CpG \
    --CHG \
    --CHH \
    --comprehensive \
    --cytosine_report \
    --genome_folder /path/to/genome/hg38/ \
    sample_R1_val_1_bismark_bt2_pe.deduplicated.bam
```python

**Output file formats:**
- `*.bismark.cov.gz`: One line per CpG — chromosome, start, end, % methylated, methylated count, unmethylated count
- `*.CpG_report.txt.gz`: Full genome coverage, including sites with zero coverage
- `CpX_context_*`: Separate files for CpG, CHG, CHH contexts

### Methylation contexts

| Context | Definition | Mammalian prevalence | Notes |
|---|---|---|---|
| CpG | CG dinucleotide | ~70% methylated | Most studied; heritable |
| CHG | Where H = A, T, or C | <1% in somatic cells | Common in plants |
| CHH | H = A, T, or C on both strands | <1% in somatic cells | Non-CpG methylation in neurons |

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

np.random.seed(42)

# Simulate a Bismark coverage file: columns are chromosome, start, end,
# pct_methylated, count_methylated, count_unmethylated
# We model three genomic contexts with realistic methylation distributions:
# CpG islands (promoters): low methylation (10)
# Gene bodies: moderate methylation (60)
# Intergenic / repeat regions: high methylation (85)

n_islands   = 3000   # CpG island CpGs
n_genebody  = 8000   # gene body CpGs
n_intergenic = 9000  # repeat / intergenic CpGs

def simulate_cpgs(n, mu_beta, sigma_beta, context_label, coverage_mean=20):
    """Simulate CpG sites with beta-distributed methylation and negative-binomial coverage."""
    # Beta parameters from mean/variance
    alpha = mu_beta * ((mu_beta * (1 - mu_beta)) / sigma_beta**2 - 1)
    beta_param = (1 - mu_beta) * ((mu_beta * (1 - mu_beta)) / sigma_beta**2 - 1)
    alpha = max(alpha, 0.5)
    beta_param = max(beta_param, 0.5)

    betas = np.random.beta(alpha, beta_param, n)
    betas = np.clip(betas, 0, 1)

    # Coverage: negative binomial (overdispersed counts typical of sequencing)
    coverage = np.random.negative_binomial(n=5, p=5/(5+coverage_mean), size=n)
    coverage = np.maximum(coverage, 1)   # at least 1 read

    count_M = np.round(betas * coverage).astype(int)
    count_U = coverage - count_M

    return pd.DataFrame({
        'beta': betas,
        'coverage': coverage,
        'count_M': count_M,
        'count_U': np.maximum(count_U, 0),
        'context': context_label
    })

df_island    = simulate_cpgs(n_islands,    mu_beta=0.08, sigma_beta=0.08, context_label='CpG Island')
df_genebody  = simulate_cpgs(n_genebody,   mu_beta=0.60, sigma_beta=0.18, context_label='Gene Body')
df_intergenic = simulate_cpgs(n_intergenic, mu_beta=0.84, sigma_beta=0.12, context_label='Intergenic/Repeat')

cpg_data = pd.concat([df_island, df_genebody, df_intergenic], ignore_index=True)

# Recompute beta from count data (as Bismark does)
cpg_data['beta_from_counts'] = cpg_data['count_M'] / (cpg_data['count_M'] + cpg_data['count_U'])

print(f"Total CpG sites simulated: {len(cpg_data):,}")
print(f"\nPer-context summary:")
print(cpg_data.groupby('context')['beta_from_counts'].describe().round(3))
print(f"\nFormula:  beta = count_M / (count_M + count_U)")
```python

## The Beta Value: Quantifying Methylation

The **beta value** (β) is the fraction of methylated reads at a CpG site:

$$\beta = \frac{M}{M + U}$$

where M = count of methylated reads, U = count of unmethylated reads. Beta values range from 0 (completely unmethylated) to 1 (fully methylated).

An alternative representation is the **M-value** (logit transform), preferred for statistical testing because it is more homoscedastic (constant variance across the [0,1] range):

$$M = \log_2\!\left(\frac{\beta + \varepsilon}{1 - \beta + \varepsilon}\right)$$

**Coverage filtering** is crucial: low-coverage CpGs produce unreliable beta values due to sampling noise. A common threshold is ≥ 10× coverage (some strict analyses require ≥ 20×).

Below we visualize the genome-wide beta distribution and show why context matters.

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

colors = {'CpG Island': '#2196F3', 'Gene Body': '#4CAF50', 'Intergenic/Repeat': '#FF5722'}

# Panel 1: Beta value distributions per context
for ctx, grp in cpg_data.groupby('context'):
    axes[0].hist(grp['beta_from_counts'], bins=50, alpha=0.65,
                 color=colors[ctx], label=ctx, density=True)
axes[0].set_xlabel('Beta value (methylation fraction)', fontsize=11)
axes[0].set_ylabel('Density', fontsize=11)
axes[0].set_title('Genome-wide Methylation Distribution\nby Genomic Context', fontsize=11)
axes[0].legend(fontsize=9)
axes[0].set_xlim(0, 1)

# Panel 2: Coverage distribution
cov_clip = np.clip(cpg_data['coverage'], 0, 60)
axes[1].hist(cov_clip, bins=40, color='#9C27B0', alpha=0.75, edgecolor='white')
axes[1].axvline(10, color='red', linestyle='--', linewidth=1.5, label='Min coverage (10×)')
axes[1].set_xlabel('Coverage (reads per CpG)', fontsize=11)
axes[1].set_ylabel('Number of CpG sites', fontsize=11)
axes[1].set_title('Read Coverage Distribution\n(negative binomial model)', fontsize=11)
axes[1].legend(fontsize=9)

# Panel 3: Beta vs M-value scatter
eps = 0.01
sample = cpg_data.sample(2000, random_state=1)
m_values = np.log2((sample['beta_from_counts'] + eps) / (1 - sample['beta_from_counts'] + eps))
scatter_colors = [colors[c] for c in sample['context']]
axes[2].scatter(sample['beta_from_counts'], m_values, c=scatter_colors, alpha=0.3, s=5)
axes[2].set_xlabel('Beta value', fontsize=11)
axes[2].set_ylabel('M-value (logit)', fontsize=11)
axes[2].set_title('Beta vs M-value Transformation\n(logit scale more homoscedastic)', fontsize=11)
# Add legend patches
patches = [mpatches.Patch(color=v, label=k) for k, v in colors.items()]
axes[2].legend(handles=patches, fontsize=8)

plt.tight_layout()
plt.savefig('wgbs_methylation_overview.png', dpi=120, bbox_inches='tight')
plt.show()
print("Figure saved.")
```python

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
