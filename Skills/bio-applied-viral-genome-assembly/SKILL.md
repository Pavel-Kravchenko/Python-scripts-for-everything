---
name: bio-applied-viral-genome-assembly
description: "**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/01_viral_genome_assembly.ipynb"
---

# Viral Genome Assembly and Variant Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/37_Virology_Bioinformatics/01_viral_genome_assembly.ipynb`*

# Viral Genome Assembly and Variant Analysis

**Tier 3 — Applied Bioinformatics | Module 37 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 17 (Genome Assembly)*

---

**By the end of this notebook you will be able to:**
1. Describe unique challenges of viral sequencing: quasispecies, high mutation rate, small genomes
2. Perform reference-guided and de novo viral genome assembly
3. Call intra-host variants and minor alleles with LoFreq or iVar
4. Annotate viral variants using SnpEff with viral reference databases
5. Construct consensus genomes for phylogenetic analysis



**Key resources:**
- [iVar documentation](https://andersen-lab.github.io/ivar/html/)
- [LoFreq documentation](https://csb5.github.io/lofreq/)
- [NCBI Virus database](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/)
- [Viral Bioinformatics Resource Center](https://virology.uvic.ca/)

## 1. Viral Sequencing Strategies

Viral genomes are tiny (5–30 kb), diverse, and evolving rapidly. This creates unique bioinformatics challenges:

| Strategy | Description | Use case |
|---|---|---|
| Amplicon (e.g., ARTIC) | Tiling PCR amplicons from known reference | SARS-CoV-2, HIV, Influenza |
| Metagenomic (shotgun) | Unbiased sequencing, host depletion needed | Unknown/novel pathogens |
| Hybrid capture | Probe-based enrichment | Low-titer samples, hepatitis |

**Viral genome types covered:**
- **ssRNA+**: SARS-CoV-2, Dengue, Zika (positive-sense RNA)
- **ssRNA−**: Influenza, Ebola (negative-sense RNA, needs anti-sense reads)
- **dsRNA**: Rotavirus (segmented — 11 segments must be assembled separately)
- **dsDNA**: Adenovirus, Herpesviruses (more stable, lower mutation rate)

**Key differences from bacterial/human assembly:**
- Quasispecies: populations of related sequences, not a single clone
- Coverage depth requirements: ≥20× for consensus, ≥200× for minority variants
- Reference genome almost always available → reference-guided preferred

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ----- Simulate a SARS-CoV-2-like coverage profile -----
np.random.seed(42)
genome_len = 29903  # SARS-CoV-2 reference length (bp)

# ARTIC V4 amplicons: alternating pools, some dropout regions
base_cov = np.random.negative_binomial(n=20, p=0.4, size=genome_len).astype(float) * 10
# Simulate amplicon dropout around primer binding sites
dropout_sites = [3000, 7500, 12000, 18500, 24000]
for site in dropout_sites:
    w = 200
    base_cov[max(0, site-w):site+w] *= 0.05

coverage = np.convolve(base_cov, np.ones(50)/50, mode='same')  # smooth

fig, axes = plt.subplots(3, 1, figsize=(14, 9))

# Panel 1: Coverage depth
axes[0].fill_between(range(genome_len), coverage, alpha=0.6, color='steelblue')
axes[0].axhline(20, color='red', linestyle='--', linewidth=1.5, label='20× threshold (consensus)')
axes[0].axhline(200, color='orange', linestyle='--', linewidth=1.5, label='200× threshold (variant calling)')
axes[0].set_ylabel('Coverage depth (×)')
axes[0].set_title('SARS-CoV-2 Amplicon Sequencing Coverage Profile')
axes[0].legend(fontsize=9)
axes[0].set_xlim(0, genome_len)

# Panel 2: Genome annotation
gene_coords = {
    'ORF1ab': (266, 21555), 'S': (21563, 25384), 'ORF3a': (25393, 26220),
    'E': (26245, 26472), 'M': (26523, 27191), 'N': (28274, 29533)
}
colors_gene = plt.cm.Set2(np.linspace(0, 1, len(gene_coords)))
for (gene, (start, end)), col in zip(gene_coords.items(), colors_gene):
    axes[1].barh(0, end - start, left=start, height=0.4, color=col, label=gene)
    axes[1].text((start + end) / 2, 0, gene, ha='center', va='center', fontsize=7, fontweight='bold')
axes[1].set_xlim(0, genome_len)
axes[1].set_yticks([])
axes[1].set_title('SARS-CoV-2 Genome Organization')
axes[1].set_xlabel('Genome position (bp)')

# Panel 3: Base quality / completeness
completeness_20x = np.sum(coverage >= 20) / genome_len * 100
completeness_200x = np.sum(coverage >= 200) / genome_len * 100
bars = axes[2].bar(['≥20× (consensus)', '≥200× (variant calling)'],
                    [completeness_20x, completeness_200x],
                    color=['steelblue', 'orange'], alpha=0.8)
axes[2].set_ylabel('Genome completeness (%)')
axes[2].set_title('Genome Coverage Completeness')
axes[2].set_ylim(0, 105)
for bar, val in zip(bars, [completeness_20x, completeness_200x]):
    axes[2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig('viral_coverage.png', dpi=100, bbox_inches='tight')
plt.show()
print(f"Genome completeness at 20×: {completeness_20x:.1f}%")
print(f"Genome completeness at 200×: {completeness_200x:.1f}%")
```

## 2. Reference-Guided Assembly Pipeline

### Step-by-step workflow for ARTIC amplicon data

```
Raw FASTQ
    │
    ▼
[1] Quality & adapter trimming (fastp or trim_galore)
    │
    ▼
[2] Primer trimming (iVar trim / Trimmomatic for ARTIC)
    │
    ▼
[3] Read alignment to reference genome (BWA-MEM2 or Minimap2)
    │
    ▼
[4] Sort and index BAM (samtools sort/index)
    │
    ▼
[5] Consensus calling (iVar consensus / samtools mpileup)
    │
    ▼
[6] Variant calling (LoFreq, iVar variants, or DeepVariant)
    │
    ▼
[7] Annotation & QC (VADR, Nextclade, coverage stats)
```

### Key tools

| Step | Tool | Notes |
|---|---|---|
| Alignment | `minimap2 -ax sr` | Short reads; `ont2d` for Nanopore |
| Primer trim | `iVar trim -p` | Soft-clips ARTIC primer regions |
| Consensus | `iVar consensus -t 0.75` | Min freq threshold 0.75 |
| Variant calling | `LoFreq call` | Poisson model for low-frequency variants |
| Annotation | `Nextclade run` | Clade assignment + quality flags |

```python
# ----- Simulate: Reference mapping statistics -----
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(7)
n_samples = 20

samples = [f'Sample_{i:02d}' for i in range(1, n_samples+1)]
# Simulate varying QC metrics
total_reads = np.random.randint(50_000, 500_000, size=n_samples)
mapped_pct = np.random.uniform(65, 99, size=n_samples)
median_depth = np.random.randint(50, 5000, size=n_samples)
genome_completeness = np.clip(0.5 + 0.5 * (np.log10(median_depth) / np.log10(5000)), 0, 1) * 100
# Add some failures
genome_completeness[[3, 11]] = np.random.uniform(20, 50, 2)

df = pd.DataFrame({
    'Sample': samples,
    'Total_reads': total_reads,
    'Mapped_pct': mapped_pct,
    'Median_depth': median_depth,
    'Completeness_pct': genome_completeness
})

# QC pass/fail
df['QC_pass'] = (df['Mapped_pct'] >= 80) & (df['Completeness_pct'] >= 90)

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Mapped reads %
colors = ['green' if p else 'red' for p in df['QC_pass']]
axes[0,0].bar(range(n_samples), df['Mapped_pct'], color=colors, alpha=0.8)
axes[0,0].axhline(80, color='black', linestyle='--', label='80% threshold')
axes[0,0].set_title('Reads Mapped to Reference (%)')
axes[0,0].set_ylabel('%'); axes[0,0].set_xticks([])
axes[0,0].legend()

# Median depth
axes[0,1].bar(range(n_samples), df['Median_depth'], color=colors, alpha=0.8)
axes[0,1].axhline(200, color='orange', linestyle='--', label='200× (variant calling)')
axes[0,1].set_title('Median Coverage Depth (×)')
axes[0,1].set_ylabel('Depth (×)'); axes[0,1].set_xticks([])
axes[0,1].legend()

# Genome completeness
axes[1,0].barh(df['Sample'], df['Completeness_pct'], color=colors, alpha=0.8)
axes[1,0].axvline(90, color='black', linestyle='--', label='90% threshold')
axes[1,0].set_title('Genome Completeness (%)')
axes[1,0].set_xlabel('%')
axes[1,0].legend()

# Scatter: depth vs completeness
axes[1,1].scatter(df['Median_depth'], df['Completeness_pct'],
                  c=['green' if p else 'red' for p in df['QC_pass']], s=60, alpha=0.8)
axes[1,1].axvline(200, color='orange', linestyle='--', alpha=0.6)
axes[1,1].axhline(90, color='black', linestyle='--', alpha=0.6)
axes[1,1].set_xscale('log')
axes[1,1].set_xlabel('Median Coverage Depth (×, log scale)')
axes[1,1].set_ylabel('Genome Completeness (%)')
axes[1,1].set_title('Depth vs Completeness')
pass_patch = mpatches.Patch(color='green', label='QC Pass')
fail_patch = mpatches.Patch(color='red', label='QC Fail')
axes[1,1].legend(handles=[pass_patch, fail_patch])

plt.tight_layout()
plt.show()
print(f"\nQC pass rate: {df['QC_pass'].sum()}/{n_samples} samples ({df['QC_pass'].mean()*100:.0f}%)")
print(df[['Sample','Mapped_pct','Median_depth','Completeness_pct','QC_pass']].to_string(index=False))
```

## 3. Intra-Host Variant Calling (Quasispecies)

Within a single infected individual, the virus exists as a **quasispecies** — a cloud of closely related sequence variants.

### Why this matters
- Drug resistance mutations may be present at 1–5% frequency
- Minority variants can become dominant after immune/drug selection
- Quasispecies structure informs transmission bottleneck analysis

### Calling minority variants

| Tool | Model | Min freq | Notes |
|---|---|---|---|
| **LoFreq** | Poisson + Bonferroni | 0.5% | Best for Illumina amplicons |
| **iVar** | Binomial exact | Configurable | ARTIC pipeline standard |
| **DeepVariant** | CNN | ~5% | High accuracy, slow |
| **VarScan2** | Fisher exact | Configurable | Somatic/viral mode |

### Filters to apply
1. Min read depth ≥ 100× at variant site
2. Min allele frequency ≥ 1% (or 3% for low-depth)
3. Strand bias filter (Fisher p-value > 0.001)
4. Base quality ≥ 25 at variant position

```python
# ----- Simulate: Intra-host viral variant calling -----
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

np.random.seed(42)

genome_len = 29903
total_variants = 45

# Simulate variant positions and allele frequencies
positions = np.sort(np.random.choice(genome_len, total_variants, replace=False))
# Mix of high-freq consensus variants and low-freq minority variants
freqs = np.concatenate([
    np.random.uniform(0.85, 0.99, 25),   # near-fixed variants
    np.random.uniform(0.01, 0.15, 12),    # minority variants (1-15%)
    np.random.uniform(0.30, 0.70, 8),     # intermediate (co-infections or drift)
])
np.random.shuffle(freqs)

# Depths
depths = np.random.randint(300, 5000, total_variants)

# Strand bias (AF on forward strand)
fwd_af = freqs + np.random.normal(0, 0.05, total_variants)
fwd_af = np.clip(fwd_af, 0.01, 0.99)
strand_bias_p = np.array([
    stats.fisher_exact([[int(d*af*ff), int(d*af*(1-ff))],
                        [int(d*(1-af)*ff), int(d*(1-af)*(1-ff))]])[1]
    for d, af, ff in zip(depths, freqs, fwd_af)
])

# Map to gene regions
gene_coords = {
    'ORF1ab': (266, 21555), 'S': (21563, 25384), 'ORF3a': (25393, 26220),
    'E': (26245, 26472), 'M': (26523, 27191), 'N': (28274, 29533)
}
def get_gene(pos):
    for gene, (s, e) in gene_coords.items():
        if s <= pos <= e:
            return gene
    return 'IGR'

df_var = pd.DataFrame({
    'Position': positions,
    'Alt_freq': freqs,
    'Depth': depths,
    'Strand_bias_p': strand_bias_p,
    'Gene': [get_gene(p) for p in positions],
})

# Apply filters
df_var['Pass_depth'] = df_var['Depth'] >= 100
df_var['Pass_freq'] = df_var['Alt_freq'] >= 0.01
df_var['Pass_strand'] = df_var['Strand_bias_p'] > 0.001
df_var['PASS'] = df_var['Pass_depth'] & df_var['Pass_freq'] & df_var['Pass_strand']

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# AF distribution
axes[0,0].hist(df_var.loc[df_var['PASS'], 'Alt_freq'], bins=25, color='steelblue', alpha=0.8)
axes[0,0].axvline(0.05, color='red', linestyle='--', label='5% threshold')
axes[0,0].set_xlabel('Alternate Allele Frequency')
axes[0,0].set_ylabel('Count')
axes[0,0].set_title('Distribution of Variant Allele Frequencies')
axes[0,0].legend()

# Genome scatter (AF vs position)
passed = df_var[df_var['PASS']]
failed = df_var[~df_var['PASS']]
axes[0,1].scatter(failed['Position'], failed['Alt_freq'], c='lightgray', s=30, label='Filtered', zorder=1)
axes[0,1].scatter(passed['Position'], passed['Alt_freq'], c='steelblue', s=50, alpha=0.8, label='PASS', zorder=2)
axes[0,1].axhline(0.05, color='red', linestyle='--', alpha=0.7)
axes[0,1].set_xlabel('Genome Position')
axes[0,1].set_ylabel('Alt Allele Freq')
axes[0,1].set_title('Variant Positions Across Genome')
axes[0,1].legend()

# Gene breakdown
gene_counts = passed.groupby('Gene').size()
axes[1,0].bar(gene_counts.index, gene_counts.values, color=plt.cm.Set2(np.linspace(0,1,len(gene_counts))))
axes[1,0].set_xlabel('Gene')
axes[1,0].set_ylabel('Number of Variants (PASS)')
axes[1,0].set_title('Variants per Gene Region')
plt.setp(axes[1,0].xaxis.get_majorticklabels(), rotation=30)

# Depth vs AF colored by strand bias
sc = axes[1,1].scatter(df_var['Depth'], df_var['Alt_freq'],
                        c=-np.log10(df_var['Strand_bias_p'] + 1e-10),
                        cmap='RdYlGn_r', s=50, alpha=0.8)
axes[1,1].axhline(0.01, color='blue', linestyle='--', label='1% AF cutoff')
axes[1,1].axvline(100, color='orange', linestyle='--', label='100× depth cutoff')
axes[1,1].set_xlabel('Coverage Depth')
axes[1,1].set_ylabel('Alt Allele Frequency')
axes[1,1].set_title('Depth vs AF (color = strand bias −log10p)')
plt.colorbar(sc, ax=axes[1,1], label='Strand bias −log10(p)')
axes[1,1].legend()

plt.tight_layout()
plt.show()
print(f"Total variants: {len(df_var)}")
print(f"PASS variants: {df_var['PASS'].sum()}")
print(f"Minority variants (AF<5%, PASS): {((df_var['Alt_freq']<0.05) & df_var['PASS']).sum()}")
```
