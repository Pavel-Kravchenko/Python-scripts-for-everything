---
name: bio-applied-wgbs-bismark
description: WGBS/RRBS Processing with Bismark with Bismark
tool_type: python
primary_tool: Bismark
---

# WGBS/RRBS Processing with Bismark

**References:** [Bismark docs](https://felixkrueger.github.io/Bismark/) | [methylKit vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/) | [ENCODE WGBS pipeline](https://www.encodeproject.org/wgbs/)

## DNA Methylation Biology

5mC added by DNMT3A/3B (de novo) and DNMT1 (maintenance). Removed by TET1/2/3 via 5hmC oxidation. Methylation occurs almost exclusively at **CpG dinucleotides** in mammals (~70–80% methylated).

| Feature | Methylation | Function |
|---|---|---|
| CpG islands (CGIs) | Mostly unmethylated | Protect promoters from silencing |
| CGI shores (±2 kb) | Variable, tissue-specific | Major DMR site |
| Gene bodies | Moderately methylated | Active transcription |
| Repeats / transposons | Heavily methylated | Silences parasitic elements |
| Imprinted regions | Hemi-methylated | Monoallelic expression |

CGIs: ≥200 bp, CpG obs/exp > 0.6, GC > 50%. ~1% of genome, overlaps ~70% of promoters.

**Disease:** Cancer → CGI hypermethylation silences TSGs (BRCA1, MLH1, CDKN2A); global hypomethylation destabilizes genome. Aging → gradual methylation drift (basis of epigenetic clocks).

## Bisulfite Chemistry

> **Unmethylated C** + bisulfite → uracil → reads as **T**
> **Methylated C** + bisulfite → protected → reads as **C**

Conversion efficiency measured via lambda spike-in or CHH/CHG sites (non-CpG methylation rare in somatic cells). Target: >99.5% conversion.

Bisulfite creates 4 distinct strands (OT, CTOT, OB, CTOB) — standard aligners cannot handle this; use Bismark.

| Method | Coverage | Cost | Best for |
|---|---|---|---|
| WGBS | All ~28M CpGs | High (>$500/sample) | Comprehensive |
| RRBS | ~5M CpGs (CpG-enriched) | Medium | Cost-effective |
| EPIC array | 850K CpGs | Low (~$200) | Large cohorts |

## Bismark Pipeline

```bash
# Step 1: Build bisulfite genome index (once per genome)
bismark_genome_preparation /path/to/genome/hg38/

# Step 2: Quality trimming (bisulfite-specific settings)
trim_galore --paired --fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# Step 3: Alignment
bismark --genome /path/to/genome/hg38/ \
    -1 sample_R1_val_1.fq.gz -2 sample_R2_val_2.fq.gz \
    --output_dir bismark_output/

# Step 4: Deduplication (essential for WGBS; less critical for RRBS)
deduplicate_bismark -p sample_R1_val_1_bismark_bt2_pe.bam

# Step 5: Methylation extraction
bismark_methylation_extractor \
    --paired-end --CpG --CHG --CHH --comprehensive \
    --cytosine_report \
    --genome_folder /path/to/genome/hg38/ \
    sample_R1_val_1_bismark_bt2_pe.deduplicated.bam
```

**Output formats:**
- `*.bismark.cov.gz`: chr, start, end, %methylated, M_count, U_count
- `*.CpG_report.txt.gz`: full genome coverage including zero-coverage sites

**Methylation contexts:**

| Context | Definition | Mammalian prevalence |
|---|---|---|
| CpG | CG dinucleotide | ~70% methylated |
| CHG | H = A, T, or C | <1% somatic |
| CHH | H = A, T, C on both strands | <1%; common in neurons |

## Beta Values and M-Values

$$\beta = \frac{M}{M + U} \quad \text{(0 = unmethylated, 1 = fully methylated)}$$

$$M\text{-value} = \log_2\!\left(\frac{\beta + \varepsilon}{1 - \beta + \varepsilon}\right) \quad \text{(more homoscedastic, preferred for testing)}$$

Coverage filtering: require ≥10× (strict analyses use ≥20×).

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

np.random.seed(42)

def simulate_cpgs(n, mu_beta, sigma_beta, context_label, coverage_mean=20):
    alpha = mu_beta * ((mu_beta * (1 - mu_beta)) / sigma_beta**2 - 1)
    beta_param = (1 - mu_beta) * ((mu_beta * (1 - mu_beta)) / sigma_beta**2 - 1)
    alpha = max(alpha, 0.5); beta_param = max(beta_param, 0.5)
    betas = np.clip(np.random.beta(alpha, beta_param, n), 0, 1)
    coverage = np.maximum(np.random.negative_binomial(n=5, p=5/(5+coverage_mean), size=n), 1)
    count_M = np.round(betas * coverage).astype(int)
    return pd.DataFrame({
        'beta': betas, 'coverage': coverage,
        'count_M': count_M, 'count_U': np.maximum(coverage - count_M, 0),
        'context': context_label
    })

df_island    = simulate_cpgs(3000,  mu_beta=0.08, sigma_beta=0.08, context_label='CpG Island')
df_genebody  = simulate_cpgs(8000,  mu_beta=0.60, sigma_beta=0.18, context_label='Gene Body')
df_intergenic = simulate_cpgs(9000, mu_beta=0.84, sigma_beta=0.12, context_label='Intergenic/Repeat')

cpg_data = pd.concat([df_island, df_genebody, df_intergenic], ignore_index=True)
cpg_data['beta_from_counts'] = cpg_data['count_M'] / (cpg_data['count_M'] + cpg_data['count_U'])

print(cpg_data.groupby('context')['beta_from_counts'].describe().round(3))

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
colors = {'CpG Island': '#2196F3', 'Gene Body': '#4CAF50', 'Intergenic/Repeat': '#FF5722'}

for ctx, grp in cpg_data.groupby('context'):
    axes[0].hist(grp['beta_from_counts'], bins=50, alpha=0.65, color=colors[ctx], label=ctx, density=True)
axes[0].set_title('Methylation Distribution by Context'); axes[0].legend(fontsize=9)

axes[1].hist(np.clip(cpg_data['coverage'], 0, 60), bins=40, color='#9C27B0', alpha=0.75)
axes[1].axvline(10, color='red', ls='--', lw=1.5, label='Min coverage (10×)')
axes[1].set_title('Coverage Distribution'); axes[1].legend()

eps = 0.01
sample = cpg_data.sample(2000, random_state=1)
m_values = np.log2((sample['beta_from_counts'] + eps) / (1 - sample['beta_from_counts'] + eps))
axes[2].scatter(sample['beta_from_counts'], m_values,
                c=[colors[c] for c in sample['context']], alpha=0.3, s=5)
axes[2].set_title('Beta vs M-value (logit)')
axes[2].legend(handles=[mpatches.Patch(color=v, label=k) for k, v in colors.items()], fontsize=8)

plt.tight_layout()
plt.savefig('wgbs_methylation_overview.png', dpi=120, bbox_inches='tight')
plt.show()
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
