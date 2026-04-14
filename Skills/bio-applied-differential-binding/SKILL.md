---
name: bio-applied-differential-binding
description: "Differential binding analysis for ChIP-seq: DiffBind workflow, consensus peaks, normalization, and MA/volcano plots. Use when comparing ChIP-seq signal between conditions."
tool_type: python
primary_tool: Matplotlib
---

# Differential Binding  Peak Annotation

- [Bioconductor DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf)
- [Bioconductor ChIPseeker vignette](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

try:
    import matplotlib_venn as venn_mod
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11
print('Setup complete.')
```python

## Differential Binding with DiffBind (R)

**DiffBind** quantifies reads in a consensus peak set across all samples and applies DESeq2 or edgeR to identify peaks that gain or lose binding between experimental conditions. The **consensus peak set** is constructed as the union of peaks called in at least 2 samples (default `minMembers = 2`), ensuring that only reproducible peaks are included in the analysis. This reduces noise from spurious peaks unique to a single replicate. The resulting matrix — samples × consensus peaks with read counts — is then processed identically to an RNA-seq count matrix.

The **DiffBind sample sheet** is a CSV file with one row per sample and columns specifying: `SampleID`, `Condition` (e.g., Control/Treatment), `Replicate` (integer), `bamReads` (path to deduplicated BAM), `ControlID` (matching Input sample ID), `bamControl` (path to Input BAM), `Peaks` (path to narrowPeak/broadPeak file), and `PeakCaller` (e.g., `macs`). The `dba()` function reads this sheet; `dba.count()` tallies reads in the consensus peak set; `dba.normalize()` applies TMM (edgeR) or RLE (DESeq2) normalization to correct for library size differences across samples.

`dba.contrast()` defines the pairwise comparison, referencing the `DBA_CONDITION` metadata column. `dba.analyze()` runs DESeq2 by default (set `method = DBA_EDGER` to switch). `dba.report()` returns a GRanges object with log₂ fold change, p-value, and FDR for each consensus peak; the `th` and `fold` arguments filter to a significance threshold. Volcano plots from `dba.plotVolcano()` and PCA plots from `dba.plotPCA()` quickly reveal whether replicates cluster by condition and whether differential binding is driven by a small number of peaks or genome-wide.


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Simulate DiffBind output (as would come from dba.report())
np.random.seed(42)
n_peaks = 45000  # consensus peak set size

log2fc = np.concatenate([
    np.random.normal(0, 0.5, 42000),   # unchanged
    np.random.normal(2.5, 0.6, 1500),  # gained binding (treatment)
    np.random.normal(-2.2, 0.6, 1500), # lost binding (treatment)
])

# p-values: most non-significant, some very significant
pval_base = np.random.exponential(0.3, n_peaks)
pval_base[:3000] = np.random.uniform(1e-30, 0.001, 3000)  # true positives
pval = np.clip(pval_base, 1e-50, 1)
fdr = np.minimum(pval * n_peaks / (np.argsort(np.argsort(pval)) + 1), 1)

df_diffbind = pd.DataFrame({
    'log2FC': log2fc[:n_peaks],
    'pvalue': pval,
    'FDR': fdr,
    'Peak': [f'peak_{i}' for i in range(n_peaks)]
})

df_diffbind['neg_log10_FDR'] = -np.log10(df_diffbind['FDR'].clip(1e-50))
df_diffbind['Significant'] = (
    (df_diffbind['FDR'] < 0.05) & (df_diffbind['log2FC'].abs() > np.log2(1.5))
)
df_diffbind['Direction'] = 'Unchanged'
df_diffbind.loc[
    df_diffbind['Significant'] & (df_diffbind['log2FC'] > 0), 'Direction'
] = 'Gained'
df_diffbind.loc[
    df_diffbind['Significant'] & (df_diffbind['log2FC'] < 0), 'Direction'
] = 'Lost'

n_gained = (df_diffbind['Direction'] == 'Gained').sum()
n_lost   = (df_diffbind['Direction'] == 'Lost').sum()
print(f'Consensus peaks: {n_peaks:,}')
print(f'Gained binding sites: {n_gained:,} ({100*n_gained/n_peaks:.1f}%)')
print(f'Lost binding sites:   {n_lost:,} ({100*n_lost/n_peaks:.1f}%)')
print(f'Unchanged:            {(df_diffbind["Direction"] == "Unchanged").sum():,}')

color_map = {'Gained': 'coral', 'Lost': 'steelblue', 'Unchanged': 'lightgray'}
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for direction, grp in df_diffbind.groupby('Direction'):
    axes[0].scatter(
        grp['log2FC'], grp['neg_log10_FDR'],
        c=color_map[direction],
        s=3 if direction == 'Unchanged' else 6,
        alpha=0.5 if direction == 'Unchanged' else 0.8,
        label=f'{direction} (n={len(grp):,})'
    )

axes[0].axhline(-np.log10(0.05), color='black', linestyle='--',
                lw=1, alpha=0.7, label='FDR = 0.05')
axes[0].axvline(np.log2(1.5), color='gray', linestyle=':', lw=1)
axes[0].axvline(-np.log2(1.5), color='gray', linestyle=':', lw=1)
axes[0].set_xlabel('log2 Fold Change (Treatment / Control)')
axes[0].set_ylabel('-log10(FDR)')
axes[0].set_title('DiffBind Volcano Plot\nDifferential ChIP-seq Binding')
axes[0].legend(markerscale=2)
axes[0].set_xlim(-8, 8)

# MA plot
sig = df_diffbind['Significant']
mean_signal = np.random.exponential(5, n_peaks) + 1
axes[1].scatter(np.log2(mean_signal[~sig]), df_diffbind['log2FC'][~sig],
                c='lightgray', s=3, alpha=0.5, label='Unchanged')
axes[1].scatter(
    np.log2(mean_signal[sig & (df_diffbind['Direction'] == 'Gained')]),
    df_diffbind['log2FC'][sig & (df_diffbind['Direction'] == 'Gained')],
    c='coral', s=6, alpha=0.8, label='Gained'
)
axes[1].scatter(
    np.log2(mean_signal[sig & (df_diffbind['Direction'] == 'Lost')]),
    df_diffbind['log2FC'][sig & (df_diffbind['Direction'] == 'Lost')],
    c='steelblue', s=6, alpha=0.8, label='Lost'
)
axes[1].axhline(0, color='black', lw=0.8)
axes[1].set_xlabel('log2 Mean Signal')
axes[1].set_ylabel('log2 Fold Change')
axes[1].set_title('MA Plot\n(Mean vs. Fold Change)')
axes[1].legend(markerscale=2)

plt.tight_layout()
plt.show()
```python

## Peak Annotation with ChIPseeker

**ChIPseeker** annotates each ChIP-seq peak to the nearest genomic feature using a **TxDb** (transcript database) annotation package for the relevant genome assembly (e.g., `TxDb.Hsapiens.UCSC.hg38.knownGene`). It assigns each peak a feature category based on its position relative to annotated transcripts: **Promoter** (within the defined window upstream/downstream of TSS), 5′ UTR, 3′ UTR, Exon, Intron, Downstream, or Distal Intergenic. A peak is assigned to the nearest TSS by default, but multiple overlapping transcripts at the same locus may result in multiple annotations.

The `annotatePeak()` function accepts a GRanges object or a BED/narrowPeak file path. The `tssRegion = c(-2000, 200)` argument defines the promoter window as 2 kb upstream to 200 bp downstream of the TSS (reflecting that TATA-box and initiator elements are typically within 200 bp downstream). Setting `annoDb = 'org.Hs.eg.db'` automatically adds gene symbols, Entrez IDs, and Ensembl IDs to the annotation table, which are required for downstream enrichment analyses with clusterProfiler.

Visualisation in ChIPseeker conveys biology directly: `plotAnnoPie()` and `plotAnnoBar()` show the proportion of peaks in each genomic category. TF ChIP-seq peaks for activating factors (CTCF, p53) are typically 25–50% promoter-associated. Enhancer marks such as H3K27ac skew strongly toward **intronic and intergenic** regions (> 60% non-promoter), reflecting their enrichment at distal enhancers. `plotDistToTSS()` plots the distribution of distances from each peak to its nearest TSS, distinguishing proximal (< 1 kb) from distal regulatory elements.


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(42)
n_peaks = 12500  # TF ChIP-seq peaks

# Simulate genomic annotation distribution (TF-like: enriched at promoters)
annotation_categories = {
    'Promoter (<=1kb)':   0.28,
    'Promoter (1-2kb)':   0.09,
    'Promoter (2-3kb)':   0.04,
    "5' UTR":             0.02,
    "3' UTR":             0.03,
    '1st Exon':           0.02,
    'Other Exon':         0.03,
    '1st Intron':         0.10,
    'Other Intron':       0.18,
    'Downstream (<=300)': 0.01,
    'Distal Intergenic':  0.20,
}

categories = list(annotation_categories.keys())
proportions = np.array(list(annotation_categories.values()))
counts = (proportions * n_peaks).astype(int)

anno_df = pd.DataFrame({
    'Feature': categories, 'Count': counts, 'Proportion': proportions
})

print('Peak Annotation Summary (TF ChIP-seq)')
print('=' * 50)
print(anno_df.to_string(index=False))
print(f'\nTotal peaks: {counts.sum():,}')
print(f'Promoter peaks (<=3kb): {counts[:3].sum():,} '
      f'({100*counts[:3].sum()/counts.sum():.1f}%)')
print(f'Intergenic peaks: {counts[-1]:,} ({100*counts[-1]/counts.sum():.1f}%)')

# Simulate distance-to-TSS distribution
dist_promoter   = np.random.normal(0, 400, counts[:3].sum())
dist_intergenic = np.concatenate([
    np.random.uniform(-50000, -3000, counts[-1]//2),
    np.random.uniform(3000, 50000, counts[-1]//2)
])
dist_intronic   = np.random.uniform(-100000, 100000, counts[7]+counts[8])
dist_all = np.concatenate([dist_promoter, dist_intronic, dist_intergenic])

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

colors_pie = plt.cm.Set3(np.linspace(0, 1, len(categories)))
wedge_props = {'edgecolor': 'white', 'linewidth': 1.5}
axes[0].pie(counts, labels=None, colors=colors_pie,
            autopct=lambda p: f'{p:.1f}%' if p > 3 else '',
            wedgeprops=wedge_props, startangle=90)
axes[0].set_title('Genomic Feature Distribution\n(ChIPseeker annotatePeak)')
axes[0].legend(categories, loc='lower left', bbox_to_anchor=(-0.3, -0.1),
               fontsize=7, ncol=2)

axes[1].hist(np.clip(dist_all, -50000, 50000), bins=80,
             color='steelblue', edgecolor='white', alpha=0.8)
axes[1].axvline(0, color='red', lw=1.5, linestyle='--', label='TSS')
axes[1].axvline(-2000, color='orange', lw=1, linestyle=':', label='+-2 kb promoter')
axes[1].axvline(2000, color='orange', lw=1, linestyle=':')
axes[1].set_xlabel('Distance to TSS (bp)')
axes[1].set_ylabel('Number of peaks')
axes[1].set_title('Distance to Nearest TSS\n(plotDistToTSS)')
axes[1].legend()
axes[1].set_xlim(-50000, 50000)

plt.tight_layout()
plt.show()
```python

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
