---
name: bio-applied-qiime2-16s
description: "16S rRNA amplicon analysis with QIIME2: DADA2 denoising, taxonomy assignment, alpha/beta diversity, and differential abundance. Use when analyzing 16S microbiome data."
tool_type: python
primary_tool: QIIME2
---

# QIIME2 16S Amplicon Workflow

All data stored as **artifacts** (`.qza`) with provenance tracking, or **visualizations** (`.qzv`) at https://view.qiime2.org.

## Import Reads

Manifest CSV maps sample IDs to absolute FASTQ paths:

```csv
sample-id,absolute-filepath,direction
Sample1,/data/Sample1_R1.fastq.gz,forward
Sample1,/data/Sample1_R2.fastq.gz,reverse
```

```bash
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.csv \
    --output-path reads.qza \
    --input-format PairedEndFastqManifestPhred33

# Check quality to determine DADA2 truncation lengths
qiime demux summarize --i-data reads.qza --o-visualization reads_summary.qzv
# Open at view.qiime2.org — truncate where median drops below Q20
```

```python
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)
read_length = 250
positions = np.arange(1, read_length + 1)

qual_mean_fwd = 37 - 0.03 * positions - 3 * np.exp(-(read_length - positions) / 30)
qual_mean_rev = 36 - 0.06 * positions - 5 * np.exp(-(read_length - positions) / 20)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, (mean_q, label) in zip(axes, [
    (qual_mean_fwd, 'Forward reads'), (qual_mean_rev, 'Reverse reads')
]):
    ax.plot(positions, np.clip(mean_q, 0, 40), color='steelblue', lw=1.5)
    ax.axhline(20, color='red', linestyle='--', lw=1, label='Q20')
    trunc_pos = positions[np.where(mean_q < 25)[0][0]] if (mean_q < 25).any() else read_length
    ax.axvline(trunc_pos, color='green', lw=2, label=f'Truncate at {trunc_pos}')
    ax.set_title(f'QIIME2 demux: {label}'); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
print('Recommended: --p-trunc-len-f 220 --p-trunc-len-r 200')
```

## Denoising with DADA2

DADA2 produces **ASVs** (exact error-corrected sequences, single-nucleotide resolution) vs. OTUs (97% similarity clusters). ASVs are reproducible across studies.

```bash
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs reads.qza \
    --p-trim-left-f 13 --p-trim-left-r 13 \
    --p-trunc-len-f 220 --p-trunc-len-r 200 \
    --p-n-threads 0 \
    --o-table feature_table.qza \
    --o-representative-sequences rep_seqs.qza \
    --o-denoising-stats stats.qza
```

**Truncation rule:** forward + reverse must overlap ≥20 bp after truncation. For V3-V4 (~460 bp amplicon): 220 + 200 = 420 bp combined gives sufficient overlap. After denoising verify ≥75% reads pass filtering and ≥60% merge.

```python
rng = np.random.default_rng(7)
n_samples = 16
input_reads = rng.integers(20_000, 80_000, size=n_samples)
filtered_pct = rng.uniform(0.82, 0.95, size=n_samples)
merged_pct   = rng.uniform(0.78, 0.92, size=n_samples)
non_chimeric = rng.uniform(0.93, 0.99, size=n_samples)

filtered = (input_reads * filtered_pct).astype(int)
merged   = (filtered * merged_pct).astype(int)
final    = (merged * non_chimeric).astype(int)
print(f'Mean retention: {(final / input_reads).mean():.1%}')
```

## Taxonomic Classification with SILVA

SILVA 138 is the standard 16S reference. Use the pre-trained Naive Bayes classifier matched to your primer pair and read length.

```bash
qiime feature-classifier classify-sklearn \
    --i-classifier silva138_nb_515_806_classifier.qza \
    --i-reads rep_seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs 1

qiime taxa barplot \
    --i-table feature_table.qza --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv --o-visualization taxa_barplot.qzv

# Remove host/contaminant sequences
qiime taxa filter-table \
    --i-table feature_table.qza --i-taxonomy taxonomy.qza \
    --p-exclude Mitochondria,Chloroplast,Eukaryota \
    --o-filtered-table feature_table_filtered.qza
```

**Confidence threshold:** default `--p-confidence 0.7` means assignments below 70% bootstrap confidence are reported as unclassified. Species-level accuracy is poor; genus-level is reliable.

```python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.cm as cm

rng = np.random.default_rng(20)
taxa = {
    'g__Bacteroides': (0.18, 0.08), 'g__Faecalibacterium': (0.10, 0.06),
    'g__Prevotella': (0.08, 0.07),  'g__Bifidobacterium': (0.07, 0.04),
    'g__Ruminococcus': (0.06, 0.04),'g__Akkermansia': (0.05, 0.04),
    'g__Blautia': (0.05, 0.03),     'g__Lachnospiraceae': (0.06, 0.04),
    'Other': (0.27, 0.05),
}
n_samples = 16
conditions = ['Healthy'] * 8 + ['IBD'] * 8
compositions = np.zeros((n_samples, len(taxa)))
for j, (genus, (mean, std)) in enumerate(taxa.items()):
    base = rng.normal(mean, std, n_samples)
    if genus == 'g__Faecalibacterium':
        base[8:] *= 0.4  # known reduction in IBD
    compositions[:, j] = np.clip(base, 0.001, 1)
compositions /= compositions.sum(axis=1, keepdims=True)
```

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
