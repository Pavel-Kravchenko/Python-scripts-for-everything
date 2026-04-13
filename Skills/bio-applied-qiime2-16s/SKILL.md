---
name: bio-applied-qiime2-16s
description: "**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 4**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/04_qiime2_16s.ipynb"
primary_tool: QIIME2
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# QIIME2 16S Amplicon Workflow

*Source: Course notebook `Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/04_qiime2_16s.ipynb`*

# QIIME2 16S Amplicon Workflow

**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 4**

*Prerequisites: Module 04 (Microbial Diversity), Notebook 1 (Taxonomic Profiling)*

---

**By the end of this notebook you will be able to:**
1. Import paired-end 16S amplicon reads into QIIME2 artifacts
2. Denoise reads with DADA2 to produce an ASV feature table
3. Assign taxonomy to ASVs using a pre-trained classifier
4. Compute alpha and beta diversity metrics
5. Identify differentially abundant taxa with ANCOM-BC



**Key resources:**
- [QIIME2 Moving Pictures tutorial](https://docs.qiime2.org/2024.10/tutorials/moving-pictures/)
- [QIIME2 documentation](https://docs.qiime2.org/)
- [Galaxy Training — 16S Metagenomics](https://training.galaxyproject.org/training-material/topics/metagenomics/)

## 1. Import Reads into QIIME2

**QIIME2** (Quantitative Insights Into Microbial Ecology 2) is a comprehensive amplicon sequencing analysis framework. All data in QIIME2 is stored as **artifacts** (`.qza` files) with provenance tracking, or **visualizations** (`.qzv` files) viewable at https://view.qiime2.org.

### Manifest Format (for paired-end reads)

The manifest CSV file maps sample IDs to absolute paths of their FASTQ files:

```csv
sample-id,absolute-filepath,direction
Sample1,/data/Sample1_R1.fastq.gz,forward
Sample1,/data/Sample1_R2.fastq.gz,reverse
Sample2,/data/Sample2_R1.fastq.gz,forward
Sample2,/data/Sample2_R2.fastq.gz,reverse
```

```bash
# Import using manifest file
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.csv \
    --output-path reads.qza \
    --input-format PairedEndFastqManifestPhred33

# Check quality to determine DADA2 truncation lengths
qiime demux summarize \
    --i-data reads.qza \
    --o-visualization reads_summary.qzv
# Open reads_summary.qzv at view.qiime2.org to see per-position quality boxplots
```

**Key decision from quality visualization:** Identify position where median quality drops below Q20 — use this as the truncation length for DADA2.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate per-position quality scores (as would be seen in QIIME2 demux summarize)
read_length = 250  # V3-V4 region, 250 bp paired-end
n_reads_vis = 5000
positions = np.arange(1, read_length + 1)

# Forward reads: quality drops in last ~30 bp
qual_mean_fwd = 37 - 0.03 * positions - 3 * np.exp(-(read_length - positions) / 30)
qual_q25_fwd  = qual_mean_fwd - 2 - 0.01 * positions
qual_q75_fwd  = qual_mean_fwd + 1.5

# Reverse reads: quality drops faster (typical paired-end behavior)
qual_mean_rev = 36 - 0.06 * positions - 5 * np.exp(-(read_length - positions) / 20)
qual_q25_rev  = qual_mean_rev - 3 - 0.02 * positions
qual_q75_rev  = qual_mean_rev + 1.5

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (mean_q, q25, q75, label) in zip(axes, [
    (qual_mean_fwd, qual_q25_fwd, qual_q75_fwd, 'Forward reads'),
    (qual_mean_rev, qual_q25_rev, qual_q75_rev, 'Reverse reads'),
]):
    ax.fill_between(positions, np.clip(q25, 0, 40), np.clip(q75, 0, 40),
                    alpha=0.3, color='steelblue', label='Q25–Q75')
    ax.plot(positions, np.clip(mean_q, 0, 40), color='steelblue', lw=1.5, label='Median Q')
    ax.axhline(20, color='red',    linestyle='--', lw=1, label='Q20 threshold')
    ax.axhline(30, color='orange', linestyle=':', lw=1, label='Q30 threshold')
    # Suggest truncation at position where median drops below Q25
    trunc_pos = positions[np.where(mean_q < 25)[0][0]] if (mean_q < 25).any() else read_length
    ax.axvline(trunc_pos, color='green', linestyle='-', lw=2, alpha=0.7,
               label=f'Truncate at pos {trunc_pos}')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Phred quality score')
    ax.set_title(f'QIIME2 demux summary: {label}')
    ax.set_xlim(1, read_length)
    ax.set_ylim(0, 42)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.show()
print('Recommended DADA2 truncation: --p-trunc-len-f 220 --p-trunc-len-r 200')
```

## 2. Denoising with DADA2

**DADA2** (Divisive Amplicon Denoising Algorithm 2) corrects sequencing errors and produces **Amplicon Sequence Variants (ASVs)** — exact sequence variants at single-nucleotide resolution — rather than the traditional OTUs (Operational Taxonomic Units) clustered at 97% similarity.

**ASV vs OTU:**
- OTU: sequences clustered at 97% similarity threshold; single representative per cluster; some diversity lost
- ASV: exact error-corrected sequences; no clustering; two ASVs can differ by a single nucleotide; reproducible across studies (same ASV = same sequence in any study)

```bash
# DADA2 denoising (paired-end)
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs reads.qza \
    --p-trim-left-f 13 \          # trim primer/adapter from 5' end of forward reads
    --p-trim-left-r 13 \          # trim primer/adapter from 5' end of reverse reads
    --p-trunc-len-f 220 \         # truncate forward reads at position 220
    --p-trunc-len-r 200 \         # truncate reverse reads at position 200
    --p-n-threads 0 \             # 0 = use all available cores
    --o-table feature_table.qza \ # ASV feature table (samples × ASVs)
    --o-representative-sequences rep_seqs.qza \  # ASV sequences
    --o-denoising-stats stats.qza
```

**Choosing truncation lengths:** Forward and reverse reads must overlap by ≥20 bp after truncation for merging. For V3-V4 (expected amplicon ~460 bp): using 220 + 200 gives 420 bp combined, which overlaps sufficiently on a 460 bp amplicon.

**Quality check:** After denoising, visualize stats to ensure ≥75% of input reads pass filtering and ≥60% are merged successfully.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)

# Simulate DADA2 denoising statistics across samples
n_samples = 16
sample_names = [f'S{i+1:02d}' for i in range(n_samples)]

input_reads    = rng.integers(20_000, 80_000, size=n_samples)
filtered_pct   = rng.uniform(0.82, 0.95, size=n_samples)   # 82-95% pass filtering
denoised_pct   = rng.uniform(0.90, 0.99, size=n_samples)   # 90-99% denoised
merged_pct     = rng.uniform(0.78, 0.92, size=n_samples)   # 78-92% merged
non_chimeric   = rng.uniform(0.93, 0.99, size=n_samples)   # 93-99% non-chimeric

filtered   = (input_reads * filtered_pct).astype(int)
denoised   = (filtered * denoised_pct).astype(int)
merged     = (denoised * merged_pct).astype(int)
final      = (merged * non_chimeric).astype(int)

stats_df = pd.DataFrame({
    'Sample':     sample_names,
    'Input':      input_reads,
    'Filtered':   filtered,
    'Denoised':   denoised,
    'Merged':     merged,
    'Non-chimeric': final,
    'Retention_%': final / input_reads * 100,
})

print('DADA2 denoising statistics (first 8 samples):')
print(stats_df.head(8)[['Sample', 'Input', 'Filtered', 'Merged', 'Non-chimeric', 'Retention_%']].to_string(index=False))
print(f'\nMean overall retention: {(final / input_reads).mean():.1%}')

# Simulate ASV counts
n_asvs_per_sample = rng.integers(200, 800, size=n_samples)
print(f'Mean ASVs per sample: {n_asvs_per_sample.mean():.0f} (range: {n_asvs_per_sample.min()}–{n_asvs_per_sample.max()})')

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Read retention waterfall
steps = ['Input', 'Filtered', 'Merged', 'Non-chimeric']
step_means = [input_reads.mean(), filtered.mean(), merged.mean(), final.mean()]
axes[0].bar(steps, step_means, color=['#4C72B0', '#55A868', '#DD8452', '#C44E52'], alpha=0.85)
axes[0].set_ylabel('Mean reads per sample')
axes[0].set_title('DADA2 read retention pipeline')
for i, (step, val) in enumerate(zip(steps, step_means)):
    axes[0].text(i, val + 500, f'{val:,.0f}', ha='center', fontsize=9)

# Retention % distribution
axes[1].hist(stats_df['Retention_%'], bins=10, color='steelblue', edgecolor='white', alpha=0.85)
axes[1].axvline(75, color='red', linestyle='--', lw=1.5, label='75% threshold')
axes[1].set_xlabel('Overall read retention (%)')
axes[1].set_ylabel('Number of samples')
axes[1].set_title('DADA2 read retention distribution')
axes[1].legend()

plt.tight_layout()
plt.show()
```

## 3. Taxonomic Classification with SILVA

**SILVA** (silva.nr_v138) is the standard reference database for 16S rRNA classification in QIIME2. Pre-trained Naive Bayes classifiers are available from the QIIME2 data resources page — use the classifier trained for your specific primer pair and read length.

```bash
# Download pre-trained SILVA 138 classifier (V3-V4, 515F/806R primers)
# Available at: https://data.qiime2.org/2024.5/common/silva-138-99-seqs-515-806.qza

# Classify ASVs
qiime feature-classifier classify-sklearn \
    --i-classifier silva138_nb_515_806_classifier.qza \
    --i-reads rep_seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs 1

# Generate taxonomy bar chart
qiime taxa barplot \
    --i-table feature_table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa_barplot.qzv

# Filter out mitochondria and chloroplast sequences (contaminants from host plant/food)
qiime taxa filter-table \
    --i-table feature_table.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude Mitochondria,Chloroplast,Eukaryota \
    --o-filtered-table feature_table_filtered.qza
```

**Confidence threshold:** The `--p-confidence 0.7` default means assignments below 70% bootstrap confidence are reported as unclassified at that taxonomic level. For species-level classification, accuracy drops substantially; genus-level is more reliable.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

rng = np.random.default_rng(20)

# Simulate taxonomy classification for 16 gut metagenome samples
taxa = {
    'g__Bacteroides':       (0.18, 0.08),
    'g__Faecalibacterium':  (0.10, 0.06),
    'g__Prevotella':        (0.08, 0.07),
    'g__Bifidobacterium':   (0.07, 0.04),
    'g__Ruminococcus':      (0.06, 0.04),
    'g__Akkermansia':       (0.05, 0.04),
    'g__Blautia':           (0.05, 0.03),
    'g__Lachnospiraceae':   (0.06, 0.04),
    'g__Streptococcus':     (0.04, 0.03),
    'g__Clostridium':       (0.04, 0.03),
    'Other':                (0.27, 0.05),
}

n_samples = 16
conditions = ['Healthy'] * 8 + ['IBD'] * 8
sample_names = [f'{c[:1]}{i+1}' for i, c in enumerate(conditions)]

n_taxa = len(taxa)
compositions = np.zeros((n_samples, n_taxa))
for j, (genus, (mean, std)) in enumerate(taxa.items()):
    base = rng.normal(mean, std, n_samples)
    # Reduce Faecalibacterium in IBD (known association)
    if genus == 'g__Faecalibacterium':
        base[8:] *= 0.4
    compositions[:, j] = np.clip(base, 0.001, 1)

# Re-normalize to sum to 1
compositions = compositions / compositions.sum(axis=1, keepdims=True)
comp_df = pd.DataFrame(compositions, columns=list(taxa.keys()), index=sample_names)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Stacked bar chart (taxonomy bar plot)
colors = cm.tab20(np.linspace(0, 1, n_taxa))
bottom = np.zeros(n_samples)
for j, (genus, col) in enumerate(zip(list(taxa.keys()), colors)):
    vals = comp_df.iloc[:, j].values
    axes[0].bar(range(n_samples), vals, bottom=bottom, label=genus, color=col, width=0.85)
    bottom += vals

axes[0].set_xticks(range(n_samples))
axes[0].set_xticklabels(sample_names, rotation=45, ha='right', fontsize=8)
axes[0].set_ylabel('Relative abundance')
axes[0].set_title('QIIME2 taxonomy bar chart (genus level)')
axes[0].legend(loc='upper right', fontsize=6, ncol=2)
# Add condition separator
axes[0].axvline(7.5, color='black', lw=2, linestyle='--')
axes[0].text(3.5, 1.02, 'Healthy', ha='center', fontsize=9, color='green')
axes[0].text(11.5, 1.02, 'IBD',    ha='center', fontsize=9, color='red')

# Faecalibacterium comparison
fpraus = comp_df['g__Faecalibacterium']
cond_colors = ['#55A868'] * 8 + ['#C44E52'] * 8
axes[1].bar(range(n_samples), fpraus.values, color=cond_colors, alpha=0.8, edgecolor='white')
axes[1].set_xticks(range(n_samples))
axes[1].set_xticklabels(sample_names, rotation=45, ha='right', fontsize=8)
axes[1].axvline(7.5, color='black', lw=1.5, linestyle='--')
axes[1].set_ylabel('Relative abundance')
axes[1].set_title('Faecalibacterium prausnitzii\n(reduced in IBD)')

from matplotlib.patches import Patch
axes[1].legend(handles=[Patch(color='#55A868', label='Healthy'),
                         Patch(color='#C44E52', label='IBD')])

plt.tight_layout()
plt.show()
```
