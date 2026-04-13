---
name: bio-applied-taxonomic-profiling
description: "**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/01_taxonomic_profiling.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Taxonomic Profiling of Shotgun Metagenomes

*Source: Course notebook `Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/01_taxonomic_profiling.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 04 (Microbial Diversity)*

---

**By the end of this notebook you will be able to:**
1. Explain the difference between 16S amplicon and shotgun metagenomics
2. Decontaminate metagenomic reads by removing host sequences
3. Classify reads taxonomically with Kraken2
4. Re-estimate species abundances with Bracken
5. Visualize community composition with Krona and bar charts



**Key resources:**
- [QIIME2 tutorials](https://docs.qiime2.org/2024.10/tutorials/)
- [Galaxy Training — Metagenomics](https://training.galaxyproject.org/training-material/topics/metagenomics/)
- [Kraken2 documentation](https://github.com/DerrickWood/kraken2)

## 1. Shotgun vs 16S Metagenomics

Two major strategies exist for characterizing microbial communities:

| Feature | 16S rRNA amplicon | Shotgun metagenomics |
|---|---|---|
| Target | 16S rRNA gene (V3-V4 region) | All DNA in the sample |
| Taxonomic resolution | Genus level typical; species sometimes | Species and strain level |
| Functional information | No | Yes (gene/pathway content) |
| Host contamination | N/A (bacteria-specific amplification) | Major issue (requires decontamination) |
| Cost per sample | Low (~$30–80) | Higher (~$200–500) |
| Bioinformatics complexity | Moderate (QIIME2, DADA2) | High (multi-step pipeline) |
| Database bias | 16S database coverage | Requires organisms in reference DB |

**When to use shotgun:** Research questions requiring species/strain resolution, functional gene content, AMR gene detection, or viral metagenomics.

### Shotgun Metagenomics Pipeline Overview

```python
Raw FASTQ → QC (fastp) → Host decontamination (Bowtie2/hg38)
  → Taxonomic profiling (Kraken2 + Bracken)
  → Functional profiling (HUMAnN3)
  → Assembly (MEGAHIT)
  → Binning (MetaBAT2)
  → MAG QC (CheckM)
```python

## 2. Host Read Removal

Shotgun metagenomic samples from human gut, skin, or blood contain significant host (human) DNA. If not removed, host reads inflate computational cost and cause false taxonomic assignments (reads misassigned to bacteria with similar sequence content).

**Tool:** Bowtie2 with `--un-conc-gz` to output read pairs that did NOT map to the host.

```bash
# Build or download hg38 Bowtie2 index (one-time setup)
# bowtie2-build hg38.fa hg38_index/hg38

# Decontaminate: align to human genome, keep unmapped reads
bowtie2 -x hg38_index/hg38 \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -p 8 \
    --very-sensitive \
    --un-conc-gz decontam_%.fastq.gz \  # unmapped pairs → decontam_1.fastq.gz, decontam_2.fastq.gz
    > /dev/null  # discard human-mapped SAM output

# Report host content
echo "Host contamination: $(cat contamination_rate.txt)%"
```python

**Typical host contamination rates:**
- Gut metagenome: 1–10% human reads
- Bronchoalveolar lavage: 40–80% human reads
- Skin: 10–30% human reads
- PBMC or blood: 60–95% human reads

**Alternative for known pathogens:** Use the **Decontam** R package to remove contaminants based on inverse correlation with DNA concentration across negative controls.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate host contamination profiles across sample types
sample_types = {
    'Gut (stool)':     rng.beta(2, 20, size=20),        # low contamination
    'Skin swab':       rng.beta(3, 10, size=20),        # moderate
    'Oral rinse':      rng.beta(3, 15, size=20),        # moderate-low
    'BAL (lung)':      rng.beta(6, 4,  size=20),        # high
    'PBMC (blood)':    rng.beta(10, 2, size=20),        # very high
}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Box plot
data = list(sample_types.values())
labels = list(sample_types.keys())
bp = axes[0].boxplot(data, labels=labels, patch_artist=True, notch=False)
colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
axes[0].set_ylabel('Host read fraction')
axes[0].set_title('Host contamination by sample type')
axes[0].tick_params(axis='x', rotation=25)
axes[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0%}'))

# Before/after decontamination: simulated total reads
n_samples = 10
total_before = rng.integers(15_000_000, 40_000_000, size=n_samples)
host_fracs   = rng.beta(2, 15, size=n_samples)
total_after  = (total_before * (1 - host_fracs)).astype(int)

x = np.arange(n_samples)
axes[1].bar(x - 0.2, total_before / 1e6, width=0.4, label='Before', color='#C44E52', alpha=0.8)
axes[1].bar(x + 0.2, total_after  / 1e6, width=0.4, label='After',  color='#4C72B0', alpha=0.8)
axes[1].set_xlabel('Sample')
axes[1].set_ylabel('Read pairs (millions)')
axes[1].set_title('Read count before/after host removal')
axes[1].legend()
axes[1].set_xticks(x)
axes[1].set_xticklabels([f'S{i+1}' for i in range(n_samples)])

plt.tight_layout()
plt.show()

print(f'Mean host fraction: {host_fracs.mean():.1%}')
print(f'Mean reads retained: {(1 - host_fracs).mean():.1%}')
```python

## 3. Taxonomic Classification with Kraken2

**Kraken2** classifies metagenomic reads by assigning each read to the lowest common ancestor (LCA) of all genomes in its reference database that share k-mer sequences with the read. Key parameters:

- **Database**: standard (bacteria + archaea + viruses + human), PlusPF (adds protozoa + fungi), or custom. Larger databases improve sensitivity but require more RAM (standard: ~55 GB; PlusPF: ~100 GB).
- `--confidence 0.1`: minimum fraction of k-mers classified to a node; helps reduce false assignments (default 0 = report everything)
- `--minimum-hit-groups 3`: minimum number of minimizer groups in the database before a read is assigned (reduces spurious assignments)

```bash
# Classify decontaminated reads
kraken2 \
    --db kraken2_db/ \
    --paired \
    --gzip-compressed \
    --threads 16 \
    --confidence 0.1 \
    --minimum-hit-groups 3 \
    --report kraken2_report.txt \
    --output kraken2_output.txt \
    decontam_1.fastq.gz decontam_2.fastq.gz
```python

**Kraken2 report columns** (`kraken2_report.txt`):
1. Percentage of reads classified to this clade
2. Reads classified at this exact taxon
3. Reads classified within this clade (includes sub-taxa)
4. Rank code (S=species, G=genus, F=family, O=order, C=class, P=phylum, K=kingdom, R=root)
5. NCBI taxonomy ID
6. Indented scientific name

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(1)

# Simulate a Kraken2 report for a gut metagenome sample
species = [
    'Bacteroides vulgatus',      'Faecalibacterium prausnitzii',
    'Prevotella copri',          'Bifidobacterium longum',
    'Ruminococcus gnavus',       'Akkermansia muciniphila',
    'Bacteroides uniformis',     'Blautia obeum',
    'Lachnospiraceae bacterium', 'Eubacterium rectale',
    'Clostridium bolteae',       'Streptococcus salivarius',
    'Collinsella aerofaciens',   'Dorea longicatena',
    'Roseburia intestinalis',
]
n_species = len(species)

# Simulate classification counts (log-normal distributed)
read_counts = rng.lognormal(mean=10.5, sigma=1.8, size=n_species).astype(int)
read_counts = np.sort(read_counts)[::-1]
total_classified = read_counts.sum()
total_reads = int(total_classified / 0.75)  # ~75% classification rate

fractions = read_counts / total_reads * 100

kraken_df = pd.DataFrame({
    'Species':   species,
    'Reads':     read_counts,
    'Percentage': fractions,
}).sort_values('Reads', ascending=False)

print('=== Kraken2 Classification Summary ===')
print(f'Total reads processed:  {total_reads:>12,}')
print(f'Total classified:       {total_classified:>12,}  ({total_classified/total_reads:.1%})')
print(f'\nTop 10 species:')
print(kraken_df[['Species', 'Reads', 'Percentage']].head(10).to_string(index=False))

fig, ax = plt.subplots(figsize=(9, 6))
top10 = kraken_df.head(10)
colors = plt.cm.tab10(np.linspace(0, 1, 10))
bars = ax.barh(range(len(top10)), top10['Percentage'][::-1], color=colors)
ax.set_yticks(range(len(top10)))
ax.set_yticklabels(top10['Species'][::-1], fontsize=9)
ax.set_xlabel('Percentage of total reads (%)')
ax.set_title('Kraken2 top species — gut metagenome')
ax.axvline(1.0, color='red', linestyle=':', lw=1, label='1% threshold')
ax.legend()
plt.tight_layout()
plt.show()
```python

## 4. Abundance Re-estimation with Bracken

**Bracken** (Bayesian Re-estimation of Abundance with KraKEN) addresses a fundamental limitation of Kraken2: reads that share k-mers with multiple species are assigned to their lowest common ancestor (LCA), inflating higher-taxon abundances. Bracken uses the k-mer distribution across all genomes in the database to **redistribute** these reads back to the most likely species level.

**Key parameters:**
- `-r 150`: read length (should match your actual sequencing read length)
- `-l S`: taxonomic level to estimate (`S`=Species, `G`=Genus, `F`=Family)
- `-t 10`: minimum reads threshold at the species level; species with fewer reads are excluded

```bash
# Species-level abundance re-estimation
bracken \
    -d kraken2_db/ \
    -i kraken2_report.txt \
    -o bracken_species.txt \
    -w bracken_species_report.txt \
    -r 150 \
    -l S \
    -t 10
```python

**Bracken output columns** (`bracken_species.txt`):
- `name`: species name
- `taxonomy_id`: NCBI taxonomy ID
- `taxonomy_lvl`: level (S)
- `kraken_assigned_reads`: reads Kraken2 assigned directly to this species
- `added_reads`: reads redistributed down from higher-level LCA nodes
- `new_est_reads`: total = kraken_assigned_reads + added_reads
- `fraction_total_reads`: new_est_reads / total_reads (the abundance estimate)

The fraction_total_reads column is what you use for downstream diversity analysis and differential abundance testing.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(2)

# Simulate Bracken output showing redistribution from LCA nodes
species = [
    'Bacteroides vulgatus',      'Faecalibacterium prausnitzii',
    'Prevotella copri',          'Bifidobacterium longum',
    'Ruminococcus gnavus',       'Akkermansia muciniphila',
    'Bacteroides uniformis',     'Blautia obeum',
    'Lachnospiraceae bacterium', 'Eubacterium rectale',
    'Clostridium bolteae',       'Streptococcus salivarius',
    'Collinsella aerofaciens',   'Dorea longicatena',
    'Roseburia intestinalis',
]

kraken_direct  = rng.lognormal(10.5, 1.8, size=len(species)).astype(int)
# Bracken adds redistributed reads (from LCA nodes)
added_reads    = (kraken_direct * rng.lognormal(-0.3, 0.8, len(species))).astype(int)
added_reads    = np.abs(added_reads)
total_est      = kraken_direct + added_reads
total_all_reads = total_est.sum() * 4  # approximate total reads

bracken_df = pd.DataFrame({
    'Species':               species,
    'Kraken2_direct':        kraken_direct,
    'Bracken_added':         added_reads,
    'Total_est_reads':       total_est,
    'Fraction_total':        total_est / total_all_reads,
}).sort_values('Total_est_reads', ascending=False)

print('=== Bracken species-level abundance ===')
print(bracken_df[['Species', 'Kraken2_direct', 'Bracken_added', 'Fraction_total']].head(10).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Stacked bar: Kraken2 direct vs Bracken added
top10 = bracken_df.head(10).reset_index(drop=True)
x = np.arange(len(top10))
axes[0].barh(x, top10['Kraken2_direct'],        label='Kraken2 direct', color='#4C72B0', alpha=0.85)
axes[0].barh(x, top10['Bracken_added'],
             left=top10['Kraken2_direct'],       label='Bracken added',  color='#DD8452', alpha=0.85)
axes[0].set_yticks(x)
axes[0].set_yticklabels(top10['Species'], fontsize=8)
axes[0].set_xlabel('Read count')
axes[0].set_title('Kraken2 vs Bracken redistribution')
axes[0].legend()

# Relative abundance pie of top 7 + other
top7 = bracken_df.head(7)
other_frac = bracken_df.iloc[7:]['Fraction_total'].sum()
pie_labels = list(top7['Species']) + ['Other']
pie_vals   = list(top7['Fraction_total']) + [other_frac]
colors_pie = plt.cm.tab10(np.linspace(0, 1, 8))
axes[1].pie(pie_vals, labels=[s.split()[-1] for s in pie_labels],
            colors=colors_pie, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 8})
axes[1].set_title('Species-level relative abundance\n(Bracken estimates)')

plt.tight_layout()
plt.show()
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
