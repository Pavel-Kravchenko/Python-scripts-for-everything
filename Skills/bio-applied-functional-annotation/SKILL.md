---
name: bio-applied-functional-annotation
description: "**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/02_functional_annotation.ipynb"
---

# Functional Annotation of Metagenomes

*Source: Course notebook `Tier_3_Applied_Bioinformatics/26_Metagenomics_Shotgun/02_functional_annotation.ipynb`*

# Functional Annotation of Metagenomes

**Tier 3 — Applied Bioinformatics | Module 26 · Notebook 2**

*Prerequisites: Notebook 1 (Taxonomic Profiling)*

---

**By the end of this notebook you will be able to:**
1. Run HUMAnN3 to profile functional pathway abundance from metagenomic reads
2. Interpret UniRef90 gene families and MetaCyc pathway tables
3. Normalize pathway abundances (copies per million, relative)
4. Perform differential pathway analysis between sample groups
5. Visualize community functional profiles with heatmaps



**Key resources:**
- [HUMAnN3 documentation](https://github.com/biobakery/humann)
- [bioBakery tutorials](https://github.com/biobakery/biobakery/wiki)
- [MetaCyc pathway database](https://metacyc.org/)

## 1. HUMAnN3 Overview

**HUMAnN3** (HMP Unified Metabolic Analysis Network 3) profiles functional gene families and metabolic pathways in metagenomic samples. Its three-stage workflow:

1. **MetaPhlAn4 taxonomic profiling**: identifies which species are present and at what abundance. Used to select species-specific pangenome databases for alignment in step 2.
2. **Pangenome alignment (Bowtie2)**: aligns reads against the nucleotide pangenomes of detected species. Fast and specific — leverages the species detected in step 1.
3. **Unaligned read translation search (DIAMOND)**: translates unaligned reads and searches them against the UniRef protein database. Catches functional genes from organisms not in step 2's pangenomes.

**Three output files produced per sample:**
- `_genefamilies.tsv`: gene family abundances in RPK (reads per kilobase), stratified by contributing species
- `_pathabundance.tsv`: MetaCyc pathway abundances in RPK, stratified by species
- `_pathcoverage.tsv`: fraction of each pathway's reactions with at least one detected gene

### HUMAnN3 Pipeline

```bash
# Merge paired reads (HUMAnN3 takes single FASTQ)
cat decontam_1.fastq.gz decontam_2.fastq.gz > merged.fastq.gz

# Run HUMAnN3 (requires ~15-50 GB RAM; runtime ~1-4 hours per sample)
humann \
    --input merged.fastq.gz \
    --output humann3_out/ \
    --threads 16 \
    --metaphlan-options '--bowtie2db metaphlan4_db' \
    --nucleotide-database chocophlan_db/ \
    --protein-database uniref90_diamond/

# Normalize by copies per million (CPM) — needed for cross-sample comparison
humann_renorm_table \
    --input humann3_out/sample_pathabundance.tsv \
    --output humann3_out/sample_pathabundance_cpm.tsv \
    --units cpm
```

## 2. Output Interpretation and Normalization

**Gene family units (RPK):** reads per kilobase of gene length. This accounts for the fact that longer genes attract more reads by chance. After normalization with `humann_renorm_table --units cpm`, values are copies per million (CPM) — proportional to the fraction of reads mapping to that gene family.

**Stratified vs unstratified output:** HUMAnN3 reports each pathway/gene family both as a total (unstratified) and broken down by contributing species (stratified). Example:
```
PYRUVATE-FERMENTATION-PWY: PWY|unclassified     152.3
PYRUVATE-FERMENTATION-PWY: Bacteroides vulgatus  89.1
PYRUVATE-FERMENTATION-PWY: Faecalibacterium      43.7
```

**Pathway coverage vs abundance:** A pathway can have high coverage (most reactions present) but low abundance (few reads). Coverage > 1 is not possible — it is capped at 1 (all reactions represented).

### Joining across samples

```bash
# After running HUMAnN3 on all samples, join into one table
humann_join_tables \
    --input humann3_out/ \
    --output all_samples_pathabundance.tsv \
    --file_name pathabundance_cpm

# Remove stratification for group-level analysis
humann_split_stratified_table \
    --input all_samples_pathabundance.tsv \
    --output humann3_stratified/
# Creates: ..._unstratified.tsv (pathways × samples) and ..._stratified.tsv
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

rng = np.random.default_rng(42)

# Simulate HUMAnN3 MetaCyc pathway abundance table (CPM, unstratified)
pathways = [
    'GLYCOLYSIS-E-D: Glycolysis (Embden-Meyerhof-Parnas)',
    'PYRUVATE-FERMENTATION-PWY: Pyruvate fermentation to acetate and lactate',
    'PWY-5667: CDP-diacylglycerol biosynthesis',
    'FASYN-ELONG-PWY: Fatty acid elongation — saturated',
    'AEROBACTINSYN-PWY: Aerobactin biosynthesis',
    'P163-PWY: L-tyrosine degradation I',
    'PWY-7111: Pyruvate fermentation to isobutanol',
    'LACTOSECAT-PWY: Lactose and galactose degradation I',
    'BUTYRATE-BIOSYNTHESIS-I: Butyrate biosynthesis I',
    'PROPIONATE-BIOSYNTHESIS-I: Propionate biosynthesis I',
    'FOLSYN-PWY: Folate biosynthesis',
    'KETOGLUCONATEUTIL-PWY: 2-Ketogluconate utilization',
]

n_samples = 12
conditions = ['Healthy'] * 4 + ['IBD'] * 4 + ['CRC'] * 4
sample_names = [f'{c[:3]}_{i%4+1}' for i, c in enumerate(conditions)]

# Base abundance (log-normal)
base = rng.lognormal(mean=6, sigma=2, size=(len(pathways), n_samples))

# IBD/CRC: reduced butyrate biosynthesis, altered fermentation
ibd_idx = [i for i, c in enumerate(conditions) if c != 'Healthy']
but_idx = next(i for i, p in enumerate(pathways) if 'Butyrate' in p)
base[but_idx, ibd_idx] *= 0.3  # butyrate biosynthesis reduced in IBD/CRC
glyc_idx = next(i for i, p in enumerate(pathways) if 'Glycolysis' in p)
base[glyc_idx, ibd_idx] *= 1.5  # glycolysis up

path_df = pd.DataFrame(base, index=pathways, columns=sample_names)

print('Pathway abundance table (CPM, first 5 pathways x 4 samples):')
print(path_df.iloc[:5, :4].round(1).to_string())

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Heatmap (log2 CPM)
log2_df = np.log2(path_df + 1)
short_labels = [p.split(':')[1].strip()[:35] if ':' in p else p[:35] for p in pathways]
sns.heatmap(log2_df, ax=axes[0], cmap='YlOrRd',
            xticklabels=sample_names, yticklabels=short_labels,
            cbar_kws={'label': 'log2(CPM + 1)'}, linewidths=0.3)
axes[0].set_title('Pathway abundances (HUMAnN3 output)')
axes[0].tick_params(axis='y', labelsize=7)
axes[0].tick_params(axis='x', labelsize=8, rotation=45)

# Butyrate biosynthesis specifically
but_vals = path_df.loc[pathways[but_idx]]
colors = {'Healthy': '#55A868', 'IBD': '#C44E52', 'CRC': '#DD8452'}
for cond, col in colors.items():
    mask = [c == cond for c in conditions]
    axes[1].bar(
        [j + {'Healthy': -0.3, 'IBD': 0, 'CRC': 0.3}[cond] for j in np.where(mask)[0]],
        but_vals[mask].values, width=0.28, color=col, label=cond, alpha=0.85
    )
axes[1].set_ylabel('CPM')
axes[1].set_title('Butyrate biosynthesis I — reduced in IBD/CRC')
axes[1].set_xticks(range(4))
axes[1].set_xticklabels([f'Rep {i+1}' for i in range(4)])
axes[1].legend()

plt.tight_layout()
plt.show()
```

## 3. Differential Pathway Analysis with MaAsLin2

**MaAsLin2** (Multivariate Association with Linear Models) is designed for microbiome compositional data. It applies log transformation and arc-sine-square-root transformation options, and fits a **generalized linear model** with random effects to each feature (pathway, gene family, or taxon) against metadata variables. It handles:
- Continuous and categorical metadata
- Random effects for repeated measurements
- Compositionality corrections

```R
library(Maaslin2)

# Load unstratified pathway abundance table
pathway_table <- read.table('humann3_stratified/all_pathabundance_unstratified.tsv',
                             sep='\t', header=TRUE, row.names=1, comment.char='#')

metadata <- read.table('metadata.tsv', sep='\t', header=TRUE, row.names=1)
# metadata columns: condition, age, BMI, batch

fit <- Maaslin2(
  input_data     = t(pathway_table),   # MaAsLin2 wants samples as rows
  input_metadata = metadata,
  output         = 'maaslin2_out/',
  fixed_effects  = c('condition'),
  random_effects = c('batch'),
  transform      = 'LOG',
  normalization  = 'TSS',              # total sum scaling (compositional correction)
  min_prevalence = 0.1,                # exclude pathways absent in > 90% of samples
  min_abundance  = 0.0001
)
```

**Output:** `significant_results.tsv` with columns:
- `feature`: pathway name
- `metadata`: metadata variable tested
- `coef`: log-fold effect size
- `pval` / `qval`: p-value and FDR-adjusted q-value

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy import stats

rng = np.random.default_rng(55)

# Simulate MaAsLin2-style differential pathway results (Healthy vs IBD)
pathways_all = [
    'Butyrate biosynthesis I',       'Propionate biosynthesis I',
    'Pyruvate fermentation',         'Glycolysis (EMP)',
    'Folate biosynthesis',           'Lactose degradation I',
    'Aerobactin biosynthesis',       'Fatty acid elongation',
    'CDP-diacylglycerol biosynthesis','L-tyrosine degradation I',
    'Isobutanol fermentation',       'Acetate biosynthesis',
    'Short-chain FA metabolism',     'Mucin degradation',
    'Bile acid transformation',
]

n_pw = len(pathways_all)
log2fc = np.concatenate([
    rng.normal(-1.8, 0.5, size=3),  # reduced in IBD (butyrate, propionate, pyruvate)
    rng.normal(1.5,  0.6, size=3),  # elevated in IBD (glycolysis, folate, lactose)
    rng.normal(0,    0.4, size=n_pw - 6),  # unchanged
])
rng.shuffle(log2fc)

se = np.abs(rng.normal(0.3, 0.1, n_pw)) + 0.1
t_stats = log2fc / se
pvals = stats.t.sf(np.abs(t_stats), df=6) * 2  # two-sided t-test, df=6 (4+4-2)
_, qvals, _, _ = multipletests(pvals, method='BH')

diff_df = pd.DataFrame({
    'Pathway':  pathways_all,
    'log2FC':   log2fc,
    'pval':     pvals,
    'qval':     qvals,
}).sort_values('qval')

sig = diff_df[diff_df['qval'] < 0.05]
print(f'Significant pathways (FDR < 5%): {len(sig)}')
print(sig[['Pathway', 'log2FC', 'pval', 'qval']].to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Volcano plot
colors_v = np.where(diff_df['qval'] < 0.05, '#C44E52', '#8172B2')
axes[0].scatter(diff_df['log2FC'], -np.log10(diff_df['pval']),
                c=colors_v, s=60, alpha=0.75)
axes[0].axhline(-np.log10(0.05), color='orange', linestyle=':', lw=1.5, label='p=0.05')
axes[0].axvline(0, color='grey', linestyle='--', lw=1)
for _, row in sig.iterrows():
    short = row['Pathway'][:20]
    axes[0].annotate(short, (row['log2FC'], -np.log10(row['pval'])),
                     textcoords='offset points', xytext=(5, 3), fontsize=7)
axes[0].set_xlabel('log2 fold change (IBD vs Healthy)')
axes[0].set_ylabel('-log10(p-value)')
axes[0].set_title('MaAsLin2: Healthy vs IBD pathway enrichment')
axes[0].legend()

# Bar plot of log2FC for significant pathways
sig_sorted = sig.sort_values('log2FC')
bar_colors = ['#C44E52' if fc > 0 else '#4C72B0' for fc in sig_sorted['log2FC']]
axes[1].barh(range(len(sig_sorted)), sig_sorted['log2FC'], color=bar_colors, alpha=0.85)
axes[1].set_yticks(range(len(sig_sorted)))
axes[1].set_yticklabels([p[:35] for p in sig_sorted['Pathway']], fontsize=8)
axes[1].axvline(0, color='black', lw=1)
axes[1].set_xlabel('log2 fold change')
axes[1].set_title('Significant differential pathways')

plt.tight_layout()
plt.show()
```

## 4. Gene Family Analysis and AMR Detection

Beyond pathways, HUMAnN3 produces **gene family** abundances (UniRef90 clusters). These can be used to:
- Detect **antimicrobial resistance (AMR) genes** by mapping to CARD, ResFinder, or ARG-ANNOT databases
- Identify **virulence factors** using VFDB
- Profile **carbohydrate-active enzymes (CAZymes)** relevant to microbiome metabolic capacity

### AMR Gene Detection with AMRFinder

```bash
# For assembled contigs (from MEGAHIT or SPAdes)
amrfinder \
    --nucleotide megahit_assembly/final.contigs.fa \
    --output amr_report.txt \
    --threads 8 \
    --plus                    # include virulence, stress resistance genes
```

**Resistance gene categories in CARD (Comprehensive Antibiotic Resistance Database):**
- **Intrinsic**: naturally present in organisms (e.g., AmpC in E. coli)
- **Acquired**: transferred via HGT (e.g., blaCTX-M, mcr-1)
- **Mutational**: point mutations conferring resistance (e.g., gyrA quinolone resistance)

AMR prevalence from shotgun data is expressed as **RPKM** (reads per kilobase per million mapped reads) or as presence/absence per sample.
