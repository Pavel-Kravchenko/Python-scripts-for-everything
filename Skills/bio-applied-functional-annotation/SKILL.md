---
name: bio-applied-functional-annotation
description: Functional Annotation of Metagenomes with NumPy
tool_type: python
primary_tool: NumPy
---

# Functional Annotation of Metagenomes

References: [HUMAnN3](https://github.com/biobakery/humann) · [bioBakery](https://github.com/biobakery/biobakery/wiki) · [MetaCyc](https://metacyc.org/)

## HUMAnN3 Overview

Three-stage workflow:
1. **MetaPhlAn4 taxonomic profiling** — identifies species and abundances; selects pangenome DBs for step 2.
2. **Pangenome alignment (Bowtie2)** — aligns reads against nucleotide pangenomes of detected species.
3. **Unaligned read translation search (DIAMOND)** — searches UniRef protein DB; catches organisms not in step 2 pangenomes.

**Output files per sample:**
- `_genefamilies.tsv` — gene family abundances (RPK), stratified by species
- `_pathabundance.tsv` — MetaCyc pathway abundances (RPK), stratified by species
- `_pathcoverage.tsv` — fraction of pathway reactions detected (capped at 1)

```bash
cat decontam_1.fastq.gz decontam_2.fastq.gz > merged.fastq.gz

humann \
    --input merged.fastq.gz \
    --output humann3_out/ \
    --threads 16 \
    --metaphlan-options '--bowtie2db metaphlan4_db' \
    --nucleotide-database chocophlan_db/ \
    --protein-database uniref90_diamond/

# Normalize to CPM for cross-sample comparison
humann_renorm_table \
    --input humann3_out/sample_pathabundance.tsv \
    --output humann3_out/sample_pathabundance_cpm.tsv \
    --units cpm
```

## Output Interpretation

**RPK** (reads per kilobase) accounts for gene length bias. After `humann_renorm_table --units cpm`, values are copies per million — proportional to fraction of reads mapping to that gene family.

**Stratified vs unstratified:** Each pathway/gene family reported as total (unstratified) and per contributing species (stratified):
```
PYRUVATE-FERMENTATION-PWY|unclassified     152.3
PYRUVATE-FERMENTATION-PWY|Bacteroides      89.1
PYRUVATE-FERMENTATION-PWY|Faecalibacterium 43.7
```

**Coverage vs abundance:** High coverage (reactions present) can coexist with low abundance (few reads). Coverage > 1 not possible.

```bash
# Join all samples
humann_join_tables --input humann3_out/ \
    --output all_samples_pathabundance.tsv --file_name pathabundance_cpm

# Split stratified/unstratified
humann_split_stratified_table \
    --input all_samples_pathabundance.tsv --output humann3_stratified/
```

```python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns

rng = np.random.default_rng(42)
pathways = [
    'GLYCOLYSIS-E-D: Glycolysis (EMP)', 'PYRUVATE-FERMENTATION-PWY: Pyruvate fermentation',
    'FASYN-ELONG-PWY: Fatty acid elongation', 'BUTYRATE-BIOSYNTHESIS-I: Butyrate biosynthesis I',
    'PROPIONATE-BIOSYNTHESIS-I: Propionate biosynthesis I', 'FOLSYN-PWY: Folate biosynthesis',
]
n_samples = 12
conditions = ['Healthy'] * 4 + ['IBD'] * 4 + ['CRC'] * 4
base = rng.lognormal(mean=6, sigma=2, size=(len(pathways), n_samples))

ibd_idx = [i for i, c in enumerate(conditions) if c != 'Healthy']
but_idx = next(i for i, p in enumerate(pathways) if 'Butyrate' in p)
base[but_idx, ibd_idx] *= 0.3   # butyrate reduced in IBD/CRC
base[0, ibd_idx] *= 1.5          # glycolysis up

path_df = pd.DataFrame(base, index=pathways,
                        columns=[f'{c[:3]}_{i%4+1}' for i, c in enumerate(conditions)])
log2_df = np.log2(path_df + 1)
short_labels = [p.split(':')[1].strip()[:35] if ':' in p else p[:35] for p in pathways]

fig, ax = plt.subplots(figsize=(10, 5))
sns.heatmap(log2_df, ax=ax, cmap='YlOrRd', yticklabels=short_labels,
            cbar_kws={'label': 'log2(CPM+1)'}, linewidths=0.3)
ax.set_title('HUMAnN3 pathway abundances'); plt.tight_layout(); plt.show()
```

## Differential Pathway Analysis with MaAsLin2

MaAsLin2 fits a GLM with optional random effects to each feature vs. metadata variables. Handles continuous/categorical metadata, repeated measures, compositionality corrections.

```R
library(Maaslin2)
pathway_table <- read.table('humann3_stratified/all_pathabundance_unstratified.tsv',
                             sep='\t', header=TRUE, row.names=1, comment.char='#')
metadata <- read.table('metadata.tsv', sep='\t', header=TRUE, row.names=1)

fit <- Maaslin2(
  input_data=t(pathway_table), input_metadata=metadata,
  output='maaslin2_out/',
  fixed_effects=c('condition'), random_effects=c('batch'),
  transform='LOG', normalization='TSS',
  min_prevalence=0.1, min_abundance=0.0001
)
# Output: significant_results.tsv — feature, metadata, coef, pval, qval
```

```python
from statsmodels.stats.multitest import multipletests
from scipy import stats

rng = np.random.default_rng(55)
pathways_all = ['Butyrate biosynthesis I', 'Propionate biosynthesis I', 'Glycolysis (EMP)',
                'Folate biosynthesis', 'Aerobactin biosynthesis', 'Fatty acid elongation',
                'Short-chain FA metabolism', 'Bile acid transformation']
n_pw = len(pathways_all)
log2fc = np.concatenate([rng.normal(-1.8, 0.5, 2), rng.normal(1.5, 0.6, 2),
                          rng.normal(0, 0.4, n_pw-4)])
se = np.abs(rng.normal(0.3, 0.1, n_pw)) + 0.1
pvals = stats.t.sf(np.abs(log2fc / se), df=6) * 2
_, qvals, _, _ = multipletests(pvals, method='BH')

diff_df = pd.DataFrame({'Pathway': pathways_all, 'log2FC': log2fc, 'qval': qvals})
print(diff_df[diff_df['qval'] < 0.05][['Pathway','log2FC','qval']].to_string(index=False))
```

## AMR Gene Detection

Gene family abundances can detect AMR genes (CARD/ResFinder/ARG-ANNOT), virulence factors (VFDB), or CAZymes.

```bash
amrfinder --nucleotide megahit_assembly/final.contigs.fa \
    --output amr_report.txt --threads 8 --plus
```

**CARD resistance categories:** Intrinsic (naturally present, e.g. AmpC), Acquired (HGT-transferred, e.g. blaCTX-M, mcr-1), Mutational (point mutations, e.g. gyrA quinolone resistance). AMR expressed as RPKM or presence/absence per sample.

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
