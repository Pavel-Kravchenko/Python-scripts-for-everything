---
name: bio-applied-mageck-gene-essentiality
description: "**Tier 3 — Applied Bioinformatics | Module 33 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/33_CRISPR_Screen_Analysis/01_mageck_gene_essentiality.ipynb"
primary_tool: MAGeCK
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# CRISPR Screen Analysis with MAGeCK

*Source: Course notebook `Tier_3_Applied_Bioinformatics/33_CRISPR_Screen_Analysis/01_mageck_gene_essentiality.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 33 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 14 (Genetic Engineering In Silico)*

---

**By the end of this notebook you will be able to:**
1. Explain pooled CRISPR screen design: positive/negative selection, sgRNA library
2. Process screen data: FASTQ → sgRNA count table with MAGeCK count
3. Run MAGeCK test for gene essentiality scoring (RRA algorithm)
4. Interpret log fold-change, p-value, and FDR for hit calling
5. Visualize screen results: volcano plot, rank plot, sgRNA-level scatter



**Key resources:**
- [MAGeCK documentation](https://sourceforge.net/p/mageck/wiki/Home/)
- [CRISPRScreenAnalyzeR](https://www.crisprscan.org/)
- [DepMap portal](https://depmap.org/portal/)

## 1. Pooled CRISPR Screen Design

### What is a pooled CRISPR screen?

A pooled CRISPR screen introduces thousands of sgRNAs simultaneously into a cell population — each cell receives (ideally) one sgRNA. The library is delivered by lentiviral transduction at low MOI (multiplicity of infection ≈ 0.3) to ensure single integration events. After selection pressure (drug treatment, growth, immune challenge), cells with different sgRNAs are enriched or depleted. Deep sequencing counts how many cells carry each sgRNA, revealing which genes, when knocked out, confer the phenotype of interest.

### Screen types

| Screen type | Selection pressure | Read-out | Example |
|---|---|---|---|
| Negative selection (viability) | None (growth over time) | sgRNA depletion | Essential genes, tumor suppressors |
| Positive selection | Drug / toxin | sgRNA enrichment | Drug resistance genes, immune evasion |
| Phenotypic | Cell sorting (FACS) | Enrichment/depletion in gate | Surface marker regulation |
| Survival | Immune kill assay | Enrichment | Immunotherapy target discovery |

### Key design parameters

**Library design:**
- **Genome-wide libraries**: Brunello (4 sgRNAs/gene, ~77K guides), GeCKO v2 (6 sgRNAs/gene)
- **Targeted libraries**: kinome, druggable genome, epigenome (~2–5K genes)
- **Non-targeting controls (NTCs)**: 500–1000 guides that don't target any gene; used for normalization and false-positive estimation

**Coverage:**
- Minimum **300× cells per sgRNA** at every timepoint (e.g., 77K guides × 300 = 23M cells minimum)
- Low MOI (0.3) means ~26% of cells receive a sgRNA (Poisson: P(k=1|λ=0.3) ≈ 0.22)
- This requires transducing 23M/0.22 ≈ 100M cells

**Replicates:**
- Minimum 2 biological replicates; 3 strongly recommended
- Correlation of sgRNA counts between replicates is a key QC metric (r > 0.95 expected)

### Why multiple sgRNAs per gene?

A single sgRNA knockout phenotype could be due to off-target cutting, incomplete editing, or chance. With 4–6 sgRNAs targeting the same gene, if 3 of 4 are depleted, this is strong evidence for an on-target effect. MAGeCK's RRA algorithm explicitly uses this redundancy.

## 2. sgRNA Counting: FASTQ → Count Table

### The MAGeCK count step

After sequencing, each read should contain exactly one sgRNA sequence. MAGeCK count:

1. **Trims** reads to remove adapter sequences (sgRNAs are typically 19–20 nt)
2. **Exact-matches** trimmed reads against the sgRNA library file (no mismatches by default, or 1 allowed)
3. **Counts** how many reads map to each sgRNA across each sample
4. Reports unmapped fraction as a quality metric (should be < 20%)

```bash
mageck count \
    --list-seq library.csv \              # sgRNA library: id, sequence, gene
    --name experiment_name \
    --sample-label "plasmid,d14_rep1,d14_rep2,d14_rep3" \
    --fastq plasmid.fastq.gz d14_rep1.fastq.gz d14_rep2.fastq.gz d14_rep3.fastq.gz \
    --sgrna-len 20 \                      # sgRNA length to extract from read
    --norm-method median                  # normalization: median, total, or none
```python

**Library file format** (CSV or space-delimited):
```python
sgRNA_ID,Sequence,Gene
sgRNA001,ATCGATCGATCGATCGATCG,GENE_A
sgRNA002,GCTAGCTAGCTAGCTAGCTA,GENE_A
...
NTC001,AAACTAGTCGATCGATCGAT,NonTargetingControl
```python

### Output: count table

The count table (`.count.txt`) has rows = sgRNAs, columns = samples:

| sgRNA | Gene | plasmid | d14_rep1 | d14_rep2 |
|---|---|---|---|---|
| sgRNA001 | GENE_A | 523 | 89 | 102 |
| sgRNA002 | GENE_A | 487 | 78 | 91 |
| NTC001 | NonTargeting | 510 | 498 | 501 |

Notice how essential gene sgRNAs (GENE_A) are depleted by day 14, while NTCs maintain counts.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import nbinom, rankdata

np.random.seed(42)

# -----------------------------------------------------------------------
# Simulate a CRISPR screen count table
# Parameters:
#   1000 genes × 4 sgRNAs each = 4000 sgRNAs
#   + 200 non-targeting controls (NTCs)
#   3 replicates (day 14 post-selection)
#   ~150 essential genes (strong negative selection)
#   ~50 resistance genes (positive selection with drug)
# -----------------------------------------------------------------------
N_GENES      = 1000
SGRNAS_GENE  = 4
N_NTC        = 200
N_ESSENTIAL  = 150   # core essential genes: depleted
N_RESISTANCE = 50    # drug-resistance genes: enriched

total_sgrnas = N_GENES * SGRNAS_GENE + N_NTC

# Build gene list
gene_ids = [f'GENE_{i:04d}' for i in range(N_GENES)]
essential_genes  = set(gene_ids[:N_ESSENTIAL])
resistance_genes = set(gene_ids[N_ESSENTIAL:N_ESSENTIAL + N_RESISTANCE])

# Build sgRNA-to-gene mapping
sgrna_ids   = []
gene_labels = []
for gene in gene_ids:
    for j in range(SGRNAS_GENE):
        sgrna_ids.append(f'{gene}_sg{j+1}')
        gene_labels.append(gene)
for k in range(N_NTC):
    sgrna_ids.append(f'NTC_{k+1:03d}')
    gene_labels.append('NonTargeting')

n_sgrna_total = len(sgrna_ids)
print(f"Library: {n_sgrna_total} sgRNAs  ({N_GENES} genes × {SGRNAS_GENE} + {N_NTC} NTCs)")

# Plasmid counts: negative binomial (overdispersed, mean 500, dispersion ~5)
plasmid_mu = 500
plasmid_size = 5   # dispersion parameter
plasmid_counts = nbinom.rvs(n=plasmid_size, p=plasmid_size/(plasmid_size+plasmid_mu),
                             size=n_sgrna_total)
plasmid_counts = np.maximum(plasmid_counts, 1)

# Day 14 counts: apply LFC based on gene category
lfc_true = np.zeros(n_sgrna_total)   # log2 fold-change

for i, (sgrna, gene) in enumerate(zip(sgrna_ids, gene_labels)):
    if gene in essential_genes:
        # Essential: strong depletion, guide-level variability
        lfc_true[i] = np.random.normal(-3.2, 0.5)
    elif gene in resistance_genes:
        # Resistance: enriched
        lfc_true[i] = np.random.normal(+2.5, 0.6)
    else:
        # Non-essential: small noise around 0
        lfc_true[i] = np.random.normal(0.0, 0.3)

# Simulate 3 replicates from the true LFC
def sim_day14(plasmid, lfc, n_reps=3, dispersion=5):
    reps = []
    for _ in range(n_reps):
        mu = plasmid * 2**lfc
        counts = nbinom.rvs(n=dispersion, p=dispersion/(dispersion+mu))
        counts = np.maximum(counts, 0)
        reps.append(counts)
    return np.array(reps).T   # shape: (n_sgrnas, n_reps)

day14_counts = sim_day14(plasmid_counts, lfc_true)  # (n_sgrnas, 3)

count_df = pd.DataFrame({
    'sgRNA': sgrna_ids,
    'Gene': gene_labels,
    'Plasmid': plasmid_counts,
    'Day14_R1': day14_counts[:, 0],
    'Day14_R2': day14_counts[:, 1],
    'Day14_R3': day14_counts[:, 2],
})

print(f"\nCount table shape: {count_df.shape}")
print(count_df.head(8).to_string(index=False))
```python

## 3. MAGeCK test: Robust Rank Aggregation

> RRA algorithm: rank sgRNAs by LFC, aggregate rankings per gene. Beta score for gene essentiality. Negative binomial model for p-value estimation.

```python
# Example: MAGeCK test
# !mageck test -k sample.count.txt -t day14_rep1,day14_rep2 -c plasmid \
#   -n mageck_output --gene-lfc-method median
```python

## 4. Hit Calling and Interpretation

> FDR threshold for positive and negative selection hits. Compare to DepMap CERES/Chronos essentiality scores. Enrichment of essential gene categories (proteasome, spliceosome, ribosome).

## 5. Visualization

> Volcano plot (LFC vs -log10 FDR). Rank plot highlighting known essential genes. sgRNA-level scatter plot for top hits. Gini index for library uniformity.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
