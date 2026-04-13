---
name: dna-methylation
description: "Bisulfite sequencing analysis: WGBS/RRBS processing, DMR calling, and epigenetic clock estimation."
tool_type: mixed
primary_tool: Pandas
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# dna-methylation

Bisulfite sequencing analysis: WGBS/RRBS processing, DMR calling, and epigenetic clock estimation.

## Quick Reference

| Step | Tool | Notes |
|------|------|-------|
| Genome prep | `bismark_genome_preparation` | Creates CT/GA converted index |
| Alignment | `bismark` | Bowtie2 backend |
| Deduplication | `deduplicate_bismark` | Remove PCR duplicates |
| Extraction | `bismark_methylation_extractor` | CpG/CHH/CHG contexts |
| Load in R | `methylKit::methRead()` | From bismark coverage files |
| DMP test | `methylKit::calculateDiffMeth()` | Logistic regression |
| DMR calling | `BSmooth::BSmooth()` | Smoothing + t-statistics |
| Age clock | `methylclock::DNAmAge()` | Horvath/GrimAge |

## Bismark Pipeline

```bash
# 1. Prepare bisulfite genome index
bismark_genome_preparation /path/to/genome/

# 2. Trim adapters (for WGBS)
trim_galore --paired --fastqc sample_R1.fastq.gz sample_R2.fastq.gz

# 3. Align
bismark --genome /path/to/genome/ -1 sample_R1_val_1.fq.gz -2 sample_R2_val_2.fq.gz

# 4. Remove duplicates
deduplicate_bismark --paired sample_bismark_bt2_pe.bam

# 5. Extract methylation
bismark_methylation_extractor --paired-end --CpG --comprehensive \
    --cytosine_report --genome_folder /path/to/genome/ \
    sample_bismark_bt2_pe.deduplicated.bam

# 6. Coverage statistics
bismark2report
```python

## Load and Visualize in Python

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load bismark coverage output (6-column format)
cpg = pd.read_csv('sample.bismark.cov.gz', sep='\t', compression='gzip',
                  names=['chrom', 'start', 'end', 'pct_meth', 'count_M', 'count_U'])
cpg['methylation'] = cpg['count_M'] / (cpg['count_M'] + cpg['count_U'])
cpg['coverage'] = cpg['count_M'] + cpg['count_U']

# Filter by coverage
cpg_filtered = cpg[cpg['coverage'] >= 10]

# Methylation distribution
plt.hist(cpg_filtered['methylation'], bins=50, edgecolor='black')
plt.xlabel('Methylation level'); plt.ylabel('Count')
plt.title('CpG methylation distribution')
plt.show()
```python

## methylKit DMR Analysis (R)

```r
library(methylKit)

# Load samples
methyl_list <- methRead(
  list('ctrl.cov', 'treated.cov'),
  sample.id = list('ctrl', 'treated'),
  assembly = 'hg38', treatment = c(0, 1),
  context = 'CpG', mincov = 10
)

# Filter and normalize
filtered <- filterByCoverage(methyl_list, lo.count=10, hi.perc=99.9)
normalized <- normalizeCoverage(filtered)

# Merge and test
united <- unite(normalized, destrand=FALSE)
dm_results <- calculateDiffMeth(united)

# Filter significant DMPs
dmps <- getMethylDiff(dm_results, difference=25, qvalue=0.05)
```python

## Horvath Epigenetic Clock (Python)

```python
import pandas as pd
import numpy as np

# Load Horvath clock coefficients (353 CpGs)
clock = pd.read_csv('horvath_clock_cpgs.csv')  # CpGmarker, CoefficientTraining

# Fetch methylation at clock CpGs from beta matrix
beta_clock = beta_matrix[clock['CpGmarker']].T

# Anti-logit transformation
def anti_trafo(x, adult_age=20):
    if x < 0:
        return (1 + adult_age) * np.exp(x) - 1
    else:
        return (1 + adult_age) * x + adult_age

dnam_age = beta_clock.T.dot(clock.set_index('CpGmarker')['CoefficientTraining'])
# Apply transformation per sample
predicted_age = dnam_age.apply(anti_trafo)
```python

## Common Pitfalls
- **Non-conversion rate > 1%**: poor bisulfite conversion; check lambda spike-in
- **Coverage asymmetry**: RRBS enriches CpG-dense regions (islands); WGBS is genome-wide
- **Strand merging**: use `destrand=TRUE` in methylKit for symmetric CpG analysis
- **Array vs WGBS**: EPIC covers ~850K CpGs; WGBS provides ~28M; overlap ~80%
- **Clock tissue bias**: Horvath is pan-tissue; Hannum is blood-specific

## Module
Tier 3 · Module 32 (DNA Methylation Analysis)
