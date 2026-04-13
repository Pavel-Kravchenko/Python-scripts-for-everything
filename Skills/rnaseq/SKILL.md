---
name: rnaseq
description: RNA-seq differential expression — DESeq2, edgeR, normalization strategies (RPKM/TPM/DESeq2 size factors), exploratory PCA, volcano plots, GSEA, STAR/Salmon/featureCounts pipelines
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, pydeseq2 0.4+, scikit-learn 1.4+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# RNA-seq Analysis & Differential Expression

## When to Use
- Differential gene expression from RNA-seq count data
- Comparing expression normalization strategies (RPKM / TPM / DESeq2 median-of-ratios)
- Building and interpreting count matrices (genes × samples)
- Exploratory analysis: PCA, hierarchical clustering of samples
- Downstream interpretation: volcano plots, MA plots, GSEA / ORA

---

## Quick Reference

### RNA-seq Workflow
```
FASTQ → QC (FastQC/fastp) → Alignment (STAR/HISAT2) → Count (featureCounts/HTSeq)
                                  OR
                          → Pseudoalignment (Salmon/kallisto) → tximport
→ Count matrix (genes × samples) → Normalization → DESeq2/edgeR → Volcano/MA plots → GSEA
```

### Normalization Formulas
| Unit | Formula | Use case |
|------|---------|----------|
| RPKM | `(C / L) / N × 10⁹` | Single-end, within-sample comparison |
| FPKM | Same as RPKM | Paired-end (counts fragments) |
| TPM | `(C/L) / Σ(Cⱼ/Lⱼ) × 10⁶` | Cross-sample comparison (sums to 1 M) |
| DESeq2 size factor | `median(sample_counts / geo_mean_across_samples)` | DE testing |

C = read count, L = gene length (bp), N = total mapped reads.  
**Use TPM for reporting/visualization; raw counts + DESeq2 for DE testing.**

### Aligner Comparison
| Aligner | Algorithm | Memory | Best for |
|---------|-----------|--------|----------|
| STAR | Suffix array (2-pass) | 30+ GB | Human/mouse; spliced alignment |
| HISAT2 | Graph FM-index | ~8 GB | Lower-memory splice-aware alignment |
| Salmon | Quasi-mapping | Low | Fast transcript-level quantification |
| kallisto | Pseudoalignment | Very low | Ultra-fast; good for exploratory work |

### Statistical Testing
| Tool | Model | Key Feature |
|------|-------|-------------|
| DESeq2 | Negative binomial, shrunk dispersion | Median-of-ratios normalisation; apeglm LFC shrinkage |
| edgeR | Negative binomial, empirical Bayes | TMM normalisation; quasi-likelihood F-test |
| limma-voom | Linear model on voom-transformed counts | Most flexible for complex designs |
| PyDESeq2 | Python port of DESeq2 | Full Python workflow; compatible result format |

---

## Key Patterns

### DESeq2 Size Factors (Python)
```python
import numpy as np
import pandas as pd

def deseq2_size_factors(count_matrix: pd.DataFrame) -> pd.Series:
    """Median-of-ratios size factors. count_matrix: genes × samples."""
    nonzero_mask = (count_matrix > 0).all(axis=1)
    filtered = count_matrix.loc[nonzero_mask]
    geo_means = np.exp(np.log(filtered).mean(axis=1))
    return pd.Series({s: np.median(filtered[s] / geo_means)
                      for s in filtered.columns})

normalized = counts_df.div(deseq2_size_factors(counts_df), axis=1)
```

### TPM from Raw Counts
```python
def counts_to_tpm(counts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """counts: 1-D array per sample; lengths: gene lengths in bp."""
    rate = counts / lengths          # normalize by gene length
    return rate / rate.sum() * 1e6  # scale to 1 M
```

### Simple DE (BH-corrected Mann-Whitney, Python)
```python
from scipy import stats

def simple_de(count_df, ctrl_samples, treat_samples):
    sf = deseq2_size_factors(count_df)
    norm = count_df.div(sf, axis=1)
    results = []
    for gene in count_df.index:
        c = norm.loc[gene, ctrl_samples]
        t = norm.loc[gene, treat_samples]
        lfc = np.log2((t.mean() + 1) / (c.mean() + 1))
        _, p = stats.mannwhitneyu(c, t, alternative='two-sided')
        results.append({'gene': gene, 'log2FC': lfc, 'pvalue': p,
                        'baseMean': (c.mean() + t.mean()) / 2})
    df = pd.DataFrame(results).set_index('gene')
    # BH correction
    p_arr = df['pvalue'].values
    n = len(p_arr)
    idx = np.argsort(p_arr)
    adj = np.minimum.accumulate((p_arr[idx] * n / (np.arange(n) + 1))[::-1])[::-1]
    df['padj'] = 0.0
    df.iloc[idx, df.columns.get_loc('padj')] = adj
    return df.sort_values('pvalue')
```

---

## Code Templates

### DESeq2 (R — production workflow)
```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData   = sample_info,
                              design    = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Treatment", "Control"))
res_shrunk <- lfcShrink(dds, coef = "condition_Treatment_vs_Control", type = "apeglm")
sig <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
```

### PyDESeq2 (Python — production workflow)
```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_matrix, metadata=sample_info,
                   design_factors="condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition", "Treatment", "Control"])
stat_res.summary()
results_df = stat_res.results_df
```

### STAR + featureCounts (shell)
```bash
# Build index
STAR --runMode genomeGenerate --genomeDir star_index/ \
     --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf

# Align
STAR --genomeDir star_index/ --readFilesIn R1.fastq R2.fastq \
     --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

# Count
featureCounts -a genes.gtf -o counts.txt -T 4 -p --countReadPairs *.bam
```

### Salmon (alignment-free quantification)
```bash
salmon index -t transcriptome.fa -i salmon_index
salmon quant -i salmon_index -l A -1 R1.fastq -2 R2.fastq \
             -o sample_quant --validateMappings
```

### Sample PCA (Python)
```python
from sklearn.decomposition import PCA

def rnaseq_pca(norm_counts: pd.DataFrame, metadata: pd.DataFrame,
               color_col: str = 'condition') -> None:
    """PCA of samples. norm_counts: genes × samples (log2 + 1 inside)."""
    import matplotlib.pyplot as plt
    log_data = np.log2(norm_counts + 1)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(log_data.T)   # samples × 2
    ev = pca.explained_variance_ratio_
    fig, ax = plt.subplots(figsize=(7, 5))
    for cond, grp in metadata.groupby(color_col):
        idx = [norm_counts.columns.get_loc(s) for s in grp.index]
        ax.scatter(coords[idx, 0], coords[idx, 1], label=cond, s=60)
    ax.set_xlabel(f'PC1 ({ev[0]:.1%})')
    ax.set_ylabel(f'PC2 ({ev[1]:.1%})')
    ax.legend(); plt.tight_layout(); plt.show()
```

### Volcano Plot (Python)
```python
import matplotlib.pyplot as plt

def volcano_plot(results_df: pd.DataFrame,
                 lfc_col='log2FC', pval_col='padj',
                 lfc_thresh=1.0, p_thresh=0.05) -> None:
    lfc = results_df[lfc_col].values
    pvals = results_df[pval_col].fillna(1).values
    neg_log10p = -np.log10(np.clip(pvals, 1e-300, 1.0))
    sig = (np.abs(lfc) > lfc_thresh) & (pvals < p_thresh)
    colors = np.where(sig, 'tomato', 'steelblue')
    plt.figure(figsize=(8, 6))
    plt.scatter(lfc, neg_log10p, c=colors, alpha=0.5, s=12)
    plt.axvline( lfc_thresh, color='grey', ls='--', lw=0.8)
    plt.axvline(-lfc_thresh, color='grey', ls='--', lw=0.8)
    plt.axhline(-np.log10(p_thresh), color='grey', ls='--', lw=0.8)
    plt.xlabel('log₂ fold change')
    plt.ylabel('-log₁₀(adj. p-value)')
    plt.title(f'Volcano — {sig.sum()} significant genes')
    plt.tight_layout(); plt.show()
```

---

## Common Pitfalls

- **Using RPKM/TPM for DE testing** — they introduce ratio-of-ratios artifacts. Always pass raw counts to DESeq2/edgeR.
- **Not accounting for library composition** — a single highly-expressed gene can suppress apparent expression of all others. DESeq2 median-of-ratios is robust to this; CPM/TPM are not.
- **Skipping dispersion shrinkage** — DESeq2 borrows information across genes to stabilize dispersion estimates for low-count genes. Manual NB tests miss this and have lower power.
- **Forgetting paired/batch design** — include batch as a covariate in `design = ~batch + condition`; omitting it inflates false positives.
- **Interpreting LFC without shrinkage** — MLE log-fold changes for low-count genes are noisy; always use `lfcShrink` (apeglm) before reporting or visualizing LFC.
- **Filtering after DE** — apply independent filtering before testing (DESeq2 does this automatically); hard count cut-offs can bias results.
- **Multiple testing in ORA** — apply BH correction across all tested pathways; individual Fisher p-values alone inflate false discovery.

---

## Related Skills
- `metagenomics-shotgun` — 16S amplicon diversity analysis, QIIME2, alpha/beta diversity
- `promoter-regulatory-analysis` — promoter scanning, CpG islands, PWM/TFBS (Module 05)
- `biostatistics-r` — hypothesis testing, BH correction, power analysis
- `numpy-pandas-wrangling` — count matrix manipulation, filtering, melting
- `cancer-transcriptomics` — tumour-specific DE, TCGA, single-cell workflows
