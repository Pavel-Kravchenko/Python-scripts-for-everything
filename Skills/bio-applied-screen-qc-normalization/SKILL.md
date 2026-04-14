---
name: bio-applied-screen-qc-normalization
description: "CRISPR screen quality control and normalization: library distribution QC, Gini index, replicate correlation, and count normalization. Use when processing CRISPR screen count matrices."
tool_type: python
primary_tool: MAGeCK
---

# CRISPR Screen QC, Normalization, and Advanced Methods

**By the end of this skill you will be able to:**
1. Assess screen quality with Gini index, read depth distribution, and sgRNA evenness
2. Normalize CRISPR screen count matrices (median ratio, total count)
3. Evaluate replicate concordance
4. Detect and handle copy-number bias

## Screen Quality Metrics

**Goal:** Assess whether a CRISPR screen has sufficient quality for downstream analysis.

**Approach:** Calculate Gini index (library evenness), read depth per guide, fraction of zero-count guides, and replicate correlation.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def gini_index(counts: np.ndarray) -> float:
    """Calculate Gini index of guide count distribution.

    Perfect evenness = 0, one guide has all reads = 1.
    Good screens: Gini < 0.1. Problematic: Gini > 0.2.
    """
    sorted_counts = np.sort(counts)
    n = len(sorted_counts)
    cumulative = np.cumsum(sorted_counts)
    return (2 * np.sum((np.arange(1, n + 1) * sorted_counts)) /
            (n * np.sum(sorted_counts)) - (n + 1) / n)


def screen_qc_report(count_df: pd.DataFrame) -> dict:
    """Generate QC metrics for a CRISPR screen count matrix.

    Args:
        count_df: DataFrame with guides as rows, samples as columns

    Returns:
        dict with QC metrics per sample
    """
    metrics = {}
    for col in count_df.columns:
        counts = count_df[col].values
        metrics[col] = {
            'total_reads': int(counts.sum()),
            'guides_detected': int((counts > 0).sum()),
            'zero_count_guides': int((counts == 0).sum()),
            'zero_fraction': float((counts == 0).mean()),
            'gini': float(gini_index(counts[counts > 0])),
            'median_count': float(np.median(counts)),
            'mean_count': float(np.mean(counts)),
        }
    return metrics


# Example usage with synthetic data
np.random.seed(42)
n_guides = 5000
count_matrix = pd.DataFrame({
    'plasmid': np.random.negative_binomial(5, 0.01, n_guides),
    'replicate_1': np.random.negative_binomial(4, 0.01, n_guides),
    'replicate_2': np.random.negative_binomial(4, 0.01, n_guides),
})

qc = screen_qc_report(count_matrix)
for sample, m in qc.items():
    print(f"{sample}: Gini={m['gini']:.3f}, zero={m['zero_fraction']:.1%}, "
          f"median={m['median_count']:.0f}")
```

## Count Normalization

**Goal:** Remove technical variation (sequencing depth differences) while preserving biological signal.

**Approach:** Median-ratio normalization (same principle as DESeq2 size factors) — divide each sample by its size factor derived from the median of ratios to a pseudo-reference.

```python
def median_ratio_normalize(count_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize counts using median-ratio method (DESeq2-style).

    1. Compute geometric mean per guide across samples (pseudo-reference)
    2. Divide each count by the pseudo-reference
    3. Size factor = median of ratios for each sample
    4. Divide raw counts by size factor
    """
    # Geometric mean per guide (exclude zeros)
    log_means = np.log(count_df + 1).mean(axis=1)
    geo_means = np.exp(log_means)

    # Ratio of each count to pseudo-reference
    ratios = count_df.div(geo_means, axis=0)

    # Size factor = median of ratios per sample
    size_factors = ratios.replace([np.inf, -np.inf], np.nan).median(axis=0)

    # Normalized counts
    normalized = count_df.div(size_factors, axis=1)
    return normalized


norm_counts = median_ratio_normalize(count_matrix)
print("Size factors:", (count_matrix.sum() / norm_counts.sum()).round(3).to_dict())
```

## Replicate Concordance

**Goal:** Verify that biological replicates agree before merging or averaging.

**Approach:** Pearson/Spearman correlation of log-fold-changes between replicates. Good screens: r > 0.7.

```python
def replicate_correlation(count_df: pd.DataFrame, control: str,
                          rep1: str, rep2: str) -> dict:
    """Calculate LFC correlation between two replicates vs control.

    Args:
        count_df: Normalized count matrix
        control: Column name for plasmid/control
        rep1, rep2: Column names for replicates

    Returns:
        dict with Pearson r, Spearman rho, and p-values
    """
    pseudo = 0.5
    lfc1 = np.log2((count_df[rep1] + pseudo) / (count_df[control] + pseudo))
    lfc2 = np.log2((count_df[rep2] + pseudo) / (count_df[control] + pseudo))

    pearson_r, pearson_p = stats.pearsonr(lfc1, lfc2)
    spearman_r, spearman_p = stats.spearmanr(lfc1, lfc2)

    return {
        'pearson_r': pearson_r, 'pearson_p': pearson_p,
        'spearman_rho': spearman_r, 'spearman_p': spearman_p,
    }


corr = replicate_correlation(norm_counts, 'plasmid', 'replicate_1', 'replicate_2')
print(f"Replicate concordance: Pearson r={corr['pearson_r']:.3f}, "
      f"Spearman rho={corr['spearman_rho']:.3f}")
```

## Copy-Number Bias Detection

**Goal:** Identify and flag guides in amplified genomic regions where high copy number inflates dropout signal.

**Approach:** Group guides by genomic segment, compare median LFC per segment to genome-wide median.

```python
def detect_cn_bias(lfc: pd.Series, guide_locations: pd.DataFrame,
                   segment_size: int = 50) -> pd.DataFrame:
    """Flag genomic segments with copy-number-biased LFC.

    Segments where median LFC deviates > 2 MADs from genome median
    are flagged as CN-biased.
    """
    guide_locations = guide_locations.copy()
    guide_locations['lfc'] = lfc.values

    # Group by chromosome and segment
    guide_locations['segment'] = guide_locations['position'] // segment_size

    seg_stats = guide_locations.groupby(['chr', 'segment'])['lfc'].agg(
        ['median', 'count']
    ).reset_index()

    genome_median = lfc.median()
    mad = np.median(np.abs(lfc - genome_median))

    seg_stats['cn_biased'] = np.abs(seg_stats['median'] - genome_median) > 2 * mad
    return seg_stats


# Synthetic guide locations
guide_locs = pd.DataFrame({
    'chr': np.random.choice(['chr1', 'chr2', 'chr3'], n_guides),
    'position': np.random.randint(0, 10000, n_guides),
})
lfc = np.log2((norm_counts['replicate_1'] + 0.5) / (norm_counts['plasmid'] + 0.5))
cn_result = detect_cn_bias(lfc, guide_locs)
print(f"CN-biased segments: {cn_result['cn_biased'].sum()} / {len(cn_result)}")
```

## Key QC Thresholds

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| Gini index | < 0.1 | 0.1–0.2 | > 0.2 |
| Zero-count guides | < 1% | 1–5% | > 5% |
| Replicate Pearson r | > 0.8 | 0.6–0.8 | < 0.6 |
| Reads per guide (median) | > 300 | 100–300 | < 100 |
| Essential gene depletion | > 5σ | 3–5σ | < 3σ |

## Pitfalls

- **Low MOI violations**: If MOI > 0.3, multiple guides per cell corrupt phenotype-guide assignment
- **Plasmid library bias**: Always compare to a plasmid control, not just between treatment groups
- **Batch effects**: Screen date and cell passage number can dominate biological signal; include batch in the model
- **Copy-number confounding**: Amplified loci show guide dropout even without essentiality — always check CN bias before calling hits
