---
name: bio-applied-copy-number-analysis
description: "DNA copy number analysis — read depth normalization, CBS segmentation, CN state calling, and genome-wide visualization"
tool_type: python
primary_tool: Pandas
---

# DNA Copy Number Analysis

## Copy Number Reference

| Type | CN | Log2 ratio | Effect |
|---|---|---|---|
| Homozygous deletion | 0 | -inf | Gene loss |
| Heterozygous deletion | 1 | -1.0 | LOH if other allele mutated |
| Normal diploid | 2 | 0.0 | Baseline |
| Single-copy gain | 3 | +0.585 | Increased gene dosage |
| Amplification | 4+ | +1.0+ | Common in oncogenes (ERBB2, MYCN) |

**Tumour purity**: at 50% purity, effective CN = 0.5*tumor_CN + 0.5*2, attenuating log2 signal. Tools like CNVkit, Sequenza, PURPLE jointly estimate purity and ploidy.

## Depth Normalization

```python
from scipy.ndimage import median_filter

def normalize_depth(depth, window=50):
    """Normalize by local rolling median to correct GC/mappability bias."""
    local_median = median_filter(depth.astype(float), size=window, mode="reflect")
    norm = depth / np.where(local_median > 0.1, local_median, 0.1)
    return norm / np.median(norm)

log2_ratio = np.log2(np.clip(normalize_depth(depth), 0.01, None))
```

## Circular Binary Segmentation (CBS)

Recursively splits chromosomes at the point maximizing absolute mean difference between resulting segments.

```python
def cbs_segment(log2_ratios, min_segment=10, alpha=0.01):
    """Simplified CBS: split at max absolute mean difference."""
    segments = []
    def _split(lo, hi):
        if hi - lo < min_segment * 2:
            segments.append((lo, hi, log2_ratios[lo:hi].mean()))
            return
        best_delta, best_k = 0, None
        for k in range(lo + min_segment, hi - min_segment):
            delta = abs(log2_ratios[lo:k].mean() - log2_ratios[k:hi].mean())
            if delta > best_delta:
                best_delta, best_k = delta, k
        if best_k is None or best_delta < 0.3:
            segments.append((lo, hi, log2_ratios[lo:hi].mean()))
            return
        _split(lo, best_k); _split(best_k, hi)
    _split(0, len(log2_ratios))
    return sorted(segments)
```

## CN State Calling

```python
def call_cn_state(mean_log2, thresholds=None):
    """Map log2 ratio to integer CN using midpoint thresholds."""
    if thresholds is None:
        thresholds = {0: float('-inf'), 1: -1.5, 2: -0.4, 3: 0.35, 4: 0.80, 5: 1.25}
    for cn in sorted(thresholds.keys(), reverse=True):
        if mean_log2 >= thresholds[cn]:
            return cn
    return 0
```

## Gene-Level Annotation

```python
import pyranges as pr
genes = pr.read_bed("hg38_genes.bed")
cnvs = pr.from_dict({"Chromosome": ["chr1"], "Start": [2500000], "End": [3750000], "CN": [1]})
overlap = cnvs.join(genes)
```

**Downstream**: amplified oncogenes (COSMIC Cancer Gene Census), deleted tumor suppressors (TP53, BRCA1/2, RB1). Focal events (<3 Mb) are more often drivers than arm-level events.

**TCGA data**: `gdc.cancer.gov` portal, `TCGAbiolinks` (R), or `xenahubs.net`. GISTIC2 identifies recurrent amplifications/deletions.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
