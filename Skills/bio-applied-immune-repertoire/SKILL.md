---
name: bio-applied-immune-repertoire
description: Immune repertoire sequencing — TRUST4/MiXCR workflows, diversity metrics, clonal tracking, Morisita-Horn overlap, VDJdb lookup, CDR3 Hamming clustering.
tool_type: python
primary_tool: Matplotlib
---

# Immune Repertoire Sequencing and Clonal Analysis

- [TRUST4](https://github.com/liulab-dfci/TRUST4)
- [MiXCR](https://docs.milaboratories.com/mixcr/)
- [VDJdb](https://vdjdb.cdr3.net/)

## Pitfalls

- **Clonotype definition matters**: Using CDR3 amino acid sequence alone merges convergent clones from different V/J genes. Use V+J+CDR3aa as the minimal clonotype key for tracking.
- **Productive/total ratio**: Expect >80% for good-quality libraries. Low ratio indicates poor library prep or contamination.
- **PCR amplification bias (MiXCR)**: Amplicon-based libraries require UMI correction before comparing clone frequencies across samples.
- **Bulk repertoire from RNA-seq (TRUST4)**: Lower sensitivity than targeted methods — small clones may not be recovered. Suitable for TIL expansion detection, not for low-frequency monitoring.
- **Morisita-Horn vs Jaccard**: Morisita-Horn accounts for clone frequency (recommended for repertoire overlap). Jaccard is presence/absence only and ignores dominant clones.
- **Public clones**: ~1–5% of T-cell repertoire is public (shared across individuals); B-cell overlap is much lower. Enrichment of known public clones in a sample can indicate antigen-driven expansion.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.

## Tools Summary

| Tool | Input | Use case |
|---|---|---|
| TRUST4 | RNA-seq BAM | TCR/BCR from bulk RNA-seq, no targeted prep |
| MiXCR | Amplicon FASTQ | Targeted BCR/TCR, gold standard for accuracy |
| Immcantation | MiXCR/IMGT output | B-cell SHM, lineage trees, isotype switching |
| VDJdb | CDR3 + V/J | Antigen-specific TCR lookup |
| GLIPH2 / TCRdist | CDR3 sequences | Group TCRs by putative specificity |

## MiXCR Command
```bash
mixcr analyze amplicon \
  --species hsa \
  --starting-material rna \
  --5-end v-primers \
  --3-end j-primers \
  sample_R1.fastq.gz sample_R2.fastq.gz result
```

## Key QC Metrics

- **Total productive reads**: reads with valid CDR3 + in-frame
- **Clonotype count**: unique V+J+CDR3 combinations
- **D50 index**: number of top clones covering 50% of reads (lower = more clonal)
- **Productive/total ratio**: >80% expected

## Diversity Metrics

```python
import numpy as np

def shannon_entropy(freqs):
    freqs = np.array(freqs); freqs = freqs[freqs > 0]
    return -np.sum(freqs * np.log(freqs))

def normalized_shannon(freqs):
    n = len([f for f in freqs if f > 0])
    h = shannon_entropy(freqs)
    return h / np.log(n) if n > 1 else 0

def simpson_d(freqs):
    return np.sum(np.array(freqs)**2)

def clonality(freqs):
    return 1 - normalized_shannon(freqs)

def d50_index(freqs):
    freqs_sorted = np.sort(freqs)[::-1]
    cumsum = np.cumsum(freqs_sorted)
    return int(np.searchsorted(cumsum, 0.5)) + 1
```

## Repertoire Overlap

```python
def morisita_horn(freq1, freq2):
    """Morisita-Horn index. Range 0 (no overlap) to 1 (identical)."""
    shared = set(freq1) & set(freq2)
    if not shared:
        return 0.0
    numerator = 2 * sum(freq1.get(c, 0) * freq2.get(c, 0) for c in shared)
    lambda1 = sum(v**2 for v in freq1.values())
    lambda2 = sum(v**2 for v in freq2.values())
    return numerator / (lambda1 + lambda2)

def jaccard(set1, set2):
    return len(set1 & set2) / len(set1 | set2)
```

## CDR3 Hamming Distance Clustering

```python
def hamming(s1, s2):
    if len(s1) != len(s2):
        return np.nan
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Pairwise distance matrix for equal-length CDR3s
cdr3_seqs = repertoire[repertoire['junction_aa'].str.len() == 14]['junction_aa'].values
dist_matrix = np.array([[hamming(a, b) for b in cdr3_seqs] for a in cdr3_seqs], dtype=float)
# Clonal families: CDR3s with Hamming ≤ 1 are likely same clone
```

## VDJdb Antigen-Specific Lookup

```python
import pandas as pd

# VDJdb entry format: CDR3aa, V-gene, J-gene, Species, Epitope, MHC, Source
vdjdb = pd.read_csv('vdjdb.tsv', sep='\t')

# Match repertoire clones to known antigen-specific sequences
matched = repertoire.merge(
    vdjdb[['junction_aa', 'epitope', 'mhc']],
    on='junction_aa', how='inner'
)
```

## Clonal Tracking (Pre/Post Treatment)

```python
def track_clones(pre_df, post_df, id_cols=('v_call', 'j_call', 'junction_aa')):
    """Classify clones as expanded, contracted, new, or lost."""
    pre  = pre_df.set_index(list(id_cols))['freq']
    post = post_df.set_index(list(id_cols))['freq']
    all_clones = pre.index.union(post.index)

    results = []
    for clone in all_clones:
        f_pre  = pre.get(clone, 0)
        f_post = post.get(clone, 0)
        if f_pre == 0:
            status = 'new'
        elif f_post == 0:
            status = 'lost'
        else:
            log2fc = np.log2(f_post / f_pre)
            status = 'expanded' if log2fc > 1 else ('contracted' if log2fc < -1 else 'persistent')
        results.append({'clone': clone, 'freq_pre': f_pre, 'freq_post': f_post, 'status': status})
    return pd.DataFrame(results)
```
