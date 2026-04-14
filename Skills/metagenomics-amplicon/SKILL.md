---
name: metagenomics-amplicon
description: 16S/amplicon metagenomics diversity and community analysis.
tool_type: mixed
primary_tool: Python
---

# 16S Amplicon Metagenomics

## When to Use
- 16S rRNA amplicon microbiome studies (OTU/ASV, alpha/beta diversity)
- Promoter scanning: TATA box detection, CpG island calling, PWM/TFBS scanning

## Quick Reference

### 16S Alpha Diversity
| Metric | Formula | Notes |
|--------|---------|-------|
| Observed species | `count(counts > 0)` | Richness only |
| Shannon H' | `-Σ pᵢ ln(pᵢ)` | Richness + evenness |
| Simpson 1-D | `1 - Σ pᵢ²` | Robust to rare taxa |
| Faith's PD | `Σ branch lengths` | Phylogenetic richness |

### 16S Beta Diversity
| Metric | Phylogenetic | Quantitative | Formula |
|--------|:-----------:|:------------:|---------|
| Jaccard | No | No | `1 - \|A∩B\| / \|A∪B\|` |
| Bray-Curtis | No | Yes | `Σ\|Xᵢⱼ - Xᵢₖ\| / Σ(Xᵢⱼ + Xᵢₖ)` |
| Unweighted UniFrac | Yes | No | unique / total observed branch length |
| Weighted UniFrac | Yes | Yes | abundance-weighted branch fractions |

### Promoter Key Elements
| Element | Position | Consensus | Notes |
|---------|----------|-----------|-------|
| TATA box | −25 to −30 bp | TATAAA (TATAWAW) | ~10–20% human genes |
| CpG island | Near TSS | GC ≥ 50%, CpG O/E ≥ 0.6, ≥ 200 bp | ~70% human promoters |
| Inr | +1 (TSS) | YYANWYY | TATA-less promoters |
| DPE | +28 to +32 | RGWYV | Pairs with Inr |

CpG O/E: `(N_CpG × L) / (N_C × N_G)`

## Key Patterns

### Shannon / Simpson
```python
def shannon_diversity(counts):
    p = counts[counts > 0] / counts[counts > 0].sum()
    return -np.sum(p * np.log(p))

def simpson_diversity(counts):
    p = counts[counts > 0] / counts[counts > 0].sum()
    return 1 - np.sum(p ** 2)
```

### Bray-Curtis Distance
```python
def bray_curtis(s1, s2):
    s1, s2 = np.array(s1, float), np.array(s2, float)
    return np.sum(np.abs(s1 - s2)) / np.sum(s1 + s2)
```

### PCoA from Distance Matrix
```python
def pcoa(dm_df):
    dm = dm_df.values; n = len(dm)
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dm ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 0
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    prop = eigvals[pos] / eigvals[pos].sum()
    return pd.DataFrame(coords[:, :2], index=dm_df.index, columns=['PC1', 'PC2']), prop
```

### CpG Island Scanner
```python
def cpg_island_scanner(seq, window=200, step=1, gc_thresh=0.5, oe_thresh=0.6):
    seq = seq.upper(); islands = []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i+window]
        nc, ng, ncpg = w.count('C'), w.count('G'), w.count('CG')
        gc = (nc + ng) / window
        oe = (ncpg * window) / (nc * ng) if nc > 0 and ng > 0 else 0
        if gc >= gc_thresh and oe >= oe_thresh:
            islands.append((i, i+window, gc, oe))
    return islands
```

### PWM Construction and Scanning
```python
def build_pfm(sites):
    pfm = {b: [0]*len(sites[0]) for b in 'ACGT'}
    for site in sites:
        for i, b in enumerate(site.upper()): pfm[b][i] += 1
    return pfm

def pfm_to_pwm(pfm, pseudocount=0.5, bg=None):
    if bg is None: bg = {b: 0.25 for b in 'ACGT'}
    n = sum(pfm[b][0] for b in 'ACGT')
    return {b: [np.log2((pfm[b][i]+pseudocount)/(n+4*pseudocount)/bg[b])
                for i in range(len(pfm['A']))] for b in 'ACGT'}

def scan_pwm(seq, pwm, threshold=0.0):
    L = len(pwm['A']); hits = []
    for i in range(len(seq) - L + 1):
        kmer = seq[i:i+L].upper()
        score = sum(pwm[b][j] for j, b in enumerate(kmer) if b in pwm)
        if score >= threshold: hits.append((i, kmer, score))
    return hits
```

## Code Templates

### PERMANOVA
```python
def permanova(dm, grouping, n_perm=999):
    dm_v = dm.values; groups = grouping.loc[dm.index].values
    unique = np.unique(groups); n = len(groups)
    def f_stat(g):
        ss_tot = np.sum(dm_v**2) / n
        ss_w = sum(np.sum(dm_v[np.ix_(g==grp,g==grp)]**2)/(2*(g==grp).sum())
                   for grp in unique if (g==grp).sum() > 1)
        df_b, df_w = len(unique)-1, n-len(unique)
        return ((ss_tot-ss_w)/df_b) / (ss_w/df_w) if df_w and ss_w else 0
    obs = f_stat(groups)
    p = (sum(f_stat(np.random.permutation(groups)) >= obs for _ in range(n_perm)) + 1) / (n_perm + 1)
    return obs, p
```

### TATA Box Detection
```python
import re
def find_tata_boxes(seq, strict=True):
    pat = 'TATAAA' if strict else r'TATA[AT]A[AT]'
    return [(m.start(), m.group()) for m in re.finditer(pat, seq.upper())]
```

## Pitfalls

- **Rarefying without curves** — rarefy only when rarefaction curves plateau; raw alpha diversity requires rarefaction (deeper samples appear more diverse).
- **PCoA axes without variance-explained** — axis 1 may explain only 20%; check percentages.
- **PERMANOVA + dispersion** — sensitive to within-group variance, not just centroids; pair with betadisper (R).
- **CpG step=1 on large sequences** — use `step=10–50` for >10 kb, then merge overlapping windows.
- **PWM threshold** — use 80–90th percentile of random-sequence scores as baseline.
- **Multiple testing in ORA** — BH-correct across all tested pathways.

## Related Skills
- `biostatistics-r` — hypothesis testing, BH correction, power analysis
- `numpy-pandas-wrangling` — count matrix manipulation
- `python-advanced-sql` — annotation lookups
