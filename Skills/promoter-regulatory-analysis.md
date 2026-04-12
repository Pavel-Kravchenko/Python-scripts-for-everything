---
name: promoter-regulatory-analysis
description: Promoter and regulatory motif analysis including CpG and PWM scanning.
---

## When to Use

Use this atomic skill for focused work on **promoter-regulatory-analysis** without bundling unrelated topics.

## Quick Reference

This skill was split from `rnaseq-metagenomics.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

---
name: rnaseq-metagenomics
description: RNA-seq differential expression (DESeq2, normalization), 16S amplicon metagenomics, and promoter/transcription factor analysis
---

# RNA-seq, Metagenomics & Regulatory Analysis

## When to Use
- Differential gene expression from RNA-seq count data
- Comparing expression normalization strategies (RPKM/TPM/DESeq2)
- 16S rRNA amplicon microbiome studies (OTU/ASV, alpha/beta diversity)
- Promoter scanning: TATA box detection, CpG island calling, PWM/TFBS scanning

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
| TPM | `(C/L) / Σ(Cⱼ/Lⱼ) × 10⁶` | Cross-sample comparison (sums to 1M) |
| DESeq2 size factor | `median(sample_counts / geo_mean_across_samples)` | DE testing |

C = read count, L = gene length (bp), N = total mapped reads. **Use TPM for reporting/visualization; raw counts + DESeq2 for DE testing.**

### 16S Alpha Diversity Metrics
| Metric | Formula | Notes |
|--------|---------|-------|
| Observed species | `count(counts > 0)` | Richness only |
| Shannon H' | `-Σ pᵢ ln(pᵢ)` | Richness + evenness; 0 = no diversity |
| Simpson 1-D | `1 - Σ pᵢ²` | Probability two reads differ; robust to rare taxa |
| Faith's PD | `Σ branch lengths (observed taxa)` | Phylogenetic richness |

### 16S Beta Diversity Metrics
| Metric | Phylogenetic | Quantitative | Formula |
|--------|:-----------:|:------------:|---------|
| Jaccard | No | No | `1 - |A∩B| / |A∪B|` |
| Bray-Curtis | No | Yes | `Σ|Xᵢⱼ - Xᵢₖ| / Σ(Xᵢⱼ + Xᵢₖ)` |
| Unweighted UniFrac | Yes | No | unique branch length / total observed branch length |
| Weighted UniFrac | Yes | Yes | abundance-weighted branch length fractions |

### Promoter Key Elements
| Element | Position | Consensus | Notes |
|---------|----------|-----------|-------|
| TATA box | −25 to −30 bp from TSS | TATAAA (TATAWAW) | ~10–20% human genes; tissue-specific |
| CpG island | Near TSS | GC ≥ 50%, CpG O/E ≥ 0.6, ≥ 200 bp | ~70% human gene promoters |
| Inr | +1 (TSS) | YYANWYY | TATA-less promoters |
| DPE | +28 to +32 | RGWYV | Pairs with Inr in TATA-less |

CpG O/E formula: `(N_CpG × L) / (N_C × N_G)`

---

## Key Patterns

### DESeq2 Size Factors (Python)
```python
def deseq2_size_factors(count_matrix):
    nonzero_mask = (count_matrix > 0).all(axis=1)
    filtered = count_matrix.loc[nonzero_mask]
    geo_means = np.exp(np.log(filtered).mean(axis=1))
    return pd.Series({s: np.median(filtered[s] / geo_means)
                      for s in filtered.columns})

normalized = counts_df.div(deseq2_size_factors(counts_df), axis=1)
```

### TPM from Raw Counts
```python
def counts_to_tpm(counts, lengths):
    rate = counts / lengths          # normalize by gene length
    return rate / rate.sum() * 1e6  # scale to 1M
```

### Shannon / Simpson (Python)
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

### PCoA (Classical MDS) from Distance Matrix
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
    return pd.DataFrame(coords[:, :2], index=dm_df.index,
                        columns=['PC1', 'PC2']), prop
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
        for i, b in enumerate(site.upper()):
            pfm[b][i] += 1
    return pfm

def pfm_to_pwm(pfm, pseudocount=0.5, bg=None):
    if bg is None: bg = {b: 0.25 for b in 'ACGT'}
    n = sum(pfm[b][0] for b in 'ACGT')
    return {b: [np.log2((pfm[b][i]+pseudocount)/(n+4*pseudocount)/bg[b])
                for i in range(len(pfm['A']))]
            for b in 'ACGT'}

def scan_pwm(seq, pwm, threshold=0.0):
    L = len(pwm['A'])
    hits = []
    for i in range(len(seq) - L + 1):
        kmer = seq[i:i+L].upper()
        score = sum(pwm[b][j] for j, b in enumerate(kmer) if b in pwm)
        if score >= threshold:
            hits.append((i, kmer, score))
    return hits
```

---

## Code Templates

### Differential Expression (Python simplified / DESeq2 R)
```python
# Python: BH-corrected Mann-Whitney on normalized counts
def simple_de(count_df, ctrl_samples, treat_samples):
    sf = deseq2_size_factors(count_df)
    norm = count_df.div(sf, axis=1)
    results = []
    for gene in count_df.index:
        c, t = norm.loc[gene, ctrl_samples], norm.loc[gene, treat_samples]
        lfc = np.log2((t.mean()+1) / (c.mean()+1))
        _, p = stats.mannwhitneyu(c, t, alternative='two-sided')
        results.append({'gene': gene, 'log2FC': lfc, 'pvalue': p,
                        'baseMean': (c.mean()+t.mean())/2})
    df = pd.DataFrame(results).set_index('gene')
    # BH correction
    p = df['pvalue'].values; n = len(p); idx = np.argsort(p)
    adj = np.minimum.accumulate((p[idx] * n / (np.arange(n)+1))[::-1])[::-1]
    df['padj'] = 0.0; df.iloc[idx, df.columns.get_loc('padj')] = adj
    return df.sort_values('pvalue')
```

```r
# R: production DESeq2 workflow
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=count_matrix,
                              colData=sample_info, design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Treatment","Control"))
res_shrunk <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")
sig <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
```

```python
# Python: pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
dds = DeseqDataSet(counts=count_matrix, metadata=sample_info, design_factors="condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition","Treatment","Control"])
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

### Salmon (alignment-free)
```bash
salmon index -t transcriptome.fa -i salmon_index
salmon quant -i salmon_index -l A -1 R1.fastq -2 R2.fastq \
             -o sample_quant --validateMappings
```

### PERMANOVA (permutation test on distance matrix)
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
    p = (sum(f_stat(np.random.permutation(groups)) >= obs
             for _ in range(n_perm)) + 1) / (n_perm + 1)
    return obs, p
```

### TATA Box Detection
```python
import re
def find_tata_boxes(seq, strict=True):
    pat = 'TATAAA' if strict else r'TATA[AT]A[AT]'
    return [(m.start(), m.group()) for m in re.finditer(pat, seq.upper())]
```

---

## Common Pitfalls

- **Using RPKM/TPM for DE testing** — they introduce ratio-of-ratios artifacts. Always pass raw counts to DESeq2/edgeR.
- **Not accounting for library composition** — a single highly-expressed gene can suppress apparent expression of all others. DESeq2 median-of-ratios is robust to this; CPM/TPM are not.
- **Skipping dispersion shrinkage** — DESeq2 borrows information across genes to stabilize dispersion estimates for low-count genes. Manual NB tests miss this and have lower power.
- **Rarefying without checking curves** — rarefy only when rarefaction curves plateau; otherwise you discard real signal.
- **Comparing raw alpha diversity without rarefaction** — deeper-sequenced samples will always appear more diverse.
- **Interpreting PCoA axes without variance-explained** — axis 1 may explain only 20% of variation; check the percentages before interpreting clustering.
- **PERMANOVA sensitivity to dispersion** — PERMANOVA is sensitive to differences in within-group variance (not just centroids). Pair it with a dispersion test (betadisper in R).
- **CpG island detection with step=1 on large sequences** — use `step=10` or `step=50` for sequences >10 kb; then merge overlapping windows with a gap tolerance.
- **PWM scanning threshold** — setting threshold too low produces many false positives. Use the 80–90th percentile of scores from random sequence as a baseline.
- **Multiple testing in ORA** — apply BH correction across all tested pathways; Fisher p-values alone inflate false discovery.

---

---

## Related Skills
- `biostatistics-r` — hypothesis testing, BH correction, power analysis
- `numpy-pandas-wrangling` — count matrix manipulation, filtering, melting
- `python-advanced-sql` — database queries for annotation lookups (gene IDs, GO terms)


## Related Skills

- `promoter-regulatory-analysis` (this file)
- `rnaseq-metagenomics` (legacy merged skill)
