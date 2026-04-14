---
name: rnaseq-analysis
description: RNA-seq differential expression analysis and normalization workflows.
tool_type: mixed
primary_tool: Python
---

# RNA-seq, Metagenomics & Regulatory Analysis

## When to Use
- Differential gene expression from RNA-seq count data
- Comparing expression normalization strategies (RPKM/TPM/DESeq2)
- 16S rRNA amplicon microbiome studies (OTU/ASV, alpha/beta diversity)
- Promoter scanning: TATA box detection, CpG island calling, PWM/TFBS scanning

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
| RPKM | `(C / L) / N × 10⁹` | Single-end, within-sample |
| FPKM | Same as RPKM | Paired-end |
| TPM | `(C/L) / Σ(Cⱼ/Lⱼ) × 10⁶` | Cross-sample (sums to 1M) |
| DESeq2 size factor | `median(sample_counts / geo_mean_across_samples)` | DE testing |

C = read count, L = gene length (bp), N = total mapped reads. **Use TPM for reporting; raw counts + DESeq2 for DE testing.**

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

### DESeq2 Size Factors (Python)
```python
def deseq2_size_factors(count_matrix):
    nonzero_mask = (count_matrix > 0).all(axis=1)
    filtered = count_matrix.loc[nonzero_mask]
    geo_means = np.exp(np.log(filtered).mean(axis=1))
    return pd.Series({s: np.median(filtered[s] / geo_means) for s in filtered.columns})

normalized = counts_df.div(deseq2_size_factors(counts_df), axis=1)
```

### TPM from Raw Counts
```python
def counts_to_tpm(counts, lengths):
    rate = counts / lengths
    return rate / rate.sum() * 1e6
```

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
    L = len(pwm['A'])
    return [(i, seq[i:i+L].upper(), sum(pwm[b][j] for j, b in enumerate(seq[i:i+L].upper()) if b in pwm))
            for i in range(len(seq) - L + 1)
            if sum(pwm[b][j] for j, b in enumerate(seq[i:i+L].upper()) if b in pwm) >= threshold]
```

## Code Templates

### Differential Expression (Python / R / pydeseq2)
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
        results.append({'gene': gene, 'log2FC': lfc, 'pvalue': p, 'baseMean': (c.mean()+t.mean())/2})
    df = pd.DataFrame(results).set_index('gene')
    p = df['pvalue'].values; n = len(p); idx = np.argsort(p)
    adj = np.minimum.accumulate((p[idx] * n / (np.arange(n)+1))[::-1])[::-1]
    df['padj'] = 0.0; df.iloc[idx, df.columns.get_loc('padj')] = adj
    return df.sort_values('pvalue')
```

```r
# R: production DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_info, design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Treatment","Control"))
res_shrunk <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")
sig <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
```

```python
# pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
dds = DeseqDataSet(counts=count_matrix, metadata=sample_info, design_factors="condition")
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["condition","Treatment","Control"])
stat_res.summary(); results_df = stat_res.results_df
```

### STAR + featureCounts
```bash
STAR --runMode genomeGenerate --genomeDir star_index/ --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
STAR --genomeDir star_index/ --readFilesIn R1.fastq R2.fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
featureCounts -a genes.gtf -o counts.txt -T 4 -p --countReadPairs *.bam
```

### Salmon (alignment-free)
```bash
salmon index -t transcriptome.fa -i salmon_index
salmon quant -i salmon_index -l A -1 R1.fastq -2 R2.fastq -o sample_quant --validateMappings
```

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

- **RPKM/TPM for DE testing** — ratio-of-ratios artifacts; always pass raw counts to DESeq2/edgeR.
- **Library composition** — one highly-expressed gene suppresses apparent expression of all others; DESeq2 median-of-ratios is robust; CPM/TPM are not.
- **Skipping dispersion shrinkage** — DESeq2 borrows across genes to stabilize low-count estimates; manual NB tests miss this.
- **Rarefying without curves** — rarefy only when rarefaction curves plateau; raw alpha diversity comparisons require rarefaction (deeper samples appear more diverse).
- **PCoA axes without variance-explained** — axis 1 may explain only 20%; check percentages before interpreting clustering.
- **PERMANOVA + dispersion** — sensitive to within-group variance differences, not just centroids; pair with betadisper (R).
- **CpG step=1 on large sequences** — use `step=10–50` for >10 kb, then merge overlapping windows.
- **PWM threshold** — use 80–90th percentile of random-sequence scores as baseline; low thresholds inflate FP.
- **Multiple testing in ORA** — BH-correct across all tested pathways.

## Related Skills
- `biostatistics-r` — hypothesis testing, BH correction, power analysis
- `numpy-pandas-wrangling` — count matrix manipulation
- `python-advanced-sql` — annotation lookups (gene IDs, GO terms)
