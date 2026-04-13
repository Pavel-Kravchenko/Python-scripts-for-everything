---
name: bio-applied-isoform-analysis
description: "**Tier 3 — Applied Bioinformatics | Module 25 · Notebook 3**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/25_Long_Read_Sequencing/03_isoform_analysis.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Isoform Analysis with Long Reads

*Source: Course notebook `Tier_3_Applied_Bioinformatics/25_Long_Read_Sequencing/03_isoform_analysis.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 25 · Notebook 3**

*Prerequisites: Notebook 1 (ONT Processing), Module 03 (RNA-seq Analysis)*

---

**By the end of this notebook you will be able to:**
1. Explain why full-length reads resolve isoforms that short-read RNA-seq cannot
2. Process ONT direct-RNA and cDNA reads for isoform analysis with FLAMES
3. Align full-length transcripts to a reference genome with Minimap2
4. Identify novel isoforms and quantify transcript-level abundance
5. Perform differential isoform usage analysis between conditions



**Key resources:**
- [FLAMES documentation](https://github.com/LuyiTian/FLAMES)
- [Iso-Seq analysis with pb-CCS/IsoSeq3](https://github.com/PacificBiosciences/IsoSeq)
- [bambu — isoform discovery for long reads](https://github.com/GoekeLab/bambu)

## 1. Long-Read Transcriptomics Overview

Short-read RNA-seq (Illumina, 75–150 bp) measures gene-level expression accurately but struggles with multi-exon transcripts because a single read rarely spans more than 1–2 exons. Isoform assembly from short reads (StringTie, STAR + Cufflinks) is **inference-based** and prone to chimeric or truncated reconstructions at complex loci.

Long-read sequencing produces **full-length cDNA reads** (500 bp – 15 kb typical) that span entire transcripts, directly resolving all splice junctions, transcription start sites, and poly-A sites in a single read.

### Technology Comparison

| Feature | Short-read (Illumina) | ONT cDNA/direct-RNA | PacBio IsoSeq (HiFi) |
|---|---|---|---|
| Read length | 75–300 bp | 1–20 kb | 1–30 kb (CCS) |
| Per-read error rate | ~0.1% | ~1–3% (R10.4.1) | ~0.1% (HiFi) |
| Isoform resolution | Inference required | Direct | Direct |
| Modification detection | No | Yes (direct-RNA only) | No |
| Throughput (Gb) | Very high | High | Medium |

### ONT Protocols for Transcriptomics

- **Direct-RNA sequencing (dRNA)**: Sequences native RNA molecules without reverse transcription. Preserves base modifications (m⁶A, pseudouridine). Median read length ~1–3 kb. Lower throughput but uniquely suited to modification studies.
- **cDNA sequencing**: Reverse-transcribed and amplified cDNA. Longer reads and higher throughput. Loses modification information. Uses the `splice` preset in Minimap2.
- **PCR-cDNA**: Additional PCR amplification; required for very low-input samples but introduces amplification bias that distorts isoform ratios.

### PacBio IsoSeq

IsoSeq uses PacBio **HiFi reads** (Circular Consensus Sequencing, CCS). After multiple passes of the insert, the consensus reaches ~99.9% per-base accuracy. The IsoSeq3 pipeline: CCS → primer removal (`lima`) → full-length non-concatemer detection → clustering → polishing, producing high-confidence isoform sequences directly. No assembly step is needed.

## 2. Aligning Transcripts to Genome

For splice-aware alignment of long cDNA or dRNA reads, **Minimap2** is used with the `splice` preset (or `splice:hq` for PacBio HiFi). This preset enables junction-aware chaining and represents introns as CIGAR `N` operations in the output BAM/SAM.

**Key flags for long-read RNA alignment:**
- `-ax splice`: long-read spliced alignment preset
- `--secondary=no`: suppress secondary alignments (one alignment per transcript read)
- `-C5`: additional cost for non-canonical splice sites (GT-AG = 0, AT-AC / GC-AG = C5 penalty)
- `--cs`: output the cs tag encoding the alignment differences

```bash
# ONT cDNA reads
minimap2 -ax splice --secondary=no -C5 --cs \
    hg38.fa cdna_reads.fastq.gz \
    | samtools sort -o cdna_aligned.bam -@ 8
samtools index cdna_aligned.bam

# PacBio IsoSeq / HiFi cDNA (splice:hq for higher-accuracy reads)
minimap2 -ax splice:hq --secondary=no -C5 \
    hg38.fa isoseq_reads.fastq.gz \
    | samtools sort -o isoseq_aligned.bam -@ 8
samtools index isoseq_aligned.bam
```python

The cell below simulates alignment QC statistics and CIGAR `N`-operation analysis for a realistic ONT cDNA run.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate ONT cDNA alignment statistics
n_reads = 8000
read_lengths = rng.lognormal(mean=8.0, sigma=0.9, size=n_reads).astype(int)
n_exons_per_read = rng.integers(1, 16, size=n_reads)
mapped = rng.random(n_reads) < 0.96
read_lengths_mapped = read_lengths[mapped]
exon_counts = n_exons_per_read[mapped]
n_introns_total = int(exon_counts.sum() - mapped.sum())
intron_lengths = rng.lognormal(mean=7.5, sigma=1.3, size=n_introns_total).astype(int)
intron_lengths = np.clip(intron_lengths, 60, 500_000)

print('=== Minimap2 cDNA alignment summary ===')
print(f'Total reads:            {n_reads:>8,}')
print(f'Mapped reads:           {mapped.sum():>8,}  ({mapped.mean():.1%})')
print(f'Median read length:     {int(np.median(read_lengths_mapped)):>8,} bp')
print(f'Multi-exon reads:       {(exon_counts > 1).sum():>8,}  ({(exon_counts > 1).mean():.1%})')
print(f'Total introns detected: {n_introns_total:>8,}')
print(f'Median intron length:   {int(np.median(intron_lengths)):>8,} bp')

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].hist(read_lengths_mapped / 1000, bins=50, color='steelblue', edgecolor='white', alpha=0.85)
axes[0].axvline(np.median(read_lengths_mapped) / 1000, color='red', linestyle='--',
                label=f'Median {np.median(read_lengths_mapped)/1000:.1f} kb')
axes[0].set_xlabel('Read length (kb)')
axes[0].set_ylabel('Count')
axes[0].set_title('cDNA read length distribution')
axes[0].legend()

axes[1].hist(exon_counts, bins=range(1, 17), color='darkorange', edgecolor='white', alpha=0.85)
axes[1].set_xlabel('Exons per read (N-count + 1)')
axes[1].set_ylabel('Count')
axes[1].set_title('Exons per mapped read')
axes[1].set_xticks(range(1, 17))

axes[2].hist(np.log10(intron_lengths + 1), bins=40, color='mediumseagreen', edgecolor='white', alpha=0.85)
axes[2].set_xlabel('log10(intron length + 1)')
axes[2].set_ylabel('Count')
axes[2].set_title('Intron length distribution (CIGAR N ops)')
axes[2].axvline(np.log10(100), color='red', linestyle=':', label='100 bp')
axes[2].axvline(np.log10(100_000), color='purple', linestyle=':', label='100 kb')
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.suptitle('Minimap2 splice alignment QC (simulated ONT cDNA)', y=1.02)
plt.show()
```python

## 3. Isoform Discovery and Quantification with bambu

**bambu** (Bioconductor) discovers and quantifies isoforms from long-read RNA-seq using a probabilistic model. Key concepts:
- Reads are assigned to isoform models based on their exon-chain compatibility
- **NDR (Novel Discovery Rate)**: controls false positive rate for novel isoforms; `NDR=0.1` requires 90% confidence that a transcript is genuinely novel
- Multi-sample mode ensures consistent isoform models across all samples (critical for differential analysis)

```R
library(bambu)

annotations <- prepareAnnotations('gencode.v44.annotation.gtf')

se <- bambu(
  reads       = c('ctrl_rep1.bam', 'ctrl_rep2.bam', 'treat_rep1.bam', 'treat_rep2.bam'),
  annotations = annotations,
  genome      = 'hg38.fa',
  NDR         = 0.1
)

writeBambuOutput(se, path = 'bambu_output/')
# counts_transcript.txt:       transcript-level counts (rows=transcripts, cols=samples)
# counts_gene.txt:             gene-level counts
# extended_annotations.gtf:   all transcripts including novel ones
```python

**Isoform categories:**
- `annotated`: exact match to a GENCODE transcript
- `novel_in_catalog`: new combination of known exons
- `novel_splice_site`: new 5' or 3' splice donor/acceptor within a known gene
- `novel_exon`: entirely new exon
- `intergenic`: in unannotated region (usually filter out)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)

N_GENES   = 300
iso_counts_arr = rng.integers(1, 8, size=N_GENES)
n_tx = iso_counts_arr.sum()

categories = rng.choice(
    ['annotated', 'novel_in_catalog', 'novel_splice_site', 'novel_exon', 'intergenic'],
    size=n_tx,
    p=[0.62, 0.18, 0.09, 0.07, 0.04]
)
gene_ids = np.repeat([f'GENE_{i:04d}' for i in range(N_GENES)], iso_counts_arr)
base_expr = rng.lognormal(mean=4, sigma=2, size=n_tx)

counts_df = pd.DataFrame({
    'tx_id':   [f'TX_{i:05d}' for i in range(n_tx)],
    'gene_id':  gene_ids,
    'category': categories,
    'Ctrl_1':   rng.poisson(base_expr).astype(int),
    'Ctrl_2':   rng.poisson(base_expr * rng.lognormal(0, 0.1, n_tx)).astype(int),
    'Treat_1':  rng.poisson(base_expr * rng.lognormal(0.3, 0.4, n_tx)).astype(int),
    'Treat_2':  rng.poisson(base_expr * rng.lognormal(0.3, 0.4, n_tx)).astype(int),
})

print(f'Total transcripts: {len(counts_df)}')
print(f'Genes:             {counts_df["gene_id"].nunique()}')
print('\nCategory breakdown:')
print(counts_df['category'].value_counts())

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
cat_counts = counts_df['category'].value_counts()
colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2']

axes[0].pie(cat_counts.values, labels=cat_counts.index, colors=colors,
            autopct='%1.1f%%', startangle=90)
axes[0].set_title('Isoform categories (bambu)')

iso_per_gene = counts_df.groupby('gene_id').size()
axes[1].hist(iso_per_gene, bins=range(1, 10), color='steelblue', edgecolor='white', alpha=0.85)
axes[1].set_xlabel('Isoforms per gene')
axes[1].set_ylabel('Number of genes')
axes[1].set_title('Isoform complexity per gene')
axes[1].set_xticks(range(1, 9))

ctrl_mean  = counts_df[['Ctrl_1', 'Ctrl_2']].mean(axis=1)
treat_mean = counts_df[['Treat_1', 'Treat_2']].mean(axis=1)
cat_color_map = dict(zip(cat_counts.index, colors))
pt_colors = counts_df['category'].map(cat_color_map)
axes[2].scatter(np.log2(ctrl_mean + 1), np.log2(treat_mean + 1),
                c=pt_colors, alpha=0.4, s=8)
lim = max(np.log2(ctrl_mean + 1).max(), np.log2(treat_mean + 1).max()) + 0.5
axes[2].plot([0, lim], [0, lim], 'k--', lw=1, alpha=0.5)
axes[2].set_xlabel('log2(Ctrl mean + 1)')
axes[2].set_ylabel('log2(Treat mean + 1)')
axes[2].set_title('Transcript-level expression')
plt.tight_layout()
plt.show()
```python

## 4. Differential Isoform Usage Analysis

Differential isoform usage (DIU) tests whether a gene changes the **relative proportion** of its isoforms between conditions — independent of total gene expression level.

### DRIMSeq: Dirichlet-Multinomial Proportions Test

DRIMSeq models isoform counts per gene as a **Dirichlet-multinomial distribution** (appropriate because isoform counts are compositional: they sum to the gene total). It estimates per-gene precision (overdispersion), then uses a likelihood ratio test.

```R
library(DRIMSeq)

counts_tx <- read.table('bambu_output/counts_transcript.txt', header=TRUE)
sample_info <- data.frame(
  sample_id = colnames(counts_tx)[-c(1,2)],
  condition = c('Control', 'Control', 'Treatment', 'Treatment')
)

d <- dmDSdata(counts = counts_tx, samples = sample_info)

# Filter: genes need >= 2 isoforms and adequate expression
d <- dmFilter(d,
  min_samps_gene_expr    = 2,
  min_samps_feature_expr = 2,
  min_gene_expr          = 10,
  min_feature_expr       = 5
)

d <- dmPrecision(d)
d <- dmFit(d)
d <- dmTest(d, coef = 'conditionTreatment')

res_gene <- results(d, level = 'gene')    # gene-level omnibus test
res_tx   <- results(d, level = 'feature') # per-transcript test
sig      <- res_gene[res_gene$adj_pvalue < 0.05, ]
```python

**Sashimi plots** visualize splice-junction read coverage as arcs between exons, with arc thickness proportional to junction read count. A changing arc thickness between conditions reveals isoform switching. Use `ggsashimi` (R) or IGV's built-in sashimi view.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

rng = np.random.default_rng(99)

n_genes = 200
pvals_raw = np.concatenate([
    rng.beta(1, 20, size=170),   # unchanged — p near 1
    rng.beta(1, 2,  size=30),    # DIU genes — enriched small p-values
])
pvals_raw = np.clip(pvals_raw + rng.exponential(0.001, n_genes), 1e-10, 1.0)

_, adj_pvals, _, _ = multipletests(pvals_raw, method='BH')

diu_df = pd.DataFrame({
    'gene':       [f'GENE_{i:04d}' for i in range(n_genes)],
    'pvalue':     pvals_raw,
    'adj_pval':   adj_pvals,
    'n_isoforms': rng.integers(2, 7, size=n_genes),
    'delta_prop': rng.uniform(-0.6, 0.6, size=n_genes),
}).sort_values('adj_pval')

sig_genes = diu_df[diu_df['adj_pval'] < 0.05]
print(f'Total genes tested:             {n_genes}')
print(f'Significant DIU (FDR < 5%):     {len(sig_genes)}')
print('\nTop 5 differentially used genes:')
print(sig_genes[['gene', 'pvalue', 'adj_pval', 'n_isoforms']].head().to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

sig_mask = diu_df['adj_pval'] < 0.05
axes[0].scatter(range(len(diu_df)), -np.log10(diu_df['pvalue']),
                c=np.where(sig_mask, '#C44E52', '#8172B2'), s=10, alpha=0.6)
axes[0].axhline(-np.log10(0.05), color='orange', linestyle=':', lw=1.2, label='p=0.05')
axes[0].set_xlabel('Gene rank')
axes[0].set_ylabel('-log10(p-value)')
axes[0].set_title(f'DRIMSeq DIU results ({len(sig_genes)} significant)')
axes[0].legend()

# Simulate isoform proportions for the top DIU gene
top = sig_genes.iloc[0]
n_iso = top['n_isoforms']
alpha_ctrl  = rng.exponential(2, size=n_iso) + 0.5
alpha_treat = alpha_ctrl.copy()
alpha_treat[0] *= 0.15   # dominant isoform drops
alpha_treat[1] *= 3.5    # another rises
ctrl_props  = rng.dirichlet(alpha_ctrl,  size=4).mean(axis=0)
treat_props = rng.dirichlet(alpha_treat, size=4).mean(axis=0)

iso_labels = [f'TX_{i+1}' for i in range(n_iso)]
x = np.arange(n_iso)
w = 0.35
axes[1].bar(x - w/2, ctrl_props,  w, label='Control',   color='#4C72B0', alpha=0.85)
axes[1].bar(x + w/2, treat_props, w, label='Treatment',  color='#DD8452', alpha=0.85)
axes[1].set_xticks(x); axes[1].set_xticklabels(iso_labels)
axes[1].set_ylabel('Isoform proportion')
axes[1].set_title(f'{top["gene"]}: isoform usage switch')
axes[1].legend(); axes[1].set_ylim(0, 1)

plt.tight_layout()
plt.show()
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
