---
name: bio-applied-dmr-analysis
description: Differentially Methylated Regions (DMRs)
tool_type: python
primary_tool: Python
---

# Differentially Methylated Regions (DMRs)

- [methylKit vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/)
- [DSS (Dispersion Shrinkage for Sequencing)](https://bioconductor.org/packages/release/bioc/html/DSS.html)
- [DMRfinder documentation](https://github.com/jmschrei/dmrfinder)

## Genome-Wide DMR Calling and Multiple Testing

When running BSmooth, DSS, or dmrseq on whole-genome data, the pipeline extends to:

### DSS workflow (R)

```r
library(DSS)
# Create BSseq objects from count data
bs_ctrl  <- makeBSseqData(list(ctrl1, ctrl2, ctrl3),  sampleNames=c("C1","C2","C3"))
bs_treat <- makeBSseqData(list(trt1,  trt2,  trt3),   sampleNames=c("T1","T2","T3"))

# Perform DML (differentially methylated loci) test
dml_test  <- DMLtest(bs_ctrl, bs_treat, smoothing=TRUE, smoothing.span=500)

# Call DMRs from DML results
dmrs <- callDMR(dml_test, p.threshold=0.001, delta=0.1, minlen=50, minCG=3)
```python

**Key parameters:**
- `smoothing.span`: bandwidth in bp for local smoothing (200–500 bp is typical)
- `p.threshold`: per-CpG Wald test p-value threshold for DMR seeding
- `delta`: minimum mean methylation difference (avoid calling DMRs with tiny effect sizes)
- `minlen`: minimum DMR length in bp
- `minCG`: minimum number of CpGs within a DMR

### Statistical model in DSS

DSS fits a beta-binomial model at each CpG:

$$Y_i \sim \text{Binomial}(n_i, p_i), \quad \text{logit}(p_i) = \mu + \epsilon_i$$

where $\epsilon_i$ captures biological overdispersion (variance > what pure binomial would predict). This overdispersion parameter is estimated by moment matching across all samples. The Wald statistic tests $H_0: \mu_{\text{ctrl}} = \mu_{\text{treat}}$.

### Multiple testing

For genome-wide DMR calling, individual p-values are not directly interpretable (~25M tests). Instead:
- **Area statistic**: the sum of test statistics across a DMR is used as the primary score (not a p-value)
- **Permutation testing**: condition labels are permuted to generate a null distribution of DMR statistics
- **FDR estimation**: by comparing observed DMR statistics to the permutation null

## DMR Annotation to Genomic Features

Calling DMRs is only the first step. To understand their functional significance, each DMR must be annotated against:

1. **Gene features**: promoters (e.g., TSS ± 2 kb), 5'UTR, exons, introns, 3'UTR, intergenic
2. **CpG island features**: CGI, shore (±2 kb), shelf (±4 kb), open sea (everything else)
3. **Regulatory elements**: ENCODE enhancers, CTCF sites, accessible chromatin (ATAC/DNase peaks)

Tools for annotation:
- **genomation** (R): `annotateWithGeneParts()`, `annotateWithFeatureFlank()`
- **annotatr** (R): flexible annotation using UCSC Genome Browser tracks
- **ChIPseeker** (R): pie charts and heatmaps of genomic feature distribution
- **Python**: `pybedtools` for overlap operations; `pyranges` for fast genomic arithmetic

### Why do hypermethylated DMRs at promoters silence genes?

Promoter methylation at CpG islands blocks transcription factor binding and recruits methyl-CpG-binding domain (MBD) proteins (MBD1, MeCP2), which in turn recruit histone deacetylases (HDACs). This creates a self-reinforcing silencing loop: methylation → deacetylation → compact chromatin → no transcription.

## Integrating DMRs with Gene Expression

One of the most powerful validations for DMR findings is to ask: do hypermethylated promoter DMRs correlate with reduced gene expression in matched RNA-seq data? The classic "methylation–expression anti-correlation" is the hallmark of epigenetic silencing.

### Workflow

1. Identify DMRs overlapping gene promoters (TSS ± 2 kb)
2. Extract the mean methylation change (Δ beta) at the promoter
3. Correlate with the RNA-seq log₂ fold-change for the same gene
4. Expect: negative correlation — genes with hypermethylated promoters are downregulated

### Caveats

- Not all promoter methylation changes affect expression — TF binding and chromatin state matter
- Gene body methylation is actually **positively** correlated with expression (counterintuitive)
- Some genes are silenced by other mechanisms (H3K27me3, mutations) and methylation is secondary
- The correlation weakens at low-methylation promoters (floor effects)

```python
np.random.seed(99)

# Simulate 80 genes with promoter DMRs and matched RNA-seq fold changes
n_genes = 80

# Promoter methylation changes (Δ beta): mix of hyper and hypo
promo_delta = np.concatenate([
    np.random.uniform(0.10, 0.65, 50),    # hypermethylated promoters
    -np.random.uniform(0.05, 0.40, 30),   # hypomethylated
])

# RNA-seq log2FC: correlated with -delta_beta (methylation silences expression)
# Add noise to reflect incomplete correlation
rna_log2fc = -2.5 * promo_delta + np.random.normal(0, 0.6, n_genes)

gene_labels = [f'Gene{i+1}' for i in range(n_genes)]
is_hyper = promo_delta > 0.25    # "strongly silenced" genes

integ_df = pd.DataFrame({
    'gene': gene_labels,
    'promo_delta_beta': promo_delta,
    'rna_log2fc': rna_log2fc,
    'strongly_silenced': is_hyper
})

# Pearson correlation
r, p_val = stats.pearsonr(promo_delta, rna_log2fc)

fig, ax = plt.subplots(figsize=(7, 5))

# Non-highlighted genes
mask_ns = ~is_hyper
ax.scatter(integ_df.loc[mask_ns, 'promo_delta_beta'],
           integ_df.loc[mask_ns, 'rna_log2fc'],
           c='#78909C', alpha=0.6, s=35, label='Other genes')

# Strongly silenced genes
ax.scatter(integ_df.loc[is_hyper, 'promo_delta_beta'],
           integ_df.loc[is_hyper, 'rna_log2fc'],
           c='#E53935', alpha=0.85, s=50, zorder=5,
           label='Hypermethylated (Δβ > 0.25)')

# Regression line
x_range = np.linspace(promo_delta.min(), promo_delta.max(), 100)
slope, intercept, _, _, _ = stats.linregress(promo_delta, rna_log2fc)
ax.plot(x_range, slope * x_range + intercept, 'k--', lw=1.5)

# Reference lines
ax.axhline(0, color='gray', lw=0.8, linestyle=':')
ax.axvline(0, color='gray', lw=0.8, linestyle=':')

ax.set_xlabel('Promoter Δ Beta (methylation change)', fontsize=11)
ax.set_ylabel('RNA-seq log₂ Fold-Change', fontsize=11)
ax.set_title(f'Promoter Methylation vs Gene Expression\nPearson r = {r:.3f}, p = {p_val:.2e}', fontsize=11)
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig('dmr_expression_integration.png', dpi=120, bbox_inches='tight')
plt.show()

print(f"Pearson r = {r:.3f}  (p = {p_val:.2e})")
print(f"Strongly silenced genes (Δβ > 0.25 AND log2FC < -1): "
      f"{((is_hyper) & (rna_log2fc < -1)).sum()}")
```python

## Methylation Heatmap Across Samples

A common visualization in DMR papers is a **heatmap of beta values** at the top DMRs across all samples. This shows: (a) which samples cluster together, (b) whether the methylation difference is consistent across replicates, and (c) which DMRs have the most distinct patterns.

```python
np.random.seed(5)

# Build a beta matrix: 50 top DMRs  6 samples (3 ctrl, 3 treat)
n_dmrs_heat = 50
n_samples   = 6
sample_names = [f'Ctrl_{i+1}' for i in range(3)] + [f'Treat_{i+1}' for i in range(3)]

# DMR base methylation
dmr_base = np.random.beta(2, 2, n_dmrs_heat)  # random baseline

# Build matrix: control gets baseline treatment gets shifted betas
ctrl_matrix  = dmr_base[:, None] + np.random.normal(0, 0.05, (n_dmrs_heat, 3))
delta_each   = np.random.uniform(-0.6, 0.6, n_dmrs_heat)   # per-DMR delta
treat_matrix = dmr_base[:, None] + delta_each[:, None] + np.random.normal(0, 0.05, (n_dmrs_heat, 3))
beta_matrix  = np.clip(np.hstack([ctrl_matrix, treat_matrix]), 0, 1)

# Sort DMRs by delta (hyper at top, hypo at bottom)
order = np.argsort(-delta_each)
beta_matrix_sorted = beta_matrix[order, :]

fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(beta_matrix_sorted, aspect='auto', cmap='RdYlBu_r',
               vmin=0, vmax=1, interpolation='nearest')

ax.set_xticks(range(n_samples))
ax.set_xticklabels(sample_names, rotation=45, ha='right', fontsize=9)
ax.set_yticks(range(0, n_dmrs_heat, 10))
ax.set_yticklabels([f'DMR {i+1}' for i in range(0, n_dmrs_heat, 10)], fontsize=8)
ax.set_title('Methylation Heatmap — Top 50 DMRs\n(sorted by Δ beta, hyper→hypo)', fontsize=11)

cbar = plt.colorbar(im, ax=ax, shrink=0.6)
cbar.set_label('Beta value', fontsize=9)

# Vertical separator between ctrl and treat
ax.axvline(2.5, color='white', lw=2.5)
ax.text(0.75, -2.5, 'Control', ha='center', va='top', fontsize=9,
        transform=ax.transData, color='#1565C0', fontweight='bold')
ax.text(3.75, -2.5, 'Treatment', ha='center', va='top', fontsize=9,
        transform=ax.transData, color='#B71C1C', fontweight='bold')

plt.tight_layout()
plt.savefig('dmr_heatmap.png', dpi=120, bbox_inches='tight')
plt.show()
```python

## Summary and Key Takeaways

1. **Why DMRs over DMPs**: Neighbouring CpGs are correlated; regional testing increases power and biological interpretability. The multiple-testing burden (~25M CpGs) makes genome-wide single-CpG FDR correction very stringent.

2. **Smoothing**: BSmooth and DSS use local averaging to reduce sampling noise at individual CpGs by borrowing information from neighbours. This is especially important at 10–30× coverage.

3. **DMR calling**: Runs of CpGs exceeding a t-statistic threshold (BSmooth) or Wald test threshold (DSS) with ≥3–5 CpGs and ≥50 bp are called as DMRs. The sum of test statistics (area statistic) is more powerful than the peak statistic.

4. **Annotation**: DMRs are annotated against CGI features, gene features, and regulatory elements. Hypermethylated promoter CGI DMRs are the most functionally interpretable (silencing of tumor suppressors in cancer).

5. **Expression integration**: Promoter Δ beta correlates negatively with RNA-seq log₂FC. This anti-correlation is the key evidence that methylation is functionally silencing gene expression rather than being a bystander mark.

**Next**: Notebook 3 covers **epigenetic clocks** — how methylation at specific CpG sites predicts biological age.

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
