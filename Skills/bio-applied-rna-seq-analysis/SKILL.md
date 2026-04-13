---
name: bio-applied-rna-seq-analysis
description: "RNA-seq (RNA sequencing) is the dominant technology for studying gene expression at the transcriptome level. This notebook covers the complete RNA-seq analysis workflow, from experimental design throu"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/03_RNA_seq_Analysis/01_rna_seq_analysis.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# RNA-seq Analysis: From Reads to Differential Expression

*Source: Course notebook `Tier_3_Applied_Bioinformatics/03_RNA_seq_Analysis/01_rna_seq_analysis.ipynb`*

# RNA-seq Analysis: From Reads to Differential Expression

RNA-seq (RNA sequencing) is the dominant technology for studying gene expression at the transcriptome level. This notebook covers the complete RNA-seq analysis workflow, from experimental design through differential expression and visualization.

---

## Learning Objectives

By the end of this notebook, you will be able to:

1. Describe the RNA-seq experimental workflow and key design decisions
2. Understand alignment-based and alignment-free quantification strategies
3. Work with count matrices and apply normalization methods
4. Perform basic differential expression analysis
5. Interpret results using volcano plots, MA plots, and gene set enrichment
6. Implement a complete mini-workflow with simulated data

## 1. Why RNA-seq?

Before RNA-seq, gene expression was measured with **microarrays**, which rely on hybridization to pre-designed probes. RNA-seq offers several advantages:

| Feature | Microarrays | RNA-seq |
|---------|-------------|--------|
| Dynamic range | ~3 orders of magnitude | ~5+ orders of magnitude |
| Prior knowledge | Requires probe design | No reference needed (de novo) |
| Novel transcripts | Cannot detect | Can discover new isoforms, fusions |
| Quantification | Relative (intensity) | Digital (read counts) |
| Background noise | Cross-hybridization | Low, especially with proper QC |
| Cost per sample | Low | Medium-high (decreasing) |

### Common Applications

- **Differential gene expression**: Compare conditions (e.g., tumor vs. normal)
- **Transcript isoform discovery**: Alternative splicing, novel transcripts
- **Allele-specific expression**: Which allele is more active?
- **Single-cell RNA-seq (scRNA-seq)**: Expression at cellular resolution
- **Spatial transcriptomics**: Gene expression in tissue context

## 2. Experimental Design

Good experimental design is critical. Mistakes here cannot be fixed computationally.

### Key Decisions

**Biological replicates**: At least 3 per condition, ideally 5+. More replicates improve statistical power far more than deeper sequencing of fewer samples.

**Sequencing depth**: For differential expression in well-annotated organisms, 10-30 million reads per sample is standard. For transcript discovery or lowly-expressed genes, 50-100M+ may be needed.

**Read length and type**:
- Single-end (SE) 50-75 bp: sufficient for gene-level DE
- Paired-end (PE) 100-150 bp: better for isoform analysis, novel transcript discovery

**Strandedness**: Strand-specific protocols (e.g., dUTP method) preserve which DNA strand was transcribed. This is now standard and helps resolve overlapping genes on opposite strands.

### Library Preparation

```
Total RNA
    |
    v
RNA selection (poly-A selection or rRNA depletion)
    |
    v
Fragmentation (chemical or enzymatic)
    |
    v
Reverse transcription -> cDNA
    |
    v
Adapter ligation
    |
    v
PCR amplification (limited cycles)
    |
    v
Sequencing (Illumina, etc.)
```

**poly-A selection** enriches for mRNA (coding genes); **rRNA depletion** retains non-coding RNA and is needed for degraded samples (e.g., FFPE tissue).

### Batch Effects

Process all samples together when possible. If batches are unavoidable, balance conditions across batches. Record batch information for statistical correction later.

## 3. RNA-seq Workflow Overview

```
FASTQ files (raw reads)
      |
      v
Quality Control (FastQC, Trim Galore, fastp)
      |
      v
  +---+---+
  |       |
  v       v
Alignment-based          Alignment-free
(STAR/HISAT2)            (Salmon/kallisto)
  |                           |
  v                           v
BAM files              Transcript quantifications
  |                           |
  v                           |
featureCounts/HTSeq           |
  |                           |
  +-------+-------------------+
          |
          v
   Count matrix (genes x samples)
          |
          v
   Normalization
          |
          v
   Differential expression (DESeq2 / edgeR)
          |
          v
   Visualization & interpretation
          |
          v
   Gene set enrichment (GSEA / ORA)
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from collections import OrderedDict

np.random.seed(42)

# Plotting defaults
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11
```

## 5. Alignment-Based Quantification

### Step 1: Alignment

Splice-aware aligners handle reads that span exon-exon junctions.

| Aligner | Speed | Memory | Notes |
|---------|-------|--------|-------|
| **STAR** | Very fast | 30+ GB RAM | Most widely used for human/mouse |
| **HISAT2** | Fast | ~8 GB RAM | Successor to TopHat2, lower memory |

```bash
# STAR: build genome index (once)
STAR --runMode genomeGenerate \
     --genomeDir star_index/ \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile genes.gtf

# STAR: align reads
STAR --genomeDir star_index/ \
     --readFilesIn sample_R1.fastq sample_R2.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts
```

### Step 2: Counting reads per gene

After alignment, we count how many reads map to each gene using the gene annotation (GTF file).

| Tool | Approach | Speed |
|------|----------|-------|
| **featureCounts** (Subread) | Assigns reads to features by overlap | Very fast |
| **HTSeq-count** | Assigns reads, strict handling of ambiguity | Slower, more conservative |

```bash
# featureCounts example
featureCounts -a genes.gtf \
              -o counts.txt \
              -T 4 \
              -p --countReadPairs \
              sample1.bam sample2.bam sample3.bam
```

The result is a **count matrix**: rows are genes, columns are samples, values are integer read counts.

## 6. Alignment-Free Quantification

Alignment-free methods skip the computationally expensive genome alignment step. Instead, they quantify transcript abundance directly from reads using lightweight algorithms (pseudo-alignment or quasi-mapping).

| Tool | Algorithm | Speed | Output |
|------|-----------|-------|--------|
| **Salmon** | Quasi-mapping + EM | ~minutes per sample | TPM per transcript |
| **kallisto** | Pseudo-alignment + EM | ~minutes per sample | TPM per transcript |

```bash
# Salmon: build transcriptome index
salmon index -t transcriptome.fa -i salmon_index

# Salmon: quantify
salmon quant -i salmon_index \
             -l A \
             -1 sample_R1.fastq \
             -2 sample_R2.fastq \
             -o sample_quant \
             --validateMappings
```

### Key Advantage: Speed and Accuracy

Salmon/kallisto are **10-100x faster** than alignment-based methods while providing comparable or better accuracy for transcript-level quantification.

### From Transcripts to Genes: tximport

Salmon and kallisto output transcript-level estimates. The R package `tximport` aggregates these to gene level for use with DESeq2/edgeR, properly handling the uncertainty in transcript assignment. In Python, this can be done with `pydeseq2` or manually summing counts per gene.

## 7. Expression Units: RPKM, FPKM, TPM

Raw counts are not directly comparable between genes (different lengths) or between samples (different sequencing depths). Several normalization units exist:

### RPKM (Reads Per Kilobase per Million mapped reads)

$$\text{RPKM}_i = \frac{C_i}{L_i \times N} \times 10^9$$

where $C_i$ = read count for gene $i$, $L_i$ = gene length in bp, $N$ = total mapped reads.

### FPKM (Fragments Per Kilobase per Million)

Same as RPKM but for paired-end data: counts *fragments* (read pairs) instead of individual reads.

### TPM (Transcripts Per Million)

$$\text{TPM}_i = \frac{C_i / L_i}{\sum_j C_j / L_j} \times 10^6$$

TPM normalizes by gene length first, then by total. **TPM values sum to 1 million in every sample**, making samples directly comparable. TPM is now the preferred unit.

### Why Not Use These for DE?

RPKM/FPKM/TPM are useful for visualization and reporting, but **not for differential expression testing**. DE tools like DESeq2 and edgeR work with raw counts and apply their own normalization that accounts for library composition bias.

### Practical: Computing TPM from Counts

```python
# Example: compute TPM from raw counts
gene_names = ['BRCA1', 'TP53', 'MYC', 'GAPDH', 'ACTB', 'IL6', 'TNF', 'EGFR']
gene_lengths = np.array([7088, 2512, 4149, 1310, 1761, 1125, 1569, 5616])  # in bp
raw_counts = np.array([320, 890, 1560, 15200, 12400, 45, 78, 520])

def counts_to_tpm(counts, lengths):
    """Convert raw counts to TPM."""
    rate = counts / lengths  # normalize by gene length
    return rate / rate.sum() * 1e6  # scale to million

tpm = counts_to_tpm(raw_counts, gene_lengths)

tpm_df = pd.DataFrame({
    'Gene': gene_names,
    'Length (bp)': gene_lengths,
    'Raw Counts': raw_counts,
    'TPM': np.round(tpm, 1)
}).set_index('Gene')

print(tpm_df)
print(f"\nTPM sum: {tpm.sum():.0f} (should be 1,000,000)")
```

## 8. The Count Matrix

The central data structure in RNA-seq analysis is the **count matrix** (or expression matrix):

- Rows = genes (or transcripts)
- Columns = samples
- Values = integer read counts

Let's simulate a realistic count matrix for our downstream analysis.

```python
def simulate_rnaseq_counts(n_genes=2000, n_samples_per_group=4, n_de_genes=200,
                           base_mean_range=(10, 5000), fold_change_range=(1.5, 4.0)):
    """
    Simulate RNA-seq count data with known differentially expressed genes.
    
    Uses a negative binomial distribution, which models the overdispersion
    typically observed in RNA-seq count data.
    """
    n_samples = n_samples_per_group * 2
    
    # Gene names
    gene_names = [f'Gene_{i:04d}' for i in range(n_genes)]
    
    # Sample names and condition labels
    sample_names = ([f'Control_{i+1}' for i in range(n_samples_per_group)] +
                    [f'Treatment_{i+1}' for i in range(n_samples_per_group)])
    conditions = (['Control'] * n_samples_per_group +
                  ['Treatment'] * n_samples_per_group)
    
    # Base mean expression for each gene (log-normal distribution)
    base_means = np.exp(np.random.uniform(
        np.log(base_mean_range[0]),
        np.log(base_mean_range[1]),
        n_genes
    ))
    
    # Dispersion parameter (inversely related to mean, as in real RNA-seq)
    dispersions = 0.1 + 1.0 / np.sqrt(base_means)
    
    # Generate counts using negative binomial
    counts = np.zeros((n_genes, n_samples), dtype=int)
    true_de = np.zeros(n_genes, dtype=bool)
    true_log2fc = np.zeros(n_genes)
    
    # Select DE genes
    de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
    true_de[de_indices] = True
    
    # Assign fold changes (half up, half down)
    fold_changes = np.ones(n_genes)
    for idx in de_indices[:n_de_genes // 2]:
        fc = np.random.uniform(fold_change_range[0], fold_change_range[1])
        fold_changes[idx] = fc
        true_log2fc[idx] = np.log2(fc)
    for idx in de_indices[n_de_genes // 2:]:
        fc = np.random.uniform(fold_change_range[0], fold_change_range[1])
        fold_changes[idx] = 1.0 / fc
        true_log2fc[idx] = -np.log2(fc)
    
    for j in range(n_samples):
        # Library size factor (slight variation between samples)
        size_factor = np.random.uniform(0.8, 1.2)
        
        for i in range(n_genes):
            mean_expr = base_means[i] * size_factor
            if conditions[j] == 'Treatment':
                mean_expr *= fold_changes[i]
            
            # Negative binomial parameterization
            p = dispersions[i] / (dispersions[i] + mean_expr)
            n_param = dispersions[i] / p if p > 0 else 1
            counts[i, j] = np.random.negative_binomial(
                max(1, int(n_param)), max(1e-10, min(1 - 1e-10, p))
            )
    
    count_df = pd.DataFrame(counts, index=gene_names, columns=sample_names)
    sample_info = pd.DataFrame({'condition': conditions}, index=sample_names)
    gene_info = pd.DataFrame({
        'true_de': true_de,
        'true_log2fc': true_log2fc
    }, index=gene_names)
    
    return count_df, sample_info, gene_info

# Generate simulated data
counts_df, sample_info, gene_info = simulate_rnaseq_counts()

print(f"Count matrix shape: {counts_df.shape} (genes x samples)")
print(f"\nSample info:")
print(sample_info)
print(f"\nDE genes in simulation: {gene_info['true_de'].sum()} out of {len(gene_info)}")
print(f"\nFirst few rows of count matrix:")
counts_df.head(10)
```
