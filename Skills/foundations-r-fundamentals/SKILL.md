---
name: foundations-r-fundamentals
description: R quick-reference for bioinformatics — syntax differences from Python, vectors, data frames, Bioconductor patterns, and statistical distributions
tool_type: python
primary_tool: Python
---

# R Fundamentals for Bioinformatics

## Python vs R Differences

| Feature | Python | R |
|---------|--------|---|
| Indexing | Starts at 0 | Starts at **1** |
| Assignment | `=` | `<-` (preferred) |
| Negative index | `x[-1]` = last element | `x[-1]` = **all except first** |
| Boolean values | `True` / `False` | `TRUE` / `FALSE` |
| Missing data | `None` | `NA` |
| Null | `None` | `NULL` |
| Auto-print | `print(x)` | Just type `x` |
| AND/OR | `&`, `\|` (vectorized) | `&`, `\|` (vectorized); `&&`, `\|\|` (scalar) |

## Core Syntax Patterns

```r
# Vectors
gene_expression <- c(2.5, 3.1, 4.2, 1.8, 5.6)
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS")
positions <- 1:10
coverage <- seq(from=0, to=100, by=10)
groups <- rep(c("control", "treatment"), each=3)

# Named vector
gc_content <- c(BRCA1=0.42, TP53=0.38, EGFR=0.55, MYC=0.61)
gc_rich <- gc_content[gc_content > 0.5]  # logical indexing

# Vectorized ops — no explicit loops needed
log2fc <- log2(counts + 1)              # pseudocount avoids log(0)
sig_up <- log2fc > 1 & pvalue < 0.05   # vectorized filter
cat("Upregulated:", genes[sig_up], "\n")
which(sig_up | sig_down)                # indices of TRUE values
```

## Data Frames

```r
gene_data <- data.frame(
    gene = c("BRCA1", "TP53", "EGFR"),
    log2fc = c(1.2, -0.5, 3.8),
    pvalue = c(0.01, 0.15, 0.001),
    stringsAsFactors = FALSE   # keep strings as character, not factor
)

# Access
gene_data$gene
gene_data[1, ]              # first row
gene_data[2, "log2fc"]      # specific cell
gene_data[gene_data$log2fc > 1 & gene_data$pvalue < 0.05, ]  # filter

# Add columns
gene_data$padj <- p.adjust(gene_data$pvalue, method="BH")   # FDR
gene_data$bonf <- p.adjust(gene_data$pvalue, method="bonferroni")

# Sort
gene_data[order(gene_data$pvalue), ]
```

## Matrices (expression matrices)

```r
expr_matrix <- matrix(
    c(5.2, 3.1, 8.5, 6.2, 4.8, 2.9, 7.1, 5.8),
    nrow=2, byrow=TRUE,
    dimnames=list(c("Sample1", "Sample2"), c("BRCA1", "TP53", "EGFR", "MYC"))
)

expr_matrix[2, "EGFR"]           # specific cell
expr_matrix[1, ]                 # full row
expr_matrix[, "BRCA1"]          # full column
apply(expr_matrix, 1, mean)     # row means (margin=1)
apply(expr_matrix, 2, sd)       # col sds (margin=2)
```

## Statistical Distributions

Consistent naming convention: `d` (density), `p` (CDF), `q` (quantile), `r` (random).

```r
# Normal
rnorm(n=100, mean=10, sd=2)     # random samples
pnorm(q=1.96, lower.tail=TRUE)  # CDF → 0.975
qnorm(p=0.975)                  # quantile → 1.96

# Negative binomial (RNA-seq count data)
rnbinom(n=100, mu=50, size=5)

# Chi-squared test (HWE, independence)
chisq.test(observed, p=expected_freqs)

set.seed(42)   # reproducibility
```

## Bioconductor Quick Patterns

```r
# Install
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")

# DESeq2 minimal workflow
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "treatment", "control"))
res_df <- as.data.frame(res)
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# Reading/writing data
df <- read.csv("data.csv", stringsAsFactors=FALSE)
df <- read.delim("data.tsv", sep="\t")
write.csv(results, "output.csv", row.names=FALSE)
```

## Pitfalls

- **R indexing starts at 1**: `x[1]` is first; `x[-1]` means "all except first" (NOT last element like Python)
- **Factors look like strings but are not**: factor stores integer codes + labels; functions expecting strings will misbehave — use `stringsAsFactors=FALSE` in `data.frame()` or `as.character()` to convert
- **`library()` vs `require()`**: use `library()` in scripts — it errors loudly if missing; `require()` returns `FALSE` silently
- **Data frames vs matrices**: Bioconductor functions often require numeric matrices; a data frame with character columns cannot be used directly — convert with `as.matrix(df[, numeric_cols])`
- **`<-` inside function arguments**: inside `f(x <- 1)`, `<-` assigns in the enclosing scope AND passes the value; use `=` for keyword arguments to avoid side effects
- **Multiple testing**: `p.adjust(..., method="BH")` (Benjamini-Hochberg FDR) is standard for genomics; Bonferroni is overly conservative for large gene sets
- **Scalar vs vectorized AND/OR**: `&&` and `\|\|` only look at first element — always use `&` and `\|` for filtering vectors
