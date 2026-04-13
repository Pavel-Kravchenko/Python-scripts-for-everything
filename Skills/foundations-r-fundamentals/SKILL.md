---
name: foundations-r-fundamentals
description: "While this course primarily uses Python, R is indispensable in bioinformatics. The **Bioconductor** ecosystem provides over 2,000 packages for genomics, transcriptomics, and proteomics. Tools like **D"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/05_R_Basics/01_r_fundamentals.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# R Fundamentals for Bioinformatics

*Source: Course notebook `Tier_0_Computational_Foundations/05_R_Basics/01_r_fundamentals.ipynb`*

# R Fundamentals for Bioinformatics

**Tier 0 -- Computational Foundations | Module 5**

---

## Why R for Bioinformatics?

While this course primarily uses Python, R is indispensable in bioinformatics. The **Bioconductor** ecosystem provides over 2,000 packages for genomics, transcriptomics, and proteomics. Tools like **DESeq2** and **edgeR** for differential expression, **Seurat** for single-cell analysis, and **ggplot2** for publication-quality figures are written in R.

This module gives you a working foundation in R so you can:
- Read and understand R code in published papers and pipelines
- Perform basic data manipulation and statistical analysis
- Use Bioconductor packages when Python alternatives are insufficient
- Create R-based visualizations

### How to Run This Notebook

This notebook uses the **R kernel** (IRkernel). You can run it in:
- **Jupyter** with IRkernel installed (`install.packages('IRkernel'); IRkernel::installspec()`)
- **RStudio** by copying code cells into an R script or R Markdown document
- **Google Colab** by changing runtime to R

---

### Learning Objectives

By the end of this module you will be able to:
- Create and manipulate vectors, matrices, and data frames
- Use vectorized operations (R's greatest strength)
- Filter and subset biological data
- Work with R's statistical distribution functions
- Create informative plots
- Read/write tabular data files
- Understand the basics of Bioconductor

## How to use this notebook
1. All code cells use the R kernel — make sure your Jupyter environment has IRkernel installed.
2. R's vectorized operations behave differently from Python loops. Pay attention to how operations apply automatically to every element of a vector.
3. When you see unfamiliar syntax, try `?function_name` in an R cell to read the documentation.
4. The bioinformatics examples throughout are intentionally realistic — the same patterns appear in real DESeq2 and Seurat workflows.

## Common stumbling points

- **R indexing starts at 1, not 0**: `x[1]` is the first element. `x[-1]` does NOT mean the last element — it means "everything except the first."
- **`<-` vs `=` for assignment**: Both work, but `<-` is the R convention. Inside function arguments, `=` has a different meaning (keyword argument), so mixing them up causes subtle bugs.
- **Factors look like strings but are not**: A factor stores an integer code plus a lookup table of labels. Functions that expect strings will misbehave if given a factor.
- **Data frames vs. matrices**: Many Bioconductor functions require specific input types. A data frame with character columns cannot be used as a numeric matrix — you must convert explicitly.
- **`library()` vs `require()`**: Use `library()` in scripts; it stops with an error if the package is missing. `require()` returns `FALSE` silently, hiding the problem.

---

## 1. Assignment, Variables, and the Vector Type

Everything in R is an object. The most fundamental object is the **vector** — even a single number is a vector of length 1. This is different from Python, where `x = 5` creates a scalar integer.

R uses `<-` for assignment (reading it as "gets") though `=` also works. The convention in R packages and published code is `<-`, so we follow that here.

### Key differences from Python

| Feature | Python | R |
|---------|--------|---|
| Indexing | Starts at 0 | Starts at **1** |
| Assignment | `=` | `<-` (preferred) or `=` |
| Negative index | `x[-1]` = last element | `x[-1]` = **all except first** |
| Boolean values | `True` / `False` | `TRUE` / `FALSE` |
| Missing data | `None` | `NA` |
| Auto-print | `print(x)` | Just type `x` |

```python
# Assignment
gene_name <- "TP53"
expression_level <- 5.7
is_oncogene <- FALSE

# Check types
cat("gene_name:", class(gene_name), "\n")
cat("expression_level:", class(expression_level), "\n")
cat("is_oncogene:", class(is_oncogene), "\n")
```

```python
# Creating vectors with c() -- "combine"
gene_expression <- c(2.5, 3.1, 4.2, 1.8, 5.6)
gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS")
is_oncogene <- c(FALSE, FALSE, TRUE, TRUE, TRUE)

# Vectors have a single type -- mixing types causes coercion
mixed <- c(1, "two", TRUE)
cat("Mixed vector:", mixed, "\n")  # Everything becomes character
cat("Type:", class(mixed), "\n")
```

```python
# Sequences -- extremely common in R
positions <- 1:10              # Integer sequence 1 to 10
cat("1:10 =", positions, "\n")

# seq() for more control
coverage <- seq(from = 0, to = 100, by = 10)
cat("seq(0, 100, 10) =", coverage, "\n")

# Evenly spaced values (useful for plotting)
gc_range <- seq(from = 0.3, to = 0.7, length.out = 5)
cat("GC content range:", gc_range, "\n")

# rep() to repeat values
groups <- rep(c("control", "treatment"), each = 3)
cat("Groups:", groups, "\n")
```

---

## 2. Vectorized Operations -- R's Superpower

In R, arithmetic and logical operations work element-by-element on entire vectors. This is both elegant and fast -- **you rarely need explicit loops in R**.

```python
x <- 1:5
y <- 6:10

cat("x =", x, "\n")
cat("y =", y, "\n\n")

# Element-wise arithmetic
cat("x + y  =", x + y, "\n")   # 7 9 11 13 15
cat("x * 2  =", x * 2, "\n")   # 2 4 6 8 10
cat("x * y  =", x * y, "\n")   # 6 14 24 36 50
cat("x^2    =", x^2, "\n")     # 1 4 9 16 25

# Common bioinformatics transformations
counts <- c(100, 500, 10, 2000, 50)
cat("\nRaw counts:  ", counts, "\n")
cat("Log2 counts: ", round(log2(counts), 2), "\n")
cat("Log2(x+1):   ", round(log2(counts + 1), 2), "\n")  # Pseudocount to avoid log(0)
```

```python
# Logical operations (the basis of all filtering)
expression <- c(2.5, 3.1, 8.2, 1.8, 5.6)

cat("expression > 3:  ", expression > 3, "\n")
cat("expression >= 5: ", expression >= 5, "\n")

# Combine conditions with & (AND) and | (OR)
in_range <- expression > 2 & expression < 6
cat("2 < expr < 6:    ", in_range, "\n")

# Summary functions on logical vectors
cat("\nHow many > 3?  ", sum(expression > 3), "\n")   # TRUE counts as 1
cat("Any > 10?      ", any(expression > 10), "\n")
cat("All > 0?       ", all(expression > 0), "\n")
```

```python
# Useful vector functions
x <- c(4, 1, 7, 2, 9, 3)

cat("length:", length(x), "\n")
cat("sum:   ", sum(x), "\n")
cat("mean:  ", mean(x), "\n")
cat("median:", median(x), "\n")
cat("sd:    ", sd(x), "\n")
cat("min:   ", min(x), "\n")
cat("max:   ", max(x), "\n")
cat("range: ", range(x), "\n")
cat("sort:  ", sort(x), "\n")
cat("order: ", order(x), "\n")  # Indices that would sort the vector
```

---

## 3. Indexing and Subsetting

R supports three powerful indexing styles:
1. **Positive integers**: Select specific positions
2. **Negative integers**: Exclude specific positions
3. **Logical vectors**: Select where `TRUE`

Remember: **R indexes from 1, not 0!**

```python
x <- c(10, 20, 30, 40, 50, 60)

# Positive indexing
cat("x[1]       =", x[1], "\n")         # First element: 10
cat("x[2:4]     =", x[2:4], "\n")       # Elements 2 to 4: 20 30 40
cat("x[c(1,3,5)]=", x[c(1, 3, 5)], "\n") # Specific positions: 10 30 50

# Negative indexing = EXCLUDE (very different from Python!)
cat("x[-1]      =", x[-1], "\n")        # All except first: 20 30 40 50 60
cat("x[-(1:3)]  =", x[-(1:3)], "\n")    # Exclude first 3: 40 50 60

# Logical indexing
cat("x[x > 25]  =", x[x > 25], "\n")    # Values > 25: 30 40 50 60
```

```python
# Named vectors -- very useful for biological data
gc_content <- c(BRCA1 = 0.42, TP53 = 0.38, EGFR = 0.55, MYC = 0.61)

cat("GC content of EGFR:", gc_content["EGFR"], "\n")
cat("Names:", names(gc_content), "\n")

# Which genes are GC-rich (> 50%)?
gc_rich <- gc_content[gc_content > 0.5]
cat("GC-rich genes:", names(gc_rich), "\n")
cat("Their GC content:", gc_rich, "\n")
```

### Bioinformatics Example: Filtering Differentially Expressed Genes

```python
# Simulated differential expression results
genes  <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PTEN", "RB1", "APC")
log2fc <- c(1.2, -0.5, 3.8, 2.1, -1.8, 0.3, -2.5, 0.1)
pvalue <- c(0.01, 0.15, 0.001, 0.005, 0.02, 0.45, 0.003, 0.72)

# Significantly upregulated: log2FC > 1 AND p < 0.05
sig_up <- log2fc > 1 & pvalue < 0.05
cat("Significantly upregulated:\n")
cat("  Genes:", genes[sig_up], "\n")
cat("  log2FC:", log2fc[sig_up], "\n")

# Significantly downregulated: log2FC < -1 AND p < 0.05
sig_down <- log2fc < -1 & pvalue < 0.05
cat("\nSignificantly downregulated:\n")
cat("  Genes:", genes[sig_down], "\n")
cat("  log2FC:", log2fc[sig_down], "\n")

# which() returns the indices of TRUE values
cat("\nIndices of significant genes:", which(sig_up | sig_down), "\n")
```

---

## 4. Matrices

Matrices are 2D arrays where all elements have the same type. They are used for expression matrices, distance matrices, and substitution matrices.

```python
# Create an expression matrix (genes x samples)
expr_matrix <- matrix(
    c(5.2, 3.1, 8.5, 6.2,
      4.8, 2.9, 7.1, 5.8,
      6.1, 3.5, 9.2, 7.0),
    nrow = 3, byrow = TRUE,
    dimnames = list(
        c("Sample1", "Sample2", "Sample3"),
        c("BRCA1", "TP53", "EGFR", "MYC")
    )
)

print(expr_matrix)
cat("\nDimensions:", dim(expr_matrix), "\n")
cat("Row names:", rownames(expr_matrix), "\n")
cat("Col names:", colnames(expr_matrix), "\n")
```

```python
# Matrix indexing: [row, column]
cat("EGFR in Sample2:", expr_matrix[2, "EGFR"], "\n")
cat("All of Sample1: ", expr_matrix[1, ], "\n")
cat("All BRCA1 values:", expr_matrix[, "BRCA1"], "\n")

# Apply functions across rows or columns
cat("\nMean per sample (rows): ", apply(expr_matrix, 1, mean), "\n")
cat("Mean per gene (columns):", apply(expr_matrix, 2, mean), "\n")
cat("SD per gene:            ", round(apply(expr_matrix, 2, sd), 2), "\n")
```

---

## 5. Data Frames -- Tabular Data

Data frames are R's equivalent of a spreadsheet or a Pandas DataFrame. Each column can have a different type. This is the primary data structure for biological datasets.

```python
# Create a data frame
gene_data <- data.frame(
    gene = c("BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PTEN"),
    expression = c(5.2, 3.1, 8.5, 6.2, 4.1, 2.8),
    log2fc = c(1.2, -0.5, 3.8, 2.1, -1.8, -2.1),
    pvalue = c(0.01, 0.15, 0.001, 0.005, 0.02, 0.008),
    is_oncogene = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE  # Keep strings as character, not factor
)

print(gene_data)
```

```python
# Inspecting data frames
cat("Dimensions:", dim(gene_data), "\n")
cat("Column names:", colnames(gene_data), "\n")
str(gene_data)   # Structure -- very useful
cat("\n")
summary(gene_data)  # Statistical summary of each column
```

```python
# Accessing data
cat("Gene column:", gene_data$gene, "\n")
cat("First row:\n")
print(gene_data[1, ])
cat("\nRow 2, column 'log2fc':", gene_data[2, "log2fc"], "\n")

# Filtering -- the most common operation
cat("\nOncogenes:\n")
print(gene_data[gene_data$is_oncogene, ])

cat("\nSignificantly upregulated (log2FC > 1, p < 0.05):\n")
sig_up <- gene_data[gene_data$log2fc > 1 & gene_data$pvalue < 0.05, ]
print(sig_up)
```

```python
# Modifying data frames

# Add a new column
gene_data$significant <- gene_data$pvalue < 0.05

# Bonferroni correction
gene_data$padj <- p.adjust(gene_data$pvalue, method = "bonferroni")

# BH/FDR correction (standard in genomics)
gene_data$fdr <- p.adjust(gene_data$pvalue, method = "BH")

# Sorting
gene_data_sorted <- gene_data[order(gene_data$pvalue), ]
print(gene_data_sorted)
```

---

## 6. Statistical Distributions

R has a consistent naming convention for distribution functions:

| Prefix | Meaning | Example |
|--------|---------|--------|
| `d` | Density (PDF value) | `dnorm(x)` |
| `p` | Cumulative probability (CDF) | `pnorm(q)` |
| `q` | Quantile (inverse CDF) | `qnorm(p)` |
| `r` | Random sampling | `rnorm(n)` |

Common distributions: `norm` (Normal), `t` (Student's t), `binom` (Binomial), `pois` (Poisson), `nbinom` (Negative binomial), `unif` (Uniform), `chisq` (Chi-squared), `f` (F).

```python
# Normal distribution
set.seed(42)  # For reproducibility

# Generate random expression values
expression_values <- rnorm(n = 100, mean = 10, sd = 2)

cat("Summary of simulated expression data:\n")
cat("  Mean:  ", round(mean(expression_values), 2), "\n")
cat("  SD:    ", round(sd(expression_values), 2), "\n")
cat("  Min:   ", round(min(expression_values), 2), "\n")
cat("  Max:   ", round(max(expression_values), 2), "\n")
cat("\n")
summary(expression_values)
```
