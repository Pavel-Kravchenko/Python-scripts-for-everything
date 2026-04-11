# Bioinformatics with Python: A Comprehensive Self-Paced Course

A modular, five-tier curriculum that takes you from computational basics to applied bioinformatics research skills. Built on materials from the Kodomo Bioinformatics Program at Moscow State University and the IAB open-source textbook, every concept is taught through hands-on Jupyter notebooks with real biological data. Whether you are a biology student learning to code, a programmer entering the life sciences, or a researcher looking to sharpen your computational toolkit, this course meets you where you are.

`107 notebooks` | `6 tiers` | `108 glossary terms` | `12 sample data files` | `30 interactive visualizations`

---

## Visual Course Map

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 0: COMPUTATIONAL FOUNDATIONS                             10 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  Skills Check │ Linux │ Git │ Bash │ Encodings │ R Basics │ Biostatistics  │
│  Probability & Statistics (Python) │ Advanced R Statistics                  │
│                                                                              │
│  Entry: No programming experience. Learn the tools every                     │
│         bioinformatician needs before writing Python.                         │
│                                                                              │
└──────────────────────────────────┬───────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 1: PYTHON FOR BIOINFORMATICS                             20 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  Skills Check │ Intro │ Variables │ Operators │ Strings │ Control Flow      │
│  Functions │ File I/O │ Lists │ Dicts │ Comprehensions │ Iterators          │
│  Regex │ OOP │ Decorators │ Error Handling │ NumPy/Pandas │ Wrangling       │
│  Visualization │ SQL                                                         │
│                                                                              │
│  Entry: Comfortable with the command line but new to Python.                 │
│         Every example uses biological data.                                  │
│                                                                              │
└──────────────────────────────────┬───────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 2: CORE BIOINFORMATICS                                   17 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  Skills Check │ Glossary │ Databases │ BioPython │ Pairwise Alignment       │
│  BLAST │ Multiple Alignment │ Phylogenetics │ Protein Structure             │
│  Nucleic Acid Structure │ Chromatograms │ Motifs & Domains                  │
│  GO & Pathways │ Comparative Genomics │ Computational Genetics              │
│  Hi-C Analysis │ Motif Discovery                                            │
│                                                                              │
│  Entry: Know Python, new to bioinformatics. The algorithmic                  │
│         and biological core of the field.                                    │
│                                                                              │
└──────────────────────────────────┬───────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 3: APPLIED BIOINFORMATICS                                28 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  Skills Check │ NGS │ Variant Calling │ RNA-seq │ Microbial Diversity       │
│  Promoters │ Statistics │ Machine Learning │ Capstone Project               │
│  Molecular Modeling │ Deep Learning │ Clinical Genomics │ Modern Workflows  │
│  Biochemistry & Enzyme Kinetics │ Genetic Engineering In Silico            │
│  Population Genetics │ Numerical Methods │ Genome Assembly                 │
│  Proteomics & Structural Methods │ GWAS │ Spatial Transcriptomics          │
│  Copy Number Analysis │ Bayesian Statistics │ TF Footprinting              │
│  Cancer Transcriptomics: Subtype Classification                             │
│                                                                              │
│  Entry: Have bioinformatics fundamentals, want real-world                    │
│         pipeline experience and advanced methods.                            │
│                                                                              │
│                            ┌──────────────┐                                  │
│                            │   CAPSTONE    │                                 │
│                            │   PROJECT     │                                 │
│                            │  (end-to-end  │                                 │
│                            │   analysis)   │                                 │
│                            └──────────────┘                                  │
│                                                                              │
└──────────────────────────────────┬───────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 4: ALGORITHMS & DATA STRUCTURES          30 notebooks + 30 interactive │
│  ──────────────────────────────────────────────────────────────────────────  │
│  Skills Check │ Complexity │ Sorting │ Searching │ Linked Lists             │
│  Stacks/Queues │ BST │ AVL │ Red-Black Trees │ Hash Tables                 │
│  Bloom Filters │ KMP │ Rabin-Karp │ Tries │ Suffix Trees │ Graphs │ DP    │
│                                                                              │
│  The CS theory behind bioinformatics tools. Study alongside                  │
│  Tiers 2--3: DP = alignment, string matching = BLAST, graphs = pathways.    │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  TIER 5: MODERN AI FOR SCIENCE                                  3 notebooks  │
│  ──────────────────────────────────────────────────────────────────────────  │
│  LLM Fine-tuning │ Vision RAG │ Diffusion & Generative Models                │
│                                                                              │
│  Entry: Tiers 1–3 complete. GPU-optional; runs on free-tier Colab.          │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘

  ENTRY POINTS ──────────────────────────────────────────────────────────────

  Take the Skills Check at the start of any tier.
  Score above 80%  ──>  skip that tier and move to the next.
  Score below 80%  ──>  work through the tier's notebooks.
  Tier 4 runs in parallel with Tiers 2--3 (not sequential).
```

---

## Who This Course Is For

| Audience | Where to Start | Tiers to Take | Estimated Time |
|----------|---------------|---------------|----------------|
| **Complete beginner** -- no programming experience, new to biology or coming from a wet-lab background | Tier 0, Module 01 (Linux Fundamentals) | All six tiers (0 through 5) | 160--215 hours |
| **Biology student** -- comfortable with a terminal and basic scripting, wants to learn Python for research | Tier 1 Skills Check, then Tier 1 Module 01 | Tiers 1, 2, 3 | 90--120 hours |
| **Experienced programmer** -- knows Python well, new to bioinformatics and biological data | Tier 2 Skills Check, then Tier 2 Module 01 | Tiers 2, 3 | 50--70 hours |
| **Bioinformatics researcher** -- has bioinformatics experience, wants applied pipelines and advanced methods | Tier 3 Skills Check, then individual Tier 3 modules | Tier 3 selectively | 25--35 hours |

---

## Complete Table of Contents

### Tier 0: Computational Foundations -- 10 notebooks

These modules cover the computing skills that every bioinformatician needs before writing Python. If you can navigate a Linux file system, use Git confidently, and write basic Bash scripts, take the Skills Check and skip ahead.

---

#### 0.00 Skills Check: Computational Foundations

[00_skills_check.ipynb](Tier_0_Computational_Foundations/00_Skills_Check/00_skills_check.ipynb) -- 52 cells

Self-assessment covering Linux commands, Git workflows, and Bash scripting. Fifteen graded questions with a detailed answer key. Score above 80% (12+ out of 15) to skip Tier 0. Includes section-specific remediation advice pointing to the exact notebook for each weak area.

`self-assessment` `linux` `git` `bash` `skip-ahead`

---

#### 0.01 Linux Fundamentals

[01_linux_fundamentals.ipynb](Tier_0_Computational_Foundations/01_Linux_Basics/01_linux_fundamentals.ipynb) -- 69 cells

File system navigation, creating/moving/copying files and directories, pipes and redirections for data-processing pipelines. Covers `grep`, `find`, `awk`, `sed`, compressed file handling (gzip, tar), process management, and SSH connections to remote servers. Every example uses bioinformatics scenarios -- FASTA files, BAM processing, FASTQ quality checks.

`file system` `pipes` `grep` `awk` `sed` `ssh` `permissions` `gzip`

---

#### 0.02 Git Version Control

[01_git_version_control.ipynb](Tier_0_Computational_Foundations/02_Git_Version_Control/01_git_version_control.ipynb) -- 59 cells

Why version control matters for reproducible research. Initializing repositories, staging and committing changes, viewing history with `log`/`diff`/`show`. Branching and merging for feature development. Collaborating with GitHub: push, pull, pull requests. Writing `.gitignore` files for bioinformatics projects (large data files, caches, OS artifacts). Undoing mistakes safely with `reset` and `revert`.

`git init` `branches` `merging` `github` `.gitignore` `reproducibility`

---

#### 0.03 Bash Scripting

[01_bash_scripting.ipynb](Tier_0_Computational_Foundations/03_Bash_Scripting/01_bash_scripting.ipynb) -- 72 cells

Writing and executing Bash scripts to automate repetitive bioinformatics tasks. Variables, conditionals, for/while loops, and functions. Batch processing of FASTQ, BAM, and VCF files. Robust scripts with `set -euo pipefail`, input validation, and error handling. Text processing with `awk`, `sed`, and `cut` inside scripts. Building reusable analysis pipeline scripts.

`bash scripts` `loops` `functions` `awk` `sed` `batch processing` `pipelines`

---

#### 0.04 Character Encodings

[01_character_encodings.ipynb](Tier_0_Computational_Foundations/04_File_Encodings/01_character_encodings.ipynb) -- 44 cells

From bits to bytes: binary, decimal, and hexadecimal representations. ASCII as the universal foundation. UTF-8 encoding and why it matters when collaborators share files across operating systems. Binary vs. text file modes in Python. Understanding encoding issues in bioinformatics data formats -- FASTA, GenBank, GFF, PDB. Practical debugging of garbled text and encoding errors.

`ASCII` `UTF-8` `binary` `text modes` `encoding errors` `cross-platform`

---

#### 0.05 R Fundamentals

[01_r_fundamentals.ipynb](Tier_0_Computational_Foundations/05_R_Basics/01_r_fundamentals.ipynb) -- 51 cells

Why R matters in bioinformatics: Bioconductor, DESeq2, edgeR, Seurat, ggplot2. R basics -- assignment, types, vectors, indexing differences from Python (1-based, negative indexing means exclusion). Vectorized operations as R's superpower. Data frames, factors, basic plotting. Enough R literacy to read R-based methods papers and run Bioconductor packages when needed.

`R syntax` `vectors` `data frames` `Bioconductor` `ggplot2` `vectorized operations`

---

#### 0.06 Biostatistics Fundamentals

[01_biostatistics_fundamentals.ipynb](Tier_0_Computational_Foundations/06_Biostatistics/01_biostatistics_fundamentals.ipynb) -- 52 cells

Populations vs. samples, parameter estimation, types of data (nominal, ordinal, discrete, continuous). Probability distributions relevant to biology: normal, Poisson, negative binomial, binomial. Hypothesis testing framework: null and alternative hypotheses, test statistics, p-values, significance levels. Multiple testing correction (Bonferroni, Benjamini-Hochberg FDR). Practical examples from differential expression, variant calling, and enrichment analysis.

`distributions` `hypothesis testing` `p-values` `multiple testing` `FDR` `Bonferroni`

---

#### 0.07 Probability and Statistics with Python

[01_probability_and_statistics_python.ipynb](Tier_0_Computational_Foundations/07_Probability_and_Statistics_Python/01_probability_and_statistics_python.ipynb) -- 68 cells

Python-native statistical analysis using scipy.stats and statsmodels. Probability distributions (normal, binomial, Poisson, negative binomial) with biological examples. Hypothesis testing: t-tests, Mann-Whitney U, chi-squared, Fisher's exact. Multiple group comparisons with ANOVA and Kruskal-Wallis. Correlation and linear regression with full diagnostics. Bootstrap confidence intervals and permutation tests. Power analysis for experimental design.

`scipy.stats` `statsmodels` `distributions` `regression` `bootstrap` `power analysis` `ANOVA`

---

#### 0.08 Advanced R Statistics

[01_r_hypothesis_testing_and_nonparametrics.ipynb](Tier_0_Computational_Foundations/08_Advanced_R_Statistics/01_r_hypothesis_testing_and_nonparametrics.ipynb) -- 57 cells
[02_r_regression_correlation_and_diagnostics.ipynb](Tier_0_Computational_Foundations/08_Advanced_R_Statistics/02_r_regression_correlation_and_diagnostics.ipynb) -- 76 cells

Hands-on R statistical computing with real biological datasets. Notebook 1 covers nonparametric methods: exact binomial tests, sign test, Wilcoxon signed-rank and rank-sum, Mann-Whitney U, Kruskal-Wallis with Dunn post-hoc, Hodges-Lehman estimation, and power analysis. Notebook 2 covers parametric methods and diagnostics: confidence intervals (Student-t, chi-squared, asymptotic), normality testing (Shapiro-Wilk, QQ plots), Pearson/Spearman/Kendall correlation with Fisher z-transform, chi-squared tests for independence and homogeneity, linear and polynomial regression with residual diagnostics, ANOVA vs Kruskal-Wallis, and Central Limit Theorem simulation. Adapted from ФББ Semester 7 biostatistics seminars by Pervushin/Muromskaya.

`R` `nonparametric tests` `Wilcoxon` `Hodges-Lehman` `regression diagnostics` `chi-squared` `CLT` `power analysis`

---

### Tier 1: Python for Bioinformatics -- 20 notebooks

Core Python programming taught entirely through biological examples. From your first program to publication-quality figures, every exercise works with sequences, gene data, or biological datasets.

---

#### 1.00 Skills Check: Python Fundamentals

[00_skills_check.ipynb](Tier_1_Python_for_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) -- 48 cells

Self-assessment covering Python basics: variables, data types, control flow, functions, file operations, and data structures. Score above 80% to skip Tier 1 and proceed directly to Tier 2 (Core Bioinformatics).

`self-assessment` `python basics` `skip-ahead`

---

#### 1.01 Python Introduction

[01_python_introduction.ipynb](Tier_1_Python_for_Bioinformatics/01_Python_Introduction/01_python_introduction.ipynb) -- 52 cells

Why Python dominates bioinformatics. Navigating Jupyter notebooks: editing cells, running code, interpreting output. Writing your first bioinformatics program. The Python ecosystem for biology -- BioPython, pandas, NumPy, matplotlib. Installation and environment setup. REPL basics and interactive exploration.

`jupyter` `python setup` `ecosystem` `first program` `REPL`

---

#### 1.02 Variables and Data Types

[01_variables_and_data_types.ipynb](Tier_1_Python_for_Bioinformatics/02_Variables_and_Data_Types/01_variables_and_data_types.ipynb) -- 64 cells

Creating variables with Python naming conventions. All basic data types: `int`, `float`, `str`, `bool`, `None`. Type conversion and type checking. Using strings to represent biological sequences. F-strings for formatted output. Real bioinformatics calculations: molecular weights, GC content, sequence length statistics.

`variables` `int` `float` `str` `bool` `type conversion` `f-strings`

---

#### 1.03 Operators and Expressions

[01_operators_and_expressions.ipynb](Tier_1_Python_for_Bioinformatics/03_Operators_and_Expressions/01_operators_and_expressions.ipynb) -- 63 cells

Arithmetic, comparison, logical, and bitwise operators. Operator precedence and associativity. Boolean expressions for filtering biological data. Chaining comparisons. Membership testing with `in`. Identity vs. equality (`is` vs. `==`). Applied exercises with sequence analysis calculations.

`arithmetic` `comparison` `logical operators` `precedence` `boolean expressions`

---

#### 1.04 Strings and Sequences

[01_strings_and_sequences.ipynb](Tier_1_Python_for_Bioinformatics/04_Strings_and_Sequences/01_strings_and_sequences.ipynb) -- 64 cells

String methods for biological sequences: indexing, slicing, and immutability. DNA/RNA/protein alphabets as strings. Reverse complement, transcription (DNA to RNA), and codon extraction. Motif finding and pattern searching. FASTA header parsing. GC content calculation. F-string formatting for scientific output.

`string methods` `slicing` `reverse complement` `transcription` `codons` `motifs` `GC content`

---

#### 1.05 Control Flow

[01_control_flow.ipynb](Tier_1_Python_for_Bioinformatics/05_Control_Flow/01_control_flow.ipynb) -- 65 cells

Conditional statements (`if`/`elif`/`else`) for classifying biological data -- purines vs. pyrimidines, amino acid properties. `for` and `while` loops for iterating over sequences and datasets. Loop control with `break`, `continue`, and the `for...else` construct. Using `enumerate()`, `zip()`, and `range()` effectively. Codon-to-amino-acid translation using a dictionary lookup inside a loop.

`if/elif/else` `for loops` `while` `break` `continue` `enumerate` `zip` `codon translation`

---

#### 1.06 Functions

[01_functions.ipynb](Tier_1_Python_for_Bioinformatics/06_Functions/01_functions.ipynb) -- 57 cells

Defining and calling functions with various parameter types. Default arguments, `*args`, and `**kwargs`. Lambda functions and `map`/`filter`/`reduce`. Variable scope and the LEGB rule. Recursive algorithms. Docstrings and function documentation. Biological applications: GC content functions, sequence validators, ORF finders.

`def` `arguments` `lambda` `map` `filter` `scope` `recursion` `*args` `**kwargs`

---

#### 1.07 File Operations

[01_file_operations.ipynb](Tier_1_Python_for_Bioinformatics/07_File_Operations/01_file_operations.ipynb) -- 66 cells

Reading and writing text files with context managers (`with` statement). File modes: read, write, append, binary. Reading methods compared (`read()`, `readline()`, `readlines()`, iteration). Parsing FASTA files from scratch -- building both a reader and a writer. Working with CSV, JSON, and Pickle formats. Understanding GenBank format. Path handling and directory operations.

`file I/O` `context managers` `FASTA parsing` `CSV` `JSON` `GenBank` `binary files`

---

#### 1.08 Lists and Tuples

[01_lists_and_tuples.ipynb](Tier_1_Python_for_Bioinformatics/08_Lists_and_Tuples/01_lists_and_tuples.ipynb) -- 50 cells

List creation, indexing, slicing, and all core methods (`append`, `extend`, `insert`, `remove`, `sort`). Tuple immutability, packing/unpacking, and named tuples. Sorting and filtering biological sequences by length and composition. Gene coordinate storage with tuples. Nested data structures for representing biological hierarchies.

`lists` `tuples` `named tuples` `sorting` `slicing` `gene coordinates`

---

#### 1.09 Dictionaries and Sets

[01_dictionaries_and_sets.ipynb](Tier_1_Python_for_Bioinformatics/09_Dictionaries_and_Sets/01_dictionaries_and_sets.ipynb) -- 48 cells

Dictionary creation, access, methods, and iteration. Nested dictionaries for complex biological data. Specialized collections: `defaultdict` and `Counter`. Set operations: union, intersection, difference, symmetric difference, and `frozenset`. Biological applications: codon tables, gene annotation lookup, k-mer counting, comparing gene sets across experiments.

`dicts` `sets` `defaultdict` `Counter` `codon tables` `k-mers` `gene sets`

---

#### 1.10 Comprehensions

[01_comprehensions.ipynb](Tier_1_Python_for_Bioinformatics/10_Comprehensions/01_comprehensions.ipynb) -- 50 cells

List comprehensions with and without conditions. Dictionary and set comprehensions. Nested comprehensions for matrix operations. Generator expressions for memory-efficient processing. Applied to biological data: GC content filtering, codon extraction, reverse complement of sequence lists, k-mer frequency tables.

`list comprehensions` `dict comprehensions` `set comprehensions` `generator expressions` `k-mers`

---

#### 1.11 Iterators and Generators

[01_iterators_and_generators.ipynb](Tier_1_Python_for_Bioinformatics/11_Iterators_and_Generators/01_iterators_and_generators.ipynb) -- 84 cells

The iteration protocol: `__iter__` and `__next__`. How `for` loops really work under the hood. Using `iter()` and `next()` manually. Creating generators with `yield` and generator expressions. The `itertools` module: `chain`, `combinations`, `permutations`, `product`, `groupby`, `islice`. Lazy evaluation for memory-efficient processing of large biological datasets (genome-scale FASTA files, streaming FASTQ).

`iterators` `generators` `yield` `itertools` `lazy evaluation` `memory efficiency`

---

#### 1.12 Regular Expressions

[01_regular_expressions.ipynb](Tier_1_Python_for_Bioinformatics/12_Regular_Expressions/01_regular_expressions.ipynb) -- 82 cells

The `re` module: `search`, `match`, `findall`, `finditer`, `sub`, `compile`. Metacharacters, character classes, quantifiers, groups, and backreferences. Lookahead and lookbehind assertions. Flags: `IGNORECASE`, `MULTILINE`, `DOTALL`, `VERBOSE`. Bioinformatics applications: restriction enzyme site patterns, NCBI gene ID extraction, PROSITE-to-regex conversion, protein motif searching, parsing structured biological text.

`re module` `metacharacters` `groups` `lookahead` `restriction sites` `PROSITE` `motif search`

---

#### 1.13 Classes and OOP

[01_classes_and_oop.ipynb](Tier_1_Python_for_Bioinformatics/13_Classes_and_OOP/01_classes_and_oop.ipynb) -- 45 cells

Why OOP matters for bioinformatics: modeling genes, sequences, and biological entities as objects. Classes, `__init__` constructor, instance attributes and methods. Dunder methods: `__str__`, `__repr__`, `__eq__`, `__lt__`. Inheritance and polymorphism. Abstract base classes. Properties, class methods, and static methods. Dataclasses for concise data modeling.

`classes` `__init__` `inheritance` `polymorphism` `dunder methods` `dataclasses` `properties`

---

#### 1.14 Decorators and Context Managers

[01_decorators_and_context_managers.ipynb](Tier_1_Python_for_Bioinformatics/14_Decorators_and_Context_Managers/01_decorators_and_context_managers.ipynb) -- 42 cells

First-class functions and closures as the foundation for decorators. Writing decorators from scratch, with and without arguments. `functools.wraps` for preserving function metadata. Context managers with `__enter__`/`__exit__` and the `contextlib` module. Applied to bioinformatics workflows: timing sequence analysis functions, caching BLAST results, managing database connections.

`decorators` `closures` `functools.wraps` `context managers` `contextlib` `with statement`

---

#### 1.15 Error Handling

[01_error_handling.ipynb](Tier_1_Python_for_Bioinformatics/15_Error_Handling/01_error_handling.ipynb) -- 46 cells

Structured error handling with `try`/`except`/`else`/`finally`. Raising exceptions and creating custom exception classes. The Python exception hierarchy. Debugging techniques: print debugging, `assert` statements, `pdb`/`breakpoint()`. The `logging` module for production bioinformatics code. Handling messy bioinformatics data: missing FASTA headers, ambiguous nucleotides, network timeouts when querying NCBI.

`try/except` `custom exceptions` `debugging` `pdb` `logging` `assert`

---

#### 1.16 NumPy and Pandas

[01_numpy_and_pandas.ipynb](Tier_1_Python_for_Bioinformatics/16_NumPy_and_Pandas/01_numpy_and_pandas.ipynb) -- 94 cells

NumPy fundamentals: creating arrays, vectorized operations vs. Python loops, broadcasting, fancy indexing, statistical functions. Pandas Series and DataFrames for tabular biological data. Reading and writing CSV, TSV, and Excel files. Indexing, selection, filtering, grouping, and aggregation. Merging and joining datasets. Applied to gene expression matrices, annotation tables, and experimental metadata.

`NumPy arrays` `vectorized operations` `broadcasting` `DataFrames` `CSV` `groupby` `merge`

---

#### 1.17 Data Wrangling

[01_data_wrangling.ipynb](Tier_1_Python_for_Bioinformatics/17_Data_Wrangling/01_data_wrangling.ipynb) -- 80 cells

Handling missing values: when to drop vs. fill, biological context for each strategy. Removing duplicates. Data type conversion and cleaning messy annotation tables. Reshaping between wide and long formats with `melt`, `pivot`, `stack`/`unstack`. String operations for extracting information from annotation columns. Custom transformations with `apply` and `applymap`. Method chaining for clean, readable pipelines.

`missing values` `duplicates` `melt` `pivot` `string operations` `apply` `method chaining`

---

#### 1.18 Data Visualization

[01_data_visualization.ipynb](Tier_1_Python_for_Bioinformatics/18_Data_Visualization/01_data_visualization.ipynb) -- 72 cells

Matplotlib fundamentals: figures, axes, the object-oriented interface. All major plot types: line, scatter, bar, histogram, box, violin, heatmap. Customization: colors, labels, legends, annotations, styles. Seaborn for statistical plots, pair plots, and clustered heatmaps. Biological plot types: volcano plots, MA plots, sequence logos, genome coverage plots. Multi-panel figures and saving publication-quality output.

`matplotlib` `seaborn` `heatmaps` `volcano plots` `publication figures` `multi-panel`

---

#### 1.19 SQL for Bioinformatics

[01_sql_for_bioinformatics.ipynb](Tier_1_Python_for_Bioinformatics/19_SQL_for_Bioinformatics/01_sql_for_bioinformatics.ipynb) -- 20 cells

Querying biological databases using SQL with Python's built-in `sqlite3` module -- no server installation required. SELECT, WHERE, ORDER BY, GROUP BY, and aggregate functions (COUNT, AVG, MAX, MIN, SUM). JOIN operations: INNER, LEFT, and self-joins for connecting gene annotation tables. Subqueries and set operations with UNION. CREATE TABLE, INSERT, UPDATE, DELETE for building local biological databases. Integrating sqlite3 with pandas: reading query results directly into DataFrames, writing DataFrames to SQL tables.

`sqlite3` `SELECT` `JOIN` `GROUP BY` `subqueries` `pandas` `biological databases`

---

### Tier 2: Core Bioinformatics -- 17 notebooks

The heart of bioinformatics: databases, algorithms, structural biology, and analysis methods. You will learn how to search databases, align sequences, build phylogenetic trees, analyze protein structures, and interpret functional annotations.

---

#### 2.00 Skills Check: Core Bioinformatics

[00_skills_check.ipynb](Tier_2_Core_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) -- 50 cells

Self-assessment covering BLAST, sequence alignment, phylogenetics, protein structure, and functional annotation. Fifteen questions test whether you can skip Tier 2 and proceed to applied bioinformatics. Includes questions on E-value interpretation, substitution matrices, tree topology, PDB format, and GO enrichment.

`self-assessment` `BLAST` `alignment` `phylogenetics` `skip-ahead`

---

#### 2.00b Bioinformatics Glossary

[01_glossary.ipynb](Tier_2_Core_Bioinformatics/00_Skills_Check/01_glossary.ipynb) -- 46 cells

A comprehensive A-Z reference of 108 bioinformatics terms used across all five tiers. Each entry includes a definition and a cross-reference to the notebook where the term is covered in depth. Use `Ctrl+F`/`Cmd+F` to search, or browse alphabetically. Covers everything from ACMG variant classification to zero-mode waveguides.

`glossary` `reference` `108 terms` `cross-references`

---

#### 2.01 Biological Databases

[01_biological_databases.ipynb](Tier_2_Core_Bioinformatics/01_Biological_Databases/01_biological_databases.ipynb) -- 52 cells

Navigating the NCBI ecosystem: GenBank, RefSeq, Gene, SRA, GEO, and the Entrez API. Understanding UniProt (Swiss-Prot vs. TrEMBL) and protein annotation. The Protein Data Bank (PDB) for structural data. Programmatic access to databases with BioPython's `Entrez` module. Cross-references between databases. Accession numbers, identifiers, and how records connect.

`NCBI` `GenBank` `UniProt` `PDB` `Entrez API` `RefSeq` `accession numbers`

---

#### 2.02 BioPython Essentials

[01_biopython_essentials.ipynb](Tier_2_Core_Bioinformatics/02_BioPython_Essentials/01_biopython_essentials.ipynb) -- 52 cells

The `Seq` object: complement, reverse complement, transcription, translation. `SeqRecord` for annotated sequences. `SeqIO` for reading and writing biological file formats (FASTA, GenBank, EMBL). Fetching sequences from NCBI programmatically. A complete gene-to-protein analysis workflow: fetch a gene, transcribe, find ORFs, translate, search for the protein.

`Seq` `SeqRecord` `SeqIO` `complement` `translate` `Entrez` `gene-to-protein`

---

#### 2.03 Pairwise Sequence Alignment

[01_pairwise_sequence_alignment.ipynb](Tier_2_Core_Bioinformatics/03_Pairwise_Sequence_Alignment/01_pairwise_sequence_alignment.ipynb) -- 64 cells

The biological motivation for aligning sequences: detecting evolutionary relationships through substitutions, insertions, and deletions. Dot plots as a visual alignment tool. Substitution matrices: PAM and BLOSUM families, how they are derived, when to use each. Gap penalty schemes: linear vs. affine. Implementing Needleman-Wunsch (global) and Smith-Waterman (local) algorithms step by step with dynamic programming. E-values and statistical significance. BioPython's `pairwise2` and `Align` modules.

`dot plots` `BLOSUM` `PAM` `Needleman-Wunsch` `Smith-Waterman` `gap penalties` `E-values`

---

#### 2.04 BLAST Searching

[01_blast_searching.ipynb](Tier_2_Core_Bioinformatics/04_BLAST_Searching/01_blast_searching.ipynb) -- 54 cells

Why BLAST is needed: the computational cost of optimal alignment at database scale. The seed-and-extend heuristic: word matching, neighborhood generation, ungapped extension, gapped extension. Choosing the right BLAST variant: `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`. Parameter tuning: word size, E-value threshold, scoring matrices. Running BLAST via NCBI web interface and programmatically with BioPython's `NCBIWWW` and `NCBIXML`. Parsing and interpreting results. Taxonomy-filtered searches.

`BLAST variants` `seed-and-extend` `E-value` `NCBIWWW` `NCBIXML` `parameter tuning`

---

#### 2.05 Multiple Sequence Alignment

[01_multiple_sequence_alignment.ipynb](Tier_2_Core_Bioinformatics/05_Multiple_Sequence_Alignment/01_multiple_sequence_alignment.ipynb) -- 56 cells

Why MSA reveals what pairwise alignment cannot: conserved regions, functional domains, evolutionary relationships across entire families. The computational challenge: NP-hardness and the need for heuristics. Progressive alignment (guide tree approach). Running ClustalW and MUSCLE from Python. Conservation scoring and consensus sequences. MSA visualization and interpretation. Identifying conserved functional residues from alignments.

`MSA` `ClustalW` `MUSCLE` `conservation` `consensus` `guide tree` `NP-hard`

---

#### 2.06 Phylogenetics

[01_phylogenetics.ipynb](Tier_2_Core_Bioinformatics/06_Phylogenetics/01_phylogenetics.ipynb) -- 61 cells

Phylogenetic tree terminology: nodes, branches, clades, rooted vs. unrooted trees. Distance-based methods: UPGMA and Neighbor-Joining. Character-based methods: maximum parsimony, maximum likelihood. Distance matrices from MSA. Building trees with BioPython's `Phylo` module. Tree visualization and annotation. Bootstrap analysis for branch support. Practical applications: gene family evolution, viral tracking, species classification.

`phylogenetic trees` `UPGMA` `Neighbor-Joining` `bootstrap` `Phylo` `distance matrix` `clades`

---

#### 2.07 Protein Structure

[01_protein_structure.ipynb](Tier_2_Core_Bioinformatics/07_Protein_Structure/01_protein_structure.ipynb) -- 50 cells

The four levels of protein structure: primary through quaternary. The PDB file format and ATOM record fields. BioPython's `Bio.PDB` module: navigating the Structure/Model/Chain/Residue/Atom hierarchy. Calculating interatomic distances, bond angles, and dihedral angles. Ramachandran plots. RMSD calculation for structural comparison. Secondary structure assignment (DSSP). 3D visualization with py3Dmol in Jupyter notebooks.

`PDB format` `Bio.PDB` `Ramachandran` `RMSD` `DSSP` `py3Dmol` `secondary structure`

---

#### 2.08 Nucleic Acid Structure

[01_nucleic_acid_structure.ipynb](Tier_2_Core_Bioinformatics/08_Nucleic_Acid_Structure/01_nucleic_acid_structure.ipynb) -- 43 cells

Chemical components of nucleotides: base, sugar, phosphate. Watson-Crick base pairing and hydrogen bonds. DNA helical forms: A-DNA, B-DNA, Z-DNA and their biological significance. Major and minor groove geometry and protein-DNA recognition. RNA secondary structure elements: stems, loops, bulges, pseudoknots. Introduction to 3DNA for nucleic acid structural analysis.

`nucleotides` `base pairing` `A/B/Z-DNA` `major/minor groove` `RNA secondary structure` `3DNA`

---

#### 2.09 Chromatogram Analysis

[01_chromatogram_analysis.ipynb](Tier_2_Core_Bioinformatics/09_Chromatogram_Analysis/01_chromatogram_analysis.ipynb) -- 46 cells

How Sanger (dideoxy chain termination) sequencing works. Reading `.ab1` chromatogram trace files with BioPython. Plotting and interpreting trace data with matplotlib -- the four fluorescence channels. Phred quality scores and base calling confidence. Identifying mixed bases, heterozygous positions, and low-quality regions. Quality trimming strategies. Building a consensus from forward and reverse reads.

`Sanger sequencing` `.ab1 files` `Phred scores` `base calling` `trace data` `consensus`

---

#### 2.10 Sequence Motifs and Domains

[01_sequence_motifs_and_domains.ipynb](Tier_2_Core_Bioinformatics/10_Sequence_Motifs_and_Domains/01_sequence_motifs_and_domains.ipynb) -- 53 cells

What sequence motifs are and why they matter biologically: transcription factor binding sites, splice signals, post-translational modification sites. Constructing Position Weight Matrices (PWMs) from aligned binding sites: PFM, PPM, PWM with pseudocounts. Information content and sequence logos. PROSITE patterns and regex conversion. Protein domains and the major databases: Pfam, InterPro, CDD. Scanning sequences for known domains.

`PWM` `sequence logos` `PROSITE` `Pfam` `InterPro` `TATA box` `binding sites`

---

#### 2.11 Gene Ontology and Pathways

[01_gene_ontology_and_pathways.ipynb](Tier_2_Core_Bioinformatics/11_Gene_Ontology_and_Pathways/01_gene_ontology_and_pathways.ipynb) -- 53 cells

The Gene Ontology system: three ontologies (Molecular Function, Biological Process, Cellular Component), DAG structure, GO terms, evidence codes. GO enrichment analysis using the hypergeometric test and Fisher's exact test. The KEGG database: pathway maps, enzymes, metabolic networks. Functional annotation workflows. Interpreting enrichment results in biological context: overrepresentation, gene set enrichment.

`Gene Ontology` `GO enrichment` `KEGG` `DAG` `hypergeometric test` `functional annotation`

---

#### 2.12 Comparative Genomics

[01_comparative_genomics.ipynb](Tier_2_Core_Bioinformatics/12_Comparative_Genomics/01_comparative_genomics.ipynb) -- 45 cells

Goals and rationale of comparative genomics: identifying conserved elements, detecting rearrangements, understanding genome evolution. Implementing dot plots from scratch and interpreting visual patterns -- diagonal lines (collinearity), breaks (rearrangements), parallel diagonals (repeats). Synteny analysis and conserved gene order. Distinguishing orthologs from paralogs. Proteome-level comparisons.

`dot plots` `synteny` `orthologs` `paralogs` `genome rearrangements` `comparative analysis`

---

#### 2.13 Computational Genetics

[01_computational_genetics.ipynb](Tier_2_Core_Bioinformatics/13_Computational_Genetics/01_computational_genetics.ipynb) -- 64 cells

Programming the genetic code: building codon tables, translation, degeneracy analysis. Codon usage bias: RSCU, Codon Adaptation Index (CAI), GC content at codon positions. Restriction enzyme analysis: recognition sites, virtual digests, gel simulation, compatible ends. Open reading frame finding across all six frames. Genetic mapping: two-point and three-point crosses, map distance, interference. Mutation analysis: transition/transversion ratios, mutation spectra, CpG deamination. Hardy-Weinberg equilibrium testing.

`genetic code` `codon usage` `CAI` `restriction enzymes` `ORF finding` `genetic mapping` `Ts/Tv` `Hardy-Weinberg`

---

#### 2.14 Hi-C Analysis

[14_hic_analysis.ipynb](Tier_2_Core_Bioinformatics/14_Hi-C_Analysis/14_hic_analysis.ipynb)

3D genome organization from Hi-C experiments using the cooler and cooltools Python stack. cooler file format: loading, inspecting metadata, and slicing contact matrices. Contact decay (expected) curves to normalize distance effects. Eigenvector decomposition for A/B compartment identification. Insulation score calculation and TAD boundary detection. Saddle plots for compartment strength visualization. Pileup (aggregate) analysis around genomic features such as CTCF sites and loop anchors. Uses public Hi-C data from the 4DN Data Portal or ENCODE.

`cooler` `cooltools` `contact matrices` `A/B compartments` `TADs` `insulation score` `saddle plots` `pileup`

---

#### 2.15 Motif Discovery

[15_motif_discovery.ipynb](Tier_2_Core_Bioinformatics/15_Motif_Discovery/15_motif_discovery.ipynb)

Quantitative motif analysis from position frequency matrices to enrichment testing. PPM/PWM construction, normalization, and pseudocount handling. Information content per position and total IC calculation. KDIC score (mean IC normalized to [0,1]). IUPAC consensus sequence generation. Score distributions: exact enumeration for short motifs (length ≤ 10), Monte Carlo with confidence intervals for longer. Motif enrichment using Fisher's exact test with Benjamini-Hochberg correction. TomTom matching concept against JASPAR 2024 and HOCOMOCO public databases. Pipeline design patterns using abstract interfaces and dataclasses; BED/FASTA I/O patterns.

`PWM` `PPM` `information content` `KDIC` `IUPAC` `Fisher enrichment` `Benjamini-Hochberg` `TomTom` `JASPAR`

---

### Tier 3: Applied Bioinformatics -- 27 notebooks

Advanced topics and real-world analysis pipelines. Each notebook covers a complete workflow from raw data to biological conclusions. Includes a capstone project integrating skills from every tier, plus specialized modules on molecular modeling, deep learning, clinical genomics, modern bioinformatics workflows, GWAS, spatial transcriptomics, copy number analysis, Bayesian statistics, and TF footprinting.

---

#### 3.00 Skills Check: Applied Bioinformatics

[00_skills_check.ipynb](Tier_3_Applied_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) -- 35 cells

Self-assessment covering NGS concepts, file formats (FASTQ, SAM/BAM, VCF), RNA-seq workflows, and statistical testing. Tests knowledge of sequencing platforms (Illumina, PacBio, Nanopore), quality control strategies, and variant interpretation. Score above 80% to selectively skip Tier 3 modules you already know.

`self-assessment` `NGS` `RNA-seq` `variant calling` `skip-ahead`

---

#### 3.01 NGS Fundamentals

[01_ngs_fundamentals.ipynb](Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb) -- 47 cells

The complete NGS pipeline from sequencing to alignment. How Illumina, PacBio SMRT, and Oxford Nanopore sequencing work -- principles, strengths, and error profiles. FASTQ format: structure, quality encoding, parsing. Quality control with FastQC: per-base quality, adapter content, duplication levels. Read trimming with Trimmomatic and fastp. SAM/BAM format: header, alignment records, FLAG fields, CIGAR strings. Alignment with BWA and HISAT2.

`FASTQ` `FastQC` `trimming` `SAM/BAM` `CIGAR` `Illumina` `PacBio` `Nanopore` `BWA`

---

#### 3.01b Bioinformatics Data Formats

[02_bio_data_formats.ipynb](Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/02_bio_data_formats.ipynb) -- 35 cells

A comprehensive guide to every major bioinformatics file format with hands-on Python parsing examples. FASTA: structure, .fai indexing, random-access fetch. FASTQ: Phred encoding, quality histograms. SAM/BAM/CRAM: FLAG decoding, CIGAR string parsing. VCF/BCF: variant records, genotype fields, INFO/FORMAT parsing. BED: 0-based intervals, overlap detection. GFF3/GTF: gene annotations, attribute parsing. WIG/BedGraph/BigWig: coverage tracks. PDB: fixed-width ATOM/HETATM records, backbone extraction, distance calculation. mmCIF/PDBx: modern structure format, AlphaFold pLDDT scores. FAST5/POD5: Oxford Nanopore raw signal. Newick/Nexus/NHX: phylogenetic tree serialization and recursive parsing.

`FASTA` `FASTQ` `SAM/BAM/CRAM` `VCF/BCF` `BED` `GFF3` `GTF` `WIG` `BedGraph` `BigWig` `PDB` `mmCIF` `FAST5` `POD5` `Newick` `Nexus` `NHX`

---

#### 3.02 Variant Calling and SNP Analysis

[01_variant_calling_and_snp_analysis.ipynb](Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/01_variant_calling_and_snp_analysis.ipynb) -- 41 cells

Classifying variant types: SNPs, indels, structural variants, CNVs. The variant calling pipeline: marking duplicates, base quality recalibration, calling with GATK HaplotypeCaller / bcftools / FreeBayes. VCF format in detail: header, INFO fields, FORMAT/genotype columns, FILTER. Hard filtering vs. VQSR. Variant annotation with SnpEff and VEP. Population genetics: allele frequencies, Hardy-Weinberg equilibrium, linkage disequilibrium. GWAS concepts and ClinVar.

`VCF` `GATK` `SNPs` `indels` `SnpEff` `VEP` `allele frequency` `GWAS` `ClinVar`

---

#### 3.03 RNA-seq Analysis

[01_rna_seq_analysis.ipynb](Tier_3_Applied_Bioinformatics/03_RNA_seq_Analysis/01_rna_seq_analysis.ipynb) -- 59 cells

Complete RNA-seq workflow from experimental design to differential expression. Advantages over microarrays. Experimental design: biological replicates, sequencing depth, batch effects. The analysis pipeline: alignment (STAR, HISAT2), quantification (featureCounts, HTSeq, salmon), count normalization (RPKM, FPKM, TPM, DESeq2 size factors). Differential expression with DESeq2 and edgeR (conceptual and R integration). Visualization: volcano plots, MA plots, heatmaps, PCA.

`RNA-seq` `DESeq2` `differential expression` `TPM` `FPKM` `volcano plots` `PCA`

---

#### 3.04 Microbial Diversity

[01_microbial_diversity.ipynb](Tier_3_Applied_Bioinformatics/04_Microbial_Diversity/01_microbial_diversity.ipynb) -- 59 cells

16S rRNA gene as a universal marker for bacteria and archaea. Amplicon sequencing workflow from environmental sample to OTU/ASV tables. Denoising with DADA2 vs. OTU clustering. Taxonomic classification against reference databases (SILVA, Greengenes, UNITE). Alpha diversity metrics: observed species, Shannon, Simpson, Chao1. Beta diversity: Bray-Curtis, UniFrac distances. Ordination and visualization: PCoA, NMDS. Statistical testing: PERMANOVA, ANOSIM. Rarefaction curves.

`16S rRNA` `OTU` `ASV` `DADA2` `alpha diversity` `beta diversity` `UniFrac` `PCoA`

---

#### 3.05 Promoter and Regulatory Analysis

[01_promoter_and_regulatory_analysis.ipynb](Tier_3_Applied_Bioinformatics/05_Promoter_and_Regulatory_Analysis/01_promoter_and_regulatory_analysis.ipynb) -- 44 cells

Gene regulation and regulatory elements: promoters, enhancers, silencers. TATA box identification and consensus sequences. CpG islands: definition, detection algorithms, role in gene regulation. Transcription factor binding site (TFBS) prediction using PWMs. Scanning genomic sequences for motif hits. DNA methylation patterns and their biological significance. Integrating regulatory predictions with expression data.

`promoters` `TATA box` `CpG islands` `TFBS` `PWM scanning` `methylation` `gene regulation`

---

#### 3.06 Statistics for Bioinformatics

[01_statistics_for_bioinformatics.ipynb](Tier_3_Applied_Bioinformatics/06_Statistics_for_Bioinformatics/01_statistics_for_bioinformatics.ipynb) -- 50 cells

Distributions in biology: normal, Poisson, negative binomial, and when each applies. Hypothesis testing review: t-tests, chi-squared, Mann-Whitney, Kruskal-Wallis. Multiple testing correction in depth: family-wise error rate, Bonferroni, Benjamini-Hochberg FDR, q-values. Survival analysis (Kaplan-Meier, log-rank test). Power analysis and sample size calculations. Resampling methods: permutation tests and bootstrapping. Practical advice for choosing the right test.

`distributions` `multiple testing` `FDR` `survival analysis` `power analysis` `bootstrap` `permutation`

---

#### 3.07 Machine Learning for Biology

[01_machine_learning_for_biology.ipynb](Tier_3_Applied_Bioinformatics/07_Machine_Learning_for_Biology/01_machine_learning_for_biology.ipynb) -- 58 cells

ML use cases in bioinformatics: variant pathogenicity, gene expression prediction, protein function classification. Classification vs. regression. The scikit-learn workflow: data splitting, feature engineering, model training, cross-validation, hyperparameter tuning. Algorithms: logistic regression, decision trees, random forests, SVM, k-nearest neighbors. Evaluation metrics: accuracy, precision, recall, F1, ROC-AUC. Feature selection for biological interpretability. Unsupervised learning: PCA, k-means, hierarchical clustering.

`scikit-learn` `random forests` `SVM` `cross-validation` `ROC-AUC` `PCA` `clustering`

---

#### 3.08 Capstone Project: From Sequence to Discovery

[01_capstone_project.ipynb](Tier_3_Applied_Bioinformatics/08_Capstone_Project/01_capstone_project.ipynb) -- 68 cells

An end-to-end integrative bioinformatics project analyzing unknown DNA sequences. You will identify sequences via BLAST, perform multiple sequence alignment, construct phylogenetic trees, analyze protein structure, and draw biological conclusions. Integrates skills from every tier: file handling (Tier 1), BioPython and alignment (Tier 2), and statistical evaluation (Tier 3). Designed to be completed in 8--12 hours across multiple sessions. The best way to verify readiness for independent bioinformatics work.

`capstone` `integrative analysis` `BLAST` `MSA` `phylogenetics` `protein structure` `end-to-end`

---

#### 3.09 Molecular Modeling and Docking

[01_molecular_modeling_and_docking.ipynb](Tier_3_Applied_Bioinformatics/09_Molecular_Modeling_and_Docking/01_molecular_modeling_and_docking.ipynb) -- 71 cells

Bridging static structural data and dynamic biological reality. Levels of molecular modeling theory: quantum mechanics, molecular mechanics, coarse-grained models. Force fields for biomolecules: AMBER, CHARMM, OPLS -- bonded and non-bonded terms. Energy minimization algorithms: steepest descent, conjugate gradient. Molecular dynamics simulation concepts: integration, thermostats, periodic boundaries. Molecular docking: rigid and flexible docking, scoring functions. Hands-on docking with AutoDock Vina. Virtual screening workflows. Adapted from Kodomo Semester 8 material.

`force fields` `energy minimization` `molecular dynamics` `docking` `AutoDock Vina` `virtual screening`

---

#### 3.10 Deep Learning for Biology

[01_deep_learning_for_biology.ipynb](Tier_3_Applied_Bioinformatics/10_Deep_Learning_for_Biology/01_deep_learning_for_biology.ipynb) -- 61 cells

When deep learning outperforms classical ML: image data, raw sequences, large-scale datasets. Neural network fundamentals: perceptrons, activation functions, backpropagation, gradient descent. Convolutional Neural Networks (CNNs) for genomic sequences and protein structure. Recurrent Neural Networks (RNNs) and LSTMs for sequential biological data. Transformer architectures and attention mechanisms. Protein language models (ESM, ProtTrans). AlphaFold and the protein structure prediction revolution. Variational autoencoders for single-cell data. Practical deep learning with PyTorch.

`neural networks` `CNNs` `transformers` `AlphaFold` `protein language models` `PyTorch` `autoencoders`

---

#### 3.11 Clinical Genomics

[01_clinical_genomics.ipynb](Tier_3_Applied_Bioinformatics/11_Clinical_Genomics/01_clinical_genomics.ipynb) -- 41 cells

Precision medicine and the journey from the Human Genome Project to routine clinical sequencing. Types of genetic testing: diagnostic, predictive, carrier, pharmacogenomic. ACMG/AMP variant classification: pathogenic, likely pathogenic, VUS, likely benign, benign. ClinVar and variant databases. Clinical sequencing approaches: gene panels, whole-exome sequencing (WES), whole-genome sequencing (WGS). Pharmacogenomics: drug-gene interactions, CYP enzymes. Ethical considerations in clinical genomics. Adapted from Kodomo Semester 9 material by V.E. Ramensky and A.A. Zharikova.

`ACMG` `ClinVar` `pharmacogenomics` `WES` `WGS` `variant classification` `precision medicine`

---

#### 3.12 Modern Bioinformatics Workflows

[01_single_cell_scanpy.ipynb](Tier_3_Applied_Bioinformatics/12_Modern_Workflows/01_single_cell_scanpy.ipynb) -- 45 cells
[02_workflow_engines.ipynb](Tier_3_Applied_Bioinformatics/12_Modern_Workflows/02_workflow_engines.ipynb) -- 23 cells
[03_testing_cicd.ipynb](Tier_3_Applied_Bioinformatics/12_Modern_Workflows/03_testing_cicd.ipynb) -- 40 cells

Contemporary tools for professional bioinformatics. Single-cell RNA-seq analysis with Scanpy: AnnData structure, quality control, normalization, dimensionality reduction (PCA/UMAP), Leiden clustering, marker gene detection. Workflow engines: Snakemake rules/wildcards/config/conda/containers, SLURM/cloud cluster execution; Nextflow DSL2 processes/channels; nf-core community pipelines (rnaseq, sarek, scrnaseq), samplesheet format, institutional configs, module installation. Testing and CI/CD: pytest fixtures and parametrized tests, GitHub Actions, code coverage, linting with ruff/black.

`single-cell` `scanpy` `Snakemake` `Nextflow` `nf-core` `DSL2` `SLURM` `cluster` `pytest` `GitHub Actions` `CI/CD`

---

#### 3.13 Biochemistry and Enzyme Kinetics

[01_biochemistry_and_enzyme_kinetics.ipynb](Tier_3_Applied_Bioinformatics/13_Biochemistry_and_Enzyme_Kinetics/01_biochemistry_and_enzyme_kinetics.ipynb) -- 54 cells

Enzyme kinetics from first principles to curve fitting. Beer-Lambert law and spectrophotometric assays. Michaelis-Menten equation: fitting Km and Vmax with scipy.optimize. Linearization methods (Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf) and their limitations. Enzyme inhibition: competitive, non-competitive, uncompetitive, and mixed — model fitting and type determination via AIC. Allosteric enzymes and the Hill equation. EC classification system. Metabolic pathway computation: stoichiometric matrices and flux balance concepts.

`Michaelis-Menten` `enzyme inhibition` `Hill equation` `EC numbers` `curve fitting` `stoichiometric matrix` `Beer-Lambert`

---

#### 3.14 Genetic Engineering In Silico

[01_genetic_engineering_in_silico.ipynb](Tier_3_Applied_Bioinformatics/14_Genetic_Engineering_In_Silico/01_genetic_engineering_in_silico.ipynb) -- 56 cells

Computational tools for molecular cloning and genome editing. In silico restriction digestion: multi-enzyme digests, fragment prediction, gel electrophoresis simulation. PCR primer design: nearest-neighbor Tm calculation, GC constraints, dimer and hairpin checking, adding restriction tails. Cloning workflows: insert/vector compatibility, reading frame verification, Gateway cloning. CRISPR guide RNA design: PAM identification for SpCas9/SaCas9/Cas12a, on-target scoring, off-target assessment. Codon optimization for E. coli, yeast, and mammalian expression. Gibson Assembly overlap design.

`restriction digests` `primer design` `Tm calculation` `CRISPR` `guide RNA` `codon optimization` `Gibson Assembly` `cloning`

---

#### 3.15 Population Genetics and Molecular Evolution

[01_population_genetics_and_molecular_evolution.ipynb](Tier_3_Applied_Bioinformatics/15_Population_Genetics_and_Molecular_Evolution/01_population_genetics_and_molecular_evolution.ipynb) -- 61 cells

Population-level processes and molecular evolution. Hardy-Weinberg equilibrium: allele frequency estimation, chi-squared testing, multi-allelic extension. Genetic drift: Wright-Fisher simulation, fixation probability, effective population size, bottleneck effects. Natural selection: directional, balancing, and purifying selection; frequency dependence; selection-drift interplay. Molecular clock: Jukes-Cantor correction, divergence time estimation. dN/dS ratio: Nei-Gojobori method, sliding window analysis, interpreting selection signals. Neutrality tests: Tajima's D, McDonald-Kreitman test. Population structure: Fst computation, linkage disequilibrium decay.

`genetic drift` `Wright-Fisher` `natural selection` `dN/dS` `molecular clock` `Tajima's D` `Fst` `linkage disequilibrium`

---

#### 3.16 Numerical Methods for Bioinformatics

[01_numerical_methods_for_bioinformatics.ipynb](Tier_3_Applied_Bioinformatics/16_Numerical_Methods_for_Bioinformatics/01_numerical_methods_for_bioinformatics.ipynb) -- 45 cells

Numerical algorithms essential for bioinformatics data analysis. Polynomial interpolation (Lagrange, Newton) for missing data points. Cubic splines for smooth curve fitting. Numerical differentiation and integration (trapezoidal, Simpson's rules) for dose-response AUC. Nonlinear least squares curve fitting with scipy.optimize.curve_fit for Michaelis-Menten and Hill equations. Gradient descent and optimization methods. Fourier transform and FFT for circadian rhythm detection and spectral analysis. Every technique demonstrated with biological applications. Adapted from ФББ Semester 7 numerical methods lectures.

`interpolation` `splines` `curve fitting` `least squares` `gradient descent` `FFT` `optimization` `AUC`

---

#### 3.17 Genome Assembly and Advanced NGS

[01_genome_assembly_and_advanced_ngs.ipynb](Tier_3_Applied_Bioinformatics/17_Genome_Assembly_and_Advanced_NGS/01_genome_assembly_and_advanced_ngs.ipynb) -- 30 cells

De novo genome assembly from short and long reads. Overlap-Layout-Consensus (OLC) approach for long reads -- implementing a greedy assembler in Python. De Bruijn graph approach for short reads -- k-mer decomposition and Eulerian path finding. Assembly quality metrics: N50, L50, NG50 calculated from scratch, BUSCO completeness, QUAST statistics. Scaffolding with paired-end, mate-pair, and Hi-C data. Read mapping algorithms in depth: Burrows-Wheeler Transform (BWT) construction and FM-index backward search implemented in Python. Advanced QC: k-mer spectra with GenomeScope, contamination screening. Long-read technologies deep dive: ONT vs PacBio HiFi, hybrid assembly strategies. Adapted from ФББ Semester 9 NGS lectures by Logacheva et al.

`de novo assembly` `de Bruijn graph` `OLC` `N50` `BUSCO` `BWT` `FM-index` `scaffolding` `long reads`

---

#### 3.18 Proteomics and Structural Methods

[01_proteomics_and_structural_methods.ipynb](Tier_3_Applied_Bioinformatics/18_Proteomics_and_Structural_Methods/01_proteomics_and_structural_methods.ipynb) -- 48 cells

Mass spectrometry-based proteomics and structural biology methods. Ionization (MALDI, ESI), mass analyzers (TOF, Orbitrap), and MS/MS fragmentation with b/y ion series -- implementing peptide mass and fragment calculators in Python. Bottom-up proteomics: in silico trypsin digestion, peptide mass fingerprinting, database searching concepts (Mascot, MaxQuant), target-decoy FDR. Quantitative proteomics: label-free (LFQ, iBAQ), isotope labeling (SILAC, TMT), DDA vs DIA. Protein engineering computational design: rational design, directed evolution library design, conservation-based mutability scoring from MSA, stability predictions (ΔΔG). Structural determination methods: X-ray crystallography (Bragg's law, R-factor, electron density maps, crystal lattices and the 7 crystal systems, space groups and asymmetric units), cryo-EM (single-particle analysis), NMR (chemical shifts, NOEs). Includes CRYST1 record parser, symmetry operation application, and Uppsala EDS density download. Adapted from ФББ Semester 9 physical-chemical methods and protein engineering lectures by Suplatov; crystallography sections ported from Kodomo archive (Lunin/IMPB RAS).

`mass spectrometry` `proteomics` `MS/MS` `peptide identification` `FDR` `TMT` `protein engineering` `X-ray crystallography` `electron density` `crystal lattice` `space group` `cryo-EM`

---

#### 3.19 Genome-Wide Association Studies (GWAS)

[19_gwas.ipynb](Tier_3_Applied_Bioinformatics/19_GWAS/19_gwas.ipynb)

GWAS study design from first principles. Case/control phenotype definition and confounder identification. Quality control: SNP and sample filtering, Hardy-Weinberg equilibrium testing, MAF thresholds. Population stratification detection via PCA on genotype data. Association testing using logistic and linear regression per SNP. Multiple testing correction with the genome-wide significance threshold (5×10⁻⁸). Manhattan and QQ plot generation from scratch. Linkage disequilibrium and clumping concepts. Downstream: fine-mapping principles and GWAS catalog lookup. Uses 1000 Genomes public data or simulated genotype matrices.

`GWAS` `population stratification` `PCA` `Manhattan plot` `QQ plot` `LD` `clumping` `fine-mapping` `5e-8`

---

#### 3.20 Spatial Transcriptomics

[20_spatial_transcriptomics.ipynb](Tier_3_Applied_Bioinformatics/20_Spatial_Transcriptomics/20_spatial_transcriptomics.ipynb)

Spatial gene expression analysis with Squidpy and Scanpy. AnnData structure with spatial coordinates; Visium and Xenium layout conventions. Quality control for spatial data: mitochondrial fraction and spot-level filtering. Normalization and dimensionality reduction in spatial context. Spatial neighborhood graph construction. Spatially variable gene detection. Cell-type deconvolution concepts (RCTD, cell2location patterns). Visualization: spatial scatter plots and expression overlays on tissue sections. Uses the public 10x Visium mouse brain dataset available through Squidpy.

`spatial transcriptomics` `Squidpy` `AnnData` `Visium` `spatially variable genes` `deconvolution` `neighborhood graph` `scanpy`

---

#### 3.21 DNA Copy Number Analysis

[21_copy_number_analysis.ipynb](Tier_3_Applied_Bioinformatics/21_Copy_Number_Analysis/21_copy_number_analysis.ipynb)

Copy number variation analysis from sequencing data. CNV concepts: gains, losses, and loss of heterozygosity. Read depth normalization approaches across genomic windows. Segmentation using the Circular Binary Segmentation (CBS) algorithm concept. Copy number state calling from segments. Genome-wide CN profile visualization. Gene-level annotation of CN events. Extends the variant calling pipeline concepts from Module 3.02.

`copy number` `CNV` `CBS segmentation` `read depth normalization` `LOH` `genome-wide profile` `somatic variants`

---

#### 3.22 Bayesian Statistics in Python

[22_bayesian_statistics_python.ipynb](Tier_3_Applied_Bioinformatics/22_Bayesian_Statistics_Python/22_bayesian_statistics_python.ipynb)

Bayesian statistical modeling in Python across seven sub-sections. Frequentist vs. Bayesian framing: posterior = likelihood × prior, credible intervals. Prior specification: informative vs. weakly informative priors, prior predictive checks. Multiple regression: collinearity diagnostics and variance inflation factor. Model comparison with WAIC and LOO-CV using ArviZ. Linear mixed-effects models with random intercepts and slopes using Bambi. GLMs in a Bayesian framework: Bernoulli, Binomial, Poisson, Negative Binomial with PyMC. Advanced: GLMM, zero-inflated models, GAM concepts, and Bayesian meta-analysis. Converted from Fränzi Korner-Nievergelt's applied statistics R course; uses the public palmerpenguins dataset.

`pymc` `bambi` `arviz` `LOO-CV` `WAIC` `credible intervals` `GLM` `mixed-effects` `prior predictive` `zero-inflated`

---

#### 3.23 TF Footprinting & Chromatin Accessibility

[23_tf_footprinting.ipynb](Tier_3_Applied_Bioinformatics/23_TF_Footprinting/23_tf_footprinting.ipynb)

Transcription factor footprinting from ATAC-seq data. ATAC-seq recap: fragment size distribution and nucleosome-free region identification. TF footprinting concept: Tn5 insertion bias around motif binding sites. Expected vs. observed cut-site profiles around motif centers. Footprint score calculation and interpretation across conditions. Genomic interval arithmetic with pybedtools: intersection, subtraction, closest-feature queries. Accumulation plots: average signal enrichment around genomic features. Extends the ngs-variant-calling and motif-discovery skills; uses public ENCODE ATAC-seq data.

`ATAC-seq` `TF footprinting` `Tn5 bias` `cut-site profiles` `pybedtools` `genomic intervals` `accumulation plots` `chromatin accessibility`

---

### Tier 4: Algorithms & Data Structures -- 30 notebooks, 927 cells

The computer science foundation behind bioinformatics tools. Study alongside Tiers 2--3 to understand *why* the algorithms work. Includes 30 interactive HTML5 visualizations and animated GIF demonstrations. See the [Tier 4 README](Tier_4_Algorithms_and_Data_Structures/README.md) for the full bioinformatics cross-reference table.

---

#### 4.00 Skills Check: Algorithms

[00_skills_check.ipynb](Tier_4_Algorithms_and_Data_Structures/00_Skills_Check/00_skills_check.ipynb) -- 24 cells

Self-assessment covering Big-O notation, sorting, searching, basic data structures, and recursion.

`self-assessment` `complexity` `data structures` `skip-ahead`

---

#### 4.01 Complexity Analysis

[01_complexity_analysis.ipynb](Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/01_complexity_analysis.ipynb) -- 34 cells
[02_basic_algorithms.ipynb](Tier_4_Algorithms_and_Data_Structures/01_Complexity_Analysis/02_basic_algorithms.ipynb) -- 52 cells

Big-O, Big-Omega, Big-Theta notation. Analyzing time and space complexity. Recurrence relations. Recursion fundamentals and basic algorithmic paradigms. **Bio connection:** choosing efficient tools for processing large genomes.

`Big-O` `time complexity` `space complexity` `recurrences` `recursion`

---

#### 4.02 Sorting Algorithms

[01_comparison_sorts.ipynb](Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/01_comparison_sorts.ipynb) -- 49 cells
[02_linear_sorts.ipynb](Tier_4_Algorithms_and_Data_Structures/02_Sorting_Algorithms/02_linear_sorts.ipynb) -- 24 cells

Bubble, Selection, Insertion, Shell, Merge, and QuickSort with step-by-step visualization. Linear-time sorts: Counting, Radix, Bucket. **Bio connection:** BAM coordinate sorting, variant prioritization.

`bubble sort` `merge sort` `quicksort` `counting sort` `radix sort`

---

#### 4.03 Searching Algorithms

[01_linear_binary_search.ipynb](Tier_4_Algorithms_and_Data_Structures/03_Searching_Algorithms/01_linear_binary_search.ipynb) -- 50 cells

Linear search, binary search and its variants, two-pointer techniques. **Bio connection:** searching sorted genomic intervals, bisecting quality thresholds.

`linear search` `binary search` `two-pointer`

---

#### 4.04 Linear Data Structures

[01_linked_lists.ipynb](Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/01_linked_lists.ipynb) -- 44 cells
[02_stacks_queues.ipynb](Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/02_stacks_queues.ipynb) -- 30 cells
[03_dynamic_arrays.ipynb](Tier_4_Algorithms_and_Data_Structures/04_Linear_Data_Structures/03_dynamic_arrays.ipynb) -- 24 cells

Singly, doubly, and circular linked lists. Stacks (LIFO), queues (FIFO), deques, priority queues. Dynamic array growth and amortized analysis.

`linked lists` `stacks` `queues` `deques` `amortized analysis`

---

#### 4.05 Tree Structures

[01_binary_search_trees.ipynb](Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/01_binary_search_trees.ipynb) -- 43 cells
[02_avl_trees.ipynb](Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/02_avl_trees.ipynb) -- 28 cells
[03_red_black_trees.ipynb](Tier_4_Algorithms_and_Data_Structures/05_Tree_Structures/03_red_black_trees.ipynb) -- 23 cells

BST insert, delete, search, and traversals (inorder, preorder, postorder, level-order). AVL self-balancing with rotations. Red-black tree properties and rebalancing. **Bio connection:** phylogenetic tree traversal, interval trees for genomic annotations.

`BST` `AVL` `red-black` `rotations` `tree traversal`

---

#### 4.06 Hash-Based Structures

[01_hash_tables_bloom.ipynb](Tier_4_Algorithms_and_Data_Structures/06_Hash_Based_Structures/01_hash_tables_bloom.ipynb) -- 20 cells

Hash functions, collision resolution (chaining, open addressing), load factors. Bloom filters for probabilistic membership testing. **Bio connection:** k-mer counting (Jellyfish), read deduplication, contamination screening.

`hash tables` `collision resolution` `Bloom filters` `k-mer counting`

---

#### 4.07 String Algorithms

[01_naive_pattern_matching.ipynb](Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/01_naive_pattern_matching.ipynb) -- 34 cells
[02_kmp_algorithm.ipynb](Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/02_kmp_algorithm.ipynb) -- 26 cells
[03_rabin_karp.ipynb](Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/03_rabin_karp.ipynb) -- 28 cells
[04_dfa_matching.ipynb](Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/04_dfa_matching.ipynb) -- 35 cells

Brute-force string matching. Knuth-Morris-Pratt with failure function. Rabin-Karp rolling hash for multiple pattern matching. DFA-based matching. **Bio connection:** the seed-and-extend strategy behind BLAST.

`KMP` `Rabin-Karp` `rolling hash` `pattern matching` `BLAST internals`

---

#### 4.08 Advanced String Structures

[01_tries.ipynb](Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/01_tries.ipynb) -- 30 cells
[02_aho_corasick.ipynb](Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/02_aho_corasick.ipynb) -- 27 cells
[03_suffix_arrays.ipynb](Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/03_suffix_arrays.ipynb) -- 36 cells
[04_suffix_trees.ipynb](Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/04_suffix_trees.ipynb) -- 34 cells

Prefix trees for autocomplete and spell-checking. Aho-Corasick multi-pattern automaton. Suffix arrays with LCP arrays. Suffix trees and Ukkonen's algorithm. **Bio connection:** genome indexing (BWT, FM-index), de Bruijn graphs for assembly, multi-motif scanning.

`tries` `Aho-Corasick` `suffix arrays` `suffix trees` `genome indexing`

---

#### 4.09 Graph Algorithms

[01_graph_representations.ipynb](Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/01_graph_representations.ipynb) -- 25 cells
[02_bfs_dfs.ipynb](Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/02_bfs_dfs.ipynb) -- 17 cells
[03_dijkstra.ipynb](Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/03_dijkstra.ipynb) -- 15 cells
[04_mst_kruskal_prim.ipynb](Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/04_mst_kruskal_prim.ipynb) -- 18 cells
[05_topological_sort.ipynb](Tier_4_Algorithms_and_Data_Structures/09_Graph_Algorithms/05_topological_sort.ipynb) -- 15 cells

Graph representations: adjacency matrix, adjacency list, and edge list with complexity trade-offs. Breadth-first and depth-first traversal with step-by-step tracing. Dijkstra's single-source shortest path with a priority queue; Bellman-Ford for negative weights; Floyd-Warshall all-pairs. Minimum spanning trees: Kruskal's (Union-Find) and Prim's (priority queue). Topological sort and cycle detection in DAGs. **Bio connection:** metabolic pathway analysis, gene regulatory networks, protein interaction graphs.

`BFS` `DFS` `Dijkstra` `MST` `topological sort` `pathways`

---

#### 4.10 Dynamic Programming

[01_intro_memoization.ipynb](Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/01_intro_memoization.ipynb) -- 24 cells
[02_tabulation.ipynb](Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/02_tabulation.ipynb) -- 19 cells
[03_knapsack.ipynb](Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/03_knapsack.ipynb) -- 17 cells
[04_sequence_alignment.ipynb](Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/04_sequence_alignment.ipynb) -- 26 cells

Memoization vs. tabulation: top-down vs. bottom-up with complexity comparison. Foundational problems: Fibonacci, climbing stairs, coin change. Sequence DP: longest common subsequence, longest increasing subsequence, edit distance. 0/1 knapsack and subset sum. Needleman-Wunsch and Smith-Waterman sequence alignment implemented as DP algorithms with traceback reconstruction. **Bio connection:** global and local sequence alignment are DP; RNA secondary structure prediction uses DP.

`memoization` `tabulation` `LCS` `edit distance` `sequence alignment` `Needleman-Wunsch` `Smith-Waterman`

---

### Tier 5: Modern AI for Science -- 3 notebooks

GPU-optional modules covering contemporary AI methods for scientific research. Each notebook is designed to run on free-tier Google Colab. Theory cells run without GPU; hands-on training cells require a T4 or better. See the [Tier 5 README](Tier_5_Modern_AI_for_Science/README.md) for GPU setup and Colab instructions.

---

#### 5.01 LLM Fine-tuning

[01_LLM_Finetuning.ipynb](Tier_5_Modern_AI_for_Science/01_LLM_Finetuning/01_LLM_Finetuning.ipynb)

Fine-tuning large language models for domain-specific instruction following. Base vs. instruction/chat models: what changes during fine-tuning and why. LoRA: low-rank adapter mathematics, rank selection, and target module identification. Quantization: 4-bit NF4 with bitsandbytes and trade-offs with output quality. Chat template formatting: system/user/assistant structure. SFTTrainer workflow: dataset preparation, training loop configuration, and evaluation. Synthetic data generation for instruction tuning. Practical tips: gradient checkpointing, batch size scheduling, learning rate warm-up. Inspired by Unsloth AI and Manuel Faysse fine-tuning patterns; uses public instruction datasets from HuggingFace Hub.

`LoRA` `quantization` `NF4` `SFTTrainer` `trl` `peft` `bitsandbytes` `chat templates` `instruction tuning` `synthetic data`

---

#### 5.02 Vision RAG

[02_Vision_RAG.ipynb](Tier_5_Modern_AI_for_Science/02_Vision_RAG/02_Vision_RAG.ipynb)

Retrieval-augmented generation with vision-language models for document understanding. VLM architecture overview: visual encoder + LLM decoder. Document understanding: page-level vs. token-level retrieval approaches. ColPali: late-interaction document retrieval concept and API pattern. RAG pipeline: retrieval → context injection → generation. Qwen2-VL inference pattern for multi-page document Q&A. Evaluation: retrieval recall and generation faithfulness metrics. Uses public PDF documents (arXiv papers, open-access reports). Inspired by Unsloth AI and Manuel Faysse patterns.

`VLM` `ColPali` `late interaction` `RAG` `document retrieval` `Qwen2-VL` `retrieval recall` `generation faithfulness`

---

#### 5.03 Diffusion & Generative Models

[03_Diffusion_Generative_Models.ipynb](Tier_5_Modern_AI_for_Science/03_Diffusion_Generative_Models/03_Diffusion_Generative_Models.ipynb)

Score-based generative models and diffusion for scientific imaging applications. Score matching intuition: score field as ∇ₓ log p(x), learned denoising. DDIM: deterministic sampling, linear and cosine noise schedules, reverse process equations. Inverse problems in imaging: denoising, inpainting, and colorization as special cases of DDRM. SVD-based degradation operators and the pseudoinverse projection step. Linear and cosine scheduler implementation patterns with explicit tensor shapes. Score field visualization with quiver plots. Scientific applications: cryo-EM denoising and medical image restoration concepts. Adapted from DDRM (Kawar et al., 2022), github.com/bahjat-kawar/ddrm; uses public MNIST/CIFAR-10 data.

`diffusion models` `DDIM` `score matching` `noise schedule` `inverse problems` `DDRM` `SVD` `cryo-EM` `image restoration`

---

## Skills Check Guide

Each of the five tiers opens with a **Skills Check** notebook -- a self-graded assessment that helps you decide whether to work through the tier or skip ahead.

### How Skills Checks Work

1. **Open the Skills Check** for the tier you think you might be ready for.
2. **Answer every question honestly** without looking at the answer key first.
3. **Grade yourself** using the detailed answer key at the end of the notebook.
4. **Interpret your score:**

| Score | Recommendation |
|-------|---------------|
| **Above 80%** | You can safely skip this tier. Proceed to the next one. |
| **60--80%** | Skim the notebooks for topics where you lost points. |
| **Below 60%** | Work through the tier's notebooks carefully. The foundation pays off later. |

### Skills Check Locations

| Tier | Notebook | Questions | Focus Areas |
|------|----------|-----------|-------------|
| **Tier 0** | [00_skills_check.ipynb](Tier_0_Computational_Foundations/00_Skills_Check/00_skills_check.ipynb) | 15 | Linux commands, Git workflows, Bash scripting |
| **Tier 1** | [00_skills_check.ipynb](Tier_1_Python_for_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) | varies | Python variables, data types, control flow, functions, data structures |
| **Tier 2** | [00_skills_check.ipynb](Tier_2_Core_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) | 15 | BLAST, alignment, phylogenetics, protein structure, GO |
| **Tier 3** | [00_skills_check.ipynb](Tier_3_Applied_Bioinformatics/00_Skills_Check/00_skills_check.ipynb) | varies | NGS, sequencing platforms, RNA-seq, variant calling |
| **Tier 4** | [00_skills_check.ipynb](Tier_4_Algorithms_and_Data_Structures/00_Skills_Check/00_skills_check.ipynb) | 24 cells | Big-O notation, sorting, searching, basic data structures, recursion |

**Recommended approach for uncertain students:** Start with the Tier 1 Skills Check. Most students with any programming background can begin at Tier 1. If you score well there, try the Tier 2 Skills Check next.

---

## Sample Data

All data files live in `Assets/data/` and are small enough to keep in the repository. Notebooks reference these files with relative paths.

### Sequence Files

| File | Size | Description | Used In |
|------|------|-------------|---------|
| `Saccharomyces_cerevisiae.fasta` | 754 B | Yeast gene nucleotide sequence | BioPython Essentials, Pairwise Alignment |
| `Bacillus_cereus.fasta` | 302 B | Bacterial gene nucleotide sequence | BioPython Essentials, Pairwise Alignment |
| `Saccharomyces_COX2.fasta` | 19 KB | Yeast COX2 gene with multiple entries for alignment | Multiple Sequence Alignment, Phylogenetics |
| `Mouse_Mus_cytb.fasta` | 277 KB | Mouse cytochrome b sequences from multiple specimens | Phylogenetics, Comparative Genomics |
| `Escherichia_coli_proteome.fasta` | 1.9 MB | Complete E. coli proteome (all protein sequences) | BLAST Searching, Sequence Motifs and Domains |

### Structure and Matrix Files

| File | Size | Description | Used In |
|------|------|-------------|---------|
| `sample_protein.pdb` | 12 KB | Small PDB structure file for analysis exercises | Protein Structure, Molecular Modeling and Docking |
| `EBLOSUM62.txt` | 2.1 KB | BLOSUM62 substitution matrix in EMBOSS format | Pairwise Alignment, BLAST Searching |

### Chromatogram Files

| File | Size | Description | Used In |
|------|------|-------------|---------|
| `WS1945_COI_F_D11_WSBS-Seq-4-08-15.ab1` | 224 KB | Forward-read Sanger chromatogram for COI barcoding gene | Chromatogram Analysis |
| `WS1945_COI_R_D12_WSBS-Seq-4-08-15.ab1` | 224 KB | Reverse-read Sanger chromatogram for COI barcoding gene | Chromatogram Analysis |

### Annotation and Variant Files

| File | Size | Description | Used In |
|------|------|-------------|---------|
| `sequence.gb` | 50 KB | GenBank-format annotated genomic sequence | Biological Databases, BioPython Essentials |
| `snp.vcf` | 31 KB | VCF file with SNP variant calls | Variant Calling and SNP Analysis, NGS Fundamentals |

---

## Learning Path Recommendations

### Path 1: The Complete Journey

**Tiers 0 --> 1 --> 2 --> 3** | ~120 hours | For complete beginners

Start from the very beginning and work through everything. This path builds the strongest foundation. Even experienced programmers occasionally find useful insights in Tier 0 (especially the Biostatistics and Bash Scripting modules).

```
Tier 0 (20h) --> Tier 1 (40h) --> Tier 2 (35h) --> Tier 3 (25h) --> Capstone
```

### Path 2: Python for Biologists

**Tiers 0 + 1** | ~50 hours | For biologists learning to code

Focus on computational fundamentals and Python programming. After completing these two tiers, you will be able to write scripts for data processing, parse biological files, create figures, and work with pandas DataFrames. A strong foundation for any bioinformatics work.

```
Tier 0 (15h) --> Tier 1 (35h)
```

### Path 3: Bioinformatics Core

**Tiers 2 + 3** | ~60 hours | For programmers entering bioinformatics

Skip the programming tiers (take the Skills Checks to confirm) and dive straight into biological databases, sequence analysis, structural biology, and applied pipelines. This is the fastest path to productive bioinformatics work for someone who already knows Python.

```
Tier 2 Skills Check --> Tier 2 (35h) --> Tier 3 (25h) --> Capstone
```

### Path 4: Applied Specialist

**Tier 3 selectively** | ~15--25 hours | For experienced bioinformaticians

Cherry-pick individual Tier 3 modules based on your needs. Each module is relatively self-contained:

| If you need... | Take... |
|----------------|---------|
| NGS pipeline skills | 3.01 NGS Fundamentals + 3.02 Variant Calling |
| Transcriptomics | 3.03 RNA-seq Analysis |
| Metagenomics | 3.04 Microbial Diversity |
| Regulatory genomics | 3.05 Promoter and Regulatory Analysis + 3.23 TF Footprinting |
| Statistical methods | 3.06 Statistics for Bioinformatics + 3.22 Bayesian Statistics |
| ML/DL for biology | 3.07 Machine Learning + 3.10 Deep Learning |
| Structural/drug design | 3.09 Molecular Modeling and Docking |
| Clinical applications | 3.11 Clinical Genomics |
| GWAS and population genetics | 3.19 GWAS + 3.15 Population Genetics |
| Single-cell and spatial omics | 3.12 Modern Workflows + 3.20 Spatial Transcriptomics |
| Cancer genomics | 3.21 Copy Number Analysis + 3.02 Variant Calling |

---

## Prerequisites and Installation

### Hardware Requirements

Any computer with 4+ GB RAM. Most notebooks run on modest hardware. Deep Learning (3.10) benefits from a GPU but includes CPU-compatible exercises.

### Software Prerequisites

- **Python 3.9 or later** (3.10+ recommended)
- **Jupyter Notebook or JupyterLab**
- A terminal/command line (all operating systems supported)

### Quick Start

```bash
# Clone the repository
git clone <repository-url>
cd Course

# Option A: pip (virtual environment)
python3 -m venv bioinfo_course
source bioinfo_course/bin/activate    # Linux / macOS
# bioinfo_course\Scripts\activate     # Windows

pip install jupyter numpy pandas matplotlib seaborn
pip install biopython
pip install scikit-learn scipy statsmodels
pip install py3Dmol requests

# Start Jupyter
jupyter notebook
```

### Conda Environment (recommended for full course)

```bash
# Create environment with all dependencies
conda create -n bioinfo python=3.11 jupyter numpy pandas matplotlib seaborn \
    biopython scikit-learn scipy statsmodels requests -y

conda activate bioinfo

# Additional packages
pip install py3Dmol

# Bioinformatics command-line tools (for Tier 2/3)
conda install -c bioconda muscle blast samtools fastqc -y
```

### Package Summary

| Package | Used In | Purpose |
|---------|---------|---------|
| `jupyter` | All tiers | Interactive notebook environment |
| `numpy` | Tier 1--3 | Numerical arrays, linear algebra |
| `pandas` | Tier 1--3 | DataFrames, tabular data manipulation |
| `matplotlib` | Tier 1--3 | Plotting and figure creation |
| `seaborn` | Tier 1--3 | Statistical visualizations, heatmaps |
| `biopython` | Tier 2--3 | Biological sequences, file parsing, NCBI access |
| `scikit-learn` | Tier 3 | Machine learning algorithms |
| `scipy` | Tier 1--3 | Statistical tests, scientific computing |
| `statsmodels` | Tier 3 | Advanced statistics, survival analysis |
| `py3Dmol` | Tier 2--3 | 3D molecular visualization in notebooks |
| `requests` | Tier 2--3 | HTTP requests to biological databases |

### External Command-Line Tools (optional)

Most Tier 2 notebooks offer both command-line and pure-Python approaches, so these are not strictly required.

| Tool | Used In | Installation |
|------|---------|-------------|
| MUSCLE | Multiple Sequence Alignment | `conda install -c bioconda muscle` |
| ClustalW | Multiple Sequence Alignment | `conda install -c bioconda clustalw` |
| BLAST+ | BLAST Searching, Capstone | `conda install -c bioconda blast` |
| FastQC | NGS Fundamentals | `conda install -c bioconda fastqc` |
| samtools | NGS, Variant Calling | `conda install -c bioconda samtools` |

---

## Sources and Acknowledgments

This course was built on the foundation of established bioinformatics education programs. We are deeply grateful to all instructors, authors, and maintainers who made high-quality bioinformatics education freely available.

### Kodomo Bioinformatics Program

**Website:** [https://kodomo.fbb.msu.ru/wiki/2017](https://kodomo.fbb.msu.ru/wiki/2017)

The Kodomo curriculum is a comprehensive 10-semester program at the **Faculty of Bioengineering and Bioinformatics, Moscow State University (MSU)**. It provided the core structure, scientific depth, and many of the exercises for this course. The original materials were in Russian and have been translated to English and adapted for self-paced learning.

#### Faculty and Instructors

| Instructor | Areas of Expertise | Semesters |
|---|---|---|
| **A.V. Golovin** | Nucleic acid structures, molecular modeling | 3, 8 |
| **S.A. Spirin** | Sequence analysis, BLAST, phylogenetics, alignment | 1--4 |
| **A.V. Alekseevsky** | Various foundational bioinformatics topics | 1--2 |
| **A. Zalevsky** | Protein structure visualization, introduction to bioinformatics | 1--2 |
| **A.S. Zlobin** | 3D bioinformatics, structural alignment, molecular recognition | 7 |
| **D. Penzar** | Python programming for bioinformatics | 1, 4 |
| **Z. Chervontseva** | Position weight matrices, PSSM, PSI-BLAST, Python | 4 |
| **I. Rusinov** | Linux, EMBOSS, data processing | 1--2 |
| **A. Zharikova** | Enzymes, KEGG, SNP analysis, transcriptomics, medical genomics | 3--4, 9 |
| **V.E. Ramensky** | Medical genomics | 9 |
| **V.Yu. Lunin** (IMPB RAS) | X-ray crystallography | 7 |
| **K.S. Mineev** (IBCh RAS) | NMR spectroscopy | 7 |
| **O.S. Sokolova** (MSU) | Cryo-electron microscopy | 7 |
| **V.D. Maslova** | Structure validation | 7 |
| **M. Khachaturyan** | De novo genome assembly | 3 |
| **D. Dibrova** | Web technologies for bioinformatics | 1 |
| **R. Kudrin** | Python programming | 4 |
| **I. Diankin** | Python programming | 4 |
| **E. Ocheredko** | Sanger sequencing | 3 |
| **A. Demkiv** | Protein contact analysis | 2 |
| **A. Ershova** | Course organization | 1 |

#### Kodomo Semester Topics

| Semester | Topic |
|---|---|
| 1 | Objects and Tools Introduction |
| 2 | Protein Structures and Sequences |
| 3 | Nucleic Acid Structures and Sequences |
| 4 | Protein Sequence Evolution + Python Programming |
| 5 | Algorithm Fundamentals + R Language and Data Analysis |
| 7 | 3D Bioinformatics |
| 8 | Biopolymer Molecular Modeling |
| 9 | Medical Genomics |
| 10 | Machine Learning in Biology |

Specific Kodomo-derived content in this course:

- **Tier 3, Module 09 (Molecular Modeling and Docking)** -- adapted from Semester 8 (A.V. Golovin)
- **Tier 3, Module 10 (Deep Learning for Biology)** -- adapted from Semester 10
- **Tier 3, Module 11 (Clinical Genomics)** -- adapted from Semester 9 (V.E. Ramensky, A. Zharikova)
- **Archive/07_Kodomo_Structural_Bioinformatics/** -- detailed structural biology exercises from Semesters 1--7

### An Introduction to Applied Bioinformatics (IAB)

**Authors:** J. Gregory Caporaso and collaborators
**Website:** [https://readiab.org/](https://readiab.org/)

IAB is an open-source interactive bioinformatics textbook from the Caporaso Lab at Northern Arizona University. Content adapted from IAB includes:

- Pairwise alignment algorithm implementations and visualizations
- Homology searching concepts and BLAST internals
- Multiple sequence alignment theory
- Phylogenetic tree construction methods
- Sequence clustering approaches
- Microbial diversity analysis concepts

The Archive directory `Archive/06_Applied_Bioinformatics_IAB/` contains adapted IAB materials used as supplementary reference.

### Summer School of Bioinformatics

Teaching methodologies and exercise design adapted from Summer School materials for:

- NGS fundamentals and quality control workflows
- Mathematical statistics for biological data
- Promoter analysis and regulatory element prediction
- Machine learning (Random Forest classifiers) applied to biological classification
- R programming for bioinformatics

### Additional Sources

- **Python lecture materials** from the MSU Bioengineering faculty -- informed the Python programming modules (Tier 1)
- **Sweet_pracks exercise collection** -- contributed practice problems and biological scenarios throughout the course
- **Term 3 practicals** -- nucleic acid and protein sequence analysis exercises
- **Bioinformatics laboratory utilities collection** -- data processing workflows and tool integration patterns

For the complete attribution details, see [CREDITS.md](CREDITS.md).

### Software Acknowledgments

This course relies on the following open-source tools and libraries:

- [BioPython](https://biopython.org/) -- biological computation in Python
- [NumPy](https://numpy.org/), [Pandas](https://pandas.pydata.org/) -- scientific computing and data manipulation
- [Matplotlib](https://matplotlib.org/), [Seaborn](https://seaborn.pydata.org/) -- visualization
- [scikit-learn](https://scikit-learn.org/) -- machine learning
- [SciPy](https://scipy.org/) -- scientific algorithms and statistical tests
- [Jupyter](https://jupyter.org/) -- interactive notebook environment
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) -- sequence similarity searching
- [MUSCLE](https://drive5.com/muscle/) / [ClustalW](http://www.clustal.org/) -- multiple sequence alignment
- [py3Dmol](https://3dmol.csb.pitt.edu/) -- molecular visualization in notebooks

---

## Repository Structure

```
Course/
├── Tier_0_Computational_Foundations/        10 notebooks
│   ├── 00_Skills_Check/
│   ├── 01_Linux_Basics/
│   ├── 02_Git_Version_Control/
│   ├── 03_Bash_Scripting/
│   ├── 04_File_Encodings/
│   ├── 05_R_Basics/
│   ├── 06_Biostatistics/
│   ├── 07_Probability_and_Statistics_Python/
│   └── 08_Advanced_R_Statistics/
│
├── Tier_1_Python_for_Bioinformatics/       20 notebooks
│   ├── 00_Skills_Check/
│   ├── 01_Python_Introduction/
│   ├── 02_Variables_and_Data_Types/
│   ├── 03_Operators_and_Expressions/
│   ├── 04_Strings_and_Sequences/
│   ├── 05_Control_Flow/
│   ├── 06_Functions/
│   ├── 07_File_Operations/
│   ├── 08_Lists_and_Tuples/
│   ├── 09_Dictionaries_and_Sets/
│   ├── 10_Comprehensions/
│   ├── 11_Iterators_and_Generators/
│   ├── 12_Regular_Expressions/
│   ├── 13_Classes_and_OOP/
│   ├── 14_Decorators_and_Context_Managers/
│   ├── 15_Error_Handling/
│   ├── 16_NumPy_and_Pandas/
│   ├── 17_Data_Wrangling/
│   ├── 18_Data_Visualization/
│   └── 19_SQL_for_Bioinformatics/
│
├── Tier_2_Core_Bioinformatics/             17 notebooks
│   ├── 00_Skills_Check/                    (includes Glossary)
│   ├── 01_Biological_Databases/
│   ├── 02_BioPython_Essentials/
│   ├── 03_Pairwise_Sequence_Alignment/
│   ├── 04_BLAST_Searching/
│   ├── 05_Multiple_Sequence_Alignment/
│   ├── 06_Phylogenetics/
│   ├── 07_Protein_Structure/
│   ├── 08_Nucleic_Acid_Structure/
│   ├── 09_Chromatogram_Analysis/
│   ├── 10_Sequence_Motifs_and_Domains/
│   ├── 11_Gene_Ontology_and_Pathways/
│   ├── 12_Comparative_Genomics/
│   ├── 13_Computational_Genetics/
│   ├── 14_Hi-C_Analysis/
│   └── 15_Motif_Discovery/
│
├── Tier_3_Applied_Bioinformatics/          26 notebooks
│   ├── 00_Skills_Check/
│   ├── 01_NGS_Fundamentals/
│   ├── 02_Variant_Calling_and_SNP_Analysis/
│   ├── 03_RNA_seq_Analysis/
│   ├── 04_Microbial_Diversity/
│   ├── 05_Promoter_and_Regulatory_Analysis/
│   ├── 06_Statistics_for_Bioinformatics/
│   ├── 07_Machine_Learning_for_Biology/
│   ├── 08_Capstone_Project/
│   ├── 09_Molecular_Modeling_and_Docking/
│   ├── 10_Deep_Learning_for_Biology/
│   ├── 11_Clinical_Genomics/
│   ├── 12_Modern_Workflows/
│   ├── 13_Biochemistry_and_Enzyme_Kinetics/
│   ├── 14_Genetic_Engineering_In_Silico/
│   ├── 15_Population_Genetics_and_Molecular_Evolution/
│   ├── 16_Numerical_Methods_for_Bioinformatics/
│   ├── 17_Genome_Assembly_and_Advanced_NGS/
│   ├── 18_Proteomics_and_Structural_Methods/
│   ├── 19_GWAS/
│   ├── 20_Spatial_Transcriptomics/
│   ├── 21_Copy_Number_Analysis/
│   ├── 22_Bayesian_Statistics_Python/
│   └── 23_TF_Footprinting/
│
├── Tier_4_Algorithms_and_Data_Structures/   30 notebooks + 30 interactive visualizations
│   ├── 01-Fundamentals/                   Complexity analysis, Big-O
│   ├── 02-Sorting/                        Comparison and linear sorts
│   ├── 03-Searching/                      Binary search, two-pointer
│   ├── 04-Linear-Data-Structures/         Linked lists, stacks, queues
│   ├── 05-Trees/                          BST, AVL, Red-Black trees
│   ├── 06-Hash-Based-Structures/          Hash tables, Bloom filters
│   ├── 07-String-Algorithms/              KMP, Rabin-Karp, DFA matching
│   ├── 08-Advanced-String-Structures/     Tries, Aho-Corasick, suffix trees
│   ├── 09-Graph-Algorithms/               BFS, DFS, Dijkstra, MST
│   ├── 10-Dynamic-Programming/            Memoization, classic DP problems
│   ├── interactive/                       30 HTML5 visualizations with controls
│   ├── assets/gifs/                       8 animated algorithm demonstrations
│   └── legacy/                            10 homework notebooks
│
├── Archive/                                Supplementary reference materials
│   ├── 06_Applied_Bioinformatics_IAB/      IAB textbook adaptations
│   └── 07_Kodomo_Structural_Bioinformatics/ Kodomo structural biology exercises
│
├── Assets/
│   └── data/                               12 sample data files
│
├── README.md                               This file
└── CREDITS.md                              Full attribution details
```

---

## How to Contribute

Contributions that improve the course are welcome:

- **Bug fixes:** Typos, broken code cells, incorrect answers in exercises
- **Clarity improvements:** Better explanations, additional diagrams, biological context
- **New exercises:** Practice problems with biological relevance and worked solutions
- **Data updates:** Newer or more representative sample datasets
- **Translations:** Adapting notebooks for other languages

To contribute, open an issue describing the change you would like to make, or submit a pull request with your improvements. Please include which notebook(s) are affected and test any code changes before submitting.

---

## Disclaimer and License

### Personal Use Only

**This repository is a personal study compilation assembled for private, non-commercial educational purposes only.** It is not intended for redistribution, resale, or commercial use of any kind.

The materials in this course were adapted, translated, and reorganized from publicly available educational resources for the sole purpose of the author's personal learning. **All intellectual property rights for the original materials remain with their respective authors and institutions.**

### Attribution of Original Works

This course draws upon and is deeply indebted to the following works and their authors:

- **Kodomo Bioinformatics Program** -- Faculty of Bioengineering and Bioinformatics, Lomonosov Moscow State University (https://kodomo.fbb.msu.ru/wiki/2017). Original curriculum developed by A.V. Golovin, S.A. Spirin, A.V. Alekseevsky, A. Zalevsky, A.S. Zlobin, D. Penzar, Z. Chervontseva, I. Rusinov, A. Zharikova, V.E. Ramensky, V.Yu. Lunin (IMPB RAS), K.S. Mineev (IBCh RAS), O.S. Sokolova (MSU), V.D. Maslova, M. Khachaturyan, D. Dibrova, R. Kudrin, I. Diankin, E. Ocheredko, A. Demkiv, A. Ershova, and other faculty members. Materials used under fair use for personal study.
- **An Introduction to Applied Bioinformatics (IAB)** -- by J. Gregory Caporaso and collaborators, Caporaso Lab, Northern Arizona University (https://readiab.org/). Used under the terms of the IAB open-source license for educational purposes.
- **Summer School of Bioinformatics** -- teaching materials used for personal study reference.
- **MSU Bioengineering Faculty Python Lectures** -- original lecture materials by faculty instructors.

### No Warranty

This material is provided "as is" without warranty of any kind, express or implied. The author makes no claims about the accuracy, completeness, or suitability of these materials for any particular purpose.

### Fair Use Statement

The use of the above-referenced materials in this personal study compilation constitutes fair use under applicable copyright law, as the materials are:
1. Used for private, non-commercial educational purposes only
2. Transformative in nature (translated, reorganized, and adapted)
3. Not a substitute for the original works
4. Not distributed for profit

**If you are an author or rights holder of any original material referenced here and wish to have your content removed, please contact the repository owner and it will be promptly removed.**

### Your Use of This Repository

If you happen to find this repository useful for your own studies, you are welcome to use it for **personal, non-commercial educational purposes only**, provided you:
1. Do not redistribute or sell these materials
2. Acknowledge the original authors listed above in any derivative work
3. Do not represent this compilation as your own original work

See [CREDITS.md](CREDITS.md) for the complete list of contributors and source materials.

---

## Getting Started

```bash
# Clone the repository
git clone <repository-url>
cd Course

# Install dependencies (see Prerequisites section above for full details)
pip install jupyter numpy pandas matplotlib seaborn biopython scikit-learn scipy

# Start Jupyter
jupyter notebook
```

**Not sure where to begin?** Open the **Tier 1 Skills Check**:

```
Tier_1_Python_for_Bioinformatics/00_Skills_Check/00_skills_check.ipynb
```

Most students with any programming background can start at Tier 1. If you score above 80%, try the Tier 2 Skills Check. If you have never used a command line before, start at Tier 0.
