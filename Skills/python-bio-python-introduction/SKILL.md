---
name: python-bio-python-introduction
description: Python Introduction for Bioinformatics
tool_type: python
primary_tool: Python
---

# Python Introduction for Bioinformatics

## Pitfalls

- **f-strings vs. `%` formatting:** f-strings (e.g., `f"GC = {gc:.2f}%"`) are the modern approach introduced in Python 3.6. Older code uses `"GC = %.2f%%" % gc` — both work, but f-strings are clearer.
- **`print()` returns `None`:** Beginners sometimes write `result = print(...)` expecting to capture the output. `print()` sends text to the screen and returns `None`. To capture a value, use a variable assignment, not `print`.
- **The kernel vs. the notebook:** The "kernel" is the running Python process. If you restart the kernel, all variables disappear — you have to re-run cells from the top.
- **`!pip install` inside a notebook:** The `!` prefix runs a shell command. This installs the package into whichever Python the kernel uses. It is the standard way to install packages from inside a Jupyter notebook.

## Why Python for Bioinformatics?

Python has become the language of choice in life-science research for several reasons:

| Strength | What it means for bioinformatics |
|----------|----------------------------------|
| **Readable syntax** | Code reads almost like English, making it accessible to biologists without a CS background |
| **Rich ecosystem** | Thousands of libraries for sequence analysis, statistics, plotting, and machine learning |
| **Community** | Large bioinformatics community = lots of tutorials, Q&A, and maintained packages |
| **Interoperability** | Easy to call external tools (BLAST, samtools, HMMER) from Python scripts |
| **Rapid prototyping** | Analyze a new dataset or test a hypothesis quickly, then scale up if needed |

### Where Python is used in biology

- **Genomics**: parsing FASTA/FASTQ files, variant calling pipelines
- **Transcriptomics**: differential expression analysis, RNA-seq workflows
- **Proteomics**: mass-spec data processing, protein structure prediction
- **Phylogenetics**: tree construction, evolutionary analysis
- **Systems biology**: network analysis, pathway modeling
- **Machine learning in biology**: drug discovery, protein function prediction

## Jupyter Notebooks

You are reading this in a **Jupyter notebook** -- an interactive document that mixes text (Markdown cells) with executable code (Code cells).

### Key concepts

- **Cell**: the basic unit of a notebook. Each cell is either *Markdown* (text) or *Code* (Python).
- **Running a cell**: press `Shift + Enter` to execute the current cell and move to the next one.
- **Kernel**: the Python process that runs your code. You can restart it via the menu if things get stuck.
- **Output**: appears directly below the code cell.

### Useful keyboard shortcuts

| Shortcut | Action |
|----------|--------|
| `Shift + Enter` | Run cell, move to next |
| `Ctrl + Enter` | Run cell, stay in place |
| `Esc` then `A` | Insert cell above |
| `Esc` then `B` | Insert cell below |
| `Esc` then `M` | Convert cell to Markdown |
| `Esc` then `Y` | Convert cell to Code |
| `Esc` then `DD` | Delete cell |


### Markdown cells

Double-click this cell to see the raw Markdown source. Press `Shift + Enter` to render it again.

Markdown lets you write:
- **Bold text**, *italic text*
- Numbered and bulleted lists
- `inline code`
- Links, images, tables, and even LaTeX math: $E = mc^2$

## Your First Bioinformatics Program

### The problem

Given a DNA string like `"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"`, count how many times each nucleotide (A, C, G, T) appears.


```python
# Count each nucleotide using the built-in .count() method
a_count = dna.count("A")
c_count = dna.count("C")
g_count = dna.count("G")
t_count = dna.count("T")

print("Nucleotide counts:")
print(f"  A: {a_count}")
print(f"  C: {c_count}")
print(f"  G: {g_count}")
print(f"  T: {t_count}")
print(f"  Total: {a_count + c_count + g_count + t_count}")
```

### Understanding the code

1. `dna = "AGCTTTTC..."` -- we stored a DNA string in a **variable** called `dna`.
2. `len(dna)` -- the built-in `len()` function returns the number of characters in a string.
3. `dna.count("A")` -- the `.count()` method counts how many times a substring appears.
4. `print(f"...")` -- **f-strings** (formatted string literals) let you embed variables directly inside text using `{variable_name}`.


### Alternative approach: using a loop

The `.count()` method is convenient, but a loop gives you more flexibility. Here is another way to solve the same problem.

```python
# Count nucleotides with a for loop
counts = {"A": 0, "C": 0, "G": 0, "T": 0}

for nucleotide in dna:
    if nucleotide in counts:
        counts[nucleotide] += 1

print("Nucleotide counts (loop method):")
for nuc, count in counts.items():
    percentage = count / len(dna) * 100
    print(f"  {nuc}: {count:>3}  ({percentage:.1f}%)")
```

### GC content -- a biologically meaningful statistic

**GC content** is the fraction of nucleotides that are guanine (G) or cytosine (C). It is one of the most important characteristics of a genome:

- High GC content correlates with higher DNA melting temperature
- Different organisms have characteristic GC content ranges
- GC-rich regions often correspond to gene-dense areas in mammalian genomes

```python
# Calculate GC content
gc_count = dna.count("G") + dna.count("C")
gc_content = gc_count / len(dna) * 100

print(f"GC content: {gc_content:.2f}%")

# Interpret the result
if gc_content < 40:
    print("This is relatively AT-rich (low GC).")
elif gc_content > 60:
    print("This is relatively GC-rich.")
else:
    print("This is in the moderate GC range.")
```

## The Python Ecosystem for Biology

One of Python's greatest strengths is its library ecosystem. Here are the key packages you will encounter in this course and in bioinformatics work:

### Core scientific stack

| Package | Purpose | Example use |
|---------|---------|-------------|
| **NumPy** | Numerical arrays and math | Store a matrix of gene expression values |
| **pandas** | Tabular data (DataFrames) | Load a CSV of patient genotypes |
| **matplotlib** | Plotting and visualization | Plot a histogram of read quality scores |
| **seaborn** | Statistical visualization | Heatmap of gene expression across tissues |
| **SciPy** | Scientific computing | Statistical tests on experimental data |

### Bioinformatics-specific

| Package | Purpose | Example use |
|---------|---------|-------------|
| **BioPython** | Sequence analysis, file parsing | Parse GenBank files, run BLAST from Python |
| **pysam** | SAM/BAM file access | Read aligned sequencing data |
| **scikit-bio** | Diversity, alignment, trees | Phylogenetic analysis |
| **scanpy** | Single-cell genomics | Analyze scRNA-seq experiments |
| **pyVCF / cyvcf2** | Variant call data | Parse and filter VCF files |

### Machine learning

| Package | Purpose | Example use |
|---------|---------|-------------|
| **scikit-learn** | Classical ML | Classify protein families |
| **TensorFlow / PyTorch** | Deep learning | Predict protein structure, variant effects |

### Quick demo: BioPython

BioPython is the most widely used Python library for bioinformatics. Here is a taste of what it can do. (If BioPython is not installed, the cell below will show an ImportError -- we will fix that in Section 5.)

```python
try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction

    # Create a BioPython Seq object
    my_seq = Seq("ATGAAACCCGGGTAA")

    print(f"Sequence:           {my_seq}")
    print(f"Complement:         {my_seq.complement()}")
    print(f"Reverse complement: {my_seq.reverse_complement()}")
    print(f"Transcribed (RNA):  {my_seq.transcribe()}")
    print(f"Translated:         {my_seq.translate()}")
    print(f"GC fraction:        {gc_fraction(my_seq):.3f}")

except ImportError:
    print("BioPython is not installed yet.")
    print("Run: pip install biopython")
    print("We will cover installation in Section 5 below.")
```

## Installing Packages and Virtual Environments

### pip -- the Python package installer

`pip` is the standard tool for installing Python packages from [PyPI](https://pypi.org/) (the Python Package Index).

```bash
# Install a single package
pip install biopython

# Install a specific version
pip install pandas==2.1.0

# Install multiple packages at once
pip install numpy pandas matplotlib seaborn

# Install from a requirements file
pip install -r requirements.txt

# See what is installed
pip list

# Upgrade a package
pip install --upgrade biopython
```

Inside a Jupyter notebook, you can run shell commands by prefixing them with `!`:


### Virtual environments -- keeping projects isolated

Different projects may need different package versions. A **virtual environment** is an isolated Python installation that keeps each project's dependencies separate.

```bash
# Create a virtual environment called 'bioenv'
python -m venv bioenv

# Activate it
#   macOS / Linux:
source bioenv/bin/activate
#   Windows:
bioenv\Scripts\activate

# Now pip installs go into this environment only
pip install biopython pandas matplotlib

# Deactivate when done
deactivate
```

### conda -- an alternative for scientific computing

Many bioinformaticians use **conda** (via Miniconda or Anaconda) because it can install non-Python tools (like BLAST, samtools) alongside Python packages.

```bash
# Create a conda environment
conda create -n bioenv python=3.11

# Activate
conda activate bioenv

# Install packages (conda or pip)
conda install -c conda-forge biopython pandas matplotlib

# Install bioinformatics tools from bioconda
conda install -c bioconda blast samtools
```

**Recommendation for this course**: use either `venv` or `conda` to create one environment for all Tier 1 notebooks.

## Python Basics: A Quick Tour

### Comments

Comments explain your code to human readers. Python ignores everything after `#`.

## Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
