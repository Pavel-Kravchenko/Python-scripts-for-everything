---
name: python-bio-python-introduction
description: Python Introduction for Bioinformatics
tool_type: python
primary_tool: Python
---

# Python Introduction for Bioinformatics

## Why Python for Bioinformatics?

| Strength | What it means for bioinformatics |
|----------|----------------------------------|
| **Readable syntax** | Accessible to biologists without a CS background |
| **Rich ecosystem** | Libraries for sequence analysis, statistics, plotting, ML |
| **Community** | Large bioinformatics community; maintained packages |
| **Interoperability** | Easy to call external tools (BLAST, samtools, HMMER) |
| **Rapid prototyping** | Test a hypothesis quickly, then scale up |

### Where Python is used in biology

- **Genomics**: parsing FASTA/FASTQ files, variant calling pipelines
- **Transcriptomics**: differential expression analysis, RNA-seq workflows
- **Proteomics**: mass-spec data processing, protein structure prediction
- **Phylogenetics**: tree construction, evolutionary analysis
- **Systems biology**: network analysis, pathway modeling
- **Machine learning in biology**: drug discovery, protein function prediction

## Your First Bioinformatics Program

Count nucleotides in a DNA string like `"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"`:

```python
# Using .count()
a_count = dna.count("A")
c_count = dna.count("C")
g_count = dna.count("G")
t_count = dna.count("T")

print(f"  A: {a_count}")
print(f"  C: {c_count}")
print(f"  G: {g_count}")
print(f"  T: {t_count}")
print(f"  Total: {a_count + c_count + g_count + t_count}")
```

Key concepts:
- `len(dna)` — number of characters in a string
- `dna.count("A")` — occurrences of substring
- `f"..."` — f-strings embed variables directly: `f"GC = {gc:.2f}%"`

### Loop-based approach

```python
counts = {"A": 0, "C": 0, "G": 0, "T": 0}

for nucleotide in dna:
    if nucleotide in counts:
        counts[nucleotide] += 1

for nuc, count in counts.items():
    percentage = count / len(dna) * 100
    print(f"  {nuc}: {count:>3}  ({percentage:.1f}%)")
```

### GC content

**GC content** = fraction of G + C nucleotides. Key properties:
- Correlates with higher DNA melting temperature
- Characteristic ranges per organism
- GC-rich regions often correspond to gene-dense areas in mammalian genomes

```python
gc_count = dna.count("G") + dna.count("C")
gc_content = gc_count / len(dna) * 100

print(f"GC content: {gc_content:.2f}%")

if gc_content < 40:
    print("AT-rich (low GC).")
elif gc_content > 60:
    print("GC-rich.")
else:
    print("Moderate GC range.")
```

## The Python Ecosystem for Biology

### Core scientific stack

| Package | Purpose | Example use |
|---------|---------|-------------|
| **NumPy** | Numerical arrays and math | Gene expression matrix |
| **pandas** | Tabular data (DataFrames) | CSV of patient genotypes |
| **matplotlib** | Plotting | Histogram of read quality scores |
| **seaborn** | Statistical visualization | Heatmap of gene expression |
| **SciPy** | Scientific computing | Statistical tests on experimental data |

### Bioinformatics-specific

| Package | Purpose | Example use |
|---------|---------|-------------|
| **BioPython** | Sequence analysis, file parsing | Parse GenBank files, run BLAST |
| **pysam** | SAM/BAM file access | Read aligned sequencing data |
| **scikit-bio** | Diversity, alignment, trees | Phylogenetic analysis |
| **scanpy** | Single-cell genomics | Analyze scRNA-seq experiments |
| **pyVCF / cyvcf2** | Variant call data | Parse and filter VCF files |

### Machine learning

| Package | Purpose | Example use |
|---------|---------|-------------|
| **scikit-learn** | Classical ML | Classify protein families |
| **TensorFlow / PyTorch** | Deep learning | Predict protein structure, variant effects |

### BioPython quick demo

```python
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

my_seq = Seq("ATGAAACCCGGGTAA")

print(f"Sequence:           {my_seq}")
print(f"Complement:         {my_seq.complement()}")
print(f"Reverse complement: {my_seq.reverse_complement()}")
print(f"Transcribed (RNA):  {my_seq.transcribe()}")
print(f"Translated:         {my_seq.translate()}")
print(f"GC fraction:        {gc_fraction(my_seq):.3f}")
```

## Installing Packages and Virtual Environments

### pip

```bash
pip install biopython
pip install pandas==2.1.0
pip install numpy pandas matplotlib seaborn
pip install -r requirements.txt
pip list
pip install --upgrade biopython
```

### Virtual environments (venv)

```bash
python -m venv bioenv
source bioenv/bin/activate   # macOS/Linux
bioenv\Scripts\activate      # Windows
pip install biopython pandas matplotlib
deactivate
```

### conda

Preferred for bioinformatics: installs non-Python tools (BLAST, samtools) alongside Python packages.

```bash
conda create -n bioenv python=3.11
conda activate bioenv
conda install -c conda-forge biopython pandas matplotlib
conda install -c bioconda blast samtools
```

## Pitfalls

- **f-strings vs. `%` formatting:** f-strings (`f"GC = {gc:.2f}%"`) are the modern approach (Python 3.6+); older code uses `"GC = %.2f%%" % gc`
- **`print()` returns `None`:** Don't write `result = print(...)`; use a variable assignment to capture values
- **Mutable default arguments:** Never use `def f(x=[])` — use `def f(x=None)` and initialize inside
- **Off-by-one errors:** Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy:** Nested structures require `copy.deepcopy()` — `list.copy()` only copies the top level
