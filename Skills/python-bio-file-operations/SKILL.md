---
name: python-bio-file-operations
description: "By the end of this module, you will be able to: - Read and write text files using context managers - Understand text vs. binary file modes - Parse FASTA files from scratch (both reader and writer) - W"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/07_File_Operations/01_file_operations.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 7: File Operations for Bioinformatics

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/07_File_Operations/01_file_operations.ipynb`*


## Reading, Writing, and Parsing Biological Data Formats

---

### Learning Objectives

By the end of this module, you will be able to:
- Read and write text files using context managers
- Understand text vs. binary file modes
- Parse FASTA files from scratch (both reader and writer)
- Work with CSV, JSON, and Pickle formats
- Understand GenBank format basics
- Process large files efficiently (line-by-line, generators)
- Handle file-related exceptions properly

---

> *According to informal estimates, up to 90% of a bioinformatician's work involves parsing files and converting between formats!*

File handling is fundamental in bioinformatics. We constantly work with FASTA (sequences), FASTQ (sequencing reads), GenBank (annotated sequences), CSV/TSV (data tables), PDB (protein structures), and many other formats.

## How to use this notebook

1. Run Section 1 first — it creates demo files in the current directory that the rest of the notebook reads.
2. The file-format sections (FASTA, CSV/TSV, JSON) can be read independently once the demo files are created.
3. When running this notebook on your own data, remember to use absolute or relative paths that actually exist on your system.
4. The `pathlib` section shows the modern way to work with file paths — prefer it over string concatenation.

## Common stumbling points

- **Always use the `with` statement:** `with open(path) as f:` automatically closes the file when the block ends, even if an exception occurs. Forgetting to close files can cause data corruption (especially when writing) and resource leaks.
- **`r`, `w`, `a` modes:** Opening with `"w"` creates the file fresh — it **destroys existing content**. Use `"a"` to append. Use `"r"` (the default) to read.
- **`strip()` when reading lines:** Every line read from a file has a trailing `\n`. Always call `.strip()` on lines you read: `for line in f: line = line.strip()`. Forgetting this causes sequences to have invisible newlines that break comparisons.
- **Binary vs text mode:** Add `"b"` for binary files (e.g., BAM files): `open(path, "rb")`. Text mode (`"r"`) decodes bytes to strings using the system encoding — this fails or gives wrong results on binary files.
- **`FileNotFoundError`:** Check that the path exists and is spelled correctly (Python is case-sensitive on Linux/macOS). Use `pathlib.Path(path).exists()` to check before opening.
- **Reading the whole file into memory:** `f.read()` loads the entire file as a string. For large files (genomes, large FASTQ files), use `for line in f:` to process one line at a time.

[← Previous: Module 6: Functions](../06_Functions/01_functions.ipynb) | [Next: Module 8: Lists and Tuples →](../08_Lists_and_Tuples/01_lists_and_tuples.ipynb)

---

## 1. Basic File Operations

### File Modes

```python
Mode   Description
----   -----------
'r'    Read (default) — file must exist
'w'    Write — creates file, overwrites if exists
'a'    Append — adds to end of existing file
'x'    Exclusive create — fails if file already exists
'r+'   Read and write
'rb'   Read in binary mode
'wb'   Write in binary mode
```python

### Run this cell first to create demo files used throughout this notebook

```python
# Run this cell first — it creates the demo files used in this notebook
sample_fasta = """>seq1 Homo sapiens BRCA1
ATGGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2 Mus musculus Brca1 ortholog
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3 Arabidopsis thaliana ATM kinase
ATGCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""

with open('demo_sequences.fasta', 'w') as f:
    f.write(sample_fasta)

print("Demo files created:")
print("  demo_sequences.fasta (FASTA with 3 sequences)")
print("\nFirst 100 characters:")
print(sample_fasta[:100])
```python

```python
# The OLD way (not recommended): manual open/close
f = open('demo_sequences.fasta', 'r')
print(f"File is open: {not f.closed}")
first_line = f.readline()
print(f"First line: {first_line.strip()}")
f.close()
print(f"File is closed: {f.closed}")
# Problem: if an error occurs before f.close(), the file stays open!
```python

```python
# The RECOMMENDED way: context manager (with statement)
# The file is automatically closed when leaving the 'with' block,
# even if an exception occurs.

with open('demo_sequences.fasta', 'r') as f:
    content = f.read()
    print(f"Inside 'with': file open = {not f.closed}")

print(f"Outside 'with': file open = {not f.closed}")
```python

---
## 2. Reading Methods

| Method | Description | Memory |
|--------|-------------|--------|
| `f.read()` | Entire file as one string | Loads all |
| `f.read(n)` | Read n characters | Partial |
| `f.readline()` | One line at a time | One line |
| `f.readlines()` | All lines into a list | Loads all |
| `for line in f:` | Iterate line by line | One line |

```python
# read() -- entire file as one string
with open('demo_sequences.fasta', 'r') as f:
    content = f.read()
print(f"read(): {len(content)} characters")
print(content[:80], "...")
```python

```python
# readline() -- one line at a time
with open('demo_sequences.fasta', 'r') as f:
    line1 = f.readline()
    line2 = f.readline()
    print(f"Line 1: {line1.strip()}")
    print(f"Line 2: {line2.strip()}")
```python

```python
# readlines() -- all lines into a list
with open('demo_sequences.fasta', 'r') as f:
    lines = f.readlines()
print(f"readlines(): {len(lines)} lines")
for i, line in enumerate(lines):
    print(f"  {i}: {line.strip()}")
```python

```python
# BEST PRACTICE: iterate directly over file object
# Memory-efficient -- only one line is in memory at a time

with open('demo_sequences.fasta', 'r') as f:
    for line_num, line in enumerate(f, 1):
        print(f"  Line {line_num}: {line.strip()}")
```python

---
## 3. Writing Files

```python
# write() -- writes a string (you must add \n yourself)
sequences = ["ATGCGATCG", "GCTAGCTAG", "TTAACCGGTT"]

with open('output.txt', 'w') as f:
    for seq in sequences:
        f.write(seq + '\n')

# Verify
with open('output.txt', 'r') as f:
    print(f.read())
```python

```python
# 'a' mode appends to the end of an existing file
with open('output.txt', 'a') as f:
    f.write('AAATTTCCC\n')  # Added at the end

with open('output.txt', 'r') as f:
    print("After appending:")
    print(f.read())
```python

```python
# writelines() -- writes a list of strings (no newlines added)
lines = ["Gene\tLength\tGC\n", "BRCA1\t7088\t42.3\n", "TP53\t2512\t51.2\n"]

with open('genes.tsv', 'w') as f:
    f.writelines(lines)

with open('genes.tsv', 'r') as f:
    print(f.read())
```python

---
## 4. Text vs. Binary Files

- **Text mode** (`'r'`, `'w'`): reads/writes strings, handles line endings
- **Binary mode** (`'rb'`, `'wb'`): reads/writes raw bytes, needed for images, compressed files, serialized data

```python
# Writing and reading binary data
data = b'\x00\x01\x02\x03ATGC'  # bytes literal

with open('binary_demo.bin', 'wb') as f:
    f.write(data)

with open('binary_demo.bin', 'rb') as f:
    raw = f.read()

print(f"Raw bytes: {raw}")
print(f"Type: {type(raw)}")
print(f"Decoded text portion: {raw[4:].decode('ascii')}")

import os
os.remove('binary_demo.bin')
```python

---
## 5. FASTA File Format

FASTA is the most common sequence format in bioinformatics.

```python
>sequence_id description text
ATGCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGA
>another_sequence more info
GGGGCCCCAAAATTTT
```python

Rules:
- Header lines start with `>`
- Sequence can span multiple lines (typically 60-80 chars per line)
- Empty lines are allowed and should be skipped

```python
# Create a more realistic FASTA file
fasta_content = """>gene1|BRCA1|Homo_sapiens DNA repair protein
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGG
TTAACCGGTTAA
>gene2|TP53|Homo_sapiens Tumor protein p53
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATT
TTCCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA
>gene3|EGFR|Homo_sapiens Epidermal growth factor receptor
ATATATATATATATATGCGCGCGCGCGCGCGCATATATATATATATATGCGCGCGCGCGCG
CGC
"""

with open('sequences.fasta', 'w') as f:
    f.write(fasta_content)

print("FASTA file created with 3 sequences.")
```python

```python
def read_fasta(filename):
    """Parse a FASTA file and return a dictionary of sequences.
    
    Args:
        filename: Path to FASTA file
    
    Returns:
        dict mapping full headers to concatenated sequences
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence!
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences

# Test
seqs = read_fasta('sequences.fasta')
for header, seq in seqs.items():
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    print(f"{header.split()[0]}")
    print(f"  Length: {len(seq)} bp, GC: {gc:.1f}%")
    print(f"  First 40 bp: {seq[:40]}...")
    print()
```python

```python
# Alternative FASTA parser using split (good for smaller files)
def read_fasta_split(filename):
    """Parse FASTA by splitting on '>' characters.
    
    Simpler but loads entire file into memory.
    
    Returns:
        list of (header, sequence) tuples
    """
    with open(filename) as f:
        entries = f.read().split('>')[1:]  # Split and drop empty first element
    
    result = []
    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])
        result.append((header, sequence))
    
    return result

for header, seq in read_fasta_split('sequences.fasta'):
    print(f"{header.split('|')[1]}: {len(seq)} bp")
```python

```python
def write_fasta(sequences, filename, line_width=60):
    """Write sequences to a FASTA file.
    
    Args:
        sequences: dict of {header: sequence} or list of (header, sequence)
        filename: Output file path
        line_width: Characters per sequence line (default: 60)
    """
    # Handle both dict and list input
    if isinstance(sequences, dict):
        items = sequences.items()
    else:
        items = sequences
    
    with open(filename, 'w') as f:
        for header, seq in items:
            f.write(f">{header}\n")
            # Write sequence in wrapped lines
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')

# Test: filter GC-rich sequences and write to new file
gc_rich = {
    header: seq for header, seq in seqs.items()
    if (seq.count('G') + seq.count('C')) / len(seq) > 0.45
}

write_fasta(gc_rich, 'gc_rich.fasta')
print(f"Wrote {len(gc_rich)} GC-rich sequences to gc_rich.fasta")

# Verify
with open('gc_rich.fasta', 'r') as f:
    print(f.read())
```python

---
## 6. CSV and TSV Files

Tabular data formats commonly used for gene expression data, annotations, and analysis results.

- **CSV** (Comma-Separated Values): fields separated by commas
- **TSV** (Tab-Separated Values): fields separated by tabs

```python
import csv

# Create sample gene expression CSV
expression_data = [
    ['gene_id', 'gene_name', 'length', 'gc_content', 'expression_sample1', 'expression_sample2'],
    ['ENSG0001', 'BRCA1', '7088', '42.3', '150.5', '175.2'],
    ['ENSG0002', 'TP53', '2512', '51.2', '89.3', '95.1'],
    ['ENSG0003', 'EGFR', '5616', '48.7', '245.8', '230.0'],
    ['ENSG0004', 'MYC', '2357', '55.1', '312.4', '298.7'],
]

# Write CSV
with open('gene_expression.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(expression_data)

print("CSV file created.")
```python

```python
# Reading CSV with csv.reader
with open('gene_expression.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # First row is the header
    print(f"Columns: {header}")
    print()
    for row in reader:
        print(f"  {row[1]:6s} | Length: {row[2]:>5s} | GC: {row[3]}% | Expr: {row[4]}")
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
