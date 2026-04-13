---
name: python-bio-data-types
description: "Split from `01_variables_and_data_types.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/02_Variables_and_Data_Types/02_data_types.ipynb"
---

# Data Types

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/02_Variables_and_Data_Types/02_data_types.ipynb`*

# Data Types

Split from `01_variables_and_data_types.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_variables_and_data_types.ipynb) · [Previous: Variables](./01_variables.ipynb)

## Why this notebook matters

The five core types — `int`, `float`, `str`, `bool`, `None` — are the building blocks of all Python programs. In bioinformatics, knowing which type to use matters: a sequence length is an `int` (never a `float`), GC content is a `float`, a DNA sequence is a `str`, and a gene annotation that has not been loaded yet should be `None`, not an empty string. Getting types right prevents silent errors when doing arithmetic or writing output files.

## How to use this notebook

1. Run cells in order. Sections 3–7 each cover one type; Section 8 covers conversions between them.
2. The string section is the longest — strings are the most important type for bioinformatics work. Do not skip it.
3. The exercises at the end (Section 10) are designed to be run and verified — try them before looking at the hints.

## Common stumbling points

- **Floating-point arithmetic is not exact:** `0.1 + 0.2` is `0.30000000000000004`, not `0.3`. This is not a Python bug — it affects every language that uses IEEE 754 doubles. Always compare floats with a tolerance (`abs(a - b) < 1e-9`), not `==`.
- **`int()` truncates, it does not round:** `int(3.9)` is `3`. Use `round(3.9)` to get `4`. This distinction matters when converting a read depth fraction to an integer.
- **String immutability:** You cannot write `dna[0] = "G"`. Strings cannot be changed in place. Every "modification" creates a new string: `mutated = "G" + dna[1:]`.
- **`None` vs. empty string:** `None` means "this value does not exist". An empty string `""` means "this value exists but has zero length". These are different and should not be used interchangeably — especially when parsing FASTA or VCF files where a missing field is semantically different from an empty field.
- **`is None` not `== None`:** Always check for None with `is`: `if result is None`. The `==` check can sometimes give surprising results if a class defines custom equality.

# Module 2: Variables and Data Types

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Create variables and follow Python naming conventions
2. Work with all basic data types: `int`, `float`, `str`, `bool`, `None`
3. Convert between types safely
4. Use strings to represent and manipulate biological sequences
5. Apply these skills to real bioinformatics calculations

---

## 1. Variables

A **variable** is a name that refers to a value stored in memory. In Python you create a variable simply by assigning a value to a name.

```
variable_name = value
```

Think of it as labeling a container:

```
gene_name -----> "BRCA1"     (the label points to the value)
seq_length ----> 5590        (another label, another value)
```

```python
# Assigning values to variables
gene_name = "BRCA1"
chromosome = 17
gc_content = 0.423
is_tumor_suppressor = True

print(gene_name)
print(chromosome)
print(gc_content)
print(is_tumor_suppressor)
```

### Naming rules and conventions

**Rules** (breaking these causes an error):
- Names can contain letters, digits, and underscores
- Names must start with a letter or underscore (not a digit)
- Names are case-sensitive (`gene` and `Gene` are different variables)
- Python keywords (`if`, `for`, `class`, ...) cannot be used as names

**Conventions** (follow these for readable code):
- Use `snake_case` for variable and function names: `gene_name`, `sequence_length`
- Use `UPPER_CASE` for constants: `AVOGADRO = 6.022e23`
- Choose descriptive names: `gc_content` is better than `gc` or `x`

```python
# Good variable names
sequence_length = 1500
gene_name = "TP53"
melting_temperature = 65.5

# Bad variable names (but technically valid)
x = 1500              # meaningless
seqLen = 1500         # not snake_case (this is camelCase)
three = 1             # misleading!

# This will cause a SyntaxError:
# 2nd_gene = "EGFR"  # cannot start with a digit
```

### Multiple assignment

Python allows assigning multiple variables in one line. This is especially useful when unpacking related values.

```python
# Assign multiple variables at once
a_count, t_count, g_count, c_count = 250, 245, 280, 275
print(f"A: {a_count}, T: {t_count}, G: {g_count}, C: {c_count}")

# Swap two variables (elegant Python idiom)
sense_strand = "ATGCGA"
antisense_strand = "TACGCT"

sense_strand, antisense_strand = antisense_strand, sense_strand
print(f"After swap -- sense: {sense_strand}, antisense: {antisense_strand}")
```

### Variables are references, not boxes

In Python, a variable does not *contain* a value -- it *points to* an object in memory. The built-in `id()` function shows the memory address of an object.

```python
sequence = "ATGCGATCG"
print(f"id of sequence: {id(sequence)}")

# Reassignment creates a new object (strings are immutable)
sequence = sequence + "AAA"
print(f"id after concatenation: {id(sequence)}")
print(f"New value: {sequence}")
```

---

## 2. Data Types Overview

Every value in Python has a **type** that determines what operations are allowed.

```
Type        Example                    Bioinformatics use
--------    -----------------------    -----------------------------------
int         42                         Sequence length, read count
float       0.487                      GC content, E-value, p-value
str         "ATGCGA"                   DNA/RNA/protein sequences, gene names
bool        True / False               Is the sequence valid? Has a stop codon?
NoneType    None                       Missing data, function with no return
```

Use `type()` to check a value's type.

```python
# Check types of bioinformatics values
read_count = 15000000
gc_fraction = 0.52
organism = "Escherichia coli"
is_model_organism = True
annotation = None

values = [read_count, gc_fraction, organism, is_model_organism, annotation]
for v in values:
    print(f"{str(v):25s} -> {type(v).__name__}")
```

---

## 3. Integers (`int`)

Integers are whole numbers with no decimal point. In bioinformatics, you use them for counts, positions, and indices.

Python integers have **unlimited precision** -- they can be as large as your memory allows.

```python
# Typical bioinformatics integers
sequence_length = 3088286401          # human genome length in bp
num_genes = 20000                     # approximate protein-coding genes
read_depth = 30                       # sequencing coverage
chromosome_number = 23                # human haploid chromosome count

print(f"Human genome:  {sequence_length:,} bp")   # comma-separated formatting
print(f"Protein-coding genes: ~{num_genes:,}")
print(f"Target coverage: {read_depth}x")
```

```python
# Integers have unlimited precision in Python
huge_number = 2 ** 100
print(f"2^100 = {huge_number}")
print(f"Number of digits: {len(str(huge_number))}")
print(f"Type: {type(huge_number)}")
```

---

## 4. Floating-Point Numbers (`float`)

Floats represent decimal numbers. They are used for measurements, percentages, scores, and probabilities.

```python
# Typical bioinformatics floats
gc_content = 0.508                    # GC fraction
melting_temp = 72.3                   # PCR primer Tm in degrees C
e_value = 1.5e-42                     # BLAST E-value (scientific notation)
p_value = 0.0031                      # statistical significance

print(f"GC content: {gc_content}")
print(f"GC content as %: {gc_content * 100:.1f}%")
print(f"Melting temp: {melting_temp} C")
print(f"E-value: {e_value}")
print(f"E-value formatted: {e_value:.2e}")
print(f"p-value: {p_value}")
```

### Floating-point precision warning

Floats are stored in binary, which means some decimal numbers cannot be represented exactly. This is a fundamental limitation of floating-point arithmetic, not a Python bug.

```python
# Floating-point surprise
print(0.1 + 0.2)          # not exactly 0.3!
print(0.1 + 0.2 == 0.3)   # False!

# In practice, compare floats with a tolerance
result = 0.1 + 0.2
expected = 0.3
tolerance = 1e-9

print(f"Close enough? {abs(result - expected) < tolerance}")  # True
```

---

## 5. Strings (`str`)

Strings are sequences of characters. They are the **most important data type in bioinformatics** because DNA, RNA, and protein sequences are all represented as strings.

### Creating strings

```python
# Three ways to create strings
dna = "ATGCGATCGATCG"          # double quotes
rna = 'AUGCGAUCGAUCG'          # single quotes (identical behavior)
protein = """MKWVTFISLLLLFSSAYS"""  # triple quotes (can span multiple lines)

# Multi-line string (useful for FASTA headers, etc.)
fasta_entry = """>sp|P04637|P53_HUMAN Cellular tumor antigen p53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLPGRNSFEV"""

print(fasta_entry[:80] + "...")
```

### String indexing

Each character in a string has an **index** (position number). Python uses zero-based indexing.

```
Index:    0   1   2   3   4   5   6   7
Seq:      A   T   G   C   G   A   T   C
Neg idx: -8  -7  -6  -5  -4  -3  -2  -1
```

```python
dna = "ATGCGATC"

# Positive indexing (from the left)
print(f"First nucleotide:  dna[0]  = {dna[0]}")
print(f"Third nucleotide:  dna[2]  = {dna[2]}")

# Negative indexing (from the right)
print(f"Last nucleotide:   dna[-1] = {dna[-1]}")
print(f"Second to last:    dna[-2] = {dna[-2]}")
```

### String slicing

Slicing extracts a substring. Syntax: `string[start:stop:step]`

- `start` is inclusive (default: 0)
- `stop` is exclusive (default: end of string)
- `step` is the increment (default: 1)

```python
dna = "ATGAAACCCGGGTAA"

# Extract the start codon (first 3 nucleotides)
start_codon = dna[0:3]
print(f"Start codon: {start_codon}")       # ATG

# Extract the stop codon (last 3 nucleotides)
stop_codon = dna[-3:]
print(f"Stop codon:  {stop_codon}")        # TAA

# Extract the coding region between start and stop
coding = dna[3:-3]
print(f"Coding region: {coding}")          # AAACCCGGG

# Every third nucleotide (first position of each codon)
first_positions = dna[0::3]
print(f"First codon positions: {first_positions}")  # AACGT

# Reverse the sequence
reversed_dna = dna[::-1]
print(f"Reversed: {reversed_dna}")         # AATGGGCCCAAAGTA
```

### String immutability

Strings in Python are **immutable** -- you cannot change individual characters. Instead, you create a new string.

```python
dna = "ATGCGA"

# This will raise a TypeError:
# dna[0] = "G"  # TypeError: 'str' object does not support item assignment

# Instead, create a new string
mutated = "G" + dna[1:]   # point mutation: A -> G at position 0
print(f"Original: {dna}")
print(f"Mutated:  {mutated}")
```

### Essential string methods for bioinformatics

Methods are functions that belong to a string. Call them with `string.method()`.

```python
dna = "atgcGATCgatc"

# Case conversion -- important when sequences come in mixed case
print(f"Upper: {dna.upper()}")       # ATGCGATCGATC
print(f"Lower: {dna.lower()}")       # atgcgatcgatc
```

```python
dna = "ATGCGATCGATCGTAGC"

# count() -- count occurrences of a substring
print(f"Number of G's: {dna.count('G')}")
print(f"Number of 'ATC' motifs: {dna.count('ATC')}")
```

```python
# find() -- find the position of a substring (-1 if not found)
dna = "ATGAAACCCGGGTAA"
print(f"Position of 'GGG': {dna.find('GGG')}")   # 9
print(f"Position of 'TTT': {dna.find('TTT')}")   # -1 (not found)
```

```python
# replace() -- replace all occurrences of a substring
dna = "ATGCGATCG"
rna = dna.replace("T", "U")    # DNA to RNA transcription
print(f"DNA: {dna}")
print(f"RNA: {rna}")
```

```python
# startswith() and endswith() -- check sequence boundaries
cds = "ATGAAACCCGGGTAA"
print(f"Starts with ATG (start codon)? {cds.startswith('ATG')}")
print(f"Ends with TAA (stop codon)?    {cds.endswith('TAA')}")

# Check for any stop codon
has_stop = cds.endswith(("TAA", "TAG", "TGA"))
print(f"Ends with any stop codon?      {has_stop}")
```

```python
# split() and join() -- essential for parsing biological file formats

# Parse a FASTA header
header = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53"
parts = header.split("|")
print(f"Database: {parts[0][1:]}")
print(f"Accession: {parts[1]}")
print(f"Entry name: {parts[2].split()[0]}")

# Join codons with a separator
codons = ["ATG", "AAA", "CCC", "GGG", "TAA"]
formatted = " - ".join(codons)
print(f"\nCodons: {formatted}")
```
