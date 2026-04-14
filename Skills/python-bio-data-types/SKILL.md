---
name: python-bio-data-types
description: "Python data types for bioinformatics: int, float, str, bool, None, and type conversions with biological examples. Use when learning how Python represents biological data."
tool_type: python
primary_tool: Python
---

# Data Types

## Pitfalls

- **Floating-point arithmetic is not exact:** `0.1 + 0.2` is `0.30000000000000004`, not `0.3`. This is not a Python bug — it affects every language that uses IEEE 754 doubles. Always compare floats with a tolerance (`abs(a - b) < 1e-9`), not `==`.
- **`int()` truncates, it does not round:** `int(3.9)` is `3`. Use `round(3.9)` to get `4`. This distinction matters when converting a read depth fraction to an integer.
- **String immutability:** You cannot write `dna[0] = "G"`. Strings cannot be changed in place. Every "modification" creates a new string: `mutated = "G" + dna[1:]`.
- **`None` vs. empty string:** `None` means "this value does not exist". An empty string `""` means "this value exists but has zero length". These are different and should not be used interchangeably — especially when parsing FASTA or VCF files where a missing field is semantically different from an empty field.
- **`is None` not `== None`:** Always check for None with `is`: `if result is None`. The `==` check can sometimes give surprising results if a class defines custom equality.

# Variables and Data Types

## Variables

```python
variable_name = value
```python


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
```python

### Variables are references, not boxes

```python
sequence = "ATGCGATCG"
print(f"id of sequence: {id(sequence)}")

# Reassignment creates a new object (strings are immutable)
sequence = sequence + "AAA"
print(f"id after concatenation: {id(sequence)}")
print(f"New value: {sequence}")
```python


## Data Types Overview

```python
Type        Example                    Bioinformatics use
--------    -----------------------    -----------------------------------
int         42                         Sequence length, read count
float       0.487                      GC content, E-value, p-value
str         "ATGCGA"                   DNA/RNA/protein sequences, gene names
bool        True / False               Is the sequence valid? Has a stop codon?
NoneType    None                       Missing data, function with no return
```python

Use `type()` to check a value's type.


## Integers (`int`)

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
```python


## Floating-Point Numbers (`float`)

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
```python

### Floating-point precision warning

Floats are stored in binary, which means some decimal numbers cannot be represented exactly. This is a fundamental limitation of floating-point arithmetic, not a Python bug.


## Strings (`str`)

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
```python

### String indexing

```python
Index:    0   1   2   3   4   5   6   7
Seq:      A   T   G   C   G   A   T   C
Neg idx: -8  -7  -6  -5  -4  -3  -2  -1
```python


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
```python

### String immutability

Strings in Python are **immutable** -- you cannot change individual characters. Instead, you create a new string.


### Essential string methods for bioinformatics

Methods are functions that belong to a string. Call them with `string.method()`.


```python
# replace() -- replace all occurrences of a substring
dna = "ATGCGATCG"
rna = dna.replace("T", "U")    # DNA to RNA transcription
print(f"DNA: {dna}")
print(f"RNA: {rna}")
```python

```python
# startswith() and endswith() -- check sequence boundaries
cds = "ATGAAACCCGGGTAA"
print(f"Starts with ATG (start codon)? {cds.startswith('ATG')}")
print(f"Ends with TAA (stop codon)?    {cds.endswith('TAA')}")

# Check for any stop codon
has_stop = cds.endswith(("TAA", "TAG", "TGA"))
print(f"Ends with any stop codon?      {has_stop}")
```python

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
```python

## Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
