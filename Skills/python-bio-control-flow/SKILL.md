---
name: python-bio-control-flow
description: "By the end of this module, you will be able to: - Use `if/elif/else` to classify biological data - Apply `for` and `while` loops to iterate over sequences - Control loop execution with `break`, `conti"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/05_Control_Flow/01_control_flow.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 5: Control Flow for Bioinformatics

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/05_Control_Flow/01_control_flow.ipynb`*


---

### Learning Objectives

By the end of this module, you will be able to:
- Use `if/elif/else` to classify biological data
- Apply `for` and `while` loops to iterate over sequences
- Control loop execution with `break`, `continue`, and `else`
- Use `enumerate()`, `zip()`, and `range()` effectively
- Translate codons to amino acids using a dictionary
- Find open reading frames (ORFs) in DNA sequences

### Prerequisites
- Python strings and string methods (Module 4)
- Basic understanding of DNA, RNA, codons, and amino acids

---

## How to use this notebook

1. Run cells in order. Section 1 covers `if`/`elif`/`else`, sections 2–3 cover `for` and `while` loops, section 4 covers `break`/`continue`/`pass`, and later sections cover advanced patterns.
2. The biological examples are not decoration — every code pattern introduced has a direct bioinformatics application shown immediately after.
3. After reading a code example, predict the output before running it. This is the fastest way to build intuition.

## Common stumbling points

- **Indentation is not optional:** Python uses indentation to define blocks. Inconsistent indentation (mixing tabs and spaces, or wrong number of spaces) causes `IndentationError`. Use 4 spaces per level — this is the universal Python convention.
- **`elif` vs separate `if` statements:** `if ... elif ... elif ...` is a single decision tree — at most one branch runs. Multiple `if` statements each run independently. For classifying a codon, you want `elif`; for applying multiple independent filters, you want separate `if`s.
- **Off-by-one in `range()`:** `range(len(seq) - 2)` stops 2 before the end, leaving room for a 3-character slice. A common mistake is `range(len(seq))` which produces an incomplete last codon.
- **`while` loops can run forever:** Every `while` loop must have a condition that eventually becomes `False`, or a `break` that exits. If your notebook hangs, the kernel is likely stuck in an infinite loop — use the Stop button in Jupyter to interrupt.
- **`break` exits only the innermost loop:** In nested loops (e.g., scanning all reading frames), `break` exits only the loop it is directly inside.
- **`pass` does nothing:** It is a syntactic placeholder when a block must exist but you want it empty. Common when sketching out code structure before filling it in.

---

[← Previous: Module 4: Strings and Biological Sequences](../04_Strings_and_Sequences/01_strings_and_sequences.ipynb) | [Next: Module 6: Functions →](../06_Functions/01_functions.ipynb)

---

## 1. Conditional Statements: `if` / `elif` / `else`

Conditionals let your program make decisions. In bioinformatics, we constantly classify, filter, and branch based on sequence properties.

```python
         condition
           / \
        True  False
         |      |
      if block  else block
```python

**Syntax:**

```python
if condition:
    # runs when condition is True
elif another_condition:
    # runs when first is False but this is True
else:
    # runs when all conditions above are False
```python

```python
# Basic example: validate a nucleotide
nucleotide = "A"

if nucleotide in "ATGC":
    print(f"{nucleotide} is a valid DNA nucleotide")
elif nucleotide == "U":
    print(f"{nucleotide} is an RNA nucleotide (uracil)")
else:
    print(f"{nucleotide} is not a standard nucleotide")

# Try changing 'nucleotide' to 'U', 'X', or 'T' to see different branches
```python

### Classifying Nucleotides: Purines vs. Pyrimidines

- **Purines** (A, G): double-ring structure, larger molecules
- **Pyrimidines** (C, T, U): single-ring structure, smaller molecules

Chargaff's rule: in double-stranded DNA, the number of purines equals the number of pyrimidines.

```python
def classify_nucleotide(nuc):
    """Classify a nucleotide as purine or pyrimidine."""
    nuc = nuc.upper()
    if nuc in "AG":
        return "Purine (double-ring)"
    elif nuc in "CTU":
        return "Pyrimidine (single-ring)"
    else:
        return "Unknown"

for base in "ATGCU":
    print(f"{base}: {classify_nucleotide(base)}")
```python

### GC Content Classification

GC content varies across genomes and has biological significance:

| GC Range | Classification | Examples |
|----------|---------------|----------|
| < 30% | AT-rich | Plasmodium falciparum (~20%) |
| 30-50% | Moderate | Homo sapiens (~41%) |
| 50-60% | High GC | Many bacterial genes |
| > 60% | Very high GC | Streptomyces (~72%) |

```python
def classify_gc_content(sequence):
    """Classify a DNA sequence by its GC content."""
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    gc_percent = (gc_count / len(seq)) * 100
    
    if gc_percent < 30:
        category = "AT-rich"
    elif gc_percent < 50:
        category = "Moderate"
    elif gc_percent < 60:
        category = "High GC"
    else:
        category = "Very high GC"
    
    return gc_percent, category

# Test with sequences of different GC content
sequences = {
    "AT-rich region":   "AAATTTAAATTTAAATTTAAATTT",
    "Human average":    "ATGCGATCAATCGTATGCATGCA",
    "Bacterial gene":   "GCGATCGCGCGATCGCGATCGCG",
    "Streptomyces-like": "GCGCGCGCGCGCGCGCGCGCGCG"
}

print(f"{'Name':<20} {'Sequence':<26} {'GC%':>5}  Category")
print("-" * 70)
for name, seq in sequences.items():
    gc, cat = classify_gc_content(seq)
    print(f"{name:<20} {seq:<26} {gc:5.1f}%  {cat}")
```python

### Detecting Sequence Type

Given an unknown biological sequence, we can determine whether it is DNA, RNA, or protein based on its character composition.

```python
def detect_sequence_type(sequence):
    """Detect whether a sequence is DNA, RNA, or protein."""
    seq = sequence.upper()
    unique_chars = set(seq)
    
    if unique_chars <= set("ATGC"):
        return "DNA"
    elif unique_chars <= set("AUGC"):
        return "RNA"
    elif unique_chars <= set("ACDEFGHIKLMNPQRSTVWY"):
        return "Protein"
    else:
        return "Unknown"

test_sequences = [
    "ATGCGATCGATCG",
    "AUGCGAUCGAUCG",
    "MKVLWAALLVLLGFANAT",
    "ATGXYZ123"
]

for seq in test_sequences:
    seq_type = detect_sequence_type(seq)
    print(f"{seq:<25} -> {seq_type}")
```python

## 2. For Loops

For loops iterate over a sequence (string, list, range, etc.). They are the workhorse of bioinformatics programming.

```python
for item in iterable:
    # process item
```python

```python
# Iterate over nucleotides in a sequence
dna = "ATGCGATC"

print("Counting nucleotides manually:")
a_count = 0
t_count = 0
g_count = 0
c_count = 0

for nucleotide in dna:
    if nucleotide == 'A':
        a_count += 1
    elif nucleotide == 'T':
        t_count += 1
    elif nucleotide == 'G':
        g_count += 1
    elif nucleotide == 'C':
        c_count += 1

print(f"A={a_count}, T={t_count}, G={g_count}, C={c_count}")
```python

### The `range()` Function

Generates a sequence of integers. Essential for position-based iteration.

```python
range(stop)              -> 0, 1, 2, ..., stop-1
range(start, stop)       -> start, start+1, ..., stop-1
range(start, stop, step) -> start, start+step, start+2*step, ...
```python

```python
# range() examples
print("range(5):        ", list(range(5)))
print("range(1, 6):     ", list(range(1, 6)))
print("range(0, 15, 3): ", list(range(0, 15, 3)))  # codon start positions
print("range(10, 0, -1):", list(range(10, 0, -1)))  # countdown
```python

```python
# Extract codons using range with step 3
dna = "ATGGCCGATCGATAGCCA"

print(f"DNA: {dna}")
print(f"Length: {len(dna)} bp")
print("\nCodons:")

codons = []
for i in range(0, len(dna) - 2, 3):
    codon = dna[i:i+3]
    if len(codon) == 3:  # only complete codons
        codons.append(codon)
        print(f"  Position {i:2d}-{i+2:2d}: {codon}")

print(f"\nTotal complete codons: {len(codons)}")
```python

### `enumerate()`: Get Index and Value

When you need both the position and the value during iteration, use `enumerate()`.

```python
# Find positions of a specific nucleotide
dna = "ATGCGATCGATCGTAG"

print(f"Sequence: {dna}")
print(f"Positions of 'G':")

g_positions = []
for i, nuc in enumerate(dna):
    if nuc == 'G':
        g_positions.append(i)

print(f"  0-based: {g_positions}")
print(f"  1-based: {[p+1 for p in g_positions]}")
```python

### `zip()`: Iterate Over Multiple Sequences in Parallel

`zip()` pairs up elements from two or more iterables. Very useful for comparing sequences, aligning positions, or combining related data.

```python
# Compare two aligned DNA sequences
seq1 = "ATGCGATCGA"
seq2 = "ATGCAATCTA"

print(f"Seq 1: {seq1}")
print(f"Seq 2: {seq2}")
print("Match: ", end="")

mismatches = 0
for nuc1, nuc2 in zip(seq1, seq2):
    if nuc1 == nuc2:
        print("|", end="")
    else:
        print("X", end="")
        mismatches += 1

print(f"\n\nMismatches: {mismatches}")
print(f"Identity: {(len(seq1) - mismatches) / len(seq1) * 100:.1f}%")
```python

```python
# zip() with enumerate() for position-aware comparison
seq1 = "ATGCGATCGA"
seq2 = "ATGCAATCTA"

print("Mutation report:")
for i, (nuc1, nuc2) in enumerate(zip(seq1, seq2)):
    if nuc1 != nuc2:
        print(f"  Position {i+1}: {nuc1} -> {nuc2}")
```python

```python
# zip() to pair gene names with their sequences
gene_names = ["BRCA1", "TP53", "EGFR"]
gene_lengths = [81189, 19149, 188307]
gc_contents = [41.8, 42.4, 48.2]

print(f"{'Gene':<10} {'Length':>10} {'GC%':>6}")
print("-" * 28)
for name, length, gc in zip(gene_names, gene_lengths, gc_contents):
    print(f"{name:<10} {length:>10,} {gc:>5.1f}%")
```python

### Nested Loops

Loops inside loops. Useful for generating combinations, searching in 2D, or analyzing all reading frames.

```python
# Generate all 64 possible codons
nucleotides = "ATGC"
all_codons = []

for first in nucleotides:
    for second in nucleotides:
        for third in nucleotides:
            all_codons.append(first + second + third)

print(f"Total codons: {len(all_codons)}")
print(f"First 8:  {all_codons[:8]}")
print(f"Last 8:   {all_codons[-8:]}")

# How many start with 'A'?
a_start = [c for c in all_codons if c.startswith('A')]
print(f"\nCodons starting with A: {len(a_start)}")
```python

## 3. While Loops

While loops keep running as long as a condition is True. Useful when you do not know in advance how many iterations you need.

```python
while condition:
    # loop body
    # must eventually make condition False!
```python

```python
# Find the first stop codon in a reading frame
dna = "ATGGCCGATCGATAGCCATAGTTAACG"
stop_codons = {"TAA", "TAG", "TGA"}

position = 0
found = False

while position <= len(dna) - 3:
    codon = dna[position:position + 3]
    if codon in stop_codons:
        print(f"Stop codon '{codon}' found at position {position}")
        found = True
        break
    position += 3  # step by codon

if not found:
    print("No stop codon found in this reading frame")
```python

```python
# Find all occurrences of a motif using while
dna = "ATGCGATGATCGATGCATG"
motif = "ATG"

positions = []
pos = dna.find(motif)
while pos != -1:
    positions.append(pos)
    pos = dna.find(motif, pos + 1)

print(f"Sequence: {dna}")
print(f"'{motif}' found at positions: {positions}")
```python

## 4. Loop Control: break, continue, and the for-else Pattern

| Statement | Effect |
|-----------|--------|
| `break` | Exit the loop immediately |
| `continue` | Skip to the next iteration |
| `else` (on loop) | Runs only if the loop completed without `break` |

```python
# BREAK: stop at the first invalid character
sequence = "ATGCXGATC"
valid_bases = set("ATGC")

print(f"Validating: {sequence}")
for i, nuc in enumerate(sequence):
    if nuc not in valid_bases:
        print(f"Invalid character '{nuc}' at position {i}. Stopping.")
        break
else:
    # This block runs only if the loop did NOT break
    print("Sequence is valid!")
```python

```python
# Validate a correct sequence -- the else clause runs
sequence = "ATGCGATC"
valid_bases = set("ATGC")

print(f"Validating: {sequence}")
for i, nuc in enumerate(sequence):
    if nuc not in valid_bases:
        print(f"Invalid character '{nuc}' at position {i}. Stopping.")
        break
else:
    print("Sequence is valid!")
```python

```python
# CONTINUE: skip non-standard amino acids
protein = "MKVXLWABLLVZLLGFANAT"
standard_aa = set("ACDEFGHIKLMNPQRSTVWY")

clean_protein = []
skipped = 0

for aa in protein:
    if aa not in standard_aa:
        skipped += 1
        continue  # skip non-standard amino acids
    clean_protein.append(aa)

print(f"Original: {protein} ({len(protein)} aa)")
print(f"Cleaned:  {''.join(clean_protein)} ({len(clean_protein)} aa)")
print(f"Skipped {skipped} non-standard residues")
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
