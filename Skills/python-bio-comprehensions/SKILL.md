---
name: python-bio-comprehensions
description: "Sequence data is almost always processed as collections: all sequences in a FASTA file, all k-mers in a read, all codons in a CDS. Comprehensions replace verbose loops with concise, readable expressio"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/10_Comprehensions/01_comprehensions.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 10: Comprehensions

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/10_Comprehensions/01_comprehensions.ipynb`*


---

### Learning Objectives
- Write list comprehensions with and without conditions
- Build dict and set comprehensions
- Understand nested comprehensions and generator expressions
- Apply comprehensions to biological data: GC filtering, codon extraction, reverse complement, k-mer analysis

---

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Order of clauses in nested comprehensions**
The loop order in a flat nested comprehension matches nested `for` loops — outer loop first:
```python
# flat: [expr for outer in outer_list for inner in inner_list]
# is the same as:
# for outer in outer_list:
#     for inner in inner_list: result.append(expr)
```python
`[[expr for inner in ...] for outer in ...]` gives a list-of-lists instead.

**2. `if` filter vs `if`/`else` transform**
- `if` after `for` is a **filter** (removes items that fail the test).
- `if`/`else` in the **expression** is a **conditional transform** (every item is kept, but mapped differently).
```python
[nt for nt in dna if nt in 'GC']              # filter: keeps only G and C
["purine" if nt in "AG" else "pyrimidine" for nt in dna]  # transform: relabels every base
```python

**3. Generator expressions are single-pass**
Once consumed (`sum(...)`, `list(...)`, etc.), the generator is exhausted. If you need to iterate the result more than once, use a list comprehension.

**4. Performance note**
Dict comprehensions with `.count()` inside recompute the count for every unique key. For k-mer counting, `collections.Counter` is faster.

```python
# Sequence lengths and GC content
sequences = ["ATATATAT", "GCGCGCGC", "ATGCATGC", "AAAGGGCCC", "TTTTAAAA"]

lengths = [len(seq) for seq in sequences]
print(f"Lengths: {lengths}")

gc_values = [
    (seq.count('G') + seq.count('C')) / len(seq) * 100
    for seq in sequences
]

for seq, gc in zip(sequences, gc_values):
    print(f"  {seq:12s}  GC = {gc:.1f}%")
```python

```python
# Extract codons from a coding sequence
cds = "ATGGCCGATCGATAGCCA"

codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]

print(f"CDS:    {cds}")
print(f"Codons: {codons}")
```python

### 1.2 List Comprehensions with Conditions

```python
# Filter sequences by GC content
sequences = ["ATATATAT", "GCGCGCGC", "ATGCATGC", "AAAGGGCCC", "TTTTAAAA"]

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

gc_rich = [seq for seq in sequences if gc_content(seq) > 50]

print("GC-rich sequences (>50%):")
for seq in gc_rich:
    print(f"  {seq}: {gc_content(seq):.1f}%")
```python

```python
# Find all positions of a motif
dna = "ATGATGCCCATGATGATG"
motif = "ATG"

positions = [i for i in range(len(dna) - len(motif) + 1) if dna[i:i+len(motif)] == motif]

print(f"Sequence: {dna}")
print(f"'{motif}' found at positions: {positions}")
```python

```python
# if-else in the expression -- conditional transformation
# Classify nucleotides as purine (A, G) or pyrimidine (C, T)
dna = "ATGCGATCG"

classification = [
    "purine" if nt in "AG" else "pyrimidine"
    for nt in dna
]

for nt, cls in zip(dna, classification):
    print(f"  {nt} -> {cls}")
```python

```python
# Filter AND transform: keep only long sequences, return (name, gc%)
fasta_data = {
    "seq1": "ATG",
    "seq2": "GCGCGCGCGCGC",
    "seq3": "ATATATATATATAT",
    "seq4": "ATGCATGCATGCATGC",
    "seq5": "GC",
}

results = [
    (name, gc_content(seq))
    for name, seq in fasta_data.items()
    if len(seq) > 10
]

print("Long sequences and their GC%:")
for name, gc in results:
    print(f"  {name}: {gc:.1f}%")
```python

---

## 2. Dictionary Comprehensions

```python
+----------------------------------------------------------------+
|              DICT COMPREHENSION ANATOMY                         |
+----------------------------------------------------------------+
|                                                                 |
|   {key_expr: val_expr  for item in iterable  if condition}      |
|                                                                 |
|   Example:                                                      |
|   {seq: len(seq) for seq in sequences}                          |
|   -> {'ATGC': 4, 'GGCAT': 5, ...}                              |
|                                                                 |
+----------------------------------------------------------------+
```python

```python
# Count nucleotides in a sequence
dna = "ATGCGATCGATCGTAGCGATCGATCG"

nt_counts = {nt: dna.count(nt) for nt in 'ATGC'}
print(f"Nucleotide counts: {nt_counts}")

# Create a complement mapping and its reverse
complement = {b: c for b, c in [('A','T'), ('T','A'), ('G','C'), ('C','G')]}
reverse_map = {v: k for k, v in complement.items()}
print(f"Complement: {complement}")
print(f"Reverse:    {reverse_map}")
```python

```python
# Map sequence IDs to GC content
fasta = {
    "gene_A": "ATGCGATCGATCG",
    "gene_B": "GCGCGCGCGCGC",
    "gene_C": "ATATATATATATAT",
    "gene_D": "ATGCATGCATGC",
}

gc_map = {name: gc_content(seq) for name, seq in fasta.items()}

print("GC content per sequence:")
for name, gc in gc_map.items():
    print(f"  {name}: {gc:.1f}%")
```python

```python
# K-mer counting with dict comprehension
dna = "ATGATGATGATG"
k = 3

kmers = [dna[i:i+k] for i in range(len(dna) - k + 1)]
kmer_counts = {kmer: kmers.count(kmer) for kmer in set(kmers)}

print(f"Sequence: {dna}")
print(f"{k}-mer counts:")
for kmer, count in sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True):
    print(f"  {kmer}: {count}")
```python

```python
# Filter a dictionary -- keep only long genes
gene_lengths = {
    "BRCA1": 81189, "TP53": 19149, "EGFR": 188307, "MYC": 4150, "KRAS": 45693,
}

long_genes = {gene: bp for gene, bp in gene_lengths.items() if bp > 50000}
print(f"Genes longer than 50 kb: {long_genes}")
```python

---

## 3. Set Comprehensions

```python
+----------------------------------------------------------------+
|               SET COMPREHENSION ANATOMY                         |
+----------------------------------------------------------------+
|                                                                 |
|   {expression  for item in iterable  if condition}              |
|                                                                 |
|   Same as list comprehension but with {} instead of []          |
|   Result: unique values only, unordered                         |
|                                                                 |
+----------------------------------------------------------------+
```python

```python
# Unique k-mers in a sequence
dna = "ATGATGATGATGCGATCG"
k = 3

unique_kmers = {dna[i:i+k] for i in range(len(dna) - k + 1)}

print(f"Sequence: {dna}")
print(f"Unique {k}-mers ({len(unique_kmers)}): {sorted(unique_kmers)}")
```python

```python
# Unique amino acids in a protein
protein = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
unique_aa = {aa for aa in protein}
print(f"Unique amino acids ({len(unique_aa)}): {sorted(unique_aa)}")

# Codons that encode a specific amino acid
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

leu_codons = {codon for codon, aa in GENETIC_CODE.items() if aa == 'L'}
ser_codons = {codon for codon, aa in GENETIC_CODE.items() if aa == 'S'}
print(f"\nLeucine codons: {sorted(leu_codons)}")
print(f"Serine codons:  {sorted(ser_codons)}")
```python

---

## 4. Nested Comprehensions

```python
+----------------------------------------------------------------+
|                NESTED COMPREHENSION                              |
+----------------------------------------------------------------+
|                                                                  |
|   Flattened (single list):                                       |
|   [expr for outer in outer_list for inner in inner_list]         |
|                                                                  |
|   Nested (list of lists):                                        |
|   [[expr for inner in inner_list] for outer in outer_list]       |
|                                                                  |
|   Note: the outer loop comes FIRST in flattened form             |
|                                                                  |
+----------------------------------------------------------------+
```python

```python
# Generate all 16 dinucleotides and all 64 codons
bases = 'ATGC'

dinucleotides = [b1 + b2 for b1 in bases for b2 in bases]
print(f"All {len(dinucleotides)} dinucleotides: {dinucleotides}")

all_codons = [b1 + b2 + b3 for b1 in bases for b2 in bases for b3 in bases]
print(f"\nTotal codons: {len(all_codons)}")
print(f"First 8: {all_codons[:8]}")
print(f"Last 8:  {all_codons[-8:]}")
```python

```python
# Three reading frames of a DNA sequence (nested -- list of lists)
dna = "ATGCGATCGATCGTAGCGATCGATCG"

reading_frames = [
    [dna[i:i+3] for i in range(frame, len(dna) - 2, 3)]
    for frame in range(3)
]

print(f"DNA: {dna}")
for idx, codons in enumerate(reading_frames):
    print(f"  Frame +{idx + 1}: {codons}")

# Flatten: all codons from all frames
all_frame_codons = [codon for frame in reading_frames for codon in frame]
print(f"\nAll codons combined: {all_frame_codons}")
```python

---

## 5. Generator Expressions

A generator expression looks like a list comprehension but uses `()` instead of `[]`. It produces items **lazily** (one at a time), saving memory when working with large datasets.

```python
+----------------------------------------------------------------+
|           GENERATOR vs LIST COMPREHENSION                       |
+----------------------------------------------------------------+
|                                                                 |
|   List:      [x**2 for x in range(1000000)]                    |
|              -> Creates 1,000,000 items in memory NOW           |
|                                                                 |
|   Generator: (x**2 for x in range(1000000))                    |
|              -> Creates items ONE AT A TIME when needed         |
|                                                                 |
+----------------------------------------------------------------+
```python

```python
import sys

n = 100_000

list_comp = [x ** 2 for x in range(n)]
gen_expr  = (x ** 2 for x in range(n))

print(f"List size:      {sys.getsizeof(list_comp):>10,} bytes")
print(f"Generator size: {sys.getsizeof(gen_expr):>10,} bytes")
print(f"Memory ratio:   {sys.getsizeof(list_comp) / sys.getsizeof(gen_expr):.0f}x")
```python

```python
# Generators work seamlessly with built-in functions
sequences = ["ATGCGATCG", "GCGCGCGCGC", "ATATATATAT", "ATGCATGCAT"]

avg_gc = sum(gc_content(s) for s in sequences) / len(sequences)
max_len = max(len(s) for s in sequences)
any_pure_gc = any(set(s) <= {'G', 'C'} for s in sequences)
all_valid = all(set(s) <= {'A', 'T', 'G', 'C'} for s in sequences)

print(f"Average GC: {avg_gc:.1f}%")
print(f"Max length: {max_len} bp")
print(f"Any pure GC: {any_pure_gc}")
print(f"All valid DNA: {all_valid}")
```python

```python
# Memory-efficient processing of a large dataset
import random
random.seed(42)

large_dataset = [
    ''.join(random.choices('ATGC', k=random.randint(50, 500)))
    for _ in range(10_000)
]

# Count sequences with GC > 55% using a generator (no intermediate list)
high_gc_count = sum(1 for seq in large_dataset if gc_content(seq) > 55)

print(f"Total sequences: {len(large_dataset):,}")
print(f"Sequences with GC > 55%: {high_gc_count:,}")
print(f"Fraction: {high_gc_count / len(large_dataset):.1%}")
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
