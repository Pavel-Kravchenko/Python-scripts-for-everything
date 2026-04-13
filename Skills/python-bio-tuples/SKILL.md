---
name: python-bio-tuples
description: "Lists are the go-to container for ordered, variable-length collections in Python. In bioinformatics: a list of codons extracted from a CDS, a list of gene names from a FASTA file, a list of quality sc"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/08_Lists_and_Tuples/02_tuples.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 8: Lists and Tuples

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/08_Lists_and_Tuples/02_tuples.ipynb`*


---

### Learning Objectives
- Master list creation, indexing, slicing, and core methods
- Understand tuple immutability, packing/unpacking, and named tuples
- Apply lists and tuples to biological data: codons, gene coordinates, sequence analysis
- Sort and filter biological sequences by length and composition

---

## How to use this notebook

1. Sections 1–2 cover lists and tuples respectively. Run them in order.
2. Section 3 covers nested structures (lists of tuples, lists of lists) which are used for expression matrices and coordinate tables.
3. The exercises in Section 5 use biological data — k-mer extraction, sliding window GC content, codon frequency tables, and gene overlap detection.

## Common stumbling points

- **Assignment does not copy a list:** `alias = original` makes both names point to the same list. `alias.append("X")` changes `original` too. To get an independent copy, use `copy = original[:]` or `copy = original.copy()`.
- **`append()` vs `extend()`:** `list.append(x)` adds `x` as a single element. `list.extend(other_list)` adds each element of `other_list` individually. `genes.append(["BRCA1", "TP53"])` adds a *list* as one element; `genes.extend(["BRCA1", "TP53"])` adds two separate elements.
- **`sort()` modifies in place, `sorted()` returns a new list:** `my_list.sort()` changes `my_list` and returns `None`. `sorted(my_list)` leaves `my_list` untouched and returns a new sorted list.
- **Index out of range:** Accessing `my_list[len(my_list)]` raises `IndexError` because the last valid index is `len(my_list) - 1`. Negative indices count from the end: `my_list[-1]` is always the last element.
- **Tuples are immutable:** `snp = ("chr7", 55259515, "T", "G")` — you cannot change `snp[2]`. This is intentional: use tuples for fixed records (coordinates, SNP data) and lists for collections you will modify.

---

[← Previous: Module 7: File Operations](../07_File_Operations/01_file_operations.ipynb) | [Next: Module 9: Dictionaries and Sets →](../09_Dictionaries_and_Sets/01_dictionaries_and_sets.ipynb)

---

## 1. Lists -- Mutable Ordered Sequences

A **list** is the workhorse data structure in Python. It holds an ordered collection of items that can be **modified** (added to, removed from, sorted).

```python
genes = ["BRCA1", "TP53", "EGFR"]   # list of strings
quality_scores = [38, 40, 37, 42, 35]  # list of ints
mixed = ["chr17", 43044295, True]      # lists can hold different types
empty = []                              # empty list
```python

Lists are **zero-indexed** — the first element is at index `0`.

```python
# Creating lists in different ways

# From literal values
genes = ["BRCA1", "TP53", "EGFR", "MYC"]
gc_values = [42.5, 48.7, 55.1, 51.2]

# Empty list
mutations = []

# From a string -- split a DNA sequence into individual nucleotides
dna = "ATGCGA"
nucleotides = list(dna)
print(f"Nucleotides: {nucleotides}")

# From range -- useful for generating positions
positions = list(range(0, 30, 3))  # every 3rd position (codon starts)
print(f"Codon start positions: {positions}")
```python

### 1.2 Indexing and Slicing

```python
# Positive and negative indexing
promoter = list("GCGCTATAAAGCGCATGC")

print(f"First nucleotide:  {promoter[0]}")
print(f"Fifth nucleotide:  {promoter[4]}")
print(f"Last nucleotide:   {promoter[-1]}")
print(f"Second-to-last:    {promoter[-2]}")
```python

```python
# Slicing: list[start:stop:step]
codons_dna = "ATGGCCGATCGATAG"
seq = list(codons_dna)

print(f"First codon:  {''.join(seq[0:3])}")
print(f"Second codon: {''.join(seq[3:6])}")
print(f"Last codon:   {''.join(seq[-3:])}")
print(f"Every other:  {''.join(seq[::2])}")
print(f"Reversed:     {''.join(seq[::-1])}")
```python

### 1.3 Storing Sequences as Lists of Codons

A common bioinformatics pattern: split a coding DNA sequence into its triplet codons.

```python
# Split a coding sequence into codons
cds = "ATGGCCGATCGATAGCCA"

codons = []
for i in range(0, len(cds) - len(cds) % 3, 3):
    codons.append(cds[i:i+3])

print(f"CDS:    {cds}")
print(f"Codons: {codons}")
print(f"Number of complete codons: {len(codons)}")
```python

### 1.4 Essential List Methods

```python
# Adding elements
discovered_genes = ["BRCA1", "TP53"]

# append() -- add a single element to the end
discovered_genes.append("EGFR")
print(f"After append: {discovered_genes}")

# extend() -- add multiple elements from an iterable
discovered_genes.extend(["MYC", "KRAS"])
print(f"After extend: {discovered_genes}")

# insert() -- add at a specific position
discovered_genes.insert(1, "RB1")
print(f"After insert: {discovered_genes}")
```python

```python
# Removing elements
samples = ["tumor", "normal", "tumor", "metastasis"]

samples.remove("tumor")  # removes first occurrence
print(f"After remove: {samples}")

last = samples.pop()     # remove and return last element
print(f"Popped: {last}, Remaining: {samples}")

# Searching
codons = ["ATG", "GCC", "GAT", "CGA", "ATG", "TAG"]
print(f"\nIndex of 'GAT':  {codons.index('GAT')}")
print(f"Count of 'ATG':  {codons.count('ATG')}")
print(f"'TAG' in codons: {'TAG' in codons}")
```python

### 1.5 Sorting Sequences by Length

```python
# sort() modifies the list in place
gene_lengths = [81189, 19149, 188307, 4150, 45693]
gene_lengths.sort()
print(f"Ascending:  {gene_lengths}")

gene_lengths.sort(reverse=True)
print(f"Descending: {gene_lengths}")

# sorted() returns a NEW list -- original is unchanged
sequences = ["ATGCGATCGATCG", "GCG", "ATATATATATATAT", "ATGCATGC"]
by_length = sorted(sequences, key=len)
print(f"\nSorted by length:")
for seq in by_length:
    print(f"  {seq} ({len(seq)} bp)")
```python

```python
# Sort by GC content
sequences = ["ATATATAT", "GCGCGCGC", "ATGCATGC", "AAAGGGCCC", "TTTTAAAA"]

def gc_content(seq):
    """Calculate GC content as a fraction."""
    return (seq.count('G') + seq.count('C')) / len(seq)

by_gc = sorted(sequences, key=gc_content)
print("Sorted by GC content:")
for seq in by_gc:
    print(f"  {seq:15s}  GC = {gc_content(seq)*100:.1f}%")
```python

### 1.6 Finding Longest and Shortest Genes

```python
gene_data = [
    ("BRCA1", 81189),
    ("TP53",  19149),
    ("EGFR",  188307),
    ("MYC",   4150),
    ("KRAS",  45693),
]

longest  = max(gene_data, key=lambda g: g[1])
shortest = min(gene_data, key=lambda g: g[1])

print(f"Longest gene:  {longest[0]} ({longest[1]:,} bp)")
print(f"Shortest gene: {shortest[0]} ({shortest[1]:,} bp)")

total = sum(length for _, length in gene_data)
avg   = total / len(gene_data)
print(f"\nTotal length:   {total:,} bp")
print(f"Average length: {avg:,.0f} bp")
```python

### 1.7 List Copying -- Avoiding a Common Trap

```python
# Assignment does NOT copy -- both names point to the same object
original = ["ATG", "GCC", "GAT"]
alias = original          # alias IS original
alias.append("TAG")
print(f"Original after alias.append: {original}")  # modified!

# Shallow copy -- creates an independent copy
original = ["ATG", "GCC", "GAT"]
copy1 = original.copy()   # method
copy2 = original[:]       # slice

copy1.append("TAG")
print(f"Original after copy1.append: {original}")  # unchanged
print(f"copy1: {copy1}")
```python

---

## 2. Tuples -- Immutable Sequences

A **tuple** is like a list that cannot be changed after creation. This makes it ideal for fixed records such as gene coordinates.

```python
+---------------------------------------------------------+
|                    TUPLE STRUCTURE                       |
+---------------------------------------------------------+
|                                                         |
|   gene_loc = ("BRCA1", "chr17", 43044295, 43125483)    |
|                 |        |         |          |         |
|                 |        |         |          +- end    |
|                 |        |         +- start             |
|                 |        +- chromosome                  |
|                 +- gene name                            |
|                                                         |
|   - Immutable (cannot add, remove, or change elements)  |
|   - Hashable (can be used as dictionary keys)           |
|   - Slightly faster and smaller than lists              |
|   - Protects data from accidental modification          |
|                                                         |
+---------------------------------------------------------+
```python

### 2.1 Creating Tuples and Immutability

```python
# Parentheses (optional but conventional)
gene_loc = ("BRCA1", "chr17", 43044295, 43125483)
print(f"Gene location: {gene_loc}")

# Without parentheses -- the comma makes it a tuple
coordinates = 43044295, 43125483
print(f"Coordinates: {coordinates}, type: {type(coordinates)}")

# Single-element tuple -- trailing comma is required!
single = ("BRCA1",)
not_a_tuple = ("BRCA1")   # this is just a string in parentheses
print(f"Single tuple: {single}, type: {type(single)}")
print(f"Not a tuple:  {not_a_tuple}, type: {type(not_a_tuple)}")
```python

```python
# Tuples cannot be modified
snp_position = ("chr7", 55259515, "T", "G")  # chromosome, position, ref, alt

# These would all raise TypeError:
# snp_position[2] = "C"
# snp_position.append("rs123")

# To "update" a tuple, create a new one
corrected = snp_position[:2] + ("C", "A") + snp_position[4:]
print(f"Original:  {snp_position}")
print(f"Corrected: {corrected}")
```python

### 2.2 Tuple Packing and Unpacking

```python
# Packing -- multiple values into one tuple
gene_record = "TP53", "chr17", 7661779, 7687538, "-"
print(f"Packed: {gene_record}")

# Unpacking -- one tuple into multiple variables
name, chrom, start, end, strand = gene_record
print(f"Gene {name} on {chrom}:{start}-{end} (strand {strand})")
print(f"Length: {end - start:,} bp")

# Star unpacking -- capture remaining items
header = ("BRCA1", "breast cancer type 1", "chr17", 43044295, 43125483, "-")
gene_id, description, *location_info = header
print(f"\nGene ID:  {gene_id}")
print(f"Location: {location_info}")

# Swap two variables
a, b = "forward", "reverse"
a, b = b, a
print(f"After swap: a={a}, b={b}")
```python

### 2.3 Tuples for Biological Records: (gene_name, start, end)

```python
# List of gene coordinate tuples
genes = [
    ("BRCA1", 43044295, 43125483),
    ("TP53",  7661779,  7687538),
    ("EGFR",  55019017, 55207337),
    ("MYC",   127735434, 127742951),
    ("KRAS",  25204789, 25250936),
]

print(f"{'Gene':<8} {'Start':>12} {'End':>12} {'Length':>10}")
print("-" * 46)
for name, start, end in genes:
    length = end - start
    print(f"{name:<8} {start:>12,} {end:>12,} {length:>9,} bp")

# Sort genes by length
sorted_genes = sorted(genes, key=lambda g: g[2] - g[1], reverse=True)
print(f"\nLongest gene: {sorted_genes[0][0]} ({sorted_genes[0][2] - sorted_genes[0][1]:,} bp)")
```python

### 2.4 Named Tuples -- Readable Records

```python
from collections import namedtuple

# Define a Gene record type
Gene = namedtuple('Gene', ['name', 'chromosome', 'start', 'end', 'strand'])

brca1 = Gene('BRCA1', 'chr17', 43044295, 43125483, '-')
tp53  = Gene('TP53',  'chr17', 7661779,  7687538,  '-')
egfr  = Gene('EGFR',  'chr7',  55019017, 55207337, '+')

# Access by name (much more readable than index)
print(f"Gene: {brca1.name}")
print(f"Location: {brca1.chromosome}:{brca1.start}-{brca1.end}")
print(f"Length: {brca1.end - brca1.start:,} bp")

# Named tuples work everywhere regular tuples do
gene_db = [brca1, tp53, egfr]
by_length = sorted(gene_db, key=lambda g: g.end - g.start, reverse=True)

print("\nGenes sorted by length:")
for gene in by_length:
    print(f"  {gene.name:6s} on {gene.chromosome}: {gene.end - gene.start:>8,} bp")
```python

### 2.5 When to Use Tuples vs. Lists

```python
+---------------------------------------------------------+
|          TUPLES vs LISTS -- WHEN TO USE                  |
+---------------------------------------------------------+
|                                                         |
|  USE TUPLES FOR:                                        |
|    - Fixed records (gene: name, chr, start, end)        |
|    - Returning multiple values from functions            |
|    - Dictionary keys (must be hashable)                  |
|    - Data that should not change (coordinates, SNPs)     |
|                                                         |
|  USE LISTS FOR:                                         |
|    - Collections that grow or shrink over time           |
|    - Sequences you need to sort or reorder               |
|    - Homogeneous data you will iterate over              |
|    - Intermediate results during computation             |
|                                                         |
+---------------------------------------------------------+
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
