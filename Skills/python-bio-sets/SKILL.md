---
name: python-bio-sets
description: "Dictionaries are the most powerful built-in data structure for bioinformatics. The standard genetic code is a dictionary (64 codons → 20 amino acids + stop). A gene annotation database is a dictionary"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/09_Dictionaries_and_Sets/02_sets.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 9: Dictionaries and Sets

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/09_Dictionaries_and_Sets/02_sets.ipynb`*

# Module 9: Dictionaries and Sets

---

### Learning Objectives
- Master dictionary creation, access, methods, and iteration
- Use nested dictionaries and specialized collections (`defaultdict`, `Counter`)
- Understand set operations: union, intersection, difference, and frozenset
- Apply these structures to biological problems: codon tables, gene annotations, k-mer analysis

---

## How to use this notebook

1. Section 1 covers dictionaries thoroughly: creation, access, modification, iteration, nested dicts, and `defaultdict`/`Counter`.
2. Section 2 covers sets: creation, set operations (union, intersection, difference), and frozensets.
3. Section 3 combines both: reverse mapping from amino acids to codons, pathway analysis.
4. The exercises (Section 4) require combining dicts and sets — work through them in order.

## Common stumbling points

- **`KeyError` vs `get()`:** `d["missing_key"]` raises a `KeyError`. `d.get("missing_key")` returns `None`. `d.get("missing_key", default)` returns `default`. In bioinformatics, unknown codons (e.g., "NNN") will cause `KeyError` if you use direct access on the codon table — use `.get()` with a default of `"X"` or `"?"`.
- **Dict keys must be hashable:** Lists and dicts cannot be dictionary keys. Strings, numbers, and tuples can. This is why genomic coordinates stored as tuples `("chr17", 43044295)` work as keys, but coordinate lists do not.
- **Iterating and modifying simultaneously:** You cannot change dictionary size while iterating over it. Create a separate dict or list to collect changes, then apply them afterward.
- **`in` checks keys, not values:** `"ATG" in codon_table` is True if `"ATG"` is a key. To check values, use `"Met" in codon_table.values()`.
- **Sets are unordered:** You cannot index into a set (`my_set[0]` raises `TypeError`). If you need the elements in order, convert to a sorted list: `sorted(my_set)`.
- **Set elements must be hashable:** Lists and dicts cannot be set members. Use frozensets when you need a set of sets.

---

[← Previous: Module 8: Lists and Tuples](../08_Lists_and_Tuples/01_lists_and_tuples.ipynb) | [Next: Module 10: Comprehensions →](../10_Comprehensions/01_comprehensions.ipynb)

---

## 1. Dictionaries -- Key-Value Mapping

A dictionary maps **keys** to **values**. Looking up a key is extremely fast (O(1)) regardless of dictionary size.

```python
codon_table = {"ATG": "Met", "TAA": "Stop", "GGG": "Gly"}   # str -> str
gene_lengths = {"BRCA1": 81189, "TP53": 19149}               # str -> int
```

**Key properties:**
- Keys must be unique and **hashable** (strings, numbers, tuples — not lists)
- Values can be anything, including other dicts or lists
- Insertion order is preserved (Python 3.7+)

```python
# Curly-brace literal
gene_lengths = {
    "BRCA1": 81189,
    "TP53":  19149,
    "EGFR":  188307,
    "MYC":   4150,
}
print(f"Gene lengths: {gene_lengths}")

# dict() constructor from keyword arguments
complement = dict(A='T', T='A', G='C', C='G')
print(f"Complement map: {complement}")

# dict() from list of (key, value) pairs
pairs = [("ATG", "Met"), ("TGG", "Trp"), ("TAA", "Stop")]
codon_names = dict(pairs)
print(f"Codon names: {codon_names}")

# Empty dictionary
annotations = {}
print(f"Empty dict: {annotations}")
```

### 1.2 Accessing Values

```python
gene_lengths = {"BRCA1": 81189, "TP53": 19149, "EGFR": 188307, "MYC": 4150}

# Direct access -- raises KeyError if key is missing
print(f"BRCA1 length: {gene_lengths['BRCA1']} bp")

# Safe access with .get() -- returns None (or a default) if missing
print(f"KRAS length:  {gene_lengths.get('KRAS')}")
print(f"KRAS length:  {gene_lengths.get('KRAS', 'not found')}")

# Membership test
print(f"\n'TP53' in dict: {'TP53' in gene_lengths}")
print(f"'KRAS' in dict: {'KRAS' in gene_lengths}")
```

### 1.3 Modifying Dictionaries

```python
gene_info = {"BRCA1": 81189, "TP53": 19149}

# Add a new entry
gene_info["EGFR"] = 188307
print(f"After add:    {gene_info}")

# Update an existing value
gene_info["TP53"] = 19200  # corrected length
print(f"After update: {gene_info}")

# Update multiple entries at once
gene_info.update({"MYC": 4150, "KRAS": 45693})
print(f"After batch:  {gene_info}")

# Remove an entry
removed_val = gene_info.pop("KRAS")
print(f"Removed KRAS ({removed_val} bp): {gene_info}")

# setdefault() -- only set if key is absent
gene_info.setdefault("TP53", 0)   # TP53 exists, so no change
gene_info.setdefault("RB1", 180388)  # RB1 is new
print(f"After setdefault: {gene_info}")
```

### 1.4 Dictionary Iteration

```python
gene_data = {
    "BRCA1": {"chr": "17", "length": 81189, "gc": 42.5},
    "TP53":  {"chr": "17", "length": 19149, "gc": 48.7},
    "EGFR":  {"chr": "7",  "length": 188307, "gc": 55.1},
}

# Iterate over keys (default)
print("Gene names:", list(gene_data.keys()))

# Iterate over values
total_length = sum(info["length"] for info in gene_data.values())
print(f"Total length: {total_length:,} bp")

# Iterate over key-value pairs
print(f"\n{'Gene':<8} {'Chr':<6} {'Length':>10} {'GC%':>6}")
print("-" * 34)
for gene, info in gene_data.items():
    print(f"{gene:<8} chr{info['chr']:<4} {info['length']:>9,} {info['gc']:>5.1f}%")
```

### 1.5 The Genetic Code as a Dictionary

The codon table is the quintessential biological dictionary: 64 codons map to 20 amino acids + stop signals.

```python
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

def translate(dna):
    """Translate a DNA coding sequence into a protein string."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)

cds = "ATGGCCGATCGATAG"
print(f"DNA:     {cds}")
print(f"Protein: {translate(cds)}")
```

### 1.6 Nucleotide Frequency Counter

```python
# Manual counting with a dictionary
sequence = "ATGCGATCGATCGTAGCGATCGATCGATGCGA"

freq = {}
for nt in sequence:
    freq[nt] = freq.get(nt, 0) + 1

print(f"Sequence ({len(sequence)} bp): {sequence}")
print("\nNucleotide frequencies:")
for nt in 'ATGC':
    count = freq.get(nt, 0)
    pct = count / len(sequence) * 100
    bar = '#' * count
    print(f"  {nt}: {count:>3} ({pct:5.1f}%)  {bar}")

gc_pct = (freq.get('G', 0) + freq.get('C', 0)) / len(sequence) * 100
print(f"\nGC content: {gc_pct:.1f}%")
```

### 1.7 Nested Dictionaries -- Gene Annotation Database

```python
gene_annotations = {
    "BRCA1": {
        "full_name": "breast cancer type 1 susceptibility protein",
        "chromosome": "17",
        "coordinates": (43044295, 43125483),
        "strand": "-",
        "go_terms": ["DNA repair", "tumor suppression", "transcription regulation"],
        "diseases": ["Breast cancer", "Ovarian cancer"],
    },
    "TP53": {
        "full_name": "tumor protein p53",
        "chromosome": "17",
        "coordinates": (7661779, 7687538),
        "strand": "-",
        "go_terms": ["apoptosis", "cell cycle arrest", "DNA repair"],
        "diseases": ["Li-Fraumeni syndrome", "Most cancers"],
    },
    "EGFR": {
        "full_name": "epidermal growth factor receptor",
        "chromosome": "7",
        "coordinates": (55019017, 55207337),
        "strand": "+",
        "go_terms": ["signal transduction", "cell proliferation"],
        "diseases": ["Non-small cell lung cancer", "Glioblastoma"],
    },
}

# Query the annotation database
for gene_id, info in gene_annotations.items():
    start, end = info["coordinates"]
    print(f"{gene_id} ({info['full_name']})")
    print(f"  chr{info['chromosome']}:{start:,}-{end:,} ({info['strand']})")
    print(f"  GO terms: {', '.join(info['go_terms'])}")
    print()
```

### 1.8 defaultdict and Counter

```python
from collections import defaultdict

# defaultdict automatically creates a default value for missing keys
# Example: group genes by chromosome
gene_locations = [
    ("BRCA1", "chr17"),
    ("TP53",  "chr17"),
    ("EGFR",  "chr7"),
    ("MYC",   "chr8"),
    ("KRAS",  "chr12"),
    ("RB1",   "chr13"),
    ("BRAF",  "chr7"),
]

# Without defaultdict you need: if chrom not in by_chrom: by_chrom[chrom] = []
by_chrom = defaultdict(list)
for gene, chrom in gene_locations:
    by_chrom[chrom].append(gene)

print("Genes grouped by chromosome:")
for chrom in sorted(by_chrom):
    print(f"  {chrom}: {by_chrom[chrom]}")
```

```python
from collections import Counter

# Counter -- specialized dict for counting
sequence = "ATGCGATCGATCGATCGATGCGATCG"

# Count nucleotides
nt_counts = Counter(sequence)
print(f"Nucleotide counts: {nt_counts}")
print(f"Most common: {nt_counts.most_common(2)}")

# Count k-mers
k = 3
kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
kmer_counts = Counter(kmers)

print(f"\n{k}-mer frequencies:")
for kmer, count in kmer_counts.most_common(5):
    print(f"  {kmer}: {count}")
```

```python
# Counter arithmetic -- compare k-mer profiles between two sequences
from collections import Counter

seq_a = "ATGATGATGCCC"
seq_b = "ATGATGGGGATG"

kmers_a = Counter(seq_a[i:i+3] for i in range(len(seq_a) - 2))
kmers_b = Counter(seq_b[i:i+3] for i in range(len(seq_b) - 2))

print(f"Seq A k-mers: {dict(kmers_a)}")
print(f"Seq B k-mers: {dict(kmers_b)}")

# Common k-mers (minimum of each count)
common = kmers_a & kmers_b
print(f"\nShared k-mers (min counts): {dict(common)}")

# Combined k-mers
combined = kmers_a + kmers_b
print(f"Combined k-mers (sum):      {dict(combined.most_common(5))}")
```

---

## 2. Sets -- Collections of Unique Elements

A set stores **unique, unordered** elements. Membership testing (`x in my_set`) is O(1), making sets perfect for checking whether something has been seen before.

```
+---------------------------------------------------------+
|                     SET STRUCTURE                        |
+---------------------------------------------------------+
|                                                         |
|   nucleotides = {'A', 'T', 'G', 'C'}                   |
|                                                         |
|   - Unordered (no index access)                         |
|   - Only unique elements (duplicates ignored)           |
|   - Fast membership testing O(1)                        |
|   - Supports union, intersection, difference            |
|                                                         |
+---------------------------------------------------------+
```

### 2.1 Creating Sets

```python
# From a literal
dna_bases = {'A', 'T', 'G', 'C'}
rna_bases = {'A', 'U', 'G', 'C'}

# From a sequence -- automatically removes duplicates
seq = "ATGCGATCGATCGATCGATCG"
unique_nt = set(seq)
print(f"Unique nucleotides in '{seq}': {unique_nt}")

# Empty set -- must use set(), NOT {} (that creates an empty dict!)
empty_set = set()
print(f"Empty set: {empty_set}, type: {type(empty_set)}")

# From a list of k-mers
k = 3
unique_kmers = set(seq[i:i+k] for i in range(len(seq) - k + 1))
print(f"Unique {k}-mers: {unique_kmers}")
```

### 2.2 Set Operations

```
+---------------------------------------------------------+
|                   SET OPERATIONS                         |
+---------------------------------------------------------+
|                                                         |
|   A = {1, 2, 3}       B = {2, 3, 4}                    |
|                                                         |
|   Union        (A | B):  {1, 2, 3, 4}                  |
|   Intersection (A & B):  {2, 3}                        |
|   Difference   (A - B):  {1}                           |
|   Sym. Diff    (A ^ B):  {1, 4}                        |
|                                                         |
|        +-----+     +-----+                              |
|        |  A  |     |  B  |                              |
|        |  1  | 2 3 |  4  |                              |
|        +-----+     +-----+                              |
|                                                         |
+---------------------------------------------------------+
```
