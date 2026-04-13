---
name: python-bio-operators
description: "Split from `01_operators_and_expressions.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/03_Operators_and_Expressions/01_operators.ipynb"
---

# Operators

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/03_Operators_and_Expressions/01_operators.ipynb`*

# Operators

Split from `01_operators_and_expressions.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_operators_and_expressions.ipynb) · [Next: Expressions](./02_expressions.ipynb)

## How to use this notebook

1. Run cells in order. Sections 1–7 each cover an operator category, followed by full worked examples in Section 8.
2. The exercises in Section 9 range from straightforward (GC content formula) to challenging (primer validation, sequence identity). Work through them in order.
3. Pay extra attention to Section 7 (operator precedence) — precedence bugs are among the most common silent errors in data-analysis scripts.

## Common stumbling points

- **`/` vs `//`:** Regular division (`/`) always returns a `float`. Floor division (`//`) truncates to the integer below. When working with codon positions you almost always want `//`: `seq_length // 3` gives the number of complete codons.
- **`%` is modulo, not percent:** `seq_length % 3` gives the number of leftover nucleotides after dividing into codons — it has nothing to do with the percentage symbol.
- **Precedence surprises:** `g + c / total * 100` computes as `g + ((c / total) * 100)` — wrong! The correct expression is `(g + c) / total * 100`. When in doubt, add parentheses.
- **`==` vs `is`:** Use `==` to check if two values are equal. Use `is` only to check identity (same object in memory) — and only for `None`, `True`, and `False`. Never write `if x == None`.
- **`and`/`or` short-circuit:** `A and B` does not evaluate `B` if `A` is falsy. `A or B` does not evaluate `B` if `A` is truthy. This is efficient but can mask bugs if `B` has side effects.
- **`in` checks keys in dicts:** `"ATG" in codon_table` checks whether `"ATG"` is a *key*, not a value. To check values, use `"ATG" in codon_table.values()`.

## 1. Arithmetic Operators

| Operator | Name | Example | Result |
|----------|------|---------|--------|
| `+` | Addition | `5 + 3` | `8` |
| `-` | Subtraction | `5 - 3` | `2` |
| `*` | Multiplication | `5 * 3` | `15` |
| `/` | Division | `5 / 2` | `2.5` |
| `//` | Floor division | `5 // 2` | `2` |
| `%` | Modulo (remainder) | `5 % 2` | `1` |
| `**` | Exponentiation | `5 ** 2` | `25` |

```python
# Basic arithmetic
a = 15
b = 4

print(f"{a} + {b}  = {a + b}")
print(f"{a} - {b}  = {a - b}")
print(f"{a} * {b}  = {a * b}")
print(f"{a} / {b}  = {a / b}")      # always returns float
print(f"{a} // {b} = {a // b}")     # truncates toward negative infinity
print(f"{a} % {b}  = {a % b}")      # remainder
print(f"{a} ** 2 = {a ** 2}")       # exponentiation
```

### Division: `/` vs `//`

This distinction matters constantly in bioinformatics:
- `/` (true division) always returns a `float`
- `//` (floor division) returns an `int` when both operands are `int`

```python
# True division vs floor division
print(f"7 / 2  = {7 / 2}")     # 3.5  (float)
print(f"7 // 2 = {7 // 2}")    # 3    (int)
print(f"7 % 2  = {7 % 2}")     # 1    (remainder)

# Verify: (a // b) * b + (a % b) == a
a, b = 7, 2
print(f"\nVerification: {a // b} * {b} + {a % b} = {(a // b) * b + (a % b)} == {a}")
```

### Bio application: Codons and reading frames

The modulo operator `%` and floor division `//` are essential for working with codons (triplets of nucleotides).

```python
# How many complete codons in a sequence?
seq_length = 1000

complete_codons = seq_length // 3
remaining_nucleotides = seq_length % 3

print(f"Sequence length: {seq_length} nt")
print(f"Complete codons: {complete_codons}")
print(f"Leftover nucleotides: {remaining_nucleotides}")
```

```python
# Which reading frame is a given position in?
# Reading frames: 0, 1, 2
positions = [0, 1, 2, 3, 6, 10, 47]

print(f"{'Position':>10} {'Reading Frame':>15}")
print("-" * 27)
for pos in positions:
    frame = pos % 3
    print(f"{pos:>10} {frame:>15}")
```

### Bio application: GC content calculation

GC content = (G + C) / total length * 100%

```python
sequence = "ATGCGATCGATCGTAGC"

g_count = sequence.count("G")
c_count = sequence.count("C")
total = len(sequence)

gc_content = (g_count + c_count) / total * 100

print(f"Sequence: {sequence}")
print(f"G: {g_count}, C: {c_count}, Total: {total}")
print(f"GC content: {gc_content:.2f}%")
```

### Bio application: Protein molecular weight

A rough estimate of protein molecular weight:
- Average amino acid mass: ~110 Da
- Water lost per peptide bond: 18 Da
- Number of peptide bonds = n - 1 (for n amino acids)

```python
num_amino_acids = 393   # human p53 protein
avg_aa_mass = 110       # Daltons
water_mass = 18         # Daltons

protein_mw = (num_amino_acids * avg_aa_mass) - ((num_amino_acids - 1) * water_mass)

print(f"Protein: p53 ({num_amino_acids} amino acids)")
print(f"Estimated MW: {protein_mw:,} Da")
print(f"Estimated MW: {protein_mw / 1000:.1f} kDa")
print(f"(Actual MW of p53: ~43.7 kDa)")
```

### Augmented assignment operators

Python provides shorthand for updating a variable:

| Shorthand | Equivalent |
|-----------|------------|
| `x += 5` | `x = x + 5` |
| `x -= 3` | `x = x - 3` |
| `x *= 2` | `x = x * 2` |
| `x /= 4` | `x = x / 4` |
| `x //= 3` | `x = x // 3` |
| `x %= 3` | `x = x % 3` |
| `x **= 2` | `x = x ** 2` |

```python
# Count GC nucleotides using augmented assignment
sequence = "ATGCGATCGATCGTAGC"
gc_count = 0

for nucleotide in sequence:
    if nucleotide in "GC":
        gc_count += 1    # same as gc_count = gc_count + 1

print(f"GC count: {gc_count}")
print(f"GC content: {gc_count / len(sequence) * 100:.1f}%")
```

---

## 2. Comparison Operators

Comparison operators return `True` or `False`. They are the foundation of all decision-making in programs.

| Operator | Meaning | Example |
|----------|---------|----------|
| `==` | Equal to | `5 == 5` -> `True` |
| `!=` | Not equal to | `5 != 3` -> `True` |
| `<` | Less than | `3 < 5` -> `True` |
| `>` | Greater than | `5 > 3` -> `True` |
| `<=` | Less than or equal | `5 <= 5` -> `True` |
| `>=` | Greater than or equal | `5 >= 3` -> `True` |

```python
# Comparing biological values
gc_content = 0.65

print(f"GC content = {gc_content}")
print(f"GC > 0.5?    {gc_content > 0.5}")      # True
print(f"GC < 0.4?    {gc_content < 0.4}")      # False
print(f"GC == 0.65?  {gc_content == 0.65}")    # True
print(f"GC != 0.5?   {gc_content != 0.5}")     # True
```

### Chained comparisons

Python allows chaining comparisons, which is more readable than using `and`.

```python
gc_content = 0.52

# Check if GC content is in the "normal" range (40-60%)
# Without chaining:
in_range_1 = gc_content >= 0.40 and gc_content <= 0.60

# With chaining (more Pythonic):
in_range_2 = 0.40 <= gc_content <= 0.60

print(f"GC = {gc_content}")
print(f"In normal range (40-60%)? {in_range_2}")
```

```python
# Bio application: classifying GC content
def classify_gc(gc_percent):
    """Classify GC content into biological categories."""
    if gc_percent < 30:
        return "Very AT-rich (e.g., Plasmodium)"
    elif 30 <= gc_percent < 40:
        return "AT-rich"
    elif 40 <= gc_percent < 60:
        return "Moderate (typical for many organisms)"
    elif 60 <= gc_percent < 70:
        return "GC-rich"
    else:
        return "Very GC-rich (e.g., Streptomyces)"

test_values = [19.4, 35.0, 41.0, 50.7, 65.0, 72.0]
for gc in test_values:
    print(f"  GC {gc:.1f}% -> {classify_gc(gc)}")
```

### Comparing strings

Strings are compared lexicographically (dictionary order), character by character.

```python
# String comparisons
print(f"'ATG' == 'ATG': {'ATG' == 'ATG'}")   # True (exact match)
print(f"'ATG' == 'atg': {'ATG' == 'atg'}")   # False (case-sensitive!)
print(f"'ATG' < 'GGG':  {'ATG' < 'GGG'}")    # True (A < G in ASCII)

# Case-insensitive comparison
seq1 = "ATGC"
seq2 = "atgc"
print(f"\nCase-insensitive match: {seq1.upper() == seq2.upper()}")
```

---

## 3. Logical Operators

Logical operators combine boolean expressions:

| Operator | Meaning | Example |
|----------|---------|--------|
| `and` | Both must be True | `True and False` -> `False` |
| `or` | At least one True | `True or False` -> `True` |
| `not` | Inverts the value | `not True` -> `False` |

### Truth tables

```
A     B     A and B    A or B    not A
True  True  True       True      False
True  False False      True      False
False True  False      True      True
False False False      False     True
```

```python
# Sequence quality filters
seq_length = 1500
gc_content = 0.48
has_ambiguous_bases = False

# A sequence passes QC if:
# - length is between 500 and 3000 bp
# - GC content is between 30% and 70%
# - no ambiguous bases

valid_length = 500 <= seq_length <= 3000
valid_gc = 0.30 <= gc_content <= 0.70
no_ambiguity = not has_ambiguous_bases

passes_qc = valid_length and valid_gc and no_ambiguity

print(f"Length ({seq_length}): valid? {valid_length}")
print(f"GC ({gc_content}):   valid? {valid_gc}")
print(f"No ambiguity:   {no_ambiguity}")
print(f"Passes QC:      {passes_qc}")
```

```python
# Bio application: checking if a codon is a stop codon
codon = "TAG"

is_stop = (codon == "TAA") or (codon == "TAG") or (codon == "TGA")
print(f"Is '{codon}' a stop codon? {is_stop}")

# More Pythonic way (using 'in'):
is_stop = codon in ["TAA", "TAG", "TGA"]
print(f"Is '{codon}' a stop codon? {is_stop}")
```

```python
# Bio application: filter sequences by multiple criteria
sequences = [
    ("seq_1", "ATGCGATCGATCGATCG", 0.53),
    ("seq_2", "AAATTTAAATTTAAATTT", 0.0),
    ("seq_3", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", 1.0),
    ("seq_4", "ATGC", 0.5),
    ("seq_5", "ATGAAACCCGGGTAA", 0.53),
]

print("Sequences passing filters (length >= 10, 0.3 <= GC <= 0.7):")
for name, seq, gc in sequences:
    if len(seq) >= 10 and 0.3 <= gc <= 0.7:
        print(f"  {name}: length={len(seq)}, GC={gc:.0%}")
```

---

## 4. Membership Operators: `in` and `not in`

The `in` operator checks whether a value exists in a collection (string, list, set, dictionary keys).

```python
# Check if a substring exists in a DNA sequence
dna = "ATGAAACCCGGGTAA"

print(f"'ATG' in sequence? {'ATG' in dna}")     # True (start codon present)
print(f"'GGG' in sequence? {'GGG' in dna}")     # True
print(f"'TTT' in sequence? {'TTT' in dna}")     # False
print(f"'TTT' not in sequence? {'TTT' not in dna}")  # True
```

```python
# Check membership in a list of stop codons
stop_codons = ["TAA", "TAG", "TGA"]

test_codons = ["ATG", "TAA", "GGG", "TGA", "CCC", "TAG"]

for codon in test_codons:
    status = "STOP" if codon in stop_codons else "coding"
    print(f"  {codon} -> {status}")
```

```python
# Validate a DNA sequence using 'in'
def is_valid_dna(sequence):
    """Check if a sequence contains only valid DNA nucleotides."""
    valid = "ATGC"
    for char in sequence.upper():
        if char not in valid:
            return False
    return True

print(is_valid_dna("ATGCGATCG"))    # True
print(is_valid_dna("ATGCXYZ"))      # False
print(is_valid_dna("augcgaucg"))    # True (after upper())
```

```python
# Membership in dictionaries checks KEYS, not values
codon_table = {
    "ATG": "Met",
    "TAA": "Stop",
    "TAG": "Stop",
    "TGA": "Stop",
    "GGG": "Gly",
}

print(f"'ATG' in codon_table? {'ATG' in codon_table}")     # True (it is a key)
print(f"'Met' in codon_table? {'Met' in codon_table}")     # False (it is a value, not a key)
```

### Performance note: `in` with sets vs lists

Checking membership in a `set` is much faster than in a `list`. For large collections, always convert to a set first.

```python
# For small collections, the difference is negligible
# For large collections, sets are dramatically faster

valid_nucleotides_list = ["A", "T", "G", "C"]
valid_nucleotides_set = {"A", "T", "G", "C"}

# Both work the same way:
print("G" in valid_nucleotides_list)  # True (scans the list)
print("G" in valid_nucleotides_set)   # True (hash lookup -- much faster)

# Rule of thumb: use sets when you only need membership testing
```

---

## 5. String Operators for Sequences

Strings support several operators that are very useful for working with biological sequences.

### Concatenation (`+`)

Join two strings end-to-end.

```python
# Concatenating sequence fragments
exon1 = "ATGAAA"
exon2 = "CCCGGG"
exon3 = "TTTTAA"

# After splicing, exons are joined
mrna = exon1 + exon2 + exon3
print(f"Exon 1: {exon1}")
print(f"Exon 2: {exon2}")
print(f"Exon 3: {exon3}")
print(f"Spliced mRNA: {mrna}")
```
