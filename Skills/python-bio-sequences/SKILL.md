---
name: python-bio-sequences
description: "By the end of this module, you will be able to: - Use Python string methods to manipulate DNA, RNA, and protein sequences - Apply slicing and indexing to extract subsequences and codons - Format scien"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/04_Strings_and_Sequences/02_sequences.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 4: Strings and Biological Sequences

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/04_Strings_and_Sequences/02_sequences.ipynb`*

# Module 4: Strings and Biological Sequences

---

### Learning Objectives

By the end of this module, you will be able to:
- Use Python string methods to manipulate DNA, RNA, and protein sequences
- Apply slicing and indexing to extract subsequences and codons
- Format scientific output with f-strings
- Perform fundamental bioinformatics operations: reverse complement, transcription, motif finding
- Parse biological sequence formats

### Prerequisites
- Basic Python variables and data types
- Familiarity with DNA/RNA/protein concepts

---

## How to use this notebook

1. Run cells in order — each section builds on techniques introduced earlier.
2. Section 5 (the `translate()` method for complements) is easy to get wrong — read the explanation carefully before running the code.
3. The exercises in Section 9 have both a blank version and a solution version. Try the blank version first; only look at the solution if you are stuck after a few minutes.
4. Exercise difficulty is marked: `*` = straightforward, `**` = requires combining techniques, `***` = challenging.

## Common stumbling points

- **Zero-based indexing:** `dna[0]` is the *first* character, not `dna[1]`. The last character is `dna[-1]` or `dna[len(dna)-1]`. Off-by-one errors in indexing are the most common beginner bug.
- **Stop index is exclusive in slices:** `dna[0:3]` gives characters at positions 0, 1, 2 — NOT 3. The stop index is always excluded. So `dna[0:3]` gives the first 3 characters.
- **Why chained `replace()` fails for reverse complement:** `dna.replace("A","T").replace("T","A")` does not work because the second `replace()` overwrites the result of the first. The `str.maketrans()` + `translate()` pattern does all substitutions simultaneously.
- **`find()` returns -1, not None:** If the motif is not found, `find()` returns `-1`, not `None` or `False`. Check with `if pos != -1:`, not `if pos:` (because `if -1:` is truthy!).
- **Strings are not lists:** `for nuc in dna` works, but `dna[0] = "G"` raises a `TypeError`. To "modify" a string, construct a new one.
- **`count()` can overlap or not:** `"ATATATAT".count("ATAT")` returns `2` (non-overlapping). If you need overlapping matches, use a loop with `find()` and `start` parameter.

Strings are **immutable** in Python -- once created, individual characters cannot be changed in place. Every string operation returns a **new** string.

```python
dna = "ATGC"
dna[0] = "C"   # TypeError! Strings cannot be modified in place.
dna = "C" + dna[1:]  # This works -- creates a new string.
```

## 2. Indexing: Accessing Individual Nucleotides

Every character in a string has a **position** (index). Python uses **0-based indexing**.

```
Sequence:   A   T   G   C   G   A   T   C
Index:      0   1   2   3   4   5   6   7
Negative:  -8  -7  -6  -5  -4  -3  -2  -1
```

Negative indices count from the end of the string, which is very convenient for accessing the last few characters.

```python
gene = "ATGCGATCGATCGTAGC"

# Access individual nucleotides
print(f"Sequence:         {gene}")
print(f"First nucleotide: {gene[0]}")
print(f"Third nucleotide: {gene[2]}")
print(f"Last nucleotide:  {gene[-1]}")
print(f"Second to last:   {gene[-2]}")
```

## 3. Slicing: Extracting Subsequences

**Syntax:** `string[start:stop:step]`

- `start` -- beginning index (inclusive, default 0)
- `stop` -- ending index (**exclusive**)
- `step` -- step size (default 1)

Think of the indices as pointing *between* characters:

```
 +---+---+---+---+---+---+
 | A | T | G | C | G | A |
 +---+---+---+---+---+---+
 0   1   2   3   4   5   6
```

```python
seq = "ATGCGATCGATCGTAGC"
print(f"Full sequence: {seq}")
print(f"Length: {len(seq)} bp")
print()

# Basic slicing
print(f"First 3 nucleotides (start codon): {seq[0:3]}")
print(f"Positions 3-8:                     {seq[3:9]}")
print(f"Last 6 nucleotides:                {seq[-6:]}")
print(f"Everything except first 3:         {seq[3:]}")
print(f"Everything except last 3:          {seq[:-3]}")
```

```python
# Using step parameter
seq = "ATGCGATCGATCGTAGC"

# Every other nucleotide
print(f"Every 2nd nucleotide:     {seq[::2]}")

# First position of each codon (every 3rd, starting at 0)
print(f"1st codon positions:      {seq[0::3]}")
print(f"2nd codon positions:      {seq[1::3]}")
print(f"3rd codon positions:      {seq[2::3]}")
```

### Extracting Codons with Slicing

Codons are triplets of nucleotides that encode amino acids. Extracting them from a sequence is one of the most common operations in bioinformatics.

```python
mrna = "AUGGCCGAUUAGCCAUAG"

# Extract individual codons by position
codon1 = mrna[0:3]
codon2 = mrna[3:6]
codon3 = mrna[6:9]

print(f"mRNA: {mrna}")
print(f"Codon 1 (start): {codon1}")
print(f"Codon 2:         {codon2}")
print(f"Codon 3:         {codon3}")

# Extract ALL codons using a list comprehension
codons = [mrna[i:i+3] for i in range(0, len(mrna), 3)]
print(f"\nAll codons: {codons}")
```

### Reversing a Sequence

A step of `-1` reverses the string. This is essential for working with the complementary DNA strand, which runs in the opposite direction (3' to 5').

```python
seq = "ATGCGATC"
reversed_seq = seq[::-1]

print(f"Original: 5'-{seq}-3'")
print(f"Reversed: 3'-{reversed_seq}-5'")

# This is also how you check if a sequence is a palindrome
test = "GAATTC"  # EcoRI recognition site
print(f"\n'{test}' reversed: '{test[::-1]}'")
print(f"Is it a string palindrome? {test == test[::-1]}")
```

## 4. Essential String Methods

Python provides many built-in methods for string manipulation. Here are the most important ones for bioinformatics.

### 4.1 Case Conversion: `upper()`, `lower()`

Sequences from databases may come in mixed case. Always normalize before analysis.

```python
# Sequences from different sources may have inconsistent case
messy_seq = "AtGcGaTcGaTaGc"

print(f"Original:  {messy_seq}")
print(f"Upper:     {messy_seq.upper()}")
print(f"Lower:     {messy_seq.lower()}")

# Always normalize input before processing
clean_seq = messy_seq.upper()
print(f"\nCleaned:   {clean_seq}")
print(f"G count:   {clean_seq.count('G')}")
```

### 4.2 Counting: `count()`

Count occurrences of a substring. Perfect for nucleotide composition and GC content.

```python
dna = "ATGCGATCGATCGTAGCATGCATGCA"

# Count individual nucleotides
a_count = dna.count('A')
t_count = dna.count('T')
g_count = dna.count('G')
c_count = dna.count('C')

print(f"Sequence: {dna}")
print(f"Length:   {len(dna)} bp")
print(f"\nNucleotide composition:")
print(f"  A: {a_count:3d} ({a_count/len(dna)*100:5.1f}%)")
print(f"  T: {t_count:3d} ({t_count/len(dna)*100:5.1f}%)")
print(f"  G: {g_count:3d} ({g_count/len(dna)*100:5.1f}%)")
print(f"  C: {c_count:3d} ({c_count/len(dna)*100:5.1f}%)")

gc_content = (g_count + c_count) / len(dna) * 100
print(f"\nGC Content: {gc_content:.1f}%")

# Count a longer motif
print(f"\nOccurrences of 'ATG': {dna.count('ATG')}")
print(f"Occurrences of 'CG' (CpG): {dna.count('CG')}")
```

### 4.3 Finding Subsequences: `find()`, `rfind()`

`find()` returns the index of the first occurrence of a substring, or `-1` if not found. `rfind()` searches from the right.

```python
dna = "ATGCGATCGATCGTAGCATGCATGCA"

# Find the first start codon
pos = dna.find("ATG")
print(f"Sequence: {dna}")
print(f"First 'ATG' at position: {pos}")

# Find starting from a specific position
pos2 = dna.find("ATG", pos + 1)
print(f"Next 'ATG' at position:  {pos2}")

# rfind: search from the right
last_pos = dna.rfind("ATG")
print(f"Last 'ATG' at position:  {last_pos}")

# Motif not found
print(f"'TATA' found at: {dna.find('TATA')}  (not found)")
```

```python
# Find ALL positions of a motif using a while loop
dna = "ATGCGATCGATCGTAGCATGCATGCA"
motif = "ATG"

positions = []
pos = dna.find(motif)
while pos != -1:
    positions.append(pos)
    pos = dna.find(motif, pos + 1)

print(f"All positions of '{motif}': {positions}")

# Visualize the matches
print(f"\n{dna}")
marker = list('.' * len(dna))
for p in positions:
    for i in range(len(motif)):
        marker[p + i] = '^'
print(''.join(marker))
```

### 4.4 Replacing: `replace()`

Replace all occurrences of a substring with another. The most direct way to perform DNA-to-RNA transcription.

```python
dna = "ATGCGATCGATCGTAG"

# Transcription: DNA -> RNA (replace T with U)
rna = dna.replace('T', 'U')
print(f"DNA: {dna}")
print(f"RNA: {rna}")

# Replace a motif
mutated = dna.replace('CGA', 'CTA')  # simulating a mutation
print(f"\nOriginal: {dna}")
print(f"Mutated:  {mutated}")

# Limit the number of replacements
single_mutation = dna.replace('CGA', 'CTA', 1)  # only first occurrence
print(f"Single:   {single_mutation}")
```

### 4.5 Splitting and Joining: `split()`, `join()`

These are essential for parsing biological data files (FASTA headers, tab-separated data, etc.).

```python
# Splitting a FASTA header
header = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens"

# Split by pipe character
parts = header.split('|')
print(f"Header: {header}")
print(f"Parts:  {parts}")
print(f"Database: {parts[0][1:]}")
print(f"Accession: {parts[1]}")
print(f"Description: {parts[2]}")

# Split by whitespace (default)
words = header.split()
print(f"\nWords: {words}")
```

```python
# Joining: assemble a sequence from codons
codons = ['ATG', 'GCC', 'GAT', 'TAG']

# Join with no separator -> continuous sequence
sequence = ''.join(codons)
print(f"Codons:   {codons}")
print(f"Sequence: {sequence}")

# Join with separator for display
formatted = '-'.join(codons)
print(f"Formatted: {formatted}")

# Build a FASTA record
header = ">seq1 Example sequence"
fasta = '\n'.join([header, sequence])
print(f"\nFASTA format:\n{fasta}")
```

### 4.6 Stripping Whitespace: `strip()`, `lstrip()`, `rstrip()`

When reading sequences from files, lines often come with trailing newlines or spaces. Always strip them before processing.

```python
# Simulating lines read from a file
line1 = "  ATGCGATCG\n"
line2 = "\tGCTAGCTA  \n"

print(f"Raw line 1: '{line1}'")
print(f"Stripped:   '{line1.strip()}'")
print(f"Left only:  '{line1.lstrip()}'")
print(f"Right only: '{line1.rstrip()}'")
print()
print(f"Raw line 2: '{line2}'")
print(f"Stripped:   '{line2.strip()}'")

# Typical file-reading pattern
raw_lines = ["ATGCGA\n", "TCGATC\n", "GTAGCA\n"]
full_seq = ''.join(line.strip() for line in raw_lines)
print(f"\nAssembled sequence: {full_seq}")
```

### 4.7 Testing String Properties: `startswith()`, `endswith()`, `in`

Useful for checking sequence features and parsing formatted data.

```python
seq = "ATGGCCGATTAGCCA"

# Does the sequence start with a start codon?
print(f"Sequence: {seq}")
print(f"Starts with ATG: {seq.startswith('ATG')}")

# Check for stop codons at the end
stop_codons = ('TAA', 'TAG', 'TGA')
last_codon = seq[-3:]
print(f"Last codon: {last_codon}")
print(f"Ends with stop codon: {last_codon in stop_codons}")

# Check membership with 'in'
print(f"\nContains 'GATTAG': {'GATTAG' in seq}")
print(f"Contains 'TATA':   {'TATA' in seq}")

# Validate a sequence contains only valid DNA characters
valid_dna = all(nuc in 'ATGC' for nuc in seq)
print(f"\nAll characters valid DNA: {valid_dna}")

bad_seq = "ATGCXGATC"
valid_bad = all(nuc in 'ATGC' for nuc in bad_seq)
print(f"'{bad_seq}' valid DNA: {valid_bad}")
```

## 5. The `translate()` Method and Complement

For the DNA complement, we need to swap A<->T and G<->C simultaneously. The `replace()` method cannot do this in one pass (replacing A->T would then be replaced back by T->A). Python's `str.translate()` with `str.maketrans()` solves this elegantly.

```python
# Why replace() fails for complement
dna = "ATGC"

# Attempt with chained replace:
wrong = dna.replace('A', 'T').replace('T', 'A')  # A->T, then ALL T->A
print(f"Wrong approach: {dna} -> {wrong}  (both A and T became A!)")

# Correct approach: str.translate()
complement_table = str.maketrans('ATGC', 'TACG')
correct = dna.translate(complement_table)
print(f"Correct:        {dna} -> {correct}")
```

```python
# The complete reverse complement operation
dna = "ATGCGATCGATCGTAG"

complement_table = str.maketrans('ATGC', 'TACG')

# Step 1: Complement
complement = dna.translate(complement_table)

# Step 2: Reverse
reverse_complement = complement[::-1]

print(f"DNA:                5'-{dna}-3'")
print(f"Complement:         3'-{complement}-5'")
print(f"Reverse complement: 5'-{reverse_complement}-3'")

# In one line:
rev_comp = dna.translate(complement_table)[::-1]
print(f"\nOne-liner result:   5'-{rev_comp}-3'")
```
