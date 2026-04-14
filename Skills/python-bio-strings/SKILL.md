---
name: python-bio-strings
description: "Python string operations for bioinformatics: reverse complement, GC content, motif finding, FASTA parsing. Bio-specific patterns and pitfalls only."
tool_type: python
primary_tool: Python
---

# Python Strings for Bioinformatics

## Core Patterns

```python
# Reverse complement — use maketrans, NOT chained replace()
complement_table = str.maketrans('ATGC', 'TACG')
rev_comp = dna.translate(complement_table)[::-1]

# Transcription
rna = dna.replace('T', 'U')

# Extract all codons
codons = [mrna[i:i+3] for i in range(0, len(mrna), 3)]

# GC content
gc = (seq.count('G') + seq.count('C')) / len(seq)

# Validate DNA
is_valid = all(c in 'ATGC' for c in seq.upper())
```

## Motif Finding

```python
# All positions of a motif (overlapping)
def find_all(seq, motif):
    positions = []
    pos = seq.find(motif)
    while pos != -1:
        positions.append(pos)
        pos = seq.find(motif, pos + 1)
    return positions

# Non-overlapping count (built-in)
"ATATATAT".count("ATAT")  # returns 2

# Regex for IUPAC patterns
import re
hits = [(m.start(), m.group()) for m in re.finditer(r'TATA[AT]A[AT]', seq)]
```

## FASTA Header Parsing

```python
# UniProt: >sp|P04637|P53_HUMAN Cellular tumor antigen...
header = ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens"
parts = header.split('|')
db, accession = parts[0][1:], parts[1]
description = parts[2].split()[0]  # P53_HUMAN
```

## Sequence Features

```python
seq = "ATGGCCGATTAGCCA"

seq.startswith('ATG')          # has start codon
seq[-3:] in ('TAA', 'TAG', 'TGA')  # ends with stop codon
'GATTAG' in seq                # motif membership
```

## Pitfalls

- **Chained `replace()` fails for complement**: `dna.replace('A','T').replace('T','A')` converts all A→T then all T→A, including the ones just created. Use `str.maketrans` + `translate()`.
- **`find()` returns -1, not None**: `if pos:` is truthy for -1. Always check `if pos != -1:`.
- **`count()` is non-overlapping**: `"ATAT".count("AT")` = 2 (correct), but `"AAAA".count("AA")` = 2, not 3. Use `find()` loop for overlapping.
- **Case sensitivity**: `.count('G')` misses lowercase 'g'. Always normalize with `.upper()` first.
- **Strings are immutable**: `seq[0] = 'C'` raises `TypeError`. Build new strings with concatenation or `''.join(...)`.
- **Bioinformatics coordinates are 1-based**: Python slices are 0-based half-open `[start, stop)`. GFF/VCF position 100 → Python index 99.
