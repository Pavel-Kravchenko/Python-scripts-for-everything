---
name: python-bio-sequences
description: "Python string operations for biological sequences: reverse complement, codon extraction, motif scanning, FASTA parsing. Bio-specific patterns and pitfalls."
tool_type: python
primary_tool: Python
---

# Biological Sequences as Python Strings

## Core Operations

```python
# Reverse complement — the only correct approach
complement_table = str.maketrans('ATGC', 'TACG')
rev_comp = dna.translate(complement_table)[::-1]

# Transcription
rna = dna.replace('T', 'U')

# All codons
codons = [mrna[i:i+3] for i in range(0, len(mrna), 3)]

# GC content
gc = (seq.count('G') + seq.count('C')) / len(seq)
```

## Motif Scanning

```python
# All positions (overlapping) — built-in find() loop
def find_all(seq, motif):
    pos, positions = seq.find(motif), []
    while pos != -1:
        positions.append(pos)
        pos = seq.find(motif, pos + 1)
    return positions

# IUPAC degenerate patterns
import re
re.findall(r'TATA[AT]A[AT]', seq)         # TATAWAW
re.finditer(r'GG[ACGT]{1,3}GG', seq)      # G-quadruplex-like
```

## FASTA Parsing

```python
# Minimal FASTA reader
def read_fasta(filepath):
    records = {}
    with open(filepath) as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    records[header] = ''.join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
        if header:
            records[header] = ''.join(seq)
    return records

# UniProt header parsing
# ">sp|P04637|P53_HUMAN ..." → db='sp', accession='P04637'
db, accession, rest = header.split('|', 2)
```

## Sequence Validation

```python
valid_dna  = set('ATGC')
valid_iupac = set('ACGTRYMKSWHBVDN')
is_valid = all(c in valid_dna for c in seq.upper())
```

## Coordinate Systems

| Format | Base | Interval type | Example: first 3 bp |
|--------|------|--------------|---------------------|
| Python | 0 | Half-open [start, stop) | `seq[0:3]` |
| BED | 0 | Half-open | `start=0, end=3` |
| VCF/GFF | 1 | Closed [start, stop] | `start=1, end=3` |
| SAM | 1 | Closed | `POS=1` |

Converting GFF→Python: `python_start = gff_start - 1`

## Pitfalls

- **Chained `replace()` for complement is wrong**: `dna.replace('A','T').replace('T','A')` corrupts — use `str.maketrans`
- **`find()` returns -1 on miss**: `if pos:` is truthy for -1; always check `if pos != -1:`
- **`count()` is non-overlapping**: `"AAAA".count("AA")` = 2, not 3
- **Case sensitivity**: always `.upper()` before scanning; lowercase = soft-masked repeats in UCSC/Ensembl
- **Strings are immutable**: `seq[0] = 'C'` raises `TypeError`; build new string with concatenation
- **Off-by-one from 1-based coordinates**: subtract 1 when converting database positions to Python indices
