---
name: python-bio-control-flow
description: Python control flow patterns specific to bioinformatics — codon iteration, frame scanning, motif search, sequence classification.
tool_type: python
primary_tool: Python
---

# Control Flow for Bioinformatics

## Bio-Specific Patterns

### Codon iteration (range with step 3)
```python
# Correct: stops before incomplete last codon
for i in range(0, len(dna) - 2, 3):
    codon = dna[i:i+3]
    if len(codon) == 3:
        process(codon)

# Wrong: range(len(dna)) produces incomplete final codon
```

### Position-aware iteration
```python
for i, nuc in enumerate(dna):
    # i is 0-based; report as i+1 for 1-based coordinates
    if nuc == 'G':
        positions.append(i)

# Aligned sequence comparison
for i, (n1, n2) in enumerate(zip(seq1, seq2)):
    if n1 != n2:
        print(f"Mismatch at position {i+1}: {n1}->{n2}")
```

### Motif scanning with while
```python
positions = []
pos = seq.find(motif)
while pos != -1:
    positions.append(pos)
    pos = seq.find(motif, pos + 1)
```

### Stop codon search
```python
stop_codons = {"TAA", "TAG", "TGA"}
pos = 0
while pos <= len(dna) - 3:
    codon = dna[pos:pos+3]
    if codon in stop_codons:
        break
    pos += 3
```

### for-else for "not found" pattern
```python
# else clause runs only when loop completes without break
for i, nuc in enumerate(sequence):
    if nuc not in valid_bases:
        print(f"Invalid '{nuc}' at {i}")
        break
else:
    print("Sequence is valid")
```

### Sequence type detection
```python
def detect_sequence_type(sequence: str) -> str:
    unique = set(sequence.upper())
    if unique <= set("ATGC"):       return "DNA"
    elif unique <= set("AUGC"):     return "RNA"
    elif unique <= set("ACDEFGHIKLMNPQRSTVWY"): return "Protein"
    return "Unknown"
```

### GC content classification
```python
GC_CLASSES = [(30, "AT-rich"), (50, "Moderate"), (60, "High GC")]

def classify_gc(sequence: str) -> tuple[float, str]:
    s = sequence.upper()
    gc = (s.count('G') + s.count('C')) / len(s) * 100
    for threshold, label in GC_CLASSES:
        if gc < threshold:
            return gc, label
    return gc, "Very high GC"
```

## Pitfalls

- **Off-by-one in codon loops**: `range(len(seq) - 2)` stops 2 before end, leaving room for a 3-char slice. `range(len(seq))` produces an incomplete last codon.
- **`elif` vs separate `if`**: use `elif` for mutually exclusive classification (codon type); use separate `if` for independent filters.
- **`break` exits only the innermost loop**: in nested frame-scanning loops, `break` does not exit the outer loop.
- **`while` without progress**: every `while` loop must modify its condition variable or use `break`; a stuck kernel means infinite loop.
- **`pass` does nothing**: a syntactic placeholder only — don't use it where you mean `continue`.
- **Mutable default arguments**: `def f(x=[])` shares the list across all calls; use `def f(x=None)` and initialize inside.
- **1-based vs 0-based reporting**: Python iteration is 0-based; bioinformatics positions are conventionally 1-based — always add 1 when reporting to users.
