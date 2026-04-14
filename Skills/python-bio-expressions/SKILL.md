---
name: python-bio-expressions
description: Python arithmetic, comparison, logical, and membership operators with bioinformatics applications — codon frames, GC content, quality filters.
tool_type: python
primary_tool: Python
---

## Arithmetic Operators

| Operator | Name | Bioinformatics use |
|----------|------|--------------------|
| `+` | Addition | Concatenate exons, count bases |
| `-` | Subtraction | Gene length = stop − start + 1 |
| `*` | Multiplication | Protein MW estimate |
| `/` | True division | GC% = (G+C)/total * 100 |
| `//` | Floor division | Complete codons = length // 3 |
| `%` | Modulo | Reading frame position = pos % 3 |
| `**` | Exponentiation | Probability calculations |

## Pitfalls

- **`/` vs `//`**: `/` always returns float; `//` truncates to int. Codon counts need `//`: `seq_length // 3` gives complete codons. `seq_length / 3` gives a float.
- **`%` is modulo, not percent**: `seq_length % 3` = leftover nucleotides after codon division.
- **Operator precedence**: `g + c / total * 100` computes as `g + ((c / total) * 100)` — wrong. Correct: `(g + c) / total * 100`. Add parentheses when combining arithmetic and division.
- **`==` vs `is`**: use `==` to compare values. Use `is` only for identity checks (`None`, `True`, `False`). Never `if x == None`.
- **`and`/`or` short-circuit**: `A and B` skips `B` when `A` is falsy. Efficient, but can mask bugs if `B` has side effects.
- **`in` checks dict keys, not values**: `"ATG" in codon_table` checks keys. To check values: `"ATG" in codon_table.values()` (slower) or invert the dict.
- **`in` on list vs set**: membership check on a list is O(n); on a set is O(1). For large stop-codon or restriction-site lookups, use a `frozenset`.

## Bio Applications

### Codons and Reading Frames
```python
seq_length = 1000
complete_codons = seq_length // 3          # 333
leftover_nt     = seq_length % 3           # 1

# Reading frame of a position (0-based)
frame = position % 3                       # 0, 1, or 2

# Which codon index contains position pos?
codon_index = pos // 3
```

### GC Content
```python
def gc_content(seq: str) -> float:
    seq = seq.upper()
    total = len(seq)
    if total == 0:
        return 0.0
    return (seq.count('G') + seq.count('C')) / total * 100

# GC categories
def classify_gc(gc_pct: float) -> str:
    if gc_pct < 30:   return "Very AT-rich (e.g. Plasmodium)"
    if gc_pct < 40:   return "AT-rich"
    if gc_pct < 60:   return "Moderate"
    if gc_pct < 70:   return "GC-rich"
    return "Very GC-rich (e.g. Streptomyces)"
```

### Stop Codon Check
```python
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

is_stop = codon in STOP_CODONS          # O(1) set lookup
# Avoid: codon in ["TAA", "TAG", "TGA"]  # O(n) list lookup
```

### QC Filter (chained comparisons)
```python
passes_qc = (
    500 <= seq_length <= 3000
    and 0.30 <= gc_content <= 0.70
    and not has_ambiguous_bases
)
```

### Protein MW Estimate
```python
# Rough estimate: ~110 Da per AA, minus 18 Da per peptide bond
def protein_mw_da(n_aa: int) -> float:
    return n_aa * 110 - (n_aa - 1) * 18
```

### Validate DNA Sequence
```python
VALID_DNA = frozenset('ACGT')

def is_valid_dna(seq: str) -> bool:
    return all(c in VALID_DNA for c in seq.upper())
    # Note: does not allow IUPAC ambiguity codes (N, R, Y, ...)
```

## Operator Precedence (high → low, bioinformatics-relevant)

1. `**`
2. Unary `-`, `+`
3. `*`, `/`, `//`, `%`
4. `+`, `-`
5. Comparison: `<`, `<=`, `>`, `>=`, `==`, `!=`
6. `not`
7. `and`
8. `or`

When in doubt, use parentheses — they cost nothing and prevent silent errors.
