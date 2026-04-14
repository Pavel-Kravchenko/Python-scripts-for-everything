---
name: python-bio-operators
description: Operator pitfalls and idioms for bioinformatics — floor division for codons, modulo for reading frames, GC content precedence, and membership testing
tool_type: python
primary_tool: Python
---

# Operators for Bioinformatics

## Pitfalls

- **`//` vs `/` for codons:** `seq_length // 3` gives complete codons (int). `/` always returns float.
- **`%` for reading frames:** `position % 3` gives the reading frame (0, 1, or 2). `seq_length % 3` gives leftover nucleotides.
- **Precedence trap in GC content:** `g + c / total * 100` is wrong. Use `(g + c) / total * 100`.
- **`==` vs `is`:** use `==` for value comparison, `is` only for `None`/`True`/`False`. Never write `if x == None`.
- **`in` checks dict keys, not values:** `"ATG" in codon_table` checks keys. Use `"Met" in codon_table.values()` for values.
- **`in` with sets is O(1), with lists is O(n).** Convert to set for repeated membership tests.

## Bio-Specific Patterns

### Codon Arithmetic

```python
seq_length = 1000
complete_codons = seq_length // 3         # 333
leftover_nucleotides = seq_length % 3     # 1
reading_frame = position % 3              # 0, 1, or 2
```

### Chained Comparisons for Range Checks

```python
# GC content in normal range?
0.40 <= gc_content <= 0.60   # more Pythonic than gc_content >= 0.40 and gc_content <= 0.60
```

### Stop Codon Check (idiomatic `in`)

```python
is_stop = codon in {"TAA", "TAG", "TGA"}   # set literal for O(1) lookup
```

### Multi-Criteria Sequence Filtering

```python
passes_qc = (500 <= seq_length <= 3000
             and 0.30 <= gc_content <= 0.70
             and not has_ambiguous_bases)
```

### Protein Molecular Weight Estimate

```python
# Rough estimate: avg AA mass ~110 Da, water lost per peptide bond = 18 Da
protein_mw = (num_amino_acids * 110) - ((num_amino_acids - 1) * 18)
```
