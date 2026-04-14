---
name: python-bio-functions
description: Bio-specific function patterns — ORF finding, sequence analysis, mutable default trap, and functional programming idioms for bioinformatics
tool_type: python
primary_tool: Python
---

# Functions for Bioinformatics

## Pitfalls

- **Mutable default arguments:** `def f(x, data=[])` is a classic bug. The default list is created once and shared across calls. Always use `def f(x, data=None): if data is None: data = []`.
- **`return` vs `print()`:** `result = my_function()` returns `None` if the function uses `print()` instead of `return`.
- **`*args`/`**kwargs` order:** required, default, `*args`, keyword-only, `**kwargs`. Getting this wrong causes `SyntaxError`.

## Bio Recipes

### ORF Finder (all 3 reading frames)

```python
def find_orfs(sequence: str, min_length: int = 30) -> list[dict]:
    """Find all ORFs (ATG to stop codon) in forward reading frames."""
    sequence = sequence.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == 'ATG':
                for j in range(i + 3, len(sequence) - 2, 3):
                    if sequence[j:j+3] in stop_codons:
                        orf_seq = sequence[i:j+3]
                        if len(orf_seq) >= min_length:
                            orfs.append({'start': i, 'end': j+3,
                                         'length': len(orf_seq), 'sequence': orf_seq})
                        i = j + 3
                        break
                else:
                    i += 3
                    continue
                continue
            i += 3
    return sorted(orfs, key=lambda x: x['start'])
```

### FASTA Header Builder with **kwargs

```python
def create_fasta_header(sequence_id, **metadata):
    """Build FASTA header from ID and arbitrary key=value metadata."""
    header = f">{sequence_id}"
    for key, value in metadata.items():
        header += f" [{key}={value}]"
    return header

# create_fasta_header("BRCA1", organism="Homo sapiens", gene="BRCA1", length=5500)
# -> ">BRCA1 [organism=Homo sapiens] [gene=BRCA1] [length=5500]"
```

### Sequence Validation with Type Hints

```python
def validate_dna_sequence(sequence: str, allow_n: bool = False) -> bool:
    """Check whether a string is valid DNA."""
    valid_chars = set('ATGCN') if allow_n else set('ATGC')
    return set(sequence.upper()) <= valid_chars
```
