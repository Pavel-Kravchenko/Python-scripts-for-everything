---
name: python-bio-classes
description: Python classes for bioinformatics — designing Sequence, Gene, and ProteinRecord classes with encapsulation, dunder methods, inheritance, ABCs, and properties.
tool_type: python
primary_tool: Python
---

# Classes for Bioinformatics

## BioSequence Hierarchy

```python
class BioSequence:
    def __init__(self, sequence: str, name: str = "unnamed"):
        self.sequence = sequence.upper()
        self.name = name

    def __len__(self):          return len(self.sequence)
    def __str__(self):          return f">{self.name}\n{self.sequence}"
    def __repr__(self):         return f"{type(self).__name__}('{self.sequence}', name='{self.name}')"
    def __contains__(self, motif): return motif.upper() in self.sequence
    def composition(self):      return {c: self.sequence.count(c) for c in sorted(set(self.sequence))}


class DNA(BioSequence):
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def gc_content(self) -> float:
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100

    def reverse_complement(self) -> "DNA":
        comp = str.maketrans('ATGCN', 'TACGN')
        return DNA(self.sequence.translate(comp)[::-1], name=f"{self.name}_revcomp")

    def transcribe(self) -> "RNA":
        return RNA(self.sequence.replace('T', 'U'), name=f"{self.name}_rna")


class RNA(BioSequence):
    def to_dna(self) -> DNA:
        return DNA(self.sequence.replace('U', 'T'), name=f"{self.name}_dna")


class Protein(BioSequence):
    AA_WEIGHTS = {
        'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146,
        'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165,
        'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117,
    }
    def molecular_weight(self) -> float:
        weight = sum(self.AA_WEIGHTS.get(aa, 110) for aa in self.sequence)
        return weight - (len(self.sequence) - 1) * 18  # subtract water per peptide bond
```

## Dunder Methods Reference

| Method | Enables |
|--------|---------|
| `__init__` | `Gene("BRCA1", "ATG...")` |
| `__str__` | `print(gene)` — human readable |
| `__repr__` | `repr(gene)` — developer view, should allow recreation |
| `__len__` | `len(seq)` |
| `__eq__` | `seq1 == seq2` |
| `__lt__` | `seq1 < seq2`, `sorted(seqs)` |
| `__contains__` | `"ATG" in seq` |
| `__enter__`/`__exit__` | `with obj as x:` |

**Note**: defining `__eq__` sets `__hash__ = None` in Python 3. Objects become unhashable (cannot be used as dict keys/set members). Define `__hash__` explicitly if you need both.

## Properties for Validation

```python
class Gene:
    def __init__(self, name: str, sequence: str, strand: str = '+'):
        self.name = name
        self.sequence = sequence   # calls setter
        self.strand = strand       # calls setter

    @property
    def sequence(self) -> str:
        return self._sequence

    @sequence.setter
    def sequence(self, value: str):
        value = value.upper()
        invalid = set(value) - set('ATGCN')
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        self._sequence = value

    @property
    def gc_content(self) -> float:   # read-only computed property
        return (self._sequence.count('G') + self._sequence.count('C')) / len(self._sequence) * 100
```

## Abstract Base Classes

```python
from abc import ABC, abstractmethod

class SequenceAnalyzer(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()

    @abstractmethod
    def validate(self) -> bool: ...

    @abstractmethod
    def summary(self) -> dict: ...


class DNAAnalyzer(SequenceAnalyzer):
    def validate(self) -> bool:
        invalid = set(self.sequence) - set('ATGCN')
        if invalid:
            raise ValueError(f"Invalid DNA bases: {invalid}")
        return True

    def summary(self) -> dict:
        gc = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100
        return {'length': len(self.sequence), 'gc_content': round(gc, 2)}
```

## classmethod and staticmethod

```python
@classmethod
def from_fasta_string(cls, fasta_text: str) -> "BioSequence":
    lines = fasta_text.strip().split('\n')
    seq_id = lines[0][1:].split()[0]
    return cls(''.join(lines[1:]), name=seq_id)

@staticmethod
def is_valid_dna(seq: str) -> bool:
    return set(seq.upper()).issubset(set('ATGCN'))
```

## Pitfalls

- **`self` is the instance, not the class**: `self.sequence` reads the instance attribute; `DNA.sequence` would be a class-level variable. Never confuse the two.
- **Mutable class attributes**: `class Gene: tags = []` — appending to `tags` affects every instance. Set mutable attributes in `__init__` with `self.tags = []`.
- **`super().__init__()` in subclasses**: forgetting it means the parent's `__init__` never runs and parent attributes are never set.
- **Properties without `_` backing store**: `self.sequence = value` inside the setter calls the setter again → infinite recursion. Use `self._sequence`.
- **`__eq__` disables `__hash__`**: after defining `__eq__`, the class becomes unhashable. Add `__hash__ = None` explicitly or implement `__hash__`.
- **Abstract class instantiation**: `SequenceAnalyzer("ATGC")` raises `TypeError` at runtime — good, that is the intent of ABCs.
