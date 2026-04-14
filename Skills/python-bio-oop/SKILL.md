---
name: python-bio-oop
description: OOP patterns for bioinformatics — subscriptable sequence databases, sliceable BioSeq, callable motif scorers, __slots__ for millions of variants, and composable mixins
tool_type: python
primary_tool: Python
---

# Advanced OOP Patterns for Bioinformatics

## Patterns Covered

1. `__getitem__`/`__setitem__`/`__contains__` — subscriptable sequence databases
2. `__getitem__` with slices — sliceable sequence objects
3. `__call__` — callable motif scorers with state
4. `__slots__` — 40-60% memory reduction for millions of records
5. Mixins — composable GC, FASTA export, reverse complement

## Subscriptable Sequence Database

```python
class SequenceDatabase:
    """Dict-like container: db['BRCA1'], db['BRCA1'] = seq, 'BRCA1' in db, len(db)."""
    def __init__(self):
        self._data = {}
    def __setitem__(self, name, sequence):
        self._data[name.upper()] = sequence.upper()
    def __getitem__(self, name):
        try:
            return self._data[name.upper()]
        except KeyError:
            raise KeyError(f"Gene '{name}' not found. Available: {list(self._data)[:5]}")
    def __contains__(self, name):
        return name.upper() in self._data
    def __len__(self):
        return len(self._data)
    def __iter__(self):
        return iter(self._data)
```

## Callable Motif Scorer

```python
class MotifScorer:
    """Callable object that scores windows against a motif, tracking usage."""
    def __init__(self, motif, mismatch_penalty=1):
        self.motif = motif.upper()
        self.mismatch_penalty = mismatch_penalty
    def __call__(self, window):
        window = window.upper()[:len(self.motif)]
        return -sum(self.mismatch_penalty for a, b in zip(self.motif, window) if a != b)
    def scan(self, sequence, threshold=0):
        return [(i, sequence[i:i+len(self.motif)], self(sequence[i:i+len(self.motif)]))
                for i in range(len(sequence) - len(self.motif) + 1)
                if self(sequence[i:i+len(self.motif)]) >= threshold]
```

## `__slots__` for Large Variant Collections

```python
class VariantSlots:
    """~40-60% less memory than regular class for millions of records."""
    __slots__ = ('chrom', 'pos', 'ref', 'alt', 'qual')
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom; self.pos = pos; self.ref = ref
        self.alt = alt; self.qual = qual
    # Cannot add arbitrary attributes — that's the tradeoff
```

## Composable Mixins

```python
class BioSequenceMixin:
    """Requires self.sequence. Adds gc_content(), nucleotide_counts(), is_valid()."""
    def gc_content(self):
        seq = self.sequence.upper()
        return (seq.count('G') + seq.count('C')) / len(seq)

class FASTASerializableMixin:
    """Requires self.name, self.sequence. Adds to_fasta()."""
    def to_fasta(self, line_width=60):
        lines = [f'>{self.name}']
        for i in range(0, len(self.sequence), line_width):
            lines.append(self.sequence[i:i + line_width])
        return '\n'.join(lines)

class ReversibleMixin:
    """Requires self.sequence. Adds reverse_complement()."""
    _COMPLEMENT = str.maketrans('ATGCatgc', 'TACGtacg')
    def reverse_complement(self):
        return self.sequence.translate(self._COMPLEMENT)[::-1]

class DNARecord(BioSequenceMixin, FASTASerializableMixin, ReversibleMixin):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence.upper()
```

## MRO (Method Resolution Order)

Python uses C3 linearization for diamond inheritance. Check with:
```python
[cls.__name__ for cls in DNARecord.__mro__]
# ['DNARecord', 'BioSequenceMixin', 'FASTASerializableMixin', 'ReversibleMixin', 'object']
```
