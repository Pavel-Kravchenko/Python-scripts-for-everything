---
name: python-bio-oop
description: "Split from `01_classes_and_oop.ipynb` for depth. Start with [Classes](./01_classes.ipynb) first."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/13_Classes_and_OOP/02_oop.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 13: Advanced OOP Patterns

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/13_Classes_and_OOP/02_oop.ipynb`*

# Module 13: Advanced OOP Patterns

Split from `01_classes_and_oop.ipynb` for depth. Start with [Classes](./01_classes.ipynb) first.

**Navigation:** [Previous: Classes](./01_classes.ipynb) · [Topic overview](./01_classes_and_oop.ipynb) · [Next: Module 14](../14_Decorators_and_Context_Managers/01_decorators_and_context_managers.ipynb)

## Why advanced OOP patterns matter in bioinformatics

Once you can write a basic class, the next challenge is designing classes that compose well, scale to real datasets, and integrate cleanly with Python's built-in protocols. A `Seq` object that supports slicing (`seq[3:10]`), a `VariantSet` that can be indexed (`vs['rs123']`), or a `GeneDatabase` that acts like a dict — all require the patterns in this notebook.

## How to use this notebook

This notebook assumes you have worked through `01_classes.ipynb`. Run cells top to bottom. The design patterns here build incrementally — each section introduces a concept then shows a complete bioinformatics example.

## Advanced OOP topics covered

**1. `__getitem__`, `__setitem__`, `__contains__`** — make objects subscriptable
**2. `__call__`** — callable objects (function-like classes)
**3. `__slots__`** — memory-efficient objects for large collections
**4. Method Resolution Order (MRO)** — how Python resolves diamond inheritance
**5. Mixins** — composable behaviors without deep inheritance
**6. Duck typing and protocols** — Pythonic interfaces
**7. `__enter__` / `__exit__`** — context manager protocol on custom classes
**8. Complete design pattern: a FASTA database object**

---

## 1. `__getitem__`, `__setitem__`, `__contains__`: Subscriptable Objects

Implement `__getitem__` to allow `obj[key]` syntax, `__setitem__` for `obj[key] = val`, and `__contains__` for `key in obj`.

```python
class SequenceDatabase:
    """A dict-like container for biological sequences.
    
    Supports:
      - db['BRCA1']           -- retrieve by gene name
      - db['BRCA1'] = seq     -- store a sequence
      - 'BRCA1' in db         -- membership test
      - len(db)               -- number of sequences
      - for name in db:       -- iteration over names
    """

    def __init__(self):
        self._data = {}

    def __setitem__(self, name, sequence):
        self._data[name.upper()] = sequence.upper()

    def __getitem__(self, name):
        try:
            return self._data[name.upper()]
        except KeyError:
            raise KeyError(f"Gene '{name}' not in database. "
                           f"Available: {list(self._data.keys())[:5]}")

    def __contains__(self, name):
        return name.upper() in self._data

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return f"SequenceDatabase({len(self)} entries)"

    def gc_filter(self, min_gc=0.5):
        """Return a new SequenceDatabase with only high-GC sequences."""
        result = SequenceDatabase()
        for name, seq in self._data.items():
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            if gc >= min_gc:
                result[name] = seq
        return result


# Build a sequence database
db = SequenceDatabase()
db['BRCA1'] = 'ATGCGATCGATCGCGATCGATCGATCGATCGATCGATCG'
db['TP53']  = 'GCGCGCGCGCGCATGATCGATCGATCGATCGATCG'
db['EGFR']  = 'ATATATATATATATATATATATCGATCGATCGATCGATCG'
db['MYC']   = 'GCGCGCGCGCGCGCGCGCGCGCGCGCATGATCGATCG'

print(db)
print(f"Has BRCA1: {'BRCA1' in db}")
print(f"Has KRAS:  {'KRAS' in db}")
print(f"BRCA1: {db['BRCA1'][:20]}...")

# Iterate
print("\nAll genes:")
for name in db:
    print(f"  {name}: {len(db[name])} bp")

# Filter
high_gc = db.gc_filter(min_gc=0.55)
print(f"\nHigh-GC sequences (>=55%): {list(high_gc)}")
```

---

## 2. Sliceable Objects: `__getitem__` with Slices

When `__getitem__` receives a `slice` object, you can make your sequence type support slicing just like strings and lists.

```python
class BioSeq:
    """A sliceable biological sequence.
    
    Supports:
      - seq[3]      -- single base by index
      - seq[3:10]   -- slice
      - seq[::3]    -- every third base (codon starts)
      - len(seq), str(seq), repr(seq)
    """

    def __init__(self, sequence, name="unnamed"):
        self.sequence = sequence.upper()
        self.name = name

    def __getitem__(self, key):
        result = self.sequence[key]
        if isinstance(key, slice):
            return BioSeq(result, name=f"{self.name}[slice]")
        return result  # single character

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"BioSeq({self.name!r}, {len(self)} bp)"

    def __add__(self, other):
        """Concatenate two BioSeq objects."""
        new_name = f"{self.name}+{other.name}"
        return BioSeq(str(self) + str(other), name=new_name)

    def gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self)

    def codons(self):
        """Yield codons from reading frame 0."""
        for i in range(0, len(self) - 2, 3):
            yield self.sequence[i:i+3]


# Test BioSeq
seq = BioSeq("ATGGCCGATCGATCGTAGCGA", name="test_gene")
print(repr(seq))
print(f"First base:  {seq[0]}")
print(f"Last base:   {seq[-1]}")
print(f"First codon: {seq[0:3]}")
print(f"GC:          {seq.gc_content():.1%}")

# Slicing returns another BioSeq
fragment = seq[3:12]
print(f"\nFragment: {repr(fragment)}")
print(f"Fragment sequence: {fragment}")

# Concatenation
seq2 = BioSeq("GCTAGCTAGC", name="exon2")
joined = seq + seq2
print(f"\nJoined: {repr(joined)}")
print(f"Codons: {list(joined.codons())}")
```

---

## 3. `__call__`: Callable Objects

Adding `__call__` to a class makes its instances work like functions. This is useful for objects that need to maintain state between calls (like a scorer or a filter).

```python
class MotifScorer:
    """A callable that scores a sequence window for a specific motif.
    
    Usage:
        scorer = MotifScorer('TATAAA', mismatch_penalty=2)
        score = scorer('TATAAAGCGT')  # called like a function
    """

    def __init__(self, motif, mismatch_penalty=1):
        self.motif = motif.upper()
        self.mismatch_penalty = mismatch_penalty
        self._calls = 0  # track usage

    def __call__(self, window):
        """Score a window: 0 = perfect match, negative = mismatches."""
        self._calls += 1
        window = window.upper()[:len(self.motif)]
        if len(window) < len(self.motif):
            return -len(self.motif) * self.mismatch_penalty
        score = 0
        for a, b in zip(self.motif, window):
            if a != b:
                score -= self.mismatch_penalty
        return score

    def scan(self, sequence, threshold=0):
        """Find all positions where the score meets the threshold."""
        sequence = sequence.upper()
        hits = []
        for i in range(len(sequence) - len(self.motif) + 1):
            window = sequence[i:i + len(self.motif)]
            score = self(window)
            if score >= threshold:
                hits.append((i, window, score))
        return hits

    def __repr__(self):
        return f"MotifScorer(motif={self.motif!r}, calls={self._calls})"


# TATA box scanner
tata_scorer = MotifScorer('TATAAA', mismatch_penalty=2)
dna = "GCGATCGTATAATGCGGTATAAAGCGATCGATATAAGCG"

hits = tata_scorer.scan(dna, threshold=-2)  # allow 1 mismatch
print(f"TATA-box candidates in {dna}:")
for pos, window, score in hits:
    print(f"  pos {pos}: {window}  score={score}")

print(f"\n{repr(tata_scorer)}")

# Can be used as a key function or passed to map/filter
sequences = ['TATAAAGCG', 'GCATCGATCG', 'ATATAAAATG', 'CGCGATCG']
scored = [(seq, tata_scorer(seq)) for seq in sequences]
best = max(scored, key=lambda x: x[1])
print(f"\nBest TATA match: {best[0]} (score={best[1]})")
```

---

## 4. `__slots__`: Memory-Efficient Objects

By default, each Python object stores its attributes in a `__dict__`. For large collections (millions of SNPs, millions of k-mers), the per-object `__dict__` overhead is significant. Declaring `__slots__` eliminates the `__dict__` and reduces memory by ~40-60%.

```python
import sys

class VariantDict:
    """A SNP/variant record without __slots__ (normal class)."""
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual


class VariantSlots:
    """A SNP/variant record with __slots__ (memory-optimized)."""
    __slots__ = ('chrom', 'pos', 'ref', 'alt', 'qual')

    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual


# Compare memory usage
n = 10_000
normal_variants = [VariantDict('chr17', i, 'A', 'G', 40.0) for i in range(n)]
slotted_variants = [VariantSlots('chr17', i, 'A', 'G', 40.0) for i in range(n)]

normal_size  = sum(sys.getsizeof(v) + sys.getsizeof(v.__dict__) for v in normal_variants)
slotted_size = sum(sys.getsizeof(v) for v in slotted_variants)

print(f"Normal (with __dict__): {normal_size / 1024:.1f} KB for {n:,} objects")
print(f"Slotted (__slots__):    {slotted_size / 1024:.1f} KB for {n:,} objects")
print(f"Memory reduction: {(1 - slotted_size/normal_size):.0%}")

# __slots__ objects work the same way
v = VariantSlots('chr17', 43045629, 'G', 'A', 99.5)
print(f"\nVariant: {v.chrom}:{v.pos} {v.ref}>{v.alt} QUAL={v.qual}")

# But you cannot add arbitrary attributes
try:
    v.annotation = "pathogenic"
except AttributeError as e:
    print(f"Cannot add new attributes: {e}")
```

---

## 5. Method Resolution Order (MRO) and Mixins

Python uses the **C3 linearization** algorithm to determine the order in which base classes are searched when a method is called. This matters for multiple inheritance.

**Mixins** are small, single-purpose classes intended to be mixed into other classes. They add a specific capability without defining a complete object on their own.

```python
from abc import ABC, abstractmethod


class BioSequenceMixin:
    """Mixin providing common sequence utilities.
    Requires self.sequence to be defined by the host class.
    """

    def gc_content(self):
        seq = self.sequence.upper()
        return (seq.count('G') + seq.count('C')) / len(seq)

    def nucleotide_counts(self):
        seq = self.sequence.upper()
        return {nt: seq.count(nt) for nt in 'ATGC'}

    def is_valid(self, valid_chars):
        return set(self.sequence.upper()) <= set(valid_chars.upper())


class FASTASerializableMixin:
    """Mixin that adds FASTA export capability.
    Requires self.name and self.sequence.
    """

    def to_fasta(self, line_width=60):
        lines = [f'>{self.name}']
        seq = self.sequence
        for i in range(0, len(seq), line_width):
            lines.append(seq[i:i + line_width])
        return '\n'.join(lines)


class ReversibleMixin:
    """Mixin for sequences that have a reverse complement."""
    _COMPLEMENT = str.maketrans('ATGCatgc', 'TACGtacg')

    def reverse_complement(self):
        return self.sequence.translate(self._COMPLEMENT)[::-1]


class DNARecord(BioSequenceMixin, FASTASerializableMixin, ReversibleMixin):
    """A concrete DNA sequence combining multiple mixins."""

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence.upper()

    def __repr__(self):
        return f"DNARecord({self.name!r}, {len(self.sequence)} bp)"


# MRO: order in which Python searches for methods
print("MRO:", [cls.__name__ for cls in DNARecord.__mro__])

# Use the class
gene = DNARecord("BRCA1_exon11", "ATGCGATCGATCGCGATCGATCGATCGATCGATCGATCGATCG")
print(f"\n{repr(gene)}")
print(f"GC content: {gene.gc_content():.1%}")
print(f"Valid DNA:  {gene.is_valid('ATGCN')}")
print(f"Rev comp:   {gene.reverse_complement()[:20]}...")
print(f"\nFASTA format:\n{gene.to_fasta(line_width=20)}")
```

---

## 6. Duck Typing and Protocols

Python does not require formal interface declarations. Any object that provides the right methods is accepted — this is **duck typing**: "if it walks like a duck and quacks like a duck, it's a duck."

For more formal type checking without ABCs, Python 3.8+ provides `typing.Protocol`.
