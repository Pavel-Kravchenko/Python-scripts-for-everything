---
name: python-bio-classes
description: "Split from `01_classes_and_oop.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/13_Classes_and_OOP/01_classes.ipynb"
---

# Classes

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/13_Classes_and_OOP/01_classes.ipynb`*

# Classes

Split from `01_classes_and_oop.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_classes_and_oop.ipynb) · [Next: Oop](./02_oop.ipynb)

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. `self` is the instance, not the class**
Every regular method receives the instance as `self`. `self.sequence` reads the instance's attribute. `Gene.sequence` would be a class-level variable. Never confuse the two.

**2. Class attributes vs instance attributes**
Class attributes are shared across all instances. Instance attributes (set with `self.x = ...` in `__init__`) are per-object. A mutable class attribute (e.g., `class Gene: sequences = []`) that you append to will affect *every* instance.

**3. `super().__init__()` in inheritance**
When a subclass defines `__init__`, you must call `super().__init__(...)` to trigger the parent's initialization. Forgetting this means the parent's attributes are never set.

**4. `__eq__` disables `__hash__`**
Defining `__eq__` automatically sets `__hash__ = None` in Python 3. Objects become unhashable (cannot be used as dict keys or in sets). Define `__hash__` explicitly if you need both.

**5. Properties are computed on access**
Accessing `gene.gc` where `gc` is a `@property` calls a function. Properties are ideal for computed values (GC%, sequence length) or validated setters.

## 2. Defining a Class: `__init__`, Attributes, Methods

```python
class Gene:
    """Represents a gene with its DNA sequence and metadata."""

    def __init__(self, name, sequence, chromosome=None):
        """
        Initialize a Gene object.

        Parameters:
            name: Gene symbol (e.g. "BRCA1")
            sequence: DNA sequence string
            chromosome: Chromosome number (optional)
        """
        self.name = name
        self.sequence = sequence.upper()
        self.chromosome = chromosome

    def length(self):
        """Return sequence length in base pairs."""
        return len(self.sequence)

    def gc_content(self):
        """Calculate GC content as a percentage."""
        gc = self.sequence.count('G') + self.sequence.count('C')
        return (gc / len(self.sequence)) * 100


# Create Gene objects (instances)
brca1 = Gene("BRCA1", "ATGGATTTCGATCGATCGTAGC", chromosome=17)
tp53 = Gene("TP53", "ATGGAGGAGCCGCAGTCAGATC", chromosome=17)

print(f"Gene: {brca1.name}")
print(f"Length: {brca1.length()} bp")
print(f"GC Content: {brca1.gc_content():.1f}%")
print(f"Chromosome: {brca1.chromosome}")
```

### How `__init__` Works

When you write `brca1 = Gene("BRCA1", "ATG...")`, Python:

1. Creates a new, empty `Gene` object
2. Calls `__init__(new_object, "BRCA1", "ATG...")` -- `self` refers to the new object
3. Returns the initialized object

`self` is not a keyword -- it is a convention. The first parameter of every instance method
is the object itself.

## 3. `__str__` and `__repr__`: Controlling How Objects Print

```python
class Sequence:
    """A biological sequence with proper string representations."""

    def __init__(self, seq_id, sequence, seq_type="DNA"):
        self.seq_id = seq_id
        self.sequence = sequence.upper()
        self.seq_type = seq_type

    def __str__(self):
        """Human-readable: called by print() and str()."""
        if len(self.sequence) > 20:
            display = f"{self.sequence[:10]}...{self.sequence[-10:]}"
        else:
            display = self.sequence
        return f"{self.seq_id} ({self.seq_type}, {len(self.sequence)} bp): {display}"

    def __repr__(self):
        """Developer-facing: unambiguous, could recreate the object."""
        return f"Sequence('{self.seq_id}', '{self.sequence}', seq_type='{self.seq_type}')"


seq = Sequence("BRCA1_exon1", "ATGGATTTCGATCGATCGTAGCGATCGATCGATCG")

# __str__ is called by print()
print(seq)

# __repr__ is called by repr() and shown in interactive prompts
print(repr(seq))
```

**Rule of thumb:**
- `__str__` -- for the user (readable output)
- `__repr__` -- for the developer (should ideally let you recreate the object)
- If only one is defined, use `__repr__` -- it serves as fallback for `__str__` too

## 4. More Dunder Methods: `__len__`, `__eq__`, `__lt__`

```python
class DNASequence:
    """DNA sequence with rich comparison support."""

    VALID_BASES = set('ATGCN')

    def __init__(self, sequence, name="unnamed"):
        self.name = name
        sequence = sequence.upper()
        invalid = set(sequence) - self.VALID_BASES
        if invalid:
            raise ValueError(f"Invalid nucleotides: {invalid}")
        self.sequence = sequence

    def __len__(self):
        """Allow len(dna_seq) to work."""
        return len(self.sequence)

    def __eq__(self, other):
        """Two DNASequence objects are equal if their sequences match."""
        if not isinstance(other, DNASequence):
            return NotImplemented
        return self.sequence == other.sequence

    def __lt__(self, other):
        """Compare by length (useful for sorting)."""
        if not isinstance(other, DNASequence):
            return NotImplemented
        return len(self.sequence) < len(other.sequence)

    def __repr__(self):
        return f"DNASequence('{self.sequence[:20]}...', name='{self.name}')"


# Equality
s1 = DNASequence("ATGCATGC", "seq1")
s2 = DNASequence("ATGCATGC", "seq2")  # different name, same sequence
s3 = DNASequence("ATGCATGCATGC", "seq3")

print(f"s1 == s2: {s1 == s2}")   # True -- same sequence
print(f"s1 == s3: {s1 == s3}")   # False
print(f"s1 < s3:  {s1 < s3}")    # True -- shorter

# Sorting works because __lt__ is defined
sequences = [s3, s1, DNASequence("AT", "tiny")]
for s in sorted(sequences):
    print(f"  {s.name}: {len(s)} bp")
```

## 5. Inheritance: Building a Sequence Hierarchy

Inheritance lets you create specialized classes from a general base class.
DNA, RNA, and protein sequences share common features (they all have a sequence
string, a name, a length) but differ in their specific operations.

```
               BioSequence (Parent)
              /        |         \
         DNA          RNA       Protein
      gc_content()  to_dna()   mol_weight()
      complement()             
      transcribe()
```

```python
class BioSequence:
    """Base class for all biological sequences."""

    def __init__(self, sequence, name="unnamed"):
        self.sequence = sequence.upper()
        self.name = name

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return f">{self.name}\n{self.sequence}"

    def __contains__(self, motif):
        """Allow 'ATG' in seq syntax."""
        return motif.upper() in self.sequence

    def composition(self):
        """Return character frequency dict."""
        return {char: self.sequence.count(char) for char in sorted(set(self.sequence))}


class DNA(BioSequence):
    """DNA sequence with DNA-specific methods."""

    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def gc_content(self):
        gc = self.sequence.count('G') + self.sequence.count('C')
        return (gc / len(self.sequence)) * 100

    def complement(self):
        """Return the complementary DNA strand."""
        comp = ''.join(self.COMPLEMENT_MAP[nt] for nt in self.sequence)
        return DNA(comp, name=f"{self.name}_comp")

    def reverse_complement(self):
        """Return the reverse complement."""
        comp = ''.join(self.COMPLEMENT_MAP[nt] for nt in reversed(self.sequence))
        return DNA(comp, name=f"{self.name}_revcomp")

    def transcribe(self):
        """Transcribe DNA to RNA (T -> U)."""
        return RNA(self.sequence.replace('T', 'U'), name=f"{self.name}_rna")


class RNA(BioSequence):
    """RNA sequence."""

    def to_dna(self):
        return DNA(self.sequence.replace('U', 'T'), name=f"{self.name}_dna")


class Protein(BioSequence):
    """Protein sequence with mass calculation."""

    # Average residue weights in Daltons
    AA_WEIGHTS = {
        'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121,
        'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131,
        'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115,
        'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117
    }

    def molecular_weight(self):
        """Calculate molecular weight in Daltons (subtract water per peptide bond)."""
        weight = sum(self.AA_WEIGHTS.get(aa, 110) for aa in self.sequence)
        weight -= (len(self.sequence) - 1) * 18  # water loss
        return weight
```

```python
# Test the hierarchy
dna = DNA("ATGCGATCGATCGTAGCGATCG", name="test_gene")

print(f"DNA: {dna.sequence}")
print(f"Length: {len(dna)} bp")
print(f"GC Content: {dna.gc_content():.1f}%")
print(f"Composition: {dna.composition()}")
print(f"Contains ATG: {'ATG' in dna}")

rev_comp = dna.reverse_complement()
print(f"\nReverse complement: {rev_comp.sequence}")

rna = dna.transcribe()
print(f"RNA: {rna.sequence}")

# isinstance checks
print(f"\ndna is DNA: {isinstance(dna, DNA)}")
print(f"dna is BioSequence: {isinstance(dna, BioSequence)}")
print(f"rna is DNA: {isinstance(rna, DNA)}")
```

## 6. Abstract Base Classes

Sometimes you want to define an interface that subclasses **must** implement.
Python's `abc` module enforces this at instantiation time.

```python
from abc import ABC, abstractmethod


class SequenceAnalyzer(ABC):
    """Abstract base class for sequence analyzers."""

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    @abstractmethod
    def validate(self):
        """Check that the sequence is valid. Must be implemented."""
        ...

    @abstractmethod
    def summary(self):
        """Return a summary dict. Must be implemented."""
        ...


# This will raise TypeError -- cannot instantiate abstract class
try:
    analyzer = SequenceAnalyzer("ATGC")
except TypeError as e:
    print(f"Cannot instantiate abstract class: {e}")


class DNAAnalyzer(SequenceAnalyzer):
    """Concrete implementation for DNA."""

    def validate(self):
        invalid = set(self.sequence) - set('ATGCN')
        if invalid:
            raise ValueError(f"Invalid DNA bases: {invalid}")
        return True

    def summary(self):
        gc = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100
        return {
            'length': len(self.sequence),
            'gc_content': round(gc, 2),
            'composition': {b: self.sequence.count(b) for b in 'ATGC'}
        }


analyzer = DNAAnalyzer("ATGCGATCGATCG")
analyzer.validate()
print(analyzer.summary())
```

## 7. Properties: Controlled Attribute Access

Properties let you add validation or computation behind attribute access,
without changing how the caller uses the object.

```python
class Gene:
    """Gene with validated attributes via properties."""

    VALID_STRANDS = {'+', '-'}

    def __init__(self, name, sequence, strand='+'):
        self.name = name
        self.sequence = sequence  # goes through the setter
        self.strand = strand      # goes through the setter

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        value = value.upper()
        invalid = set(value) - set('ATGCN')
        if invalid:
            raise ValueError(f"Invalid nucleotides in sequence: {invalid}")
        self._sequence = value

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value):
        if value not in self.VALID_STRANDS:
            raise ValueError(f"Strand must be '+' or '-', got '{value}'")
        self._strand = value

    @property
    def gc_content(self):
        """Read-only computed property."""
        gc = self._sequence.count('G') + self._sequence.count('C')
        return (gc / len(self._sequence)) * 100


gene = Gene("TP53", "ATGGAGGAGCCGCAGTCAGATC", strand='+')
print(f"{gene.name}: GC = {gene.gc_content:.1f}%")

# Validation in action
try:
    gene.strand = 'X'
except ValueError as e:
    print(f"Caught: {e}")

try:
    gene.sequence = "ATGXYZ"
except ValueError as e:
    print(f"Caught: {e}")
```
