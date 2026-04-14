---
name: python-bio-iterators
description: Custom iterators for bioinformatics — codon iteration, k-mer generation, streaming FASTA/FASTQ, and lazy sequence processing patterns
tool_type: python
primary_tool: Python
---

# Iterators for Bioinformatics

## Key Concepts

- **Class-based iterator vs generator function:** both do the same job, but generators are far less code. Use a class only when you need extra state (resettable position, `peek()` method).
- **Single-pass exhaustion:** once an iterator is exhausted, `list(it)` returns `[]`. Re-create if needed.
- **`itertools` returns lazy iterators.** Wrap in `list()` to materialize, or consume in a for loop.

## Class-Based Codon Iterator

Use when you need resettable state or custom methods beyond what a generator provides.

```python
class CodonIterator:
    """Iterate over a DNA sequence in codons (triplets)."""
    def __init__(self, sequence):
        self.sequence = sequence
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index + 3 > len(self.sequence):
            raise StopIteration
        codon = self.sequence[self.index:self.index + 3]
        self.index += 3
        return codon
```

## Generator Equivalents (Preferred)

```python
def codon_generator(sequence):
    for i in range(0, len(sequence) - 2, 3):
        yield sequence[i:i+3]

def kmer_generator(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]
```

## Translation Pipeline (lazy)

```python
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate_generator(sequence):
    """Lazily translate DNA to amino acids, stopping at stop codon."""
    for codon in codon_generator(sequence):
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            return
        yield aa

# Pipeline: codons -> amino acids -> charged residues (no intermediate lists)
codons = codon_generator(dna)
amino_acids = (CODON_TABLE.get(c, 'X') for c in codons)
charged = (aa for aa in amino_acids if aa in 'DEKRH')
charged_count = sum(1 for _ in charged)
```

## Memory: List vs Generator

```python
import sys
dna = "A" * 1_000_000
kmers_list = [dna[i:i+3] for i in range(len(dna) - 2)]  # ~8 MB
kmers_gen  = kmer_generator(dna, 3)                        # ~200 bytes
```
