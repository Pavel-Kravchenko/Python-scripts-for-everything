---
name: python-bio-iterators
description: "Split from `01_iterators_and_generators.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/11_Iterators_and_Generators/01_iterators.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Iterators

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/11_Iterators_and_Generators/01_iterators.ipynb`*

# Iterators

Split from `01_iterators_and_generators.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_iterators_and_generators.ipynb) · [Next: Generators](./02_generators.ipynb)

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. Iterables vs iterators**
A list is *iterable* but not an *iterator*. `iter(my_list)` returns an iterator. Most objects you loop over are iterables that produce a fresh iterator each time; an iterator is always also an iterable (returns `self` from `__iter__`), but is single-pass.

**2. Single-pass exhaustion**
Once you call `list(my_iterator)` or drain it with a for loop, the iterator is empty. Calling `list(my_iterator)` again returns `[]`. This is intentional — it keeps memory constant.

**3. Class-based iterator vs generator function**
Both do the same job, but a generator function is far less code. Use a class-based iterator only when you need additional state (e.g., a position counter that can be reset, or a `peek()` method).

**4. `itertools` returns lazy iterators**
Every `itertools` function returns a lazy iterator. Wrap in `list()` to see the values, or consume in a for loop.

```python
# Any iterable can be turned into an iterator with iter()
nucleotides = ["A", "T", "G", "C"]
nuc_iter = iter(nucleotides)

print(type(nucleotides))   # list -- iterable
print(type(nuc_iter))      # list_iterator -- iterator
```

```python
# Manually advance the iterator with next()
print(next(nuc_iter))  # A
print(next(nuc_iter))  # T
print(next(nuc_iter))  # G
print(next(nuc_iter))  # C
```

```python
# When an iterator is exhausted, StopIteration is raised
try:
    print(next(nuc_iter))
except StopIteration:
    print("Iterator exhausted! No more nucleotides.")
```

```python
# Strings are iterable too -- iterate over a DNA sequence
dna = "ATGCGA"
seq_iter = iter(dna)

print(next(seq_iter))  # A
print(next(seq_iter))  # T
print(next(seq_iter))  # G
```

### How `for` Loops Really Work

A `for` loop is syntactic sugar around `iter()` and `next()`:

```python
# This for loop:
for nuc in "ATGC":
    print(nuc, end=" ")
print()

# Is equivalent to this:
iterator = iter("ATGC")
while True:
    try:
        nuc = next(iterator)
        print(nuc, end=" ")
    except StopIteration:
        break
print()
```

### Building a Custom Iterator Class

Let's build a **CodonIterator** that walks through a DNA sequence in triplets:

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


dna = "ATGAAACCCGGGTTTAAA"
print(f"Sequence: {dna}")
print("Codons:")
for codon in CodonIterator(dna):
    print(f"  {codon}")
```

```python
# Iterators are single-pass -- once exhausted, they are done
ci = CodonIterator("ATGAAA")

print("First pass:", list(ci))
print("Second pass:", list(ci))  # empty -- already exhausted
```

---
## 2. Generators with `yield`

Writing a full iterator class is verbose. **Generators** are functions that use `yield` instead of `return`. Each call to `next()` resumes execution right after the last `yield`.

Generators automatically implement the iterator protocol.

```python
def codon_generator(sequence):
    """Yield codons (triplets) from a DNA sequence."""
    for i in range(0, len(sequence) - 2, 3):
        yield sequence[i:i+3]


gen = codon_generator("ATGAAACCCGGGTTT")
print(f"Type: {type(gen)}")
print(f"First codon: {next(gen)}")
print(f"Second codon: {next(gen)}")
print(f"Remaining: {list(gen)}")
```

```python
# Generator for k-mers (subsequences of length k)
def kmer_generator(sequence, k):
    """Yield all k-mers from a sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]


dna = "ATGCGATCG"
print(f"All 3-mers from {dna}:")
print(list(kmer_generator(dna, 3)))
```

```python
# Generators can be infinite -- Fibonacci example for modeling populations
def fibonacci():
    """Infinite Fibonacci generator (e.g., rabbit population growth model)."""
    a, b = 0, 1
    while True:
        yield a
        a, b = b, a + b


fib = fibonacci()
print("First 12 Fibonacci numbers:")
for _ in range(12):
    print(next(fib), end=" ")
print()
```

```python
# Translation generator -- yield amino acids one at a time
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def translate_generator(sequence):
    """Lazily translate a DNA sequence into amino acids, stopping at stop codon."""
    for codon in codon_generator(sequence):
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            return  # StopIteration
        yield aa


coding_seq = "ATGAAAGCCTTTGGGTGA"
protein = ''.join(translate_generator(coding_seq))
print(f"DNA:     {coding_seq}")
print(f"Protein: {protein}")
```

---
## 3. Generator Expressions

Generator expressions look like list comprehensions but use `()` instead of `[]`. They produce values lazily.

```python
dna = "ATGCGATCGATCGATCG"

# List comprehension -- all values computed and stored in memory
gc_list = [1 for nuc in dna if nuc in 'GC']

# Generator expression -- values computed on demand
gc_gen = (1 for nuc in dna if nuc in 'GC')

print(f"List: {gc_list}")
print(f"Generator: {gc_gen}")
print(f"GC count via generator: {sum(1 for nuc in dna if nuc in 'GC')}")
```

```python
# GC content with a generator expression
def gc_content(seq):
    """Calculate GC content using a generator expression."""
    return sum(1 for nuc in seq if nuc in 'GC') / len(seq)


sequences = {
    'AT-rich promoter': 'AAATTTATTTAAAGCTA',
    'CpG island':       'GCGCGGCGCGCGCCGCG',
    'Random':           'ATGCGATCGATCGATCG'
}

for name, seq in sequences.items():
    print(f"{name:20s}: GC = {gc_content(seq):.1%}")
```

---
## 4. Memory Efficiency

Generators are critical when working with data that does not fit in memory.

```python
import sys

# Compare memory usage: list vs generator
dna = "A" * 1_000_000  # 1 million bases

# All 3-mers as a list
kmers_list = [dna[i:i+3] for i in range(len(dna) - 2)]

# All 3-mers as a generator
kmers_gen = kmer_generator(dna, 3)

print(f"List of {len(kmers_list):,} k-mers: {sys.getsizeof(kmers_list):,} bytes")
print(f"Generator object:              {sys.getsizeof(kmers_gen)} bytes")
print(f"\nThe generator uses ~{sys.getsizeof(kmers_list) / sys.getsizeof(kmers_gen):,.0f}x less memory!")
```

```python
# Chaining generators for pipeline processing (no intermediate lists)
dna = "ATGAAAGCCTTTGGGTGATCGATCG" * 100

# Pipeline: codons -> amino acids -> only charged residues
codons = codon_generator(dna)
amino_acids = (CODON_TABLE.get(c, 'X') for c in codons)
charged = (aa for aa in amino_acids if aa in 'DEKRH')

# Nothing has been computed yet! Only when we consume:
charged_count = sum(1 for _ in charged)
print(f"Charged residues in {len(dna)} bp sequence: {charged_count}")
```

---
## 5. Bioinformatics Application: Streaming FASTQ Reader

FASTQ files can be tens of gigabytes. A generator reads one record at a time.

```python
# First, create a sample FASTQ file
fastq_content = """@SEQ_ID_001 length=30
ATGCGATCGATCGATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID_002 length=30
GCTAGCTAGCTAGCTAGCTAGCTAGCTANN
+
IIIIIIIIIIIIIIIIIIIIIIIIIII!!!
@SEQ_ID_003 length=30
TTAACCGGTTAACCGGTTAACCGGTTAACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID_004 length=30
GGCCAATTGGCCAATTGGCCAATTGGCCAA
+
IIIIIIIIIIIIIIIIIII!!!!!!!!!!!
"""

with open('sample.fastq', 'w') as f:
    f.write(fastq_content)

print("Created sample.fastq")
```

```python
def read_fastq(filename):
    """Streaming FASTQ reader -- yields one record at a time.
    
    Each record is a dict with keys: 'id', 'sequence', 'quality'.
    Memory usage is constant regardless of file size.
    """
    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()  # '+' line
            quality = f.readline().strip()
            
            yield {
                'id': header[1:].split()[0],
                'sequence': sequence,
                'quality': quality
            }


# Process without loading entire file
print("Streaming FASTQ records:")
for record in read_fastq('sample.fastq'):
    avg_qual = sum(ord(c) - 33 for c in record['quality']) / len(record['quality'])
    gc = gc_content(record['sequence'])
    print(f"  {record['id']}: GC={gc:.1%}, mean_qual={avg_qual:.1f}")
```

```python
# Generator pipeline: filter low-quality reads
def quality_filter(records, min_avg_quality=30):
    """Yield only records with average quality >= threshold."""
    for record in records:
        avg_qual = sum(ord(c) - 33 for c in record['quality']) / len(record['quality'])
        if avg_qual >= min_avg_quality:
            yield record


def trim_n_bases(records):
    """Trim trailing N bases from sequences."""
    for record in records:
        seq = record['sequence'].rstrip('N')
        qual = record['quality'][:len(seq)]
        yield {**record, 'sequence': seq, 'quality': qual}


# Chain generators into a processing pipeline
raw_reads = read_fastq('sample.fastq')
trimmed = trim_n_bases(raw_reads)
filtered = quality_filter(trimmed, min_avg_quality=25)

print("Filtered and trimmed reads:")
for record in filtered:
    print(f"  {record['id']}: {record['sequence'][:40]}... (len={len(record['sequence'])})")
```

---
## 6. Bioinformatics Application: Streaming FASTA Reader

```python
# Create a sample FASTA file
fasta_content = """>gene1 BRCA1 tumor suppressor
ATGAAAGCCTTTGGGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGTAA
>gene2 TP53 guardian of the genome
ATGGCCCCCGGGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCTGA
>gene3 EGFR growth factor receptor
ATGCCCAAATTTGGGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGTAG
"""

with open('sample.fasta', 'w') as f:
    f.write(fasta_content)

print("Created sample.fasta")
```

```python
def read_fasta(filename):
    """Memory-efficient FASTA parser. Yields (header, sequence) tuples."""
    header = None
    seq_parts = []
    
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield (header, ''.join(seq_parts))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    
    if header is not None:
        yield (header, ''.join(seq_parts))


print("Parsed FASTA records:")
for header, seq in read_fasta('sample.fasta'):
    gene_id = header.split()[0]
    print(f"  {gene_id}: {len(seq)} bp, GC={gc_content(seq):.1%}")
```

---
## 7. Bioinformatics Application: All Possible k-mers

Generate all possible DNA k-mers of a given length (4^k total).

```python
import itertools


def all_possible_kmers(k):
    """Generator for all possible DNA k-mers of length k."""
    for combo in itertools.product('ATGC', repeat=k):
        yield ''.join(combo)


# All possible dimers (4^2 = 16)
print("All 16 possible 2-mers:")
print(list(all_possible_kmers(2)))

# How many 5-mers exist? (4^5 = 1024)
count = sum(1 for _ in all_possible_kmers(5))
print(f"\nTotal possible 5-mers: {count}")
```
