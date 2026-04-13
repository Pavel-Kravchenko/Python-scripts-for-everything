---
name: python-bio-generators
description: "Split from `01_iterators_and_generators.ipynb` for depth. Start with [Iterators](./01_iterators.ipynb) first."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/11_Iterators_and_Generators/02_generators.ipynb"
---

# Module 11: Generators

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/11_Iterators_and_Generators/02_generators.ipynb`*

# Module 11: Generators

Split from `01_iterators_and_generators.ipynb` for depth. Start with [Iterators](./01_iterators.ipynb) first.

**Navigation:** [Previous: Iterators](./01_iterators.ipynb) · [Topic overview](./01_iterators_and_generators.ipynb) · [Next: Regular Expressions →](../12_Regular_Expressions/01_regular_expressions.ipynb)

## Why generators matter in bioinformatics

A streaming FASTQ reader that `yield`s one record at a time uses the same memory regardless of whether the file has 1,000 or 100,000,000 reads. Generator pipelines — filter → trim → quality-score → translate — produce results on demand, never storing the full dataset in RAM. This is how production bioinformatics tools like `pysam` and streaming FASTA parsers work internally.

## How to use this notebook

This notebook assumes you have worked through `01_iterators.ipynb`. Run all cells top to bottom; later sections depend on functions defined earlier.

## Complicated moments explained

**1. `yield` vs `return`**
`return` ends the function. `yield` pauses it. A function with any `yield` statement is a generator function; calling it returns a generator object without executing any code yet.

**2. Generator expressions vs generator functions**
`(x for x in seq if condition)` is a generator expression — concise, single expression. A `def` function with `yield` supports loops, conditionals, and multiple yields — use it for anything more complex.

**3. Chaining generators: no intermediate lists**
```python
codons  = codon_generator(dna)                 # lazy
aas     = (CODON_TABLE[c] for c in codons)     # lazy
charged = (aa for aa in aas if aa in 'RKHDE')  # lazy
list(charged)  # computation happens only here
```

**4. Generators are single-pass**
Once exhausted, a generator is empty. If you need to iterate the result more than once, convert to a list first.

**5. `send()` for coroutines**
Advanced generators can receive values via `gen.send(value)`, enabling coroutines. Rarely needed in typical bioinformatics scripts.

---

## 1. Generator Functions: `yield`

A generator function uses `yield` instead of `return`. Each call to `next()` resumes execution right after the last `yield`.

```
GENERATOR EXECUTION FLOW
+----------------------------------------------------------+
|  def gen_func():      # calling gen_func() returns obj   |
|      yield 1          # execution pauses here            |
|      yield 2          # resumes here on next next()      |
|      yield 3          # resumes here on next next()      |
|                       # StopIteration raised after 3     |
+----------------------------------------------------------+
```

```python
# The simplest possible generator
def count_up(start, stop):
    i = start
    while i <= stop:
        yield i
        i += 1  # execution resumes HERE on next next()

# Calling the function returns a generator object -- no code runs yet
gen = count_up(1, 5)
print(f"Generator object: {gen}")
print(f"Type: {type(gen)}")

# Advance manually with next()
print(f"\nnext(): {next(gen)}")
print(f"next(): {next(gen)}")

# Drain the rest with a for loop
print("Remaining:", end=" ")
for val in gen:
    print(val, end=" ")
print()

# In bioinformatics: yield positions of a motif
def find_motif_positions(sequence, motif):
    seq_len = len(sequence)
    motif_len = len(motif)
    for i in range(seq_len - motif_len + 1):
        if sequence[i:i + motif_len] == motif:
            yield i

dna = "ATGATGCCCATGATGATG"
print(f"\nATG positions in {dna}:")
for pos in find_motif_positions(dna, "ATG"):
    print(f"  position {pos}")
```

---

## 2. Generator Pipelines: Processing Without Intermediate Lists

Each stage of a generator pipeline produces one value at a time. No stage holds the full dataset in memory simultaneously.

```python
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def codons(sequence):
    for i in range(0, len(sequence) - 2, 3):
        yield sequence[i:i+3]

def translate(codon_gen):
    for codon in codon_gen:
        aa = CODON_TABLE.get(codon.upper(), 'X')
        if aa == '*':
            return
        yield aa

def only_charged(aa_gen):
    charged = set('RKHDE')
    for aa in aa_gen:
        if aa in charged:
            yield aa

dna = "ATGAAGCGCGATGAAATCGATGAAGTGGTTGAAATCGAA" * 3

# Pipeline: each stage is lazy -- no intermediate list is created
pipeline = only_charged(translate(codons(dna)))
charged_residues = list(pipeline)
print(f"DNA: {dna[:40]}...")
print(f"Charged residues: {''.join(charged_residues)}")
print(f"Count: {len(charged_residues)}")
```

---

## 3. Sliding Window Generator

```python
def sliding_window(sequence, window_size, step=1):
    for i in range(0, len(sequence) - window_size + 1, step):
        yield i, sequence[i:i + window_size]

def gc_content(seq):
    s = seq.upper()
    return (s.count('G') + s.count('C')) / len(s) if s else 0.0

dna = "GCGCGCATATATATGCGCGCATATATATGCGCGCGC"
print(f"Sequence: {dna}")
print(f"Sliding window GC content (window=8, step=4):")

for pos, window in sliding_window(dna, window_size=8, step=4):
    gc = gc_content(window)
    bar = '#' * int(gc * 20)
    print(f"  pos {pos:3d}: {window}  GC={gc:.0%}  {bar}")
```

---

## 4. Memory Efficiency: List vs Generator

```python
import sys

n = 1_000_000
list_result = [i**2 for i in range(n)]
gen_result  = (i**2 for i in range(n))

print(f"List:      {sys.getsizeof(list_result):>12,} bytes ({sys.getsizeof(list_result)/1e6:.1f} MB)")
print(f"Generator: {sys.getsizeof(gen_result):>12,} bytes")
print(f"Ratio:     {sys.getsizeof(list_result) / sys.getsizeof(gen_result):.0f}x")

# For bioinformatics: count GC-rich windows in a large sequence without storing them
import random
random.seed(42)
genome = ''.join(random.choices('ATGC', weights=[0.3,0.3,0.2,0.2], k=500_000))

window_size = 500
high_gc_count = sum(
    1 for _, window in sliding_window(genome, window_size, step=100)
    if gc_content(window) > 0.55
)
total_windows = (len(genome) - window_size) // 100 + 1

print(f"\nGenome: {len(genome):,} bp")
print(f"Windows scanned:  {total_windows:,}")
print(f"GC>55% windows:   {high_gc_count:,} ({high_gc_count/total_windows:.1%})")
print("(No intermediate list of windows was ever created)")
```

---

## 5. Streaming FASTQ Reader

FASTQ files can be 50-100 GB. A generator-based reader uses constant memory.

```python
import os, tempfile

# Create a sample FASTQ file
records_data = [
    ("SEQ001", "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
     "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    ("SEQ002", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
     "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    ("SEQ003", "ATGCGATCGATCGATCGNNCGATCGATCGATCGATCGATCGATCGATCGAT",
     "IIIIIIIIIIIIIIIII##IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    ("SEQ004", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
     "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"),
    ("SEQ005", "ATATATATATATATATATATATATATATATATATATATATATATATATATA",
     "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"),
]
fastq_lines = []
for seq_id, seq, qual in records_data:
    fastq_lines.extend([f"@{seq_id}", seq, "+", qual])
fastq_content = "\n".join(fastq_lines)

fastq_path = tempfile.mktemp(suffix='.fastq')
with open(fastq_path, 'w') as f:
    f.write(fastq_content)

def read_fastq(filename):
    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq  = f.readline().strip()
            f.readline()  # '+' separator
            qual = f.readline().strip()
            if header.startswith('@'):
                yield {
                    'id':       header[1:],
                    'sequence': seq,
                    'quality':  qual,
                    'avg_qual': sum(ord(c) - 33 for c in qual) / len(qual) if qual else 0,
                }

def filter_quality(records, min_avg_qual=30):
    for rec in records:
        if rec['avg_qual'] >= min_avg_qual:
            yield rec

def filter_no_n(records):
    for rec in records:
        if 'N' not in rec['sequence'].upper():
            yield rec

# Generator pipeline: quality filter -> N filter
raw     = read_fastq(fastq_path)
q_filt  = filter_quality(raw, min_avg_qual=30)
n_filt  = filter_no_n(q_filt)

print("Reads passing QC (avg_qual>=30, no N bases):")
for rec in n_filt:
    print(f"  {rec['id']}: {len(rec['sequence'])} bp, avg_qual={rec['avg_qual']:.1f}")

os.remove(fastq_path)
```

---

## 6. Streaming FASTA Reader

```python
def read_fasta(text):
    header = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq_parts)
            header = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line)
    if header is not None:
        yield header, ''.join(seq_parts)


fasta_text = (
    ">BRCA1 Breast cancer susceptibility protein 1\n"
    "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA\n"
    ">TP53 Tumor protein p53\n"
    "GCGCGCGCGCGCGCGCGCGCGCGCGCGATCGATCGATCGATCGATCGATCGT\n"
    "ATCGATCGATCGATCGATCGATCGATCG\n"
    ">EGFR Epidermal growth factor receptor\n"
    "ATATATATATATATATATATATATATATATATATATATATATCGATCGATCGAT\n"
    "CGATCGATCGATCGATCGATCGATCG\n"
)

print("FASTA records:")
for header, seq in read_fasta(fasta_text):
    gene = header.split()[0]
    print(f"  {gene}: {len(seq)} bp, GC={gc_content(seq):.1%}")
```

---

## 7. `yield from`: Delegating to Sub-generators

`yield from` delegates to another iterable or generator, forwarding all values.

```python
import itertools

def all_kmers_multiK(sequence, k_values):
    for k in k_values:
        yield from _kmers(sequence, k)  # delegate

def _kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield (k, sequence[i:i+k])

dna = "ATGCGATCG"
print(f"All 2-mers and 3-mers of '{dna}':")
for k, kmer in all_kmers_multiK(dna, [2, 3]):
    print(f"  {k}-mer: {kmer}")
```

---

## 8. `itertools` with Generators

```python
# itertools.takewhile: consume until condition fails
quality_scores = [40, 39, 38, 37, 35, 30, 20, 10, 5, 2]
high_quality = list(itertools.takewhile(lambda q: q >= 30, quality_scores))
print(f"High-quality scores (>=30): {high_quality}")
print(f"Quality drop at position: {len(high_quality)}")

# itertools.groupby: group homopolymer runs
dna_run = "AAATTTGGGGCCAATTGCCCCCC"
print(f"\nHomopolymer runs in '{dna_run}':")
for base, group in itertools.groupby(dna_run):
    run = list(group)
    print(f"  {base} x {len(run)}")

# itertools.chain: concatenate multiple generators
exon1_seq = "ATGAAAGCC"
exon2_seq  = "TTTGGGTGA"
mrna = ''.join(itertools.chain(exon1_seq, exon2_seq))
print(f"\nSpliced mRNA (exon1+exon2): {mrna}")
```

---

## Exercises

### Exercise 1: Codon Counter Generator

Write `codon_counter(fasta_text)` that yields `(gene_id, codon_counts_dict)` tuples without storing all sequences in memory.

```python
# Exercise 1: Your code here
from collections import Counter

def codon_counter(fasta_text):
    for header, seq in read_fasta(fasta_text):
        gene_id = header.split()[0]
        # YOUR CODE: count codons in seq
        pass

# test_fasta = ">gene1\nATGGCCGATCGAGCC\n>gene2\nATGAAACCCGGGTTTTGA\n"
# for gene_id, counts in codon_counter(test_fasta):
#     print(f"{gene_id}: {dict(counts)}")
```
