---
name: python-bio-generators
description: Generator patterns for streaming bioinformatics data — FASTA/FASTQ readers, sliding windows, translation pipelines, and memory-efficient sequence processing
tool_type: python
primary_tool: Python
---

# Generators for Bioinformatics

## Key Concepts

- A streaming FASTQ reader that `yield`s one record at a time uses constant memory regardless of file size.
- Generator pipelines (filter -> trim -> translate) produce results on demand, never storing full datasets in RAM.
- **Generators are single-pass.** Once exhausted, they are empty. Convert to list if you need multiple iterations.

## Generator Pipelines: No Intermediate Lists

```python
codons  = codon_generator(dna)                 # lazy
aas     = (CODON_TABLE[c] for c in codons)     # lazy
charged = (aa for aa in aas if aa in 'RKHDE')  # lazy
list(charged)  # computation happens only here
```

## Sliding Window Generator

```python
def sliding_window(sequence, window_size, step=1):
    for i in range(0, len(sequence) - window_size + 1, step):
        yield i, sequence[i:i + window_size]
```

## Streaming FASTQ Reader

```python
def read_fastq(filename):
    """Yield one record at a time. Memory: O(single record)."""
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
                    'id':       header[1:].split()[0],
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

# Pipeline: quality filter -> N filter (nothing computed until consumed)
passing = list(filter_no_n(filter_quality(read_fastq('reads.fastq'))))
```

## Streaming FASTA Reader

```python
def read_fasta(filename):
    """Yield (header, sequence) tuples. Memory: O(single record)."""
    header, parts = None, []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)
```

## `yield from` for Multi-K K-mer Generation

```python
def all_kmers_multiK(sequence, k_values):
    for k in k_values:
        yield from (sequence[i:i+k] for i in range(len(sequence) - k + 1))
```

## Useful `itertools` Patterns

```python
import itertools

# takewhile: trim quality scores until they drop
high_quality = list(itertools.takewhile(lambda q: q >= 30, quality_scores))

# groupby: detect homopolymer runs
for base, group in itertools.groupby(dna):
    run_length = sum(1 for _ in group)

# product: generate all possible k-mers of length k (4^k)
all_kmers = (''.join(c) for c in itertools.product('ATGC', repeat=k))
```
