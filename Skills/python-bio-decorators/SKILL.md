---
name: python-bio-decorators
description: "Split from `01_decorators_and_context_managers.ipynb` to keep this topic self-contained."
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/14_Decorators_and_Context_Managers/01_decorators.ipynb"
---

# Decorators

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/14_Decorators_and_Context_Managers/01_decorators.ipynb`*

# Decorators

Split from `01_decorators_and_context_managers.ipynb` to keep this topic self-contained.

**Navigation:** [Topic overview](./01_decorators_and_context_managers.ipynb) · [Next: Context Managers](./02_context_managers.ipynb)

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. A decorator is just a callable that takes a callable**
`@timer` above `def analyze(seq):` is exactly `analyze = timer(analyze)`. It replaces the function with its decorated version at definition time.

**2. Always use `functools.wraps`**
Without `@functools.wraps(func)`, the wrapper has `__name__ = 'wrapper'` and no docstring. This breaks `help()`, logging, and any tool that inspects function metadata.

**3. Decorator factories need three levels**
When the decorator takes arguments (`@validate_sequence('ATGC')`), add one extra level:
```python
def validate_sequence(valid_chars):   # factory — takes the argument
    def decorator(func):              # actual decorator — takes the function
        @functools.wraps(func)
        def wrapper(*args, **kwargs): # called at runtime
            ...
        return wrapper
    return decorator
```

**4. Stacking decorators applies bottom-up**
`@A @B def f():` is `f = A(B(f))`. The decorator closest to the function (`B`) is applied first.

**5. `lru_cache` requires hashable arguments**
`@lru_cache` caches by argument value. All arguments must be hashable — lists and dicts cannot be cached directly.

```python
def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def at_content(seq):
    """Calculate AT content of a DNA sequence."""
    seq = seq.upper()
    return (seq.count('A') + seq.count('T')) / len(seq) * 100

# Functions are objects -- assign to variable, put in list
analyzers = [gc_content, at_content]

test_seq = "ATGCGATCGATCGTAGC"
for func in analyzers:
    print(f"{func.__name__}: {func(test_seq):.1f}%")
```

```python
# Higher-order function: takes a function as argument
def apply_to_sequences(sequences, analysis_func):
    """Apply an analysis function to multiple sequences."""
    return {name: analysis_func(seq) for name, seq in sequences.items()}

seqs = {
    "BRCA1": "ATGGATTTCGATCGATCGTAGC",
    "TP53":  "ATGGAGGAGCCGCAGTCAGATC",
    "MYC":   "GGCCAATTGGCCAATTGGCC",
}

gc_results = apply_to_sequences(seqs, gc_content)
for gene, gc in gc_results.items():
    print(f"{gene}: {gc:.1f}% GC")
```

## 2. Closures: Functions That Remember

A closure is an inner function that captures variables from its enclosing scope.
This lets you create specialized functions at runtime.

```python
def make_motif_counter(motif):
    """Return a function that counts a specific motif in any sequence."""
    motif = motif.upper()

    def counter(sequence):
        return sequence.upper().count(motif)

    return counter

# Create specialized counters
count_cpg = make_motif_counter("CG")
count_tata = make_motif_counter("TATA")

promoter = "GCTATAAAAGGCGCGCGTATAATCGCG"
print(f"CG dinucleotides: {count_cpg(promoter)}")
print(f"TATA boxes: {count_tata(promoter)}")

# The inner function 'remembers' the motif even after make_motif_counter returns
print(f"\nClosure variables: {count_cpg.__closure__[0].cell_contents}")
```

## 3. Decorators: The Basics

A decorator is a function that takes a function and returns a modified version of it.

```
DECORATOR PATTERN
+-------------------------------------------+
|  @decorator        is equivalent to        |
|  def func():  -->  func = decorator(func)  |
+-------------------------------------------+
```

```python
import functools


def log_call(func):
    """Decorator that logs when a function is called."""
    @functools.wraps(func)  # preserves __name__, __doc__, etc.
    def wrapper(*args, **kwargs):
        print(f"Calling {func.__name__}...")
        result = func(*args, **kwargs)
        print(f"{func.__name__} returned {result}")
        return result
    return wrapper


@log_call
def count_nucleotides(seq):
    """Count each nucleotide in a DNA sequence."""
    seq = seq.upper()
    return {base: seq.count(base) for base in 'ATGC'}


result = count_nucleotides("ATGCGATCGATCG")
print(f"\nResult: {result}")

# functools.wraps preserves the original function's identity
print(f"Function name: {count_nucleotides.__name__}")
print(f"Docstring: {count_nucleotides.__doc__}")
```

### Why `functools.wraps` Matters

Without `@functools.wraps(func)`, the decorated function would have:
- `__name__` = `"wrapper"` instead of the original name
- `__doc__` = the wrapper's docstring instead of the original

This breaks introspection, help(), and debugging. Always use `functools.wraps`.

## 4. Timing Decorator: Benchmarking Bioinformatics Algorithms

```python
import time


def timer(func):
    """Measure and print execution time of a function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"[timer] {func.__name__}: {elapsed:.4f} seconds")
        return result
    return wrapper
```

```python
import random


@timer
def gc_content_loop(seq):
    """GC content via explicit loop."""
    gc_count = 0
    for base in seq:
        if base in 'GCgc':
            gc_count += 1
    return gc_count / len(seq) * 100


@timer
def gc_content_count(seq):
    """GC content via str.count()."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100


# Generate a 1 million bp test sequence
big_seq = ''.join(random.choices('ATGC', k=1_000_000))

print("Comparing GC content methods on 1M bp:")
r1 = gc_content_loop(big_seq)
r2 = gc_content_count(big_seq)
print(f"Both give: {r1:.2f}% vs {r2:.2f}%")
```

## 5. Sequence Validation Decorator

```python
def validate_sequence(valid_chars, seq_type="DNA"):
    """Decorator factory: validate that the first argument is a valid sequence."""
    valid_set = set(valid_chars.upper())

    def decorator(func):
        @functools.wraps(func)
        def wrapper(seq, *args, **kwargs):
            invalid = set(seq.upper()) - valid_set
            if invalid:
                raise ValueError(
                    f"Invalid {seq_type} characters {invalid} in input to {func.__name__}()"
                )
            return func(seq, *args, **kwargs)
        return wrapper
    return decorator


@validate_sequence('ATGC', seq_type='DNA')
def complement(seq):
    """Return the DNA complement."""
    comp_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp_map[b] for b in seq.upper())


@validate_sequence('ACDEFGHIKLMNPQRSTVWY', seq_type='protein')
def protein_mass(seq):
    """Estimate protein mass (rough, 110 Da per residue)."""
    return len(seq) * 110


# Valid input
print(complement("ATGCGATC"))

# Invalid DNA
try:
    complement("ATGXYZ")
except ValueError as e:
    print(f"Caught: {e}")

# Invalid protein
try:
    protein_mass("MVLS123")
except ValueError as e:
    print(f"Caught: {e}")
```

## 6. Decorators with Arguments

When a decorator needs parameters (like `@validate_sequence('ATGC')`),
you need a three-level nesting: a factory that returns the actual decorator.

```python
def repeat(n=3):
    """Decorator factory: run the function n times and return all results."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return [func(*args, **kwargs) for _ in range(n)]
        return wrapper
    return decorator


@repeat(n=5)
def random_gc_sample(seq, sample_size=100):
    """GC content of a random subsequence."""
    start = random.randint(0, len(seq) - sample_size)
    sample = seq[start:start + sample_size].upper()
    gc = (sample.count('G') + sample.count('C')) / len(sample) * 100
    return round(gc, 1)


gc_samples = random_gc_sample(big_seq, sample_size=500)
print(f"GC content from 5 random 500-bp windows: {gc_samples}")
print(f"Mean: {sum(gc_samples) / len(gc_samples):.1f}%")
```

## 7. Memoization: Caching Expensive Results

```python
def memoize(func):
    """Cache function results for repeated calls with the same arguments."""
    cache = {}

    @functools.wraps(func)
    def wrapper(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]

    wrapper.cache = cache
    wrapper.clear_cache = cache.clear
    return wrapper


# Recursive Needleman-Wunsch alignment score (simplified)
@memoize
def align_score(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """Recursive pairwise alignment score with memoization."""
    if not seq1:
        return len(seq2) * gap
    if not seq2:
        return len(seq1) * gap

    score = match if seq1[-1] == seq2[-1] else mismatch

    return max(
        align_score(seq1[:-1], seq2[:-1]) + score,
        align_score(seq1[:-1], seq2) + gap,
        align_score(seq1, seq2[:-1]) + gap,
    )


s1 = "ACGTACGT"
s2 = "ACGTAGCT"
print(f"Alignment score for '{s1}' vs '{s2}': {align_score(s1, s2)}")
print(f"Cache entries: {len(align_score.cache)}")
```

```python
# Python's built-in lru_cache is more efficient and has a size limit
from functools import lru_cache


@lru_cache(maxsize=128)
def translate_codon(codon):
    """Translate a single codon to amino acid (cached)."""
    table = {
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
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    return table.get(codon.upper(), 'X')


@timer
def translate_sequence(dna):
    """Translate a DNA sequence using the cached codon lookup."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aa = translate_codon(dna[i:i+3])
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


# Translate a long repeated sequence to see caching in action
long_dna = "ATGGCTGCTTAG" * 1000
protein = translate_sequence(long_dna)
print(f"Protein (first 30 aa): {protein[:30]}...")
print(f"Cache info: {translate_codon.cache_info()}")
```

## 8. Stacking Multiple Decorators

Decorators are applied bottom-up (closest to the function first).

```python
@timer
@validate_sequence('ATGC', seq_type='DNA')
def analyze_sequence(seq):
    """Full analysis of a DNA sequence."""
    seq = seq.upper()
    return {
        'length': len(seq),
        'gc_content': (seq.count('G') + seq.count('C')) / len(seq) * 100,
        'starts_with_atg': seq.startswith('ATG'),
        'a': seq.count('A'),
        't': seq.count('T'),
        'g': seq.count('G'),
        'c': seq.count('C'),
    }


# Valid sequence -- gets validated then timed
result = analyze_sequence("ATGCGATCGATCGATCGATCG")
print(f"Analysis: {result}")

# Invalid -- validation catches it before timing even starts
try:
    analyze_sequence("ATGXYZ")
except ValueError as e:
    print(f"\nCaught: {e}")
```

---

## 9. Context Managers: Safe Resource Handling

Context managers ensure resources (files, connections, locks) are properly
acquired and released, even when exceptions occur.

```
CONTEXT MANAGER LIFECYCLE
+-----------------------------------+
|  with ctx as resource:            |
|      |                            |
|      v                            |
|  __enter__() -> acquire resource  |
|      |                            |
|      v                            |
|  (your code runs here)            |
|      |                            |
|      v                            |
|  __exit__() -> release resource   |
+-----------------------------------+
```
