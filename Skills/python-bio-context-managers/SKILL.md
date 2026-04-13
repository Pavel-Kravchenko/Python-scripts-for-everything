---
name: python-bio-context-managers
description: "Bioinformatics pipelines open many files — FASTAs, BAMs, VCFs, temporary alignments — and may crash mid-pipeline. Without a context manager, file handles leak, temporary files persist, and databases a"
tool_type: python
source_notebook: "Tier_1_Python_for_Bioinformatics/14_Decorators_and_Context_Managers/02_context_managers.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 14: Decorators and Context Managers

*Source: Course notebook `Tier_1_Python_for_Bioinformatics/14_Decorators_and_Context_Managers/02_context_managers.ipynb`*


## Elegant Python Patterns for Bioinformatics Workflows

---

### Learning Objectives
- Understand first-class functions and closures
- Write decorators from scratch, with and without arguments
- Use `functools.wraps` to preserve function metadata
- Build context managers with `__enter__`/`__exit__` and `contextlib`
- Apply these patterns to bioinformatics: timing, validation, and safe file handling

---

## How to use this notebook

Run cells from top to bottom on the first pass — later cells depend on functions and variables defined earlier. Once you have run through the notebook, feel free to modify parameters and re-run individual sections.

Each section has runnable examples first, followed by exercises. Try the exercise before looking at the solution cell below it.

## Complicated moments explained

**1. `__exit__` receives the exception, not `__enter__`**
If the `with` block raises, Python calls `__exit__(exc_type, exc_val, tb)`. Returning `True` from `__exit__` suppresses the exception; returning `None` or `False` re-raises it.

**2. `@contextmanager` yield position**
In a `@contextmanager` function, code before `yield` is setup (`__enter__`), code after `yield` is teardown (`__exit__`). Always wrap teardown in `try/finally`:
```python
@contextmanager
def managed_resource():
    resource = acquire()
    try:
        yield resource
    finally:
        release(resource)  # always runs, even on exception
```python

**3. Context managers for SQLite**
`with conn:` automatically commits on success or rolls back on exception. Essential for multi-step genomic database inserts.

**4. Nesting with a single `with` statement**
`with open(f1) as fasta, open(f2) as vcf:` is idiomatic — equivalent to two nested `with` statements.

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
```python

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
```python

## 3. Decorators: The Basics

A decorator is a function that takes a function and returns a modified version of it.

```python
DECORATOR PATTERN
+-------------------------------------------+
|  @decorator        is equivalent to        |
|  def func():  -->  func = decorator(func)  |
+-------------------------------------------+
```python

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
```python

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
```python

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
```python

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
```python

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
```python

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
```python

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
```python

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
```python

---

## 9. Context Managers: Safe Resource Handling

Context managers ensure resources (files, connections, locks) are properly
acquired and released, even when exceptions occur.

```python
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
```python

```python
class FastaWriter:
    """Context manager for writing FASTA files."""

    def __init__(self, filename, line_width=80):
        self.filename = filename
        self.line_width = line_width
        self.file = None
        self.record_count = 0

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
        if exc_type is not None:
            print(f"Error occurred while writing: {exc_val}")
        print(f"Wrote {self.record_count} records to {self.filename}")
        return False  # do not suppress exceptions

    def write_record(self, seq_id, sequence, description=""):
        """Write a single FASTA record with line wrapping."""
        header = f">{seq_id}" + (f" {description}" if description else "")
        self.file.write(header + "\n")
        for i in range(0, len(sequence), self.line_width):
            self.file.write(sequence[i:i + self.line_width] + "\n")
        self.record_count += 1


# Write sequences safely -- file is always closed
with FastaWriter("output_test.fasta") as writer:
    writer.write_record("BRCA1", "ATGGATTTCGATCG" * 10, "breast cancer gene")
    writer.write_record("TP53", "ATGGAGGAGCCGCAG" * 8, "tumor protein p53")
    writer.write_record("EGFR", "ATGCGACCCTCCGGG" * 12, "epidermal growth factor receptor")

# Verify
with open("output_test.fasta") as f:
    print(f.read()[:300] + "...")
```python

## Common Pitfalls

- **Mutable default arguments**: Never use `def f(x=[])` — use `def f(x=None)` and set inside the function
- **Off-by-one errors**: Python ranges are half-open `[start, stop)` — bioinformatics coordinates are often 1-based
- **Deep vs shallow copy**: Nested data structures require `copy.deepcopy()` — `list.copy()` only copies the top level
