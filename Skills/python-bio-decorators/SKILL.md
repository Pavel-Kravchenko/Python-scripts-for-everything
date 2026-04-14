---
name: python-bio-decorators
description: "Python decorator patterns for bioinformatics: timing, validation, memoization, and stacking."
tool_type: python
primary_tool: Python
---

# Decorators for Bioinformatics

`@decorator` above `def f():` is exactly `f = decorator(f)` at definition time.

## Core Pattern

```python
import functools

def timer(func):
    @functools.wraps(func)   # preserves __name__, __doc__
    def wrapper(*args, **kwargs):
        import time
        start = time.perf_counter()
        result = func(*args, **kwargs)
        print(f"[timer] {func.__name__}: {time.perf_counter() - start:.4f}s")
        return result
    return wrapper

@timer
def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100
```

## Decorator Factory (with arguments)

Three levels: factory → decorator → wrapper.

```python
def validate_sequence(valid_chars: str, seq_type: str = "DNA"):
    valid_set = set(valid_chars.upper())

    def decorator(func):
        @functools.wraps(func)
        def wrapper(seq, *args, **kwargs):
            invalid = set(seq.upper()) - valid_set
            if invalid:
                raise ValueError(
                    f"Invalid {seq_type} characters {invalid} in {func.__name__}()"
                )
            return func(seq, *args, **kwargs)
        return wrapper
    return decorator


@validate_sequence('ATGC', seq_type='DNA')
def complement(seq: str) -> str:
    return seq.upper().translate(str.maketrans('ATGC', 'TACG'))
```

## Memoization

```python
# Manual memoize (cache by args tuple)
def memoize(func):
    cache = {}
    @functools.wraps(func)
    def wrapper(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]
    wrapper.cache = cache
    return wrapper

# Built-in (prefer this)
from functools import lru_cache

@lru_cache(maxsize=128)
def translate_codon(codon: str) -> str:
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    return table.get(codon.upper(), 'X')
```

## Stacking Decorators

Applied bottom-up (closest to the function first).

```python
@timer                              # applied second: f = timer(validate_sequence(...)(f))
@validate_sequence('ATGC')         # applied first
def analyze(seq: str) -> dict:
    seq = seq.upper()
    return {'length': len(seq), 'gc': (seq.count('G') + seq.count('C')) / len(seq) * 100}
```

## Closure Pattern

```python
def make_motif_counter(motif: str):
    motif = motif.upper()
    def counter(sequence: str) -> int:
        return sequence.upper().count(motif)
    return counter

count_cpg = make_motif_counter("CG")   # counter "remembers" motif
count_cpg("GCGCGCATCG")               # 3
```

## Pitfalls

- **Always use `functools.wraps`**: without it, `func.__name__` becomes `'wrapper'`, breaking logging, `help()`, and stack traces.
- **Decorator factories need 3 levels**: `@validate_sequence('ATGC')` requires factory → decorator → wrapper; a 2-level decorator receives the *argument* as the function, causing a confusing `TypeError`.
- **Stacking order matters**: `@A @B def f` = `A(B(f))`; validation should be inner (runs first), timing outer (measures total including validation).
- **`lru_cache` requires hashable arguments**: lists, dicts, and numpy arrays cannot be cached; convert to `tuple` or `bytes` before passing.
- **`lru_cache` holds strong references**: cached results are never GC'd until cache is cleared; use `maxsize` and call `.cache_clear()` when needed.
