---
name: python-advanced-sql
description: OOP for bioinformatics classes, decorators, context managers, error handling, and SQL queries for biological databases.
tool_type: python
primary_tool: Python
---

# Advanced Python & SQL for Bioinformatics

## When to Use
- Modeling sequences, genes, variants, alignments as Python objects
- Adding caching, validation, timing, or retry logic to pipeline functions
- Safely handling files, DB connections, and temp files via context managers
- Querying biological databases (SQLite, Ensembl MySQL, UCSC) with SQL

## Quick Reference

### OOP Dunders
| Method | Purpose |
|--------|---------|
| `__init__` | Constructor |
| `__str__` / `__repr__` | User / debug string |
| `__len__` | `len(obj)` |
| `__eq__`, `__lt__` | Comparison / sorting |
| `__contains__` | `"ATG" in seq` syntax |
| `__enter__` / `__exit__` | Context manager |

### SQL Clauses
| Clause | Use |
|--------|-----|
| `WHERE biotype = 'protein_coding'` | Filter rows |
| `GROUP BY tissue, condition` | Aggregate groups |
| `HAVING AVG(tpm) > 50` | Filter after grouping |
| `INNER JOIN` | Matching rows only |
| `LEFT JOIN` | All left rows, NULLs for no match |
| Subquery with `IN (SELECT ...)` | Multi-condition filter |


## Key Patterns

### Decorators
```python
import functools, time

def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t = time.perf_counter()
        result = func(*args, **kwargs)
        print(f"[timer] {func.__name__}: {time.perf_counter()-t:.4f}s")
        return result
    return wrapper

def memoize(func):
    cache = {}
    @functools.wraps(func)
    def wrapper(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]
    wrapper.cache = cache
    return wrapper

def validate_sequence(valid_chars: str, seq_type: str = "DNA"):
    valid_set = set(valid_chars.upper())
    def decorator(func):
        @functools.wraps(func)
        def wrapper(seq, *args, **kwargs):
            invalid = set(seq.upper()) - valid_set
            if invalid:
                raise ValueError(f"Invalid {seq_type} characters {invalid} in {func.__name__}()")
            return func(seq, *args, **kwargs)
        return wrapper
    return decorator

def retry(max_attempts: int = 3, delay: float = 0.5):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_err = None
            for attempt in range(1, max_attempts + 1):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_err = e
                    if attempt < max_attempts:
                        time.sleep(delay)
            raise last_err
        return wrapper
    return decorator

# Stacking: applied bottom-up (validate runs first, then timer wraps it)
@timer
@validate_sequence('ATGC', seq_type='DNA')
def analyze(seq: str): ...
```

### Context Managers
```python
class FastaWriter:
    def __init__(self, filename: str, line_width: int = 80):
        self.filename, self.line_width = filename, line_width
        self.file = None

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()
        return False  # never suppress exceptions

    def write_record(self, seq_id: str, seq: str, desc: str = ""):
        header = f">{seq_id}" + (f" {desc}" if desc else "")
        self.file.write(header + "\n")
        for i in range(0, len(seq), self.line_width):
            self.file.write(seq[i:i+self.line_width] + "\n")


# contextlib for simple cases
from contextlib import contextmanager
import os, tempfile

@contextmanager
def temp_fasta(sequences: dict[str, str]):
    fd, path = tempfile.mkstemp(suffix='.fasta')
    try:
        with os.fdopen(fd, 'w') as f:
            for sid, seq in sequences.items():
                f.write(f">{sid}\n{seq}\n")
        yield path
    finally:
        os.unlink(path)
```

### GeneAnnotation Dataclass
```python
from dataclasses import dataclass, field

@dataclass(order=True)
class GeneAnnotation:
    chromosome: str
    start: int
    end: int
    name: str      = field(compare=False, default="")
    strand: str    = field(compare=False, default='+')
    gene_type: str = field(compare=False, default="protein_coding")

    @property
    def length(self): return self.end - self.start

    def overlaps(self, other: "GeneAnnotation") -> bool:
        return (self.chromosome == other.chromosome
                and self.start < other.end
                and other.start < self.end)
```

### SQL — SQLite Setup
```python
import sqlite3, pandas as pd

conn = sqlite3.connect(':memory:')
conn.executescript('''
    CREATE TABLE genes (
        gene_id INTEGER PRIMARY KEY, symbol TEXT, chromosome TEXT,
        start_pos INTEGER, end_pos INTEGER, strand TEXT, biotype TEXT
    );
    CREATE TABLE expression (
        expr_id INTEGER PRIMARY KEY, gene_id INTEGER REFERENCES genes(gene_id),
        sample_id TEXT, tissue TEXT, tpm REAL, condition TEXT
    );
    CREATE TABLE variants (
        variant_id INTEGER PRIMARY KEY, gene_id INTEGER REFERENCES genes(gene_id),
        position INTEGER, ref_allele TEXT, alt_allele TEXT, clinical_significance TEXT
    );
''')
```

### SQL — Common Bio Queries
```python
# Gene lengths > 100 kb
pd.read_sql_query("""
    SELECT symbol, (end_pos - start_pos) AS length
    FROM genes WHERE (end_pos - start_pos) > 100000
    ORDER BY length DESC
""", conn)

# Avg expression per tissue/condition
pd.read_sql_query("""
    SELECT tissue, condition, ROUND(AVG(tpm),2) AS avg_tpm, COUNT(*) AS n
    FROM expression GROUP BY tissue, condition
""", conn)

# Top tumor-expressed genes (HAVING filters aggregated groups)
pd.read_sql_query("""
    SELECT g.symbol, ROUND(AVG(e.tpm),2) AS avg_tumor_tpm
    FROM genes g JOIN expression e ON g.gene_id = e.gene_id
    WHERE e.condition = 'tumor'
    GROUP BY g.symbol HAVING AVG(e.tpm) > 50
    ORDER BY avg_tumor_tpm DESC
""", conn)

# All genes with variant count including 0s (LEFT JOIN)
pd.read_sql_query("""
    SELECT g.symbol, COUNT(v.variant_id) AS n_variants
    FROM genes g LEFT JOIN variants v ON g.gene_id = v.gene_id
    GROUP BY g.symbol ORDER BY n_variants DESC
""", conn)

# Genes highly expressed in tumors AND carrying pathogenic variants
pd.read_sql_query("""
    SELECT symbol FROM genes
    WHERE gene_id IN (
        SELECT gene_id FROM expression WHERE condition='tumor'
        GROUP BY gene_id HAVING AVG(tpm) > 50
    ) AND gene_id IN (
        SELECT gene_id FROM variants WHERE clinical_significance='pathogenic'
    )
""", conn)
```


## Pitfalls

- **Missing `@functools.wraps`**: decorated function loses `__name__` and `__doc__`, breaking introspection and logging.
- **Bare `except:`**: catches `SystemExit` and `KeyboardInterrupt`; always catch specific types.
- **`__exit__` returning `True`**: suppresses all exceptions silently — only do this intentionally.
- **`@lru_cache` on methods**: caches `self`, leaking instances; use on module-level or static functions only.
- **Stacking decorators**: applied bottom-up — `@timer` wraps `@validate_sequence`, so validation runs first.
- **SQL `HAVING` vs `WHERE`**: `WHERE` filters rows before grouping; `HAVING` filters after aggregation.
- **`LEFT JOIN` count**: use `COUNT(v.variant_id)` not `COUNT(*)` to count non-NULL matches.
- **`raise ... from e`**: preserves original traceback; omitting `from e` hides the root cause.
- **Properties without `_` backing attribute**: `self.sequence = value` inside a `sequence` setter causes infinite recursion; use `self._sequence`.
