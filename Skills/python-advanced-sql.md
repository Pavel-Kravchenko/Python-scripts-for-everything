---
name: python-advanced-sql
description: OOP for bioinformatics classes, decorators, context managers, error handling, and SQL queries for biological databases
---

# Advanced Python & SQL for Bioinformatics

## When to Use
- Modeling sequences, genes, variants, alignments as Python objects
- Adding caching, validation, timing, or retry logic to pipeline functions
- Safely handling files, DB connections, and temp files via context managers
- Catching and propagating errors with domain-specific exception types
- Querying biological databases (SQLite, Ensembl MySQL, UCSC) with SQL

---

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

### Decorator Patterns
| Decorator | Use case |
|-----------|---------|
| `@timer` | Benchmark algorithms |
| `@memoize` / `@lru_cache` | Cache expensive calls |
| `@validate_sequence(chars)` | Guard invalid input |
| `@retry(n, delay)` | Unreliable network calls |
| `@contextmanager` | Simple context managers |

### SQL Clauses
| Clause | Use |
|--------|-----|
| `WHERE biotype = 'protein_coding'` | Filter rows |
| `GROUP BY tissue, condition` | Aggregate groups |
| `HAVING AVG(tpm) > 50` | Filter after grouping |
| `INNER JOIN` | Matching rows only |
| `LEFT JOIN` | All left rows, NULLs for no match |
| Subquery with `IN (SELECT ...)` | Multi-condition filter |

---

## Key Patterns

### BioSequence Hierarchy
```python
class BioSequence:
    def __init__(self, sequence, name="unnamed"):
        self.sequence = sequence.upper()
        self.name = name
    def __len__(self):      return len(self.sequence)
    def __str__(self):      return f">{self.name}\n{self.sequence}"
    def __contains__(self, motif): return motif.upper() in self.sequence
    def composition(self):
        return {c: self.sequence.count(c) for c in sorted(set(self.sequence))}

class DNA(BioSequence):
    VALID = set('ATGCN')
    def gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100
    def reverse_complement(self):
        table = str.maketrans('ATGC', 'TACG')
        return DNA(self.sequence.translate(table)[::-1], self.name)

class Protein(BioSequence):
    VALID = set('ACDEFGHIKLMNPQRSTVWY*')
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

    def overlaps(self, other):
        return self.chromosome == other.chromosome and self.start < other.end and other.start < self.end
```

### Properties for Validation
```python
@property
def sequence(self): return self._sequence

@sequence.setter
def sequence(self, value):
    value = value.upper()
    invalid = set(value) - set('ATGCN')
    if invalid:
        raise ValueError(f"Invalid nucleotides: {invalid}")
    self._sequence = value
```

### classmethod / staticmethod
```python
@classmethod
def from_fasta_string(cls, fasta_text):
    lines = fasta_text.strip().split('\n')
    parts = lines[0][1:].split(None, 1)
    return cls(parts[0], ''.join(lines[1:]), parts[1] if len(parts) > 1 else "")

@staticmethod
def is_valid_dna(seq):
    return set(seq.upper()).issubset(set('ATGCN'))
```

---

## Code Templates

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

# Decorator factory (parameterized)
def validate_sequence(valid_chars, seq_type="DNA"):
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

def retry(max_attempts=3, delay=0.5):
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

# Stack: applied bottom-up (validate first, then time)
@timer
@validate_sequence('ATGC', seq_type='DNA')
def analyze(seq): ...
```

### Context Managers
```python
# Class-based — for complex state
class FastaWriter:
    def __init__(self, filename, line_width=80):
        self.filename, self.line_width = filename, line_width
        self.file, self.record_count = None, 0

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()
        return False  # never suppress exceptions

    def write_record(self, seq_id, seq, desc=""):
        header = f">{seq_id}" + (f" {desc}" if desc else "")
        self.file.write(header + "\n")
        for i in range(0, len(seq), self.line_width):
            self.file.write(seq[i:i+self.line_width] + "\n")
        self.record_count += 1

# contextlib — for simple cases
from contextlib import contextmanager
import os, tempfile

@contextmanager
def timed_section(name):
    t = time.perf_counter()
    try:
        yield
    finally:
        print(f"[{name}] {time.perf_counter()-t:.4f}s")

@contextmanager
def temp_fasta(sequences):
    fd, path = tempfile.mkstemp(suffix='.fasta')
    try:
        with os.fdopen(fd, 'w') as f:
            for sid, seq in sequences.items():
                f.write(f">{sid}\n{seq}\n")
        yield path
    finally:
        os.unlink(path)
```

### Error Handling Hierarchy
```python
class BioinformaticsError(Exception):
    """Base for all pipeline errors."""

class InvalidSequenceError(BioinformaticsError):
    def __init__(self, invalid_chars, seq_type="DNA"):
        super().__init__(f"Invalid {seq_type} characters: {invalid_chars}")
        self.invalid_chars = invalid_chars

class SequenceLengthError(BioinformaticsError):
    def __init__(self, actual, minimum=None, maximum=None):
        msg = f"Sequence length {actual}"
        if minimum and actual < minimum: msg += f" < minimum {minimum}"
        if maximum and actual > maximum: msg += f" > maximum {maximum}"
        super().__init__(msg)

class FastaParseError(BioinformaticsError):
    def __init__(self, message, line_number=None, line_content=None):
        super().__init__(message)
        self.line_number, self.line_content = line_number, line_content

class TranslationError(BioinformaticsError): pass
```

### Graceful Batch Processing
```python
def batch_gc_analysis(sequences, strict=False):
    results, errors = {}, []
    for seq_id, seq in sequences.items():
        try:
            if not seq: raise SequenceLengthError(0, minimum=1)
            invalid = set(seq.upper()) - set('ATGCN')
            if invalid: raise InvalidSequenceError(invalid)
            s = seq.upper()
            results[seq_id] = round((s.count('G') + s.count('C')) / len(s) * 100, 2)
        except BioinformaticsError as e:
            if strict: raise
            errors.append((seq_id, str(e)))
    return results, errors
```

### Exception Chaining
```python
try:
    value = int(fields[3])
except ValueError as e:
    raise FastaParseError(f"Non-integer start position", line_number=n) from e
```

### Logging
```python
import logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger('bioinfo_pipeline')
# Use DEBUG for detail, INFO for progress, WARNING for bad-but-recoverable, ERROR for failures
```

### SQL — SQLite setup
```python
import sqlite3, pandas as pd

conn = sqlite3.connect(':memory:')  # or 'my_db.sqlite'
cursor = conn.cursor()
cursor.executescript('''
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
conn.commit()
```

### SQL — Common Bio Queries
```python
# Genes on a chromosome
pd.read_sql_query(
    "SELECT symbol, start_pos, end_pos FROM genes WHERE chromosome = 'chr17'", conn)

# Gene lengths > 100 kb
pd.read_sql_query("""
    SELECT symbol, (end_pos - start_pos) AS length
    FROM genes WHERE (end_pos - start_pos) > 100000
    ORDER BY length DESC
""", conn)

# Avg expression per tissue/condition
pd.read_sql_query("""
    SELECT tissue, condition, ROUND(AVG(tpm),2) AS avg_tpm, COUNT(*) AS n
    FROM expression GROUP BY tissue, condition ORDER BY tissue, condition
""", conn)

# Top tumor-expressed genes (HAVING filters aggregated groups)
pd.read_sql_query("""
    SELECT g.symbol, ROUND(AVG(e.tpm),2) AS avg_tumor_tpm
    FROM genes g JOIN expression e ON g.gene_id = e.gene_id
    WHERE e.condition = 'tumor'
    GROUP BY g.symbol HAVING AVG(e.tpm) > 50
    ORDER BY avg_tumor_tpm DESC
""", conn)

# Pathogenic variants with gene names (INNER JOIN)
pd.read_sql_query("""
    SELECT g.symbol, v.position, v.ref_allele, v.alt_allele
    FROM variants v INNER JOIN genes g ON v.gene_id = g.gene_id
    WHERE v.clinical_significance = 'pathogenic'
""", conn)

# All genes with variant count including 0s (LEFT JOIN)
pd.read_sql_query("""
    SELECT g.symbol, COUNT(v.variant_id) AS n_variants
    FROM genes g LEFT JOIN variants v ON g.gene_id = v.gene_id
    GROUP BY g.symbol ORDER BY n_variants DESC
""", conn)

# Genes highly expressed in tumors AND carrying pathogenic variants (subquery)
pd.read_sql_query("""
    SELECT symbol FROM genes
    WHERE gene_id IN (
        SELECT gene_id FROM expression WHERE condition='tumor'
        GROUP BY gene_id HAVING AVG(tpm) > 50
    ) AND gene_id IN (
        SELECT gene_id FROM variants WHERE clinical_significance='pathogenic'
    )
""", conn)

# Store DE results and query with CASE
pd.read_sql_query("""
    SELECT g.symbol, d.log2fc, d.padj,
           CASE WHEN d.significant=1 THEN 'Yes' ELSE 'No' END AS significant
    FROM de_results d JOIN genes g ON d.gene_id = g.gene_id
    ORDER BY d.padj
""", conn)
```

---

## Common Pitfalls

- **Missing `@functools.wraps`**: decorated function loses `__name__` and `__doc__`, breaking introspection.
- **Bare `except:`**: catches `SystemExit` and `KeyboardInterrupt`; always catch specific types.
- **`__exit__` returning `True`**: suppresses all exceptions silently — only do this intentionally.
- **`@lru_cache` on methods**: caches `self`, leaking instances; use on module-level or static functions.
- **Stacking decorators**: applied bottom-up — `@timer` wraps `@validate_sequence`, so validation runs first.
- **SQL `HAVING` vs `WHERE`**: `WHERE` filters rows before grouping; `HAVING` filters after aggregation.
- **`LEFT JOIN` count includes NULL**: use `COUNT(v.variant_id)` not `COUNT(*)` to count non-NULL matches.
- **`raise ... from e`**: preserve original traceback with exception chaining for pipeline debugging.
- **Properties without `_` backing attribute**: infinite recursion if setter assigns to `self.sequence` (not `self._sequence`).

---

## Related Skills
- `numpy-pandas-wrangling` — pandas merge/groupby mirror SQL JOIN/GROUP BY
- `biopython-databases` — SeqRecord, SeqIO replace the BioSequence hierarchy for production use
- `ml-deep-learning-bio` — sklearn pipelines benefit from the same decorator/context manager patterns
