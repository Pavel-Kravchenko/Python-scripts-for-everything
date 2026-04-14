---
name: python-bio-context-managers
description: "Context managers for safe file, database, and resource handling in bioinformatics pipelines."
tool_type: python
primary_tool: Python
---

# Context Managers for Bioinformatics

## Lifecycle

```
with ctx as resource:
    __enter__() → acquire
    (body)
    __exit__() → release (always runs, even on exception)
```

`__exit__(exc_type, exc_val, tb)` — return `True` to suppress the exception, `False`/`None` to re-raise.

## Class-Based Context Manager

```python
class FastaWriter:
    def __init__(self, filename: str, line_width: int = 80):
        self.filename = filename
        self.line_width = line_width
        self.file = None
        self.record_count = 0

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, tb):
        if self.file:
            self.file.close()
        return False  # never suppress exceptions

    def write_record(self, seq_id: str, sequence: str, description: str = "") -> None:
        header = f">{seq_id}" + (f" {description}" if description else "")
        self.file.write(header + "\n")
        for i in range(0, len(sequence), self.line_width):
            self.file.write(sequence[i:i + self.line_width] + "\n")
        self.record_count += 1


with FastaWriter("output.fasta") as writer:
    writer.write_record("BRCA1", "ATGGATTTCGATCG" * 10, "breast cancer gene")
```

## `@contextmanager` (generator style)

Code before `yield` = `__enter__`. Code after `yield` = `__exit__`. Always use `try/finally` so teardown runs on exception.

```python
from contextlib import contextmanager
import tempfile, os

@contextmanager
def temp_fasta(sequences: dict[str, str]):
    """Write sequences to a temp file, yield path, auto-delete."""
    fd, path = tempfile.mkstemp(suffix=".fasta")
    try:
        with os.fdopen(fd, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n{seq}\n")
        yield path
    finally:
        os.unlink(path)


with temp_fasta({"seq1": "ATGC", "seq2": "TTAA"}) as path:
    # path is valid here; deleted after the block
    pass
```

## SQLite Transactions

```python
import sqlite3

with sqlite3.connect("variants.db") as conn:
    # auto-commits on success, rolls back on exception
    conn.execute("INSERT INTO variants VALUES (?, ?, ?)", ("chr1", 100, "A"))
```

## Multi-resource `with`

```python
with open("genome.fasta") as fasta, open("variants.vcf") as vcf:
    # both files open; both closed on exit
    ...
```

## Pitfalls

- **`yield` position matters**: code after `yield` is teardown — wrap it in `try/finally` or an exception in the body skips cleanup.
- **`__exit__` receives the exception, not `__enter__`**: if the `with` block raises, Python calls `__exit__` with exception info; forgetting to return `False` can accidentally suppress exceptions.
- **`with conn:` on SQLite commits, not closes**: `conn.close()` is separate; use `with closing(conn):` if you also want auto-close.
- **Nesting vs stacking**: `with A() as a, B() as b:` is equivalent to two nested `with` statements; `B.__exit__` runs before `A.__exit__`.
