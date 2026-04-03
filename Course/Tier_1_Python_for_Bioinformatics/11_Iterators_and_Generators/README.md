# 1.11 Iterators and Generators

**Tier 1: Python for Bioinformatics**

The iteration protocol under the hood, creating generators with `yield`, and the `itertools` module for lazy evaluation. Enables memory-efficient processing of genome-scale FASTA files and streaming FASTQ data.

## Topics Covered

- The iteration protocol: `__iter__` and `__next__`
- How `for` loops work under the hood
- Using `iter()` and `next()` manually
- Creating generators with `yield` and generator expressions
- `itertools`: `chain`, `combinations`, `permutations`, `product`, `groupby`, `islice`
- Lazy evaluation for large biological datasets

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_iterators_and_generators.ipynb](01_iterators_and_generators.ipynb) | Lazy iteration for memory-efficient bioinformatics data processing |

## Prerequisites

- [1.10 Comprehensions](../10_Comprehensions/)

---

[<- Previous Module](../10_Comprehensions/) | [Back to Course Overview](../../README.md) | [Next Module ->](../12_Regular_Expressions/)
