---
name: bio-applied-testing-cicd
description: "Testing and CI/CD for bioinformatics: pytest patterns, fixtures, GitHub Actions workflows, bioinformatics-specific test strategies"
tool_type: python
primary_tool: Python
---

# Testing and CI/CD for Bioinformatics

## Bioinformatics Test Module Example

```python
# bio_utils.py
def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    seq = sequence.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def reverse_complement(sequence: str) -> str:
    complement = {'A':'T','T':'A','G':'C','C':'G','N':'N',
                  'a':'t','t':'a','g':'c','c':'g','n':'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def find_motif(sequence: str, motif: str) -> list:
    positions, start = [], 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1: break
        positions.append(pos + 1)  # 1-based
        start = pos + 1
    return positions
```

## pytest Patterns

```python
import pytest

class TestGCContent:
    def test_balanced(self):       assert gc_content("ATGC") == 50.0
    def test_all_gc(self):         assert gc_content("GGCC") == 100.0
    def test_empty(self):          assert gc_content("") == 0.0
    def test_lowercase(self):      assert gc_content("atgc") == 50.0
    def test_with_n(self):         assert gc_content("GCNN") == 50.0

@pytest.mark.parametrize("seq,expected", [
    ("GGGG", 100.0), ("AAAA", 0.0), ("ATGC", 50.0),
])
def test_gc_parametrized(seq, expected):
    assert gc_content(seq) == expected
```

### pytest Commands

```bash
pytest -v                              # verbose
pytest test_file.py::TestClass::test   # specific test
pytest -x                              # stop at first failure
pytest --cov=bio_utils --cov-report=html  # coverage
pytest -k "gc or reverse"             # keyword filter
```

## Fixtures (conftest.py)

```python
@pytest.fixture
def sample_fasta_file(tmp_path):
    content = ">gene1\nATGCGATCGATCG\n>gene2\nGCTAGCTAGCTAG\n"
    path = tmp_path / "test.fasta"
    path.write_text(content)
    return path
```

## GitHub Actions CI

```yaml
# .github/workflows/tests.yml
name: Tests
on:
  push: { branches: [main] }
  pull_request: { branches: [main] }
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.9", "3.10", "3.11", "3.12"]
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v5
      with: { python-version: "${{ matrix.python-version }}" }
    - run: pip install pytest pytest-cov && pip install -r requirements.txt
    - run: pytest --cov=src --cov-report=xml --cov-report=term-missing
```

## Bioinformatics Testing Checklist

| Category | What to test |
|----------|-------------|
| Edge cases | Empty, single-base, very long sequences |
| Case handling | Lowercase, uppercase, mixed |
| Ambiguous bases | N, R, Y, IUPAC codes |
| Coordinates | 0-based vs 1-based, inclusive vs exclusive |
| Strand | Forward, reverse, reverse complement |
| File formats | Malformed, empty, compressed |
| Numeric | Float comparisons with tolerance |

## Project Structure

```
my_tool/
├── .github/workflows/{tests.yml, lint.yml}
├── src/my_tool/{__init__.py, io.py, analysis.py, utils.py}
├── tests/{conftest.py, test_io.py, test_analysis.py, data/}
├── pyproject.toml
└── requirements.txt
```

## Pitfalls

- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive — off-by-one errors
- **Float comparison**: use `pytest.approx()` for p-values and scores
- **Test data**: include small test files in `tests/data/`, not full datasets
