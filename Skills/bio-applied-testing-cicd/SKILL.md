---
name: bio-applied-testing-cicd
description: "1. Explain **why** testing is critical for bioinformatics code reliability 2. Write **unit tests** with pytest for bioinformatics functions 3. Use **fixtures** and **parametrized tests** for efficient"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/12_Modern_Workflows/03_testing_cicd.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: biopython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Testing and CI/CD for Bioinformatics

*Source: Course notebook `Tier_3_Applied_Bioinformatics/12_Modern_Workflows/03_testing_cicd.ipynb`*

# Testing and CI/CD for Bioinformatics

**Tier 3 -- Applied Bioinformatics**

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Explain **why** testing is critical for bioinformatics code reliability
2. Write **unit tests** with pytest for bioinformatics functions
3. Use **fixtures** and **parametrized tests** for efficient test design
4. Set up **GitHub Actions** for continuous integration
5. Configure **code coverage** and quality checks
6. Apply best practices for testing biological data transformations

---

## 1. Why Test Bioinformatics Code?

### The Stakes Are High

Bioinformatics code directly influences scientific conclusions and, in clinical settings, patient care. Bugs can lead to:

| Bug Type | Example | Consequence |
|----------|---------|-------------|
| Off-by-one | 0-based vs 1-based coordinates | Wrong variants called |
| Strand error | Reverse complement not applied | Protein sequence wrong |
| Missing data | N bases or gaps not handled | Silent failures |
| Edge case | Empty sequence input | Crash during analysis |
| Rounding | Float precision in p-values | Wrong significance calls |

### The Benefits of Testing

- **Confidence**: Know your code works before analyzing real data
- **Regression prevention**: Changes don't break existing functionality
- **Documentation**: Tests show how functions should be used
- **Refactoring**: Safely improve code structure

---

## 2. Writing Tests with pytest

**pytest** is the most popular Python testing framework due to its:
- Simple syntax (plain `assert` statements)
- Powerful fixtures for test setup
- Rich plugin ecosystem

### A Bioinformatics Module to Test

Let's write tests for common bioinformatics functions:

```python
# bio_utils.py - Our bioinformatics module

def gc_content(sequence: str) -> float:
    """Calculate GC content as a percentage.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC percentage (0-100)
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence (A, T, G, C, N)
        
    Returns:
        Reverse complement string
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))


def find_motif(sequence: str, motif: str) -> list:
    """Find all positions of a motif in a sequence (1-based).
    
    Args:
        sequence: DNA/protein sequence to search
        motif: Pattern to find
        
    Returns:
        List of 1-based positions where motif starts
    """
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # Convert to 1-based
        start = pos + 1
    return positions


# Quick verification
print(f"GC content of 'ATGC': {gc_content('ATGC')}%")
print(f"Reverse complement of 'ATGC': {reverse_complement('ATGC')}")
print(f"'AT' positions in 'ATATGAT': {find_motif('ATATGAT', 'AT')}")
```

```python
# test_bio_utils.py - Our test file

test_code = '''
import pytest
from bio_utils import gc_content, reverse_complement, find_motif


class TestGCContent:
    """Tests for the gc_content function."""
    
    def test_balanced_sequence(self):
        """Test sequence with equal ATGC."""
        assert gc_content("ATGC") == 50.0
    
    def test_all_gc(self):
        """Test sequence with only G and C."""
        assert gc_content("GGCC") == 100.0
    
    def test_no_gc(self):
        """Test sequence with only A and T."""
        assert gc_content("AATT") == 0.0
    
    def test_empty_sequence(self):
        """Empty sequence should return 0, not raise error."""
        assert gc_content("") == 0.0
    
    def test_lowercase(self):
        """Should handle lowercase input."""
        assert gc_content("atgc") == 50.0
    
    def test_with_n_bases(self):
        """N bases count toward length but not GC."""
        assert gc_content("GCNN") == 50.0  # 2 GC out of 4 total


class TestReverseComplement:
    """Tests for the reverse_complement function."""
    
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"
    
    def test_palindromic_sequence(self):
        """Restriction site GAATTC (EcoRI) is palindromic."""
        assert reverse_complement("GAATTC") == "GAATTC"
    
    def test_handles_n(self):
        """N should remain N."""
        assert reverse_complement("ATNG") == "CNAT"
    
    def test_preserves_case(self):
        """Lowercase input produces lowercase output."""
        assert reverse_complement("atgc") == "gcat"
    
    def test_single_base(self):
        assert reverse_complement("A") == "T"


class TestFindMotif:
    """Tests for the find_motif function."""
    
    def test_single_occurrence(self):
        assert find_motif("ATGCGATCGA", "GAT") == [5]
    
    def test_multiple_occurrences(self):
        assert find_motif("ATGATGATG", "ATG") == [1, 4, 7]
    
    def test_overlapping_occurrences(self):
        """Overlapping matches should all be found."""
        assert find_motif("AAAA", "AA") == [1, 2, 3]
    
    def test_not_found(self):
        assert find_motif("ATGC", "GGG") == []
    
    def test_empty_sequence(self):
        assert find_motif("", "ATG") == []


# Parametrized test: run same test logic with multiple inputs
@pytest.mark.parametrize("seq,expected", [
    ("GGGG", 100.0),
    ("CCCC", 100.0),
    ("AAAA", 0.0),
    ("TTTT", 0.0),
    ("ATGC", 50.0),
    ("ATGCATGC", 50.0),
    ("GC", 100.0),
    ("AT", 0.0),
])
def test_gc_content_parametrized(seq, expected):
    """Test GC content with various sequences."""
    assert gc_content(seq) == expected
'''

print(test_code)
```

```python
# pytest command examples
commands = '''
# Run all tests in current directory
pytest

# Verbose output (show each test name)
pytest -v

# Run specific test file
pytest test_bio_utils.py

# Run specific test class
pytest test_bio_utils.py::TestGCContent

# Run specific test method
pytest test_bio_utils.py::TestGCContent::test_empty_sequence

# Stop at first failure
pytest -x

# Show print statements (not captured)
pytest -s

# Run with coverage report
pytest --cov=bio_utils --cov-report=html

# Run tests matching a keyword
pytest -k "gc or reverse"
'''
print(commands)
```

---

## 3. Fixtures and Test Organization

**Fixtures** provide reusable test data and setup. They're essential for:
- Loading test files
- Creating complex test objects
- Setting up and tearing down resources

```python
# conftest.py - Shared fixtures for all tests
fixtures_code = '''
import pytest
from io import StringIO
from pathlib import Path


@pytest.fixture
def sample_fasta_content():
    """Provide sample FASTA data as a string."""
    return """>gene1 Homo sapiens beta-globin
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC
GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG
>gene2 Homo sapiens alpha-globin
ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCG
CACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACC
"""


@pytest.fixture
def sample_fasta_file(tmp_path, sample_fasta_content):
    """Create a temporary FASTA file."""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(sample_fasta_content)
    return fasta_path


@pytest.fixture
def human_gene_sequences():
    """Dictionary of human gene sequences for testing."""
    return {
        "BRCA1_exon1": "ATGCGATCGATCGATCGATCG",
        "TP53_exon1": "GCTAGCTAGCTAGCTAGCTAG",
        "EGFR_exon1": "ATATATATATGCGCGCGCGC",
    }
'''

print(fixtures_code)
```

```python
# Using fixtures in tests
fixture_usage = '''
from Bio import SeqIO
from io import StringIO

def test_parse_fasta(sample_fasta_content):
    """Test that we can parse FASTA content."""
    records = list(SeqIO.parse(StringIO(sample_fasta_content), "fasta"))
    assert len(records) == 2
    assert records[0].id == "gene1"


def test_fasta_file_exists(sample_fasta_file):
    """Test that fixture creates a real file."""
    assert sample_fasta_file.exists()
    assert sample_fasta_file.suffix == ".fasta"


def test_gene_gc_content(human_gene_sequences):
    """Test GC content of fixture gene sequences."""
    from bio_utils import gc_content
    
    for gene, seq in human_gene_sequences.items():
        gc = gc_content(seq)
        assert 0 <= gc <= 100, f"Invalid GC for {gene}"
'''

print(fixture_usage)
```

---

## 4. GitHub Actions for CI/CD

**Continuous Integration (CI)** automatically runs tests on every push and pull request. This catches bugs before they reach the main branch.

**GitHub Actions** is GitHub's built-in CI/CD platform.

```python
# .github/workflows/tests.yml
github_actions_yaml = '''
name: Tests

on:
  push:
    branches: [main, master]
  pull_request:
    branches: [main, master]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.9", "3.10", "3.11", "3.12"]
    
    steps:
    - name: Checkout code
      uses: actions/checkout@v4
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v5
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install pytest pytest-cov
        pip install -r requirements.txt
    
    - name: Run tests with coverage
      run: |
        pytest --cov=src --cov-report=xml --cov-report=term-missing
    
    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v4
      with:
        file: coverage.xml
        fail_ci_if_error: false
'''

print(github_actions_yaml)
```

```python
# .github/workflows/lint.yml - Code quality checks
lint_yaml = '''
name: Lint and Format

on: [push, pull_request]

jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v4
    
    - name: Set up Python
      uses: actions/setup-python@v5
      with:
        python-version: "3.11"
    
    - name: Install linters
      run: pip install ruff black isort mypy
    
    - name: Check formatting with black
      run: black --check src/ tests/
    
    - name: Check import sorting with isort
      run: isort --check-only src/ tests/
    
    - name: Lint with ruff
      run: ruff check src/ tests/
    
    - name: Type check with mypy
      run: mypy src/ --ignore-missing-imports
'''

print(lint_yaml)
```

---

## 5. Best Practices

### Project Structure for Testable Code

```python
# Recommended project layout
project_structure = '''
my_bioinformatics_tool/
├── .github/
│   └── workflows/
│       ├── tests.yml
│       └── lint.yml
├── src/
│   └── my_tool/
│       ├── __init__.py
│       ├── io.py           # File parsing functions
│       ├── analysis.py     # Core algorithms
│       └── utils.py        # Helper functions
├── tests/
│   ├── __init__.py
│   ├── conftest.py         # Shared fixtures
│   ├── test_io.py
│   ├── test_analysis.py
│   └── data/               # Test data files
│       ├── sample.fasta
│       └── expected_output.txt
├── pyproject.toml          # Modern Python packaging
├── requirements.txt
└── README.md
'''

print(project_structure)
```

### Testing Checklist for Bioinformatics

| Category | What to Test |
|----------|-------------|
| **Edge cases** | Empty sequences, single-base sequences, very long sequences |
| **Case handling** | Lowercase, uppercase, mixed case |
| **Ambiguous bases** | N, R, Y, and other IUPAC codes |
| **Coordinate systems** | 0-based vs 1-based, inclusive vs exclusive |
| **Strand handling** | Forward, reverse, reverse complement |
| **File formats** | Malformed files, empty files, compressed files |
| **Numeric precision** | Float comparisons with tolerance |

---

## Exercises

### Exercise 1: Write Tests for Translation

Write pytest tests for a `translate()` function that converts DNA to protein.
