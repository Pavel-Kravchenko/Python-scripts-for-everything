# Rosalind - Assignments

Python wrapper and learning tools for working with [Rosalind](https://rosalind.info) bioinformatics problems.

## Files

### `rosalind.py`
A Python wrapper for fetching Rosalind problems, testing solutions, and generating notebook templates.
Includes a `RosalindTutor` for AI-assisted (Socratic) learning without giving away answers.

## How to Use

```python
from rosalind import Rosalind

ros = Rosalind()

# List first 10 problems
for p in ros.list_problems()[:10]:
    print(f"{p['id']} - {p['title']}")

# Get problem details and test your solution
def solve_dna(data):
    s = data.strip()
    return f"{s.count('A')} {s.count('C')} {s.count('G')} {s.count('T')}"

ros.test_solution("DNA", solve_dna)
# ✓ Solution passed for DNA!

# Generate a notebook scaffold for a problem
ros.create_notebook("RNA", output_dir="./problems")
```

## Learning Mode (AI-Assisted)

```python
from rosalind import RosalindTutor, copilot_help, learning_mode

tutor = RosalindTutor()
print(tutor.get_hint("GC", level=1))   # Vague hint
print(tutor.get_hint("GC", level=2))   # More specific
print(tutor.get_hint("GC", level=3))   # Most direct (still no answer!)
```

## Resources

- **Rosalind**: https://rosalind.info
- **BioPython**: https://biopython.org
