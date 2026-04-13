---
name: algo-naive-pattern-matching
description: "**Problem Definition:** Given a text `T[0..n-1]` and a pattern `P[0..m-1]`, find all occurrences of pattern `P` in text `T`."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/01_naive_pattern_matching.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Naive Pattern Matching Algorithm

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/01_naive_pattern_matching.ipynb`*


## The Pattern Matching Problem

**Problem Definition:** Given a text `T[0..n-1]` and a pattern `P[0..m-1]`, find all occurrences of pattern `P` in text `T`.

This is one of the fundamental problems in computer science with applications in:
- Text editors (find/replace)
- Search engines
- DNA sequence analysis
- Plagiarism detection
- Network intrusion detection

## Naive (Brute Force) Approach

The naive algorithm uses a **sliding window** approach:
1. Place the pattern at each position `i` in the text (where `0 ≤ i ≤ n - m`)
2. Compare the pattern character by character with the text starting at position `i`
3. If all characters match, record position `i` as a match
4. Move to the next position regardless of match result

## ASCII Art Visualization

### Successful Match Example

```python
Text:    A B C A B C A B D
Pattern: A B C A B D

Attempt 1 (i=0):
Text:    A B C A B C A B D
         ↓ ↓ ↓ ↓ ↓ ✗
Pattern: A B C A B D
Match fails at position 5 (C ≠ D)

Attempt 2 (i=1):
Text:    A B C A B C A B D
           ✗
Pattern:   A B C A B D
Match fails at position 0 (B ≠ A)

Attempt 3 (i=2):
Text:    A B C A B C A B D
             ✗
Pattern:     A B C A B D
Match fails at position 0 (C ≠ A)

Attempt 4 (i=3):
Text:    A B C A B C A B D
               ↓ ↓ ↓ ↓ ↓ ↓
Pattern:       A B C A B D
✓ MATCH FOUND at position 3!
```python

### Why Naive is Inefficient

Consider the worst-case scenario:

```python
Text:    A A A A A A A A A B
Pattern: A A A A B

Attempt 1: A A A A ✗ (compare 5 chars, fail at last)
Attempt 2:   A A A A ✗ (compare 5 chars again!)
Attempt 3:     A A A A ✗
Attempt 4:       A A A A ✗
Attempt 5:         A A A A ✗
Attempt 6:           A A A A ↓ (finally matches!)

We keep re-comparing characters we already know match!
Total: O(n × m) comparisons
```python

The key insight: **We throw away useful information after each failed match.**

When the match fails at attempt 1, we already know that `T[1:5] = "AAAA"`. 
But we start fresh at attempt 2 and compare everything again!

## Complexity Analysis

| Aspect | Complexity | Explanation |
|--------|------------|-------------|
| **Time (best)** | O(n) | Pattern not in text, first char mismatches quickly |
| **Time (worst)** | O(n × m) | Many partial matches (e.g., `T = "AAA...AAB"`, `P = "AA...AB"`) |
| **Time (average)** | O(n × m) | Depends on alphabet size and pattern structure |
| **Space** | O(1) | Only a few index variables (O(k) for result storage with k matches) |

Where:
- `n` = length of text
- `m` = length of pattern

---

## Implementation

```python
def naive_search(text: str, pattern: str) -> list[int]:
    """
    Find all occurrences of pattern in text using naive approach.
    
    Args:
        text: The text to search in
        pattern: The pattern to find
        
    Returns:
        List of starting indices where pattern was found
        
    Time Complexity: O(n × m) where n=len(text), m=len(pattern)
    Space Complexity: O(1) auxiliary (O(k) for result with k matches)
    """
    n = len(text)
    m = len(pattern)
    positions = []
    
    # Edge cases
    if m == 0:
        return list(range(n + 1))  # Empty pattern matches everywhere
    if m > n:
        return []  # Pattern longer than text
    
    # Check each possible starting position
    for i in range(n - m + 1):
        # Compare pattern with text starting at position i
        match = True
        for j in range(m):
            if text[i + j] != pattern[j]:
                match = False
                break
        
        if match:
            positions.append(i)
    
    return positions
```python

### Alternative Implementation (Pythonic using slicing)

This version is simpler but uses string slicing internally.

```python
def naive_search_simple(text: str, pattern: str) -> list[int]:
    """
    Simpler implementation using Python string slicing.
    
    Note: String slicing creates a new string object each time,
    so this has O(m) space per comparison instead of O(1).
    """
    n, m = len(text), len(pattern)
    positions = []
    
    for i in range(n - m + 1):
        if text[i:i + m] == pattern:
            positions.append(i)
    
    return positions
```python

```python
# Basic examples
text = "ABCABCABD"
pattern = "ABCABD"

result = naive_search(text, pattern)
print(f"Text:    '{text}'")
print(f"Pattern: '{pattern}'")
print(f"Found at positions: {result}")
```python

```python
# Multiple occurrences
text = "ABABABA"
pattern = "ABA"

result = naive_search(text, pattern)
print(f"Text:    '{text}'")
print(f"Pattern: '{pattern}'")
print(f"Found at positions: {result}")
print("\nNote: Overlapping matches are found!")
```python

### Step-by-Step Execution with Comparison Counting

```python
def naive_search_verbose(text: str, pattern: str) -> tuple[list[int], int]:
    """
    Verbose version that shows each step and counts comparisons.
    
    Returns:
        Tuple of (positions, total_comparisons)
    """
    n, m = len(text), len(pattern)
    positions = []
    total_comparisons = 0
    
    print(f"Text:    '{text}' (length {n})")
    print(f"Pattern: '{pattern}' (length {m})")
    print(f"Need to check {n - m + 1} positions\n")
    print("=" * 50)
    
    for i in range(n - m + 1):
        print(f"\nAttempt at position i={i}:")
        print(f"  Text:    {' ' * i}{text[i:i+m]} (substring)")
        print(f"  Pattern: {' ' * i}{pattern}")
        
        match = True
        comparisons_this_attempt = 0
        
        for j in range(m):
            total_comparisons += 1
            comparisons_this_attempt += 1
            
            if text[i + j] == pattern[j]:
                print(f"  Compare T[{i+j}]='{text[i+j]}' with P[{j}]='{pattern[j]}': ✓ Match")
            else:
                print(f"  Compare T[{i+j}]='{text[i+j]}' with P[{j}]='{pattern[j]}': ✗ Mismatch")
                match = False
                break
        
        print(f"  → Comparisons in this attempt: {comparisons_this_attempt}")
        
        if match:
            positions.append(i)
            print(f"  ★ MATCH FOUND at position {i}!")
    
    print("\n" + "=" * 50)
    print(f"\nTotal comparisons: {total_comparisons}")
    print(f"Matches found at: {positions}")
    
    return positions, total_comparisons
```python

```python
# Run verbose example
positions, comparisons = naive_search_verbose("ABCABCABD", "ABCABD")
```python

### Best Case Example: Pattern Not Found Quickly

```python
def count_comparisons(text: str, pattern: str) -> int:
    """Count total character comparisons in naive search."""
    n, m = len(text), len(pattern)
    total = 0
    
    for i in range(n - m + 1):
        for j in range(m):
            total += 1
            if text[i + j] != pattern[j]:
                break
    return total


# Best case: pattern starts with unique character not in most of text
text_best = "AAAAAAAAAAAAAAAAAAAAB"
pattern_best = "BAAAA"

comparisons = count_comparisons(text_best, pattern_best)
n, m = len(text_best), len(pattern_best)

print("BEST CASE EXAMPLE")
print(f"Text:    '{text_best}'")
print(f"Pattern: '{pattern_best}'")
print(f"\nText length (n): {n}")
print(f"Pattern length (m): {m}")
print(f"Positions to check: {n - m + 1}")
print(f"Total comparisons: {comparisons}")
print(f"Maximum possible (n × m): {n * m}")
print(f"\nOnly ~1 comparison per position because 'B' ≠ 'A' immediately!")
```python

### Worst Case Example: Many Partial Matches

```python
# Worst case: many partial matches
text_worst = "AAAAAAAAAB"
pattern_worst = "AAAAB"

comparisons = count_comparisons(text_worst, pattern_worst)
n, m = len(text_worst), len(pattern_worst)

print("WORST CASE EXAMPLE")
print(f"Text:    '{text_worst}'")
print(f"Pattern: '{pattern_worst}'")
print(f"\nText length (n): {n}")
print(f"Pattern length (m): {m}")
print(f"Positions to check: {n - m + 1}")
print(f"Total comparisons: {comparisons}")
print(f"Maximum possible (n × m): {(n - m + 1) * m}")
print(f"\nAlmost all comparisons match until the last character!")
```python

```python
# Detailed worst case walkthrough
positions, comparisons = naive_search_verbose("AAAAAAAAAB", "AAAAB")
```python

---

## Performance Measurement

```python
import time
import random
import matplotlib.pyplot as plt

def generate_random_string(length: int, alphabet: str = "AB") -> str:
    """Generate a random string of given length."""
    return "".join(random.choice(alphabet) for _ in range(length))


def benchmark_naive_search(sizes: list[int], pattern_ratio: float = 0.1) -> tuple[list[int], list[float]]:
    """Benchmark naive search for various text sizes."""
    times = []
    
    for n in sizes:
        m = max(1, int(n * pattern_ratio))
        text = generate_random_string(n)
        pattern = generate_random_string(m)
        
        start = time.perf_counter()
        naive_search(text, pattern)
        end = time.perf_counter()
        
        times.append(end - start)
    
    return sizes, times
```python

```python
# Run benchmark
sizes = list(range(100, 10001, 500))
x, y = benchmark_naive_search(sizes, pattern_ratio=0.5)

plt.figure(figsize=(10, 6))
plt.plot(x, y, 'bo-', markersize=4)
plt.xlabel('Text Length (n)')
plt.ylabel('Time (seconds)')
plt.title('Naive Pattern Matching: Time vs Text Length\n(pattern length = n/2)')
plt.grid(True, alpha=0.3)
plt.show()

print(f"\nNote: With pattern length = n/2, we see O(n²/2) ≈ O(n²) behavior")
```python

---

## When is Naive Actually Fine?

Despite its worst-case O(n × m) complexity, the naive algorithm is often acceptable:

### 1. **Small Patterns**
When `m` is small (e.g., searching for a word in a document), O(n × m) ≈ O(n).

```python
# Small pattern example - effectively O(n)
large_text = generate_random_string(100000, "ABCDEFGHIJ")
small_pattern = "ABCD"

start = time.perf_counter()
result = naive_search(large_text, small_pattern)
end = time.perf_counter()

print(f"Text length: {len(large_text):,}")
print(f"Pattern length: {len(small_pattern)}")
print(f"Matches found: {len(result)}")
print(f"Time: {(end-start)*1000:.2f} ms")
print("\n→ Small patterns make naive search very fast!")
```python

### 2. **Large Alphabet**
With more characters, mismatches are found quickly.

```python
import string

# Compare binary alphabet vs full alphabet
n = 50000
m = 100

# Binary alphabet (worst case more likely)
text_binary = generate_random_string(n, "AB")
pattern_binary = generate_random_string(m, "AB")

start = time.perf_counter()
naive_search(text_binary, pattern_binary)
time_binary = time.perf_counter() - start

# Full ASCII alphabet (quick mismatches)
text_ascii = generate_random_string(n, string.ascii_letters)
pattern_ascii = generate_random_string(m, string.ascii_letters)

start = time.perf_counter()
naive_search(text_ascii, pattern_ascii)
time_ascii = time.perf_counter() - start

print(f"Binary alphabet (2 chars):  {time_binary*1000:.2f} ms")
print(f"Full alphabet (52 chars):   {time_ascii*1000:.2f} ms")
print(f"\nSpeedup: {time_binary/time_ascii:.1f}x faster with larger alphabet")
```python

### 3. **Few Expected Occurrences**
If the pattern rarely matches, comparisons are short.

```python
# Pattern with unique prefix
text = generate_random_string(100000, "ABCD")
unique_pattern = "XYZA"  # X and Y not in alphabet

comparisons = count_comparisons(text, unique_pattern)
print(f"Text length: {len(text):,}")
print(f"Pattern: '{unique_pattern}'")
print(f"Total comparisons: {comparisons:,}")
print(f"Positions checked: {len(text) - len(unique_pattern) + 1:,}")
print(f"Average comparisons per position: {comparisons / (len(text) - len(unique_pattern) + 1):.2f}")
print("\n→ Only 1 comparison per position (first char mismatch)!")
```python

---

## Summary & Motivation for Better Algorithms

### The Problem with Naive Search

The naive algorithm **forgets** what it learned during partial matches:

```python
Text:    A B C A B C A B D
         ↓ ↓ ↓ ↓ ↓ ✗
Pattern: A B C A B D

After this failure, we know:
- T[0:5] = "ABCAB"
- The pattern "ABCABD" has "AB" at both start and middle

Smart insight: We could skip ahead knowing T[3:5] = "AB" = P[0:2]!
But naive search starts fresh at position 1...
```python

### Better Algorithms

| Algorithm | Time Complexity | Key Idea |
|-----------|-----------------|----------|
| **KMP** | O(n + m) | Precompute failure function to skip redundant comparisons |
| **Rabin-Karp** | O(n + m) average | Use hashing to compare substrings in O(1) |
| **Boyer-Moore** | O(n/m) best | Skip large portions of text using bad character rule |

**Coming up next:** The KMP algorithm uses the **prefix function** to remember partial match information!

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
