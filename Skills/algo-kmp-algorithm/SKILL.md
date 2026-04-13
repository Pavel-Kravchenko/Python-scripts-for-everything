---
name: algo-kmp-algorithm
description: "The KMP algorithm is one of the most elegant string matching algorithms, achieving **O(n + m)** time complexity by avoiding redundant comparisons. This notebook covers:"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/02_kmp_algorithm.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Knuth-Morris-Pratt (KMP) Algorithm

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/02_kmp_algorithm.ipynb`*

# Knuth-Morris-Pratt (KMP) Algorithm

The KMP algorithm is one of the most elegant string matching algorithms, achieving **O(n + m)** time complexity by avoiding redundant comparisons. This notebook covers:

1. **The Key Insight** - Why naive string matching wastes work
2. **Prefix Function (Failure Function)** - The heart of KMP
3. **Building the Prefix Table** - Naive O(n²) vs Optimized O(n)
4. **KMP Algorithm** - Using the prefix table for efficient search

## 1. The Key Insight: Don't Restart from Scratch

### Why Does the Naive Algorithm Waste Work?

In naive string matching, after a mismatch, we shift the pattern by 1 and restart comparison from the beginning. But this throws away valuable information!

```
Text:    A B A B A B A C A B A
Pattern: A B A B A C A
         ✓ ✓ ✓ ✓ ✓ ✗     <- mismatch at position 5

Naive approach: Shift pattern by 1, restart from scratch
Text:    A B A B A B A C A B A
Pattern:   A B A B A C A
           ✗               <- immediately fails!

We already knew that text[2..4] = "ABA" because we just matched it!
Why compare again?
```

### The Smart Observation

When we have a mismatch after matching some characters, we know:
- What characters we just matched in the text
- These characters are part of the pattern

**Key Question:** After matching `P[0..j-1]`, can we find a shorter prefix of the pattern that would still match the end of what we just saw?

```
Pattern matched so far: A B A B A
                        ↑───↑       This prefix "ABA"...
                            ↑───↑   ...equals this suffix!

So after mismatch, we can continue from position 3 in the pattern!
We don't need to re-match "ABA" - we know it's already there.
```

## 2. The Prefix Function (Failure Function)

### Definition

For a pattern `P` of length `m`, the **prefix function** π is defined as:

$$\pi[i] = \text{length of the longest proper prefix of } P[0..i] \text{ that is also a suffix}$$

**Proper prefix:** A prefix that is not the entire string.

### Example: Pattern "ABABACA"

```
Pattern: A B A B A C A

Index:   0   1   2   3   4   5   6
Char:    A   B   A   B   A   C   A
π[i]:    0   0   1   2   3   0   1
```

Let's verify each value:

```
i=0: P[0..0] = "A"
     No proper prefix exists → π[0] = 0

i=1: P[0..1] = "AB"
     Proper prefixes: "A"
     Suffixes: "B"
     No match → π[1] = 0

i=2: P[0..2] = "ABA"
     Proper prefixes: "A", "AB"
     Suffixes: "A", "BA"
     Match: "A" (length 1) → π[2] = 1

i=3: P[0..3] = "ABAB"
     Proper prefixes: "A", "AB", "ABA"
     Suffixes: "B", "AB", "BAB"
     Match: "AB" (length 2) → π[3] = 2

i=4: P[0..4] = "ABABA"
     Proper prefixes: "A", "AB", "ABA", "ABAB"
     Suffixes: "A", "BA", "ABA", "BABA"
     Match: "ABA" (length 3) → π[4] = 3

i=5: P[0..5] = "ABABAC"
     No proper prefix equals any suffix → π[5] = 0

i=6: P[0..6] = "ABABACA"
     Match: "A" (length 1) → π[6] = 1
```

### Visual Explanation for π[4] = 3

```
P[0..4] = "ABABA"
Longest proper prefix that is also suffix = "ABA" (length 3)

A B A B A
↑─────↑      prefix  P[0..2]
    ↑─────↑  suffix  P[2..4]
    
Both are "ABA" → π[4] = 3
```

## 3. Building the Prefix Table

### 3.1 Naive Approach - O(n²)

For each position, check all possible prefix lengths:

```python
def build_prefix_naive(pattern: str) -> list[int]:
    """
    Build prefix function using naive O(n²) approach.
    
    For each position i, checks all possible prefix lengths from longest
    to shortest until finding one that matches the suffix.
    
    Args:
        pattern: The pattern string to build prefix table for
        
    Returns:
        List where pi[i] = length of longest proper prefix of pattern[0..i]
        that is also a suffix
        
    Time Complexity: O(n²) where n = len(pattern)
    Space Complexity: O(n) for the prefix table
    """
    n = len(pattern)
    pi = [0] * n
    
    for i in range(1, n):
        # Try all possible prefix lengths from longest to shortest
        for length in range(i, 0, -1):
            prefix = pattern[0:length]
            suffix = pattern[i - length + 1:i + 1]
            if prefix == suffix:
                pi[i] = length
                break
    
    return pi


# Test
pattern = "ABABACA"
print(f"Pattern: {pattern}")
print(f"Prefix table (naive): {build_prefix_naive(pattern)}")
```

### 3.2 Optimized Approach - O(n)

The key insight: **use previously computed values to avoid redundant work!**

If we know π[i-1] = k, then for position i:
- We already have P[0..k-1] = P[i-k..i-1]
- Just check if P[k] == P[i]
- If yes: π[i] = k + 1
- If no: fall back to π[k-1] and try again

### Step-by-Step Example

```
Pattern: A B A B A C A
         0 1 2 3 4 5 6

Initialize: pi[0] = 0, j = 0

i=1: Compare P[1]='B' with P[0]='A'
     B ≠ A, j=0 can't fall back
     → π[1] = 0

i=2: Compare P[2]='A' with P[0]='A'
     A = A → j becomes 1
     → π[2] = 1

i=3: Compare P[3]='B' with P[1]='B'
     B = B → j becomes 2
     → π[3] = 2

i=4: Compare P[4]='A' with P[2]='A'
     A = A → j becomes 3
     → π[4] = 3

i=5: Compare P[5]='C' with P[3]='B'
     C ≠ B → fall back: j = π[2] = 1
     Compare P[5]='C' with P[1]='B'
     C ≠ B → fall back: j = π[0] = 0
     Compare P[5]='C' with P[0]='A'
     C ≠ A, j=0 can't fall back
     → π[5] = 0

i=6: Compare P[6]='A' with P[0]='A'
     A = A → j becomes 1
     → π[6] = 1

Final: [0, 0, 1, 2, 3, 0, 1]
```

```python
def build_prefix(pattern: str) -> list[int]:
    """
    Build prefix function using optimized O(n) approach.
    
    Uses the key insight that if π[i-1] = k, then we only need to check
    if pattern[k] == pattern[i]. If not, we fall back using previously
    computed π values.
    
    Args:
        pattern: The pattern string to build prefix table for
        
    Returns:
        List where pi[i] = length of longest proper prefix of pattern[0..i]
        that is also a suffix
        
    Time Complexity: O(n) where n = len(pattern)
    Space Complexity: O(n) for the prefix table
    """
    n = len(pattern)
    if n == 0:
        return []
    
    pi = [0] * n
    j = 0  # length of previous longest prefix-suffix
    
    for i in range(1, n):
        # Fall back until we find a match or reach the beginning
        while j > 0 and pattern[j] != pattern[i]:
            j = pi[j - 1]
        
        # If characters match, extend the prefix-suffix
        if pattern[j] == pattern[i]:
            j += 1
        
        pi[i] = j
    
    return pi


# Test and verify
pattern = "ABABACA"
print(f"Pattern: {pattern}")
print(f"Prefix table (optimized): {build_prefix(pattern)}")
print(f"Prefix table (naive):     {build_prefix_naive(pattern)}")
assert build_prefix(pattern) == build_prefix_naive(pattern), "Results should match!"
```

### Why is the Optimized Version O(n)?

It looks like there's a nested loop (while inside for), suggesting O(n²). But consider:

- `j` can only increase by 1 per iteration of the outer loop (at most n increases)
- `j` decreases in the while loop, but can never go below 0
- Total decreases ≤ total increases ≤ n

Therefore, the total number of operations is **O(n)**.

This is called **amortized analysis** - while individual iterations may do more work, the total work is bounded.

```python
def build_prefix_verbose(pattern: str) -> list[int]:
    """
    Build prefix function with step-by-step output for learning.
    
    Args:
        pattern: The pattern string
        
    Returns:
        The prefix table
    """
    n = len(pattern)
    if n == 0:
        return []
    
    pi = [0] * n
    j = 0
    
    print(f"Pattern: {pattern}")
    print(f"         {''.join(str(i) for i in range(n))}")
    print(f"\nBuilding prefix table:\n")
    print(f"pi[0] = 0 (by definition)\n")
    
    for i in range(1, n):
        print(f"i={i}: Checking P[{i}]='{pattern[i]}'")
        
        while j > 0 and pattern[j] != pattern[i]:
            print(f"      P[{j}]='{pattern[j]}' ≠ P[{i}]='{pattern[i]}' → fall back to j=pi[{j-1}]={pi[j-1]}")
            j = pi[j - 1]
        
        if pattern[j] == pattern[i]:
            j += 1
            print(f"      P[{j-1}]='{pattern[j-1]}' = P[{i}]='{pattern[i]}' → extend to j={j}")
        else:
            print(f"      P[{j}]='{pattern[j]}' ≠ P[{i}]='{pattern[i]}', j=0 → no extension")
        
        pi[i] = j
        print(f"      → pi[{i}] = {j}\n")
    
    print(f"Final prefix table: {pi}")
    return pi


# Demonstrate with verbose output
_ = build_prefix_verbose("ABABACA")
```

## 4. The KMP Algorithm

### How KMP Uses the Prefix Function

When searching for pattern `P` in text `T`:
1. If `T[i] == P[j]`: both advance (match)
2. If `T[i] != P[j]` and `j > 0`: use prefix function to shift pattern
   - Set `j = π[j-1]` (the pattern slides so the prefix matches what we've seen)
3. If `T[i] != P[j]` and `j == 0`: move to next character in text

### Visual Example

```
Text:    A B A B A B A C A B A
Pattern: A B A B A C A
π:       [0,0,1,2,3,0,1]

Step 1-5: Match first 5 characters (positions 0-4)
Text:    A B A B A B A C A B A
         ↓ ↓ ↓ ↓ ↓ ✗
Pattern: A B A B A C A
         0 1 2 3 4 5
                   ^ mismatch at j=5 (P[5]='C' ≠ T[5]='B')

Instead of restarting, use π[4] = 3:
- We know P[0..2] = P[2..4] = "ABA"
- The text characters T[2..4] = "ABA" (we just matched them!)
- So P[0..2] already matches T[2..4]
- Jump: j = π[j-1] = π[4] = 3

Text:    A B A B A B A C A B A
               ↓ ↓ ↓ ↓ ↓ ↓ ↓
Pattern:       A B A B A C A
               0 1 2 3 4 5 6
                     ^ continue from j=3, i stays at 5

Text index never goes backward! Always O(n) for text traversal.
```

```python
def kmp_search(text: str, pattern: str) -> list[int]:
    """
    Find all occurrences of pattern in text using KMP algorithm.
    
    The KMP algorithm preprocesses the pattern to build a prefix table,
    then uses it to avoid redundant comparisons during search. The text
    pointer never moves backward.
    
    Args:
        text: The text to search in
        pattern: The pattern to search for
        
    Returns:
        List of starting indices where pattern occurs in text
        
    Time Complexity: O(n + m) where n = len(text), m = len(pattern)
    Space Complexity: O(m) for the prefix table
    """
    if not pattern:
        return []
    if len(pattern) > len(text):
        return []
    
    # Build prefix table
    pi = build_prefix(pattern)
    
    matches = []
    j = 0  # position in pattern
    
    for i in range(len(text)):
        # On mismatch, use prefix table to shift pattern
        while j > 0 and text[i] != pattern[j]:
            j = pi[j - 1]
        
        # If characters match, advance in pattern
        if text[i] == pattern[j]:
            j += 1
        
        # If complete pattern matched
        if j == len(pattern):
            matches.append(i - j + 1)
            j = pi[j - 1]  # Look for next occurrence
    
    return matches


# Test examples
print("Example 1:")
text1 = "ABABABABACA"
pattern1 = "ABABACA"
result1 = kmp_search(text1, pattern1)
print(f"Text:    '{text1}'")
print(f"Pattern: '{pattern1}'")
print(f"Found at positions: {result1}")

print("\nExample 2 (overlapping matches):")
text2 = "AAAAAA"
pattern2 = "AAA"
result2 = kmp_search(text2, pattern2)
print(f"Text:    '{text2}'")
print(f"Pattern: '{pattern2}'")
print(f"Found at positions: {result2}")

print("\nExample 3:")
text3 = "ABCABCABCAB"
pattern3 = "ABCABCAB"
result3 = kmp_search(text3, pattern3)
print(f"Text:    '{text3}'")
print(f"Pattern: '{pattern3}'")
print(f"Found at positions: {result3}")
```

```python
def kmp_search_verbose(text: str, pattern: str) -> list[int]:
    """
    KMP search with detailed step-by-step output.
    
    Args:
        text: The text to search in
        pattern: The pattern to search for
        
    Returns:
        List of starting indices where pattern occurs
    """
    if not pattern or len(pattern) > len(text):
        return []
    
    pi = build_prefix(pattern)
    print(f"Text:    '{text}'")
    print(f"Pattern: '{pattern}'")
    print(f"Prefix:  {pi}")
    print("\n" + "="*50 + "\n")
    
    matches = []
    j = 0
    comparisons = 0
    
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            comparisons += 1
            print(f"i={i}, j={j}: T[{i}]='{text[i]}' ≠ P[{j}]='{pattern[j]}' → shift: j=pi[{j-1}]={pi[j-1]}")
            j = pi[j - 1]
        
        comparisons += 1
        if text[i] == pattern[j]:
            print(f"i={i}, j={j}: T[{i}]='{text[i]}' = P[{j}]='{pattern[j]}' → match, j++")
            j += 1
        else:
            print(f"i={i}, j={j}: T[{i}]='{text[i]}' ≠ P[{j}]='{pattern[j]}' → no match, move to next i")
        
        if j == len(pattern):
            pos = i - j + 1
            matches.append(pos)
            print(f"  *** FOUND at position {pos}! ***")
            j = pi[j - 1]
    
    print(f"\nTotal comparisons: {comparisons}")
    return matches


# Demonstrate
print("KMP Search Step-by-Step:")
print()
_ = kmp_search_verbose("ABABABABACA", "ABABACA")
```
