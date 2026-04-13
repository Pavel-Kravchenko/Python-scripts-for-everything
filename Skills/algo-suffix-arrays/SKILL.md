---
name: algo-suffix-arrays
description: "A comprehensive guide to suffix arrays - a space-efficient alternative to suffix trees for string processing."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/03_suffix_arrays.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Suffix Arrays

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/08_Advanced_String_Structures/03_suffix_arrays.ipynb`*

# Suffix Arrays

A comprehensive guide to suffix arrays - a space-efficient alternative to suffix trees for string processing.

## 1. What is a Suffix Array?

A **suffix array** is a sorted array of all suffixes of a string. More precisely:

- Given a string `S` of length `n`
- The suffix array `SA` is an array of integers from `0` to `n-1`
- `SA[i]` = starting position of the `i`-th smallest suffix (in lexicographic order)

### Key Properties

1. **Space-efficient**: Uses O(n) space (just an array of integers), compared to O(n) with larger constants for suffix trees
2. **Enables fast pattern matching**: Binary search in O(m log n) time
3. **Foundation for many algorithms**: Used in bioinformatics, data compression, and text indexing

### Suffix Array Definition

```
SA[i] = starting position of the i-th smallest suffix in lexicographic order
```

This means if we sort all suffixes alphabetically, `SA[i]` tells us where the `i`-th suffix in that sorted order originally started in the string.

## ASCII Art Visualization: Suffix Array for "banana$"

```
All suffixes:
Index  Suffix
  0    banana$
  1    anana$
  2    nana$
  3    ana$
  4    na$
  5    a$
  6    $

Sorted suffixes:
Rank   Suffix      Original position (SA value)
  0    $           6
  1    a$          5
  2    ana$        3
  3    anana$      1
  4    banana$     0
  5    na$         4
  6    nana$       2

Suffix Array SA = [6, 5, 3, 1, 0, 4, 2]

Reading the array:
  SA[0] = 6  →  Smallest suffix starts at position 6: "$"
  SA[1] = 5  →  2nd smallest starts at position 5: "a$"
  SA[2] = 3  →  3rd smallest starts at position 3: "ana$"
  ...
```

## 2. Naive Construction - O(n² log n)

The simplest approach: generate all suffixes and sort them.

### Algorithm
1. Generate all n suffixes with their starting positions
2. Sort the suffixes lexicographically
3. Extract the starting positions

### Complexity Analysis
- Generating suffixes: O(n)
- Sorting: O(n log n) comparisons
- Each comparison: O(n) in worst case (comparing strings)
- **Total: O(n² log n)**

```python
def build_suffix_array_naive(text):
    """
    Build suffix array using naive O(n² log n) approach.
    
    Creates all suffixes, sorts them lexicographically, and returns
    an array of starting positions in sorted order.
    
    Parameters
    ----------
    text : str
        Input string (should end with a unique terminal character like '$')
    
    Returns
    -------
    list[int]
        Suffix array - SA[i] is the starting position of the i-th smallest suffix
    
    Examples
    --------
    >>> build_suffix_array_naive('banana$')
    [6, 5, 3, 1, 0, 4, 2]
    
    Time Complexity: O(n² log n)
    Space Complexity: O(n²) for storing all suffixes
    """
    n = len(text)
    # Create list of (suffix, starting_position) pairs
    suffixes = [(text[i:], i) for i in range(n)]
    # Sort by suffix string
    suffixes.sort(key=lambda x: x[0])
    # Extract starting positions
    return [pos for suffix, pos in suffixes]
```

```python
# Test naive construction
text = "banana$"
sa = build_suffix_array_naive(text)
print(f"Text: {text}")
print(f"Suffix Array: {sa}")
print()
print("Verification - suffixes in sorted order:")
for i, pos in enumerate(sa):
    print(f"  SA[{i}] = {pos}  →  '{text[pos:]}'")
```

### Step-by-Step Naive Construction

```python
def build_suffix_array_naive_verbose(text):
    """
    Build suffix array with detailed step-by-step output.
    
    Parameters
    ----------
    text : str
        Input string
    
    Returns
    -------
    list[int]
        Suffix array
    """
    n = len(text)
    
    print("Step 1: Generate all suffixes")
    print("-" * 40)
    suffixes = []
    for i in range(n):
        suffix = text[i:]
        suffixes.append((suffix, i))
        print(f"  Position {i}: '{suffix}'")
    
    print()
    print("Step 2: Sort suffixes lexicographically")
    print("-" * 40)
    suffixes.sort(key=lambda x: x[0])
    for i, (suffix, pos) in enumerate(suffixes):
        print(f"  Rank {i}: '{suffix}' (originally at position {pos})")
    
    print()
    print("Step 3: Extract starting positions")
    print("-" * 40)
    sa = [pos for suffix, pos in suffixes]
    print(f"  Suffix Array SA = {sa}")
    
    return sa

sa = build_suffix_array_naive_verbose("banana$")
```

## 3. Manber-Myers Algorithm - O(n log n)

The Manber-Myers algorithm uses a clever **doubling technique** to achieve O(n log n) construction time.

### Key Idea
Instead of comparing entire suffixes (which can be O(n) per comparison), we:
1. First sort by the first character (k=1)
2. Then sort by the first 2 characters (k=2)
3. Then sort by the first 4 characters (k=4)
4. Continue doubling until k ≥ n

The trick: When sorting by 2k characters, we use the **ranks from the k-character sort** as keys. This makes each comparison O(1)!

### ASCII Art: Manber-Myers Doubling Process

```
Text: "banana$"

═══════════════════════════════════════════════════════════════════
Round 1 (k=1): Sort by first 1 character
═══════════════════════════════════════════════════════════════════
  Position:  0   1   2   3   4   5   6
  Character: b   a   n   a   n   a   $
  Rank:      2   1   3   1   3   1   0
                │   │   │   │   │   │
                └───┴───┴───┘   │   │
                  'a' gets rank 1   │
                                    │
                              '$' gets rank 0 (smallest)

═══════════════════════════════════════════════════════════════════
Round 2 (k=2): Sort by first 2 characters
═══════════════════════════════════════════════════════════════════
  Use pairs: (rank[i], rank[i+k]) where k=1
  
  Position 0: (rank[0], rank[1]) = (2, 1)  → "ba"
  Position 1: (rank[1], rank[2]) = (1, 3)  → "an"
  Position 2: (rank[2], rank[3]) = (3, 1)  → "na"
  Position 3: (rank[3], rank[4]) = (1, 3)  → "an"
  Position 4: (rank[4], rank[5]) = (3, 1)  → "na"
  Position 5: (rank[5], rank[6]) = (1, 0)  → "a$"
  Position 6: (rank[6], -)       = (0, -)  → "$"

═══════════════════════════════════════════════════════════════════
Round 3 (k=4): Sort by first 4 characters
═══════════════════════════════════════════════════════════════════
  Use pairs: (rank[i], rank[i+k]) where k=2
  Continue until all suffixes have unique ranks or k ≥ n

After log(n) rounds → Fully sorted suffix array!
```

```python
import collections

def build_suffix_array_manber_myers(text):
    """
    Build suffix array using Manber-Myers doubling algorithm.
    
    Uses recursive bucket sorting with doubling to achieve O(n log n)
    time complexity. Each round doubles the number of characters used
    for comparison.
    
    Parameters
    ----------
    text : str
        Input string (should end with a unique terminal character like '$')
    
    Returns
    -------
    list[int]
        Suffix array - SA[i] is the starting position of the i-th smallest suffix
    
    Examples
    --------
    >>> build_suffix_array_manber_myers('banana$')
    [6, 5, 3, 1, 0, 4, 2]
    
    Time Complexity: O(n log n)
    Space Complexity: O(n)
    """
    def sort_bucket(s, bucket, order):
        """
        Recursively sort a bucket of suffix indices.
        
        Parameters
        ----------
        s : str
            The input string
        bucket : iterable
            Indices to sort
        order : int
            Current prefix length for comparison (doubles each round)
        
        Returns
        -------
        list[int]
            Sorted indices
        """
        # Group indices by their prefix of length 'order'
        d = collections.defaultdict(list)
        for i in bucket:
            key = s[i:i + order]
            d[key].append(i)
        
        result = []
        # Process groups in sorted order
        for k, v in sorted(d.items()):
            if len(v) > 1:
                # Multiple suffixes share this prefix - recurse with doubled order
                result += sort_bucket(s, v, order * 2)
            else:
                # Single suffix - it's in its final position
                result.append(v[0])
        return result
    
    return sort_bucket(text, range(len(text)), 1)
```

```python
# Test Manber-Myers construction
text = "banana$"
sa_mm = build_suffix_array_manber_myers(text)
print(f"Text: {text}")
print(f"Suffix Array (Manber-Myers): {sa_mm}")
print()
print("Verification:")
for i, pos in enumerate(sa_mm):
    print(f"  SA[{i}] = {pos}  →  '{text[pos:]}'")
```

```python
def build_suffix_array_manber_myers_verbose(text):
    """
    Build suffix array with Manber-Myers algorithm showing each round.
    
    Parameters
    ----------
    text : str
        Input string
    
    Returns
    -------
    list[int]
        Suffix array
    """
    def sort_bucket(s, bucket, order, depth=0):
        indent = "  " * depth
        bucket = list(bucket)
        
        if depth == 0:
            print(f"\n{'='*60}")
            print(f"Round: Sorting by first {order} character(s)")
            print(f"{'='*60}")
        
        # Group by prefix
        d = collections.defaultdict(list)
        for i in bucket:
            key = s[i:i + order]
            d[key].append(i)
        
        print(f"{indent}Bucket {bucket} grouped by prefix (len={order}):")
        for k, v in sorted(d.items()):
            print(f"{indent}  '{k}': positions {v}")
        
        result = []
        for k, v in sorted(d.items()):
            if len(v) > 1:
                print(f"{indent}  → '{k}' has {len(v)} suffixes, need to recurse with order={order*2}")
                result += sort_bucket(s, v, order * 2, depth + 1)
            else:
                print(f"{indent}  → '{k}' unique, position {v[0]} is final")
                result.append(v[0])
        
        return result
    
    print(f"Text: '{text}'")
    print(f"Length: {len(text)}")
    sa = sort_bucket(text, range(len(text)), 1)
    print(f"\nFinal Suffix Array: {sa}")
    return sa

sa = build_suffix_array_manber_myers_verbose("banana$")
```

### Alternative Manber-Myers Implementation (Iterative with Ranks)

```python
def build_suffix_array_manber_myers_iterative(text):
    """
    Iterative Manber-Myers suffix array construction.
    
    Uses rank arrays and pair comparisons for efficient sorting.
    
    Parameters
    ----------
    text : str
        Input string
    
    Returns
    -------
    list[int]
        Suffix array
    
    Time Complexity: O(n log n)
    Space Complexity: O(n)
    """
    n = len(text)
    if n == 0:
        return []
    if n == 1:
        return [0]
    
    # Initialize suffix array and rank array
    sa = list(range(n))
    rank = [ord(c) for c in text]
    tmp = [0] * n
    
    k = 1
    while k < n:
        # Sort by (rank[i], rank[i+k]) pairs
        def compare_key(i):
            return (rank[i], rank[i + k] if i + k < n else -1)
        
        sa.sort(key=compare_key)
        
        # Compute new ranks
        tmp[sa[0]] = 0
        for i in range(1, n):
            prev_key = compare_key(sa[i - 1])
            curr_key = compare_key(sa[i])
            tmp[sa[i]] = tmp[sa[i - 1]] + (1 if curr_key > prev_key else 0)
        
        rank, tmp = tmp, rank
        
        # Early termination if all ranks are unique
        if rank[sa[n - 1]] == n - 1:
            break
        
        k *= 2
    
    return sa

# Test
text = "banana$"
sa = build_suffix_array_manber_myers_iterative(text)
print(f"Text: {text}")
print(f"Suffix Array: {sa}")
```

## 4. Pattern Search using Suffix Array - O(m log n)

Once we have a suffix array, we can search for patterns using **binary search**.

### Key Insight
- If pattern P occurs in text T, it must be a **prefix** of some suffix
- Suffixes are sorted in the suffix array
- Binary search finds where P would fit in the sorted order

### ASCII Art: Pattern Search with Binary Search

```
Text = "banana$"
SA = [6, 5, 3, 1, 0, 4, 2]
Find: "ana"

Sorted suffixes (what SA represents):
  SA[0]=6 → "$"
  SA[1]=5 → "a$"
  SA[2]=3 → "ana$"      ← "ana" is prefix! ✓
  SA[3]=1 → "anana$"    ← "ana" is prefix! ✓
  SA[4]=0 → "banana$"
  SA[5]=4 → "na$"
  SA[6]=2 → "nana$"

Binary search for lower bound of "ana":
  mid = 3 → "anana$"
  Compare: "ana" vs "anana$"[:3] = "ana" → equal, but need lower bound
  Search left...
  
Binary search for upper bound of "ana":
  Find first suffix where "ana" is NOT a prefix
  
Result: All occurrences in range [2, 4) of SA
  → Positions SA[2]=3 and SA[3]=1
  → Pattern "ana" found at positions: 1, 3
```
