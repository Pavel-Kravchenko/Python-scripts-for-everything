---
name: algo-rabin-karp
description: "A hash-based string matching algorithm that uses rolling hash for efficient pattern searching."
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/03_rabin_karp.ipynb"
---

# Rabin-Karp Algorithm

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/03_rabin_karp.ipynb`*

# Rabin-Karp Algorithm

A hash-based string matching algorithm that uses rolling hash for efficient pattern searching.

## Key Insight

**If `hash(window) ≠ hash(pattern)`, there is definitely no match.**

This allows us to skip expensive character-by-character comparisons for most positions. We only need to verify when hashes match (which might be a collision).

## Topics Covered

1. **Hash-based Pattern Matching** - compare hashes instead of character strings
2. **Rolling Hash** - update hash in O(1) when sliding window
3. **Handling Collisions** - hash match doesn't guarantee string match
4. **Multiple Pattern Search** - where Rabin-Karp excels

---

## 1. Hash-based Pattern Matching

### The Basic Idea

Instead of comparing strings character by character, we compare their hash values:

```
Naive approach:     Compare m characters → O(m) per position
Hash-based:         Compare 2 integers  → O(1) per position
```

### Why It Works

- **Different strings usually have different hashes** → Quick rejection
- **Same hash might mean collision** → Need verification
- **Net effect**: Dramatically fewer comparisons on average

---

## 2. Rolling Hash Concept

The magic of Rabin-Karp is the **rolling hash** - updating the hash in O(1) as we slide the window.

```
Text: A B C D E F G
      └───┴───┘
      window = "ABCD"
      hash = h("ABCD")
      
Slide window:
Text: A B C D E F G
        └───┴───┘
        window = "BCDE"
        
Instead of recomputing h("BCDE") from scratch:
  new_hash = roll(old_hash, remove='A', add='E')

Rolling: O(1) instead of O(m)!
```

### Rolling Hash Formula

We treat a string as a number in base `d` (typically 256 for ASCII):

```
h(s[i..i+m-1]) = s[i]*d^(m-1) + s[i+1]*d^(m-2) + ... + s[i+m-1]*d^0
```

To roll from position `i` to `i+1`:

```
h(s[i+1..i+m]) = d * (h(s[i..i+m-1]) - s[i]*d^(m-1)) + s[i+m]
```

**In plain English:**
1. Subtract contribution of outgoing character (scaled by d^(m-1))
2. Multiply by d (shift all positions left)
3. Add incoming character

```
Base d = 256 (ASCII), prime q = 101

Pattern: "ABC"  →  h = (65*256² + 66*256 + 67) mod 101

Text: "XYZABC"
Window "XYZ": h = (88*256² + 89*256 + 90) mod 101

Roll to "YZA":
  1. Remove 'X': (h - 88*256²) mod 101
  2. Shift left: result * 256
  3. Add 'A': result + 65
  4. Mod: result mod 101

Formula: h_new = (d * (h_old - char_out * d^(m-1)) + char_in) mod q
```

### Why Use Modular Arithmetic?

Without modulo, hash values would grow exponentially:

```
Pattern length m = 10, base d = 256
Hash could be as large as: 256^10 = 1,208,925,819,614,629,174,706,176
```

**Problems without modulo:**
- Integer overflow
- Arithmetic becomes expensive (arbitrary precision)
- Comparison becomes O(m) for big integers

**Benefits of modular arithmetic:**
- Hash values stay bounded (0 to q-1)
- All arithmetic is O(1)
- Choose prime q to minimize collisions

**Why prime q?**
- Distributes hash values more uniformly
- Reduces patterns in the hash function
- Common choice: large prime like 101, 1000000007

---

## 3. Handling Collisions (Spurious Hits)

A **hash collision** occurs when different strings have the same hash.

```
Text:    "AABAACAADAAB"
Pattern: "AAB"
Pattern hash = 65*256² + 65*256 + 66 = 4276546 mod 101 = 80

Position 0: "AAB" → hash = 80 → MATCH! Verify: "AAB" = "AAB" ✓
Position 1: "ABA" → hash = 82 → no match
Position 2: "BAA" → hash = 95 → no match
Position 3: "AAC" → hash = 81 → no match
Position 4: "ACA" → hash = 72 → no match
...

Hash collision example:
Position 5: "CAA" → hash = 80 → HASH MATCH! Verify: "CAA" ≠ "AAB" ✗
                                 (spurious hit - false positive)
```

### The Rule

**Always verify on hash match!**

```python
if hash_window == hash_pattern:
    if text[i:i+m] == pattern:  # Actual comparison
        found_match(i)
```

This verification costs O(m) but happens rarely with a good hash function.

---

## 4. Search Process Visualization

```
Text:    [A][B][R][A][C][A][D][A][B][R][A]
Pattern: [A][B][R][A]

Step 1: Compute pattern hash
        h("ABRA") = (65·256³ + 66·256² + 82·256 + 65) mod q

Step 2: Compute first window hash
        h(text[0:4]) = h("ABRA") 

Step 3: Compare and slide

  pos=0: [A][B][R][A] C  A  D  A  B  R  A
         └────────┘
         hash = h_p → MATCH! Verify ✓ → Found at 0

  pos=1:  A [B][R][A][C] A  D  A  B  R  A
            └────────┘
         hash = roll(h, -'A', +'C') → ≠ h_p

  pos=2:  A  B [R][A][C][A] D  A  B  R  A
               └────────┘
         hash = roll(...) → ≠ h_p

  ...continue rolling...

  pos=7:  A  B  R  A  C  A  D [A][B][R][A]
                              └────────┘
         hash = h_p → MATCH! Verify ✓ → Found at 7
```

---

## Complexity Analysis

| Case | Time Complexity | Notes |
|------|-----------------|-------|
| **Average/Expected** | O(n + m) | Few hash collisions |
| **Worst case** (many collisions) | O(n × m) | Every position needs verification |
| **Space** | O(1) | Just storing hash values |

### When Does Worst Case Occur?

```
Text:    "AAAAAAAAAA"
Pattern: "AAA"

Every window has the same hash → verify at every position!
```

**With a good hash function and prime modulus, worst case is rare in practice.**

---

## Implementation

```python
from typing import Optional, List

def rabin_karp_search(text: str, pattern: str, d: int = 256, q: int = 101) -> list[int]:
    """
    Find all occurrences of pattern in text using Rabin-Karp algorithm.
    
    Uses rolling hash to achieve O(n+m) average time complexity.
    
    Parameters
    ----------
    text : str
        The text to search in.
    pattern : str
        The pattern to search for.
    d : int, optional
        Base for the hash function (default 256 for ASCII).
    q : int, optional
        Prime modulus to prevent overflow (default 101).
        Larger primes reduce collision probability.
    
    Returns
    -------
    list[int]
        List of starting indices where pattern occurs in text.
    
    Examples
    --------
    >>> rabin_karp_search("AABAACAADAABAABA", "AABA")
    [0, 9, 12]
    >>> rabin_karp_search("abcdef", "xyz")
    []
    
    Notes
    -----
    Rolling hash formula:
        h(s[i+1..i+m]) = (d * (h(s[i..i+m-1]) - s[i] * d^(m-1)) + s[i+m]) mod q
    
    The algorithm verifies actual string equality on hash matches to handle
    collisions (spurious hits).
    """
    n = len(text)
    m = len(pattern)
    matches = []
    
    # Edge cases
    if m == 0 or m > n:
        return matches
    
    # Precompute d^(m-1) mod q for rolling hash
    # This is the multiplier for the leftmost character
    h = 1
    for _ in range(m - 1):
        h = (h * d) % q
    
    # Compute initial hash values
    pattern_hash = 0
    window_hash = 0
    
    for i in range(m):
        pattern_hash = (d * pattern_hash + ord(pattern[i])) % q
        window_hash = (d * window_hash + ord(text[i])) % q
    
    # Slide pattern over text
    for i in range(n - m + 1):
        # Check if hashes match
        if pattern_hash == window_hash:
            # Verify actual string match (handle collisions)
            if text[i:i+m] == pattern:
                matches.append(i)
        
        # Compute hash for next window using rolling hash
        if i < n - m:
            # Remove leading character, add trailing character
            window_hash = (d * (window_hash - ord(text[i]) * h) + ord(text[i + m])) % q
            
            # Handle negative modulo (Python handles this, but be explicit)
            if window_hash < 0:
                window_hash += q
    
    return matches
```

```python
# Test the basic implementation
print("Basic tests:")
print(f"  'AABAACAADAABAABA' contains 'AABA' at: {rabin_karp_search('AABAACAADAABAABA', 'AABA')}")
print(f"  'abcdef' contains 'xyz' at: {rabin_karp_search('abcdef', 'xyz')}")
print(f"  'ABRACADABRA' contains 'ABRA' at: {rabin_karp_search('ABRACADABRA', 'ABRA')}")
```

---

## Step-by-Step Hash Computation Example

```python
def demonstrate_hash_computation(text: str, pattern: str, d: int = 256, q: int = 101):
    """
    Demonstrate step-by-step hash computation in Rabin-Karp.
    """
    n, m = len(text), len(pattern)
    
    print(f"Text:    '{text}'")
    print(f"Pattern: '{pattern}' (length m={m})")
    print(f"Base d={d}, Prime q={q}")
    print("=" * 60)
    
    # Compute d^(m-1) mod q
    h = pow(d, m-1, q)
    print(f"\nd^(m-1) mod q = {d}^{m-1} mod {q} = {h}")
    
    # Compute pattern hash step by step
    print(f"\n--- Pattern Hash Computation ---")
    pattern_hash = 0
    for i, char in enumerate(pattern):
        old_hash = pattern_hash
        pattern_hash = (d * pattern_hash + ord(char)) % q
        print(f"  char '{char}' (ASCII {ord(char)}): ({d} × {old_hash} + {ord(char)}) mod {q} = {pattern_hash}")
    print(f"  Pattern hash = {pattern_hash}")
    
    # Compute first window hash
    print(f"\n--- Initial Window Hash ---")
    window_hash = 0
    for i in range(m):
        old_hash = window_hash
        window_hash = (d * window_hash + ord(text[i])) % q
        print(f"  char '{text[i]}' (ASCII {ord(text[i])}): ({d} × {old_hash} + {ord(text[i])}) mod {q} = {window_hash}")
    print(f"  Window[0:{m}] = '{text[0:m]}', hash = {window_hash}")
    
    # Rolling hash demonstration
    print(f"\n--- Rolling Hash in Action ---")
    for i in range(min(5, n - m)):  # Show first 5 rolls
        old_hash = window_hash
        char_out = text[i]
        char_in = text[i + m]
        
        # Rolling hash formula
        window_hash = (d * (window_hash - ord(char_out) * h) + ord(char_in)) % q
        if window_hash < 0:
            window_hash += q
        
        match_status = "HASH MATCH!" if window_hash == pattern_hash else ""
        print(f"  Roll: remove '{char_out}', add '{char_in}'")
        print(f"        Window[{i+1}:{i+1+m}] = '{text[i+1:i+1+m]}', hash = {window_hash} {match_status}")
        
        if match_status:
            actual_match = text[i+1:i+1+m] == pattern
            print(f"        Verify: '{text[i+1:i+1+m]}' == '{pattern}' → {actual_match}")

# Demonstrate with a simple example
demonstrate_hash_computation("ABRACADABRA", "ABR")
```

---

## Demonstrating Spurious Hits (Collisions)

```python
def find_spurious_hits(text: str, pattern: str, d: int = 256, q: int = 11):
    """
    Find and display spurious hits (hash collisions that aren't actual matches).
    Uses small prime q to increase collision probability for demonstration.
    """
    n, m = len(text), len(pattern)
    h = pow(d, m-1, q)
    
    # Compute pattern hash
    pattern_hash = 0
    for char in pattern:
        pattern_hash = (d * pattern_hash + ord(char)) % q
    
    # Compute first window hash
    window_hash = 0
    for i in range(m):
        window_hash = (d * window_hash + ord(text[i])) % q
    
    print(f"Text: '{text}'")
    print(f"Pattern: '{pattern}' (hash = {pattern_hash} with small prime q={q})")
    print("=" * 60)
    
    true_matches = []
    spurious_hits = []
    
    for i in range(n - m + 1):
        window = text[i:i+m]
        
        if window_hash == pattern_hash:
            if window == pattern:
                true_matches.append((i, window))
                print(f"  pos {i}: '{window}' hash={window_hash} → TRUE MATCH ✓")
            else:
                spurious_hits.append((i, window))
                print(f"  pos {i}: '{window}' hash={window_hash} → SPURIOUS HIT ✗ (collision!)")
        
        # Roll hash
        if i < n - m:
            window_hash = (d * (window_hash - ord(text[i]) * h) + ord(text[i + m])) % q
            if window_hash < 0:
                window_hash += q
    
    print("\nSummary:")
    print(f"  True matches: {len(true_matches)}")
    print(f"  Spurious hits: {len(spurious_hits)}")
    return true_matches, spurious_hits

# Use small prime to force collisions
find_spurious_hits("ABCDEFGHIJKLMNOPQRSTUVWXYZABC", "ABC", q=11)
```

```python
# Another collision example
print("\nWith very small prime q=7:")
find_spurious_hits("AABAACAADAABAABCAAB", "AAB", q=7)
```

---

## Multiple Pattern Search

**This is where Rabin-Karp really shines!**

For searching k patterns of the same length:
- Naive: O(k × n × m)
- Rabin-Karp: O(n + k × m) average

We compute one rolling hash over the text, and check against all pattern hashes.
