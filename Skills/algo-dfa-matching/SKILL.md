---
name: algo-dfa-matching
description: "This notebook covers **Deterministic Finite Automaton (DFA)** based string pattern matching - a technique that preprocesses the pattern to build a state machine enabling O(n) search time with no backt"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/04_dfa_matching.ipynb"
---

# DFA-Based Pattern Matching

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/07_String_Algorithms/04_dfa_matching.ipynb`*

# DFA-Based Pattern Matching

This notebook covers **Deterministic Finite Automaton (DFA)** based string pattern matching - a technique that preprocesses the pattern to build a state machine enabling O(n) search time with no backtracking.

## Table of Contents
1. [Finite Automaton Basics](#1-finite-automaton-basics)
2. [DFA for Pattern Matching](#2-dfa-for-pattern-matching)
3. [Building the Transition Table](#3-building-the-transition-table)
4. [Searching with DFA](#4-searching-with-dfa)
5. [Implementation](#5-implementation)
6. [Examples](#6-examples)
7. [DFA vs KMP Comparison](#7-dfa-vs-kmp-comparison)

---
## 1. Finite Automaton Basics

A **Deterministic Finite Automaton (DFA)** is formally defined as a 5-tuple:

$$\text{DFA} = (Q, \Sigma, \delta, q_0, F)$$

where:
- **Q** = finite set of states
- **Σ** (Sigma) = input alphabet (finite set of symbols)
- **δ** (delta) = transition function: Q × Σ → Q
- **q₀** = start state (q₀ ∈ Q)
- **F** = set of accepting (final) states (F ⊆ Q)

### Key Properties

1. **Deterministic**: For each state and input symbol, there is exactly ONE next state
2. **Complete**: Every state has a transition defined for every symbol in Σ
3. **No ε-transitions**: Only transitions on actual input symbols

### DFA Operation

```
Given input string w = a₁a₂...aₙ:

1. Start at state q₀
2. For each character aᵢ:
   current_state = δ(current_state, aᵢ)
3. Accept if final state ∈ F
```

---
## 2. DFA for Pattern Matching

For pattern matching, we construct a DFA where:

- **State i** means "we have matched the first i characters of the pattern"
- **State 0** = no characters matched (start state)
- **State m** = all m characters matched (accepting state)

### Formal Definition for Pattern P[0..m-1]

$$\text{DFA} = (Q, \Sigma, \delta, 0, \{m\})$$

- **Q** = {0, 1, 2, ..., m} — (m + 1) states
- **Σ** = alphabet of the text
- **q₀** = 0 (no match)
- **F** = {m} (complete match)

### State Meaning

```
Pattern: "ABABC" (m = 5)

State 0: No characters matched     → matched ""
State 1: Matched P[0]              → matched "A"
State 2: Matched P[0..1]           → matched "AB"
State 3: Matched P[0..2]           → matched "ABA"
State 4: Matched P[0..3]           → matched "ABAB"
State 5: Matched P[0..4]           → matched "ABABC" ✓ ACCEPTING
```

### ASCII Art: DFA for Pattern "ABABC"

```
States: 0, 1, 2, 3, 4, 5 (5 = accepting state)

      A       B       A       B       C
→(0)────→(1)────→(2)────→(3)────→(4)────→((5))
   ↑           │       │
   │     B     │   A   │
   └───────────┘       │
         │             │
         └─────────────┘
              B

Legend:
  →(0)   = Start state
  ((5))  = Accepting (double circle)
  ────→  = Transition on matching character
```

### Complete Transition Table

```
┌─────────────────────────────────────────────┐
│ State │  A  │  B  │  C  │ (other)           │
├───────┼─────┼─────┼─────┼───────────────────┤
│   0   │  1  │  0  │  0  │    0              │
│   1   │  1  │  2  │  0  │    0              │
│   2   │  3  │  0  │  0  │    0              │
│   3   │  1  │  4  │  0  │    0              │
│   4   │  3  │  0  │  5  │    0              │
│   5   │  -  │  -  │  -  │    - (accepting)  │
└───────┴─────┴─────┴─────┴───────────────────┘

Reading the table:
- Row = current state
- Column = input character
- Cell value = next state
```

---
## 3. Building the Transition Table

The key insight is computing δ(state, char) for ALL combinations:

### Transition Function Logic

For state `q` and character `c`:

```
δ(q, c) = length of LONGEST prefix of P that is also
          a suffix of (P[0..q-1] + c)
```

In other words:
- We've matched P[0..q-1]
- We see character c
- What's the longest prefix of P that ends with this new state?

### Two Cases

**Case 1: Character matches next pattern character**
```
If c == P[q], then δ(q, c) = q + 1
(We extend the match by one character)
```

**Case 2: Character doesn't match (mismatch)**
```
Find the longest proper prefix of P[0..q-1] + c
that is also a suffix of the pattern.

This is where we use the PREFIX FUNCTION!
```

### Building Transitions Using Prefix Function

```
Pattern: "ABABC"
Prefix function π: [0, 0, 1, 2, 0]

═══════════════════════════════════════════════════════
Example 1: δ(2, 'A')  — state=2 (matched "AB"), char='A'
═══════════════════════════════════════════════════════

  Does 'A' match P[2]? P[2]='A' → YES!
  δ(2, 'A') = 2 + 1 = 3

  Meaning: "AB" + "A" = "ABA" = P[0..2] ✓

═══════════════════════════════════════════════════════
Example 2: δ(2, 'C')  — state=2 (matched "AB"), char='C'
═══════════════════════════════════════════════════════

  Does 'C' match P[2]? P[2]='A' → NO
  
  We need to find: longest prefix of P that is suffix of "ABC"
  
  Method: Use prefix function recursively
  - We had matched "AB" (state 2)
  - Fall back to π[1] = 0 (longest proper prefix-suffix of "A")
  - From state 0, check: does 'C' match P[0]? P[0]='A' → NO
  - δ(2, 'C') = 0

═══════════════════════════════════════════════════════
Example 3: δ(4, 'A')  — state=4 (matched "ABAB"), char='A'
═══════════════════════════════════════════════════════

  Does 'A' match P[4]? P[4]='C' → NO
  
  Fall back using prefix function:
  - π[3] = 2 ("ABAB" has prefix-suffix "AB" of length 2)
  - From state 2, check: does 'A' match P[2]? P[2]='A' → YES!
  - δ(4, 'A') = 2 + 1 = 3

  Meaning: "ABAB" + "A" → keep suffix "ABA" = P[0..2]
```

### Naive vs Optimized Construction

**Naive Approach: O(m³|Σ|)**
```python
for state in range(m + 1):
    for char in alphabet:
        # Try all possible prefix lengths
        for k in range(min(m, state + 1), 0, -1):
            if pattern[:k] == (pattern[:state] + char)[-k:]:
                delta[state][char] = k
                break
```

**Optimized Approach: O(m|Σ|)** using prefix function
```python
for state in range(m + 1):
    for char in alphabet:
        # Use prefix function for O(1) amortized lookup
        k = state
        while k > 0 and pattern[k] != char:
            k = prefix[k - 1]
        if pattern[k] == char:
            k += 1
        delta[state][char] = k
```

---
## 4. Searching with DFA

Once the DFA is built, searching is simple and fast:

```
DFA_SEARCH(text, dfa):
    state = 0
    for i = 0 to len(text) - 1:
        state = dfa[state][text[i]]
        if state == m:           # Accepting state
            report match at position (i - m + 1)
```

### Search Example

```
Pattern: "ABABC" (m = 5)
Text:    "ABABABC"

Step-by-step state transitions:

Position:  0    1    2    3    4    5    6
Character: A    B    A    B    A    B    C
           ↓    ↓    ↓    ↓    ↓    ↓    ↓
State:  0→ 1 → 2 → 3 → 4 → 3 → 4 → 5
                          │         ↑
                          │    FOUND! ✓
                          │
                    Mismatch at state 4:
                    Expected 'C', got 'A'
                    δ(4,'A') = 3 (fall back)

Match found at position: 6 - 5 + 1 = 2

Verification: text[2:7] = "ABABC" ✓
```

### Finding All Occurrences

For overlapping matches, continue searching after finding a match:

```
Pattern: "ABA"
Text:    "ABABABA"

Char:  A  B  A  B  A  B  A
       ↓  ↓  ↓  ↓  ↓  ↓  ↓
State: 0→1→2→3→2→3→2→3
             ↑     ↑     ↑
          Match  Match  Match
          pos=0  pos=2  pos=4

Note: After reaching accepting state 3,
      we transition using dfa[3][next_char]
      to continue finding overlapping matches.
```

---
## 5. Implementation

```python
def compute_prefix_function(pattern: str) -> list:
    """
    Compute the prefix function (failure function) for the pattern.
    
    The prefix function π[i] is the length of the longest proper prefix
    of pattern[0..i] that is also a suffix of pattern[0..i].
    
    Parameters
    ----------
    pattern : str
        The pattern to compute the prefix function for.
    
    Returns
    -------
    list
        Prefix function array of length m.
    
    Time Complexity: O(m) where m = len(pattern)
    Space Complexity: O(m)
    
    Example
    -------
    >>> compute_prefix_function("ABABC")
    [0, 0, 1, 2, 0]
    """
    m = len(pattern)
    if m == 0:
        return []
    
    prefix = [0] * m
    k = 0  # length of previous longest prefix-suffix
    
    for i in range(1, m):
        # Fall back using prefix function until we find a match or reach start
        while k > 0 and pattern[i] != pattern[k]:
            k = prefix[k - 1]
        
        if pattern[i] == pattern[k]:
            k += 1
        
        prefix[i] = k
    
    return prefix


# Test the prefix function
print("Prefix function examples:")
for pattern in ["ABABC", "AAAA", "ABCABC", "ATTCTGATTT"]:
    print(f"  π('{pattern}') = {compute_prefix_function(pattern)}")
```

```python
def build_dfa(pattern: str, alphabet: str) -> dict:
    """
    Build a DFA transition table for pattern matching.
    
    Constructs a deterministic finite automaton where state i represents
    having matched the first i characters of the pattern.
    
    Parameters
    ----------
    pattern : str
        The pattern to search for.
    alphabet : str
        String containing all possible characters in the text.
    
    Returns
    -------
    dict
        Nested dictionary: dfa[state][char] -> next_state
        States range from 0 to m (inclusive).
    
    Time Complexity: O(m × |Σ|) where m = len(pattern), |Σ| = len(alphabet)
    Space Complexity: O(m × |Σ|)
    
    Example
    -------
    >>> dfa = build_dfa("AB", "ABC")
    >>> dfa[0]['A']  # From state 0, reading 'A'
    1
    >>> dfa[1]['B']  # From state 1, reading 'B'
    2  # Accepting state
    """
    m = len(pattern)
    
    # Handle empty pattern
    if m == 0:
        return {0: {c: 0 for c in alphabet}}
    
    # Compute prefix function for the pattern
    prefix = compute_prefix_function(pattern)
    
    # Build transition table for states 0 to m
    dfa = {}
    
    for state in range(m + 1):
        dfa[state] = {}
        
        for char in alphabet:
            if state < m and char == pattern[state]:
                # Character matches: advance to next state
                dfa[state][char] = state + 1
            else:
                # Mismatch: use prefix function to find fallback state
                # Find longest prefix of pattern that is suffix of (pattern[0..state-1] + char)
                k = state
                while k > 0:
                    k = prefix[k - 1]
                    if k < m and char == pattern[k]:
                        k += 1
                        break
                    if k == 0:
                        if char == pattern[0]:
                            k = 1
                        break
                dfa[state][char] = k
    
    return dfa


# Test DFA construction
test_dfa = build_dfa("AB", "ABC")
print("DFA for pattern 'AB':")
for state in test_dfa:
    print(f"  State {state}: {test_dfa[state]}")
```

```python
def build_dfa_via_suffix_check(pattern: str, alphabet: str) -> dict:
    """
    Build DFA using direct suffix matching (alternative implementation).
    
    For each state and character, computes the longest prefix of pattern
    that matches a suffix of (pattern[0:state] + char).
    
    Parameters
    ----------
    pattern : str
        The pattern to search for.
    alphabet : str
        String containing all possible characters in the text.
    
    Returns
    -------
    dict
        Nested dictionary: dfa[state][char] -> next_state
    
    Time Complexity: O(m² × |Σ|) - less efficient but more intuitive
    Space Complexity: O(m × |Σ|)
    """
    m = len(pattern)
    dfa = {}
    
    for state in range(m + 1):
        dfa[state] = {}
        current_matched = pattern[:state]  # What we've matched so far
        
        for char in alphabet:
            # After reading char, what's the longest prefix of pattern
            # that is a suffix of (current_matched + char)?
            test_string = current_matched + char
            
            # Try all possible prefix lengths, from longest to shortest
            for k in range(min(m, len(test_string)), -1, -1):
                if test_string.endswith(pattern[:k]):
                    dfa[state][char] = k
                    break
    
    return dfa


# Verify both implementations give same results
pattern = "ABABC"
alphabet = "ABC"

dfa1 = build_dfa(pattern, alphabet)
dfa2 = build_dfa_via_suffix_check(pattern, alphabet)

print(f"Both implementations match: {dfa1 == dfa2}")
```
