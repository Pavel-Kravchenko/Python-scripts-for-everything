---
name: algo-dfa-matching
description: "DFA-based exact pattern matching: O(m|Σ|) build via prefix function, O(n) search with no backtracking."
tool_type: python
primary_tool: Python
---

# DFA-Based Pattern Matching

## How It Works

Build a DFA with m+1 states where **state i** = "matched first i chars of pattern".
Accepting state = m. Transition function: `δ(state, char)` = longest prefix of pattern
that is a suffix of `pattern[:state] + char`.

## vs KMP

| | DFA | KMP |
|--|-----|-----|
| Build | O(m \| Σ \|) | O(m) |
| Search | O(n) | O(n) |
| Space | O(m \| Σ \|) | O(m) |
| Backtracking | None (table lookup) | Via failure links |
| Multi-text reuse | Ideal (precomputed) | Rebuild failure links per pattern |

Use DFA when you search one pattern against many texts. Use KMP to save memory.

## Build Transition Table

```python
def compute_prefix_function(pattern: str) -> list[int]:
    m = len(pattern)
    pi = [0] * m
    k = 0
    for i in range(1, m):
        while k > 0 and pattern[i] != pattern[k]:
            k = pi[k - 1]
        if pattern[i] == pattern[k]:
            k += 1
        pi[i] = k
    return pi


def build_dfa(pattern: str, alphabet: str) -> dict[int, dict[str, int]]:
    """Build DFA transition table. O(m * |alphabet|)."""
    m = len(pattern)
    pi = compute_prefix_function(pattern)
    dfa: dict[int, dict[str, int]] = {}

    for state in range(m + 1):
        dfa[state] = {}
        for ch in alphabet:
            if state < m and ch == pattern[state]:
                dfa[state][ch] = state + 1
            else:
                # fall back using prefix function
                k = state
                while k > 0:
                    k = pi[k - 1]
                    if k < m and ch == pattern[k]:
                        k += 1
                        break
                    if k == 0:
                        k = 1 if ch == pattern[0] else 0
                        break
                dfa[state][ch] = k

    return dfa
```

## Search

```python
def dfa_search(text: str, pattern: str, alphabet: str) -> list[int]:
    """Return list of start positions (0-based) of all occurrences."""
    if not pattern:
        return []
    m = len(pattern)
    dfa = build_dfa(pattern, alphabet)
    state = 0
    results = []
    for i, ch in enumerate(text):
        state = dfa[state].get(ch, 0)   # unknown chars go to state 0
        if state == m:
            results.append(i - m + 1)
    return results
```

## Example

```
Pattern: "ABABC"   Text: "ABABABC"

State: 0→1→2→3→4→3→4→5
                  ↑
          δ(4,'A')=3 (fell back; "ABAB"+"A" → suffix "ABA" = pattern[:3])

Match at position 2: text[2:7] = "ABABC"
```

## Pitfalls

- **Alphabet must cover all text characters**: characters not in `alphabet` get no entry in `dfa[state]`; handle via `.get(ch, 0)` or enumerate the actual text to build the alphabet.
- **Build is O(m × |Σ|), not O(m)**: for large alphabets (Unicode), this is expensive — prefer KMP or Aho-Corasick.
- **State m is accepting but search continues**: after reaching state m, the DFA naturally transitions to the correct fallback state to find overlapping matches — do not reset to 0.
- **Empty pattern edge case**: `m=0` means every position matches; handle before building the DFA.
- **Fallback loop in build**: the while loop mirrors KMP's failure-link traversal; an off-by-one in the `pi` index (`pi[k-1]` vs `pi[k]`) silently produces wrong transitions.
