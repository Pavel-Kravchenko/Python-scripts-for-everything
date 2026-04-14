---
name: algo-tabulation
description: "Bottom-up DP with tabulation — edit distance, LCS, space optimization with rolling arrays"
tool_type: python
primary_tool: Python
---

# Tabulation (Bottom-Up DP)

## Memoization vs Tabulation

| Aspect | Memoization (Top-Down) | Tabulation (Bottom-Up) |
|--------|----------------------|----------------------|
| Direction | Recurse down, cache | Iterate up, fill table |
| Stack overflow | Risk with deep recursion | No risk |
| Subproblems | Only needed ones | All subproblems |
| Space optimization | Harder | Easy (rolling arrays) |

## Fibonacci — O(1) Space

```python
def fib_optimized(n):
    if n <= 1: return n
    prev2, prev1 = 0, 1
    for _ in range(2, n + 1):
        prev2, prev1 = prev1, prev2 + prev1
    return prev1
```

## Edit Distance (Levenshtein) — O(mn) time, O(mn) space

```python
def edit_distance(s1, s2):
    """dp[i][j] = min edits to transform s1[:i] into s2[:j]."""
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1): dp[i][0] = i
    for j in range(n + 1): dp[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
    return dp[m][n]
```

## Longest Common Subsequence (LCS)

```python
def lcs_with_traceback(s1, s2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    # Traceback
    result = []
    i, j = m, n
    while i > 0 and j > 0:
        if s1[i-1] == s2[j-1]:
            result.append(s1[i-1]); i -= 1; j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            i -= 1
        else:
            j -= 1
    return ''.join(reversed(result))
```

LCS similarity score: `lcs_length / max(len(s1), len(s2))` — quick conservation metric.

## Space Optimization — O(min(m,n))

Only the previous row is needed for the current computation:

```python
def lcs_space_optimized(s1, s2):
    if len(s2) > len(s1): s1, s2 = s2, s1
    m, n = len(s1), len(s2)
    prev = [0] * (n + 1)
    curr = [0] * (n + 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                curr[j] = prev[j-1] + 1
            else:
                curr[j] = max(prev[j], curr[j-1])
        prev, curr = curr, [0] * (n + 1)
    return prev[n]
```

## Pitfalls

- Space-optimized DP loses traceback — you can only recover the optimal value, not the actual alignment/subsequence
- Edit distance traceback: walk backwards from dp[m][n] checking which predecessor was used
- LCS and edit distance are related: `edit_distance >= len(s1) + len(s2) - 2 * lcs_length`
