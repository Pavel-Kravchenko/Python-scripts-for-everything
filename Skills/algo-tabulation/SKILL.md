---
name: algo-tabulation
description: "1. Understand bottom-up DP (tabulation) vs top-down (memoization) 2. Build DP tables iteratively 3. Optimize space complexity 4. Apply to sequence comparison problems"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/02_tabulation.ipynb"
---

# 📊 Tabulation: Bottom-Up Dynamic Programming

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/02_tabulation.ipynb`*

# 📊 Tabulation: Bottom-Up Dynamic Programming

## Learning Objectives

1. Understand bottom-up DP (tabulation) vs top-down (memoization)
2. Build DP tables iteratively
3. Optimize space complexity
4. Apply to sequence comparison problems

---

## 1. Memoization vs Tabulation

| Aspect | Memoization (Top-Down) | Tabulation (Bottom-Up) |
|--------|----------------------|----------------------|
| Direction | Start big, recurse down | Start small, build up |
| Implementation | Recursion + cache | Iteration + table |
| Subproblems solved | Only needed ones | All subproblems |
| Stack overflow risk | Yes (deep recursion) | No |
| Space optimization | Harder | Easier |

```python
# Fibonacci: Tabulation approach

def fib_tabulation(n):
    """Bottom-up Fibonacci with O(n) space."""
    if n <= 1:
        return n
    
    dp = [0] * (n + 1)
    dp[0], dp[1] = 0, 1
    
    for i in range(2, n + 1):
        dp[i] = dp[i-1] + dp[i-2]
    
    return dp[n]

def fib_optimized(n):
    """Space-optimized Fibonacci with O(1) space."""
    if n <= 1:
        return n
    
    prev2, prev1 = 0, 1
    for _ in range(2, n + 1):
        prev2, prev1 = prev1, prev2 + prev1
    
    return prev1

print(f"F(50) = {fib_optimized(50)}")
```

## 2. 🧬 Edit Distance (Levenshtein Distance)

Classic DP problem directly relevant to sequence alignment. Find minimum operations (insert, delete, substitute) to transform one string into another.

```python
def edit_distance(s1, s2):
    """
    Compute edit distance using tabulation.
    dp[i][j] = edit distance between s1[:i] and s2[:j]
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Base cases: transforming to/from empty string
    for i in range(m + 1):
        dp[i][0] = i  # Delete all characters
    for j in range(n + 1):
        dp[0][j] = j  # Insert all characters
    
    # Fill table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1]  # No operation needed
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],    # Delete from s1
                    dp[i][j-1],    # Insert into s1
                    dp[i-1][j-1]   # Substitute
                )
    
    return dp[m][n]

# Compare DNA sequences
seq1 = "ATCGATCG"
seq2 = "ATGGATGG"
print(f"Edit distance between {seq1} and {seq2}: {edit_distance(seq1, seq2)}")
```

## Visualizing the Edit Distance Table

Printing the DP table with row and column labels makes it easy to trace how each cell was computed and to see the minimum-cost transformation path from the top-left to the bottom-right corner.

```python
def edit_distance_verbose(s1, s2):
    """Edit distance with table visualization."""
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

    # Print with labels
    header = "     " + "  -  " + "".join(f"  {c}  " for c in s2)
    print(header)
    for i in range(m + 1):
        label = "-" if i == 0 else s1[i-1]
        row = "".join(f"{dp[i][j]:>5}" for j in range(n + 1))
        print(f"  {label} {row}")
    return dp[m][n]

edit_distance_verbose("ATCG", "ATGG")
```

## 3. Longest Common Subsequence (LCS)

Foundation of sequence alignment algorithms.

```python
def lcs(s1, s2):
    """Find length of longest common subsequence."""
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    return dp[m][n]

def lcs_with_traceback(s1, s2):
    """Find LCS with actual subsequence."""
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
            result.append(s1[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            i -= 1
        else:
            j -= 1
    
    return ''.join(reversed(result))

seq1 = "ATGCAATGC"
seq2 = "ATCCAATCC"
print(f"LCS of {seq1} and {seq2}: {lcs_with_traceback(seq1, seq2)}")
print(f"LCS length: {lcs(seq1, seq2)}")
```

## Biological Significance of LCS

LCS length relative to sequence length gives a quick measure of conservation between two sequences.

A simple similarity score: `LCS_length / max(len(s1), len(s2))`

- Score near 1.0: highly conserved sequences (e.g., orthologous genes across closely related species)
- Score near 0.5: moderate conservation (e.g., paralogous genes after divergence)
- Score near 0.0: little shared structure

This is simpler and faster than full pairwise alignment (Needleman-Wunsch / Smith-Waterman) and useful for quick filtering — e.g., pre-screening candidate sequences before running a computationally expensive aligner.

## 4. Space Optimization

LCS uses O(mn) space, but we only need the previous row → O(min(m,n)) space.

```python
def lcs_space_optimized(s1, s2):
    """LCS with O(min(m,n)) space."""
    if len(s2) > len(s1):
        s1, s2 = s2, s1
    
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

print(f"LCS (space-optimized): {lcs_space_optimized(seq1, seq2)}")
```

## 🧬 Exercise: Longest Increasing Subsequence

Implement LIS using tabulation. Given expression values over time, find longest increasing trend.

```python
def lis_tabulation(arr):
    """Find length of LIS using O(n²) DP."""
    n = len(arr)
    if n == 0:
        return 0
    
    # dp[i] = length of LIS ending at index i
    dp = [1] * n
    
    for i in range(1, n):
        for j in range(i):
            if arr[j] < arr[i]:
                dp[i] = max(dp[i], dp[j] + 1)
    
    return max(dp)

expression = [2.1, 3.5, 1.2, 4.8, 3.9, 5.2, 6.1, 4.5, 7.0]
print(f"Gene expression: {expression}")
print(f"Longest increasing trend: {lis_tabulation(expression)} time points")
```

### Exercise 2 (**)

**Edit Distance Operations**

Extend `edit_distance` to also return the sequence of edit operations (Insert, Delete, Substitute) via traceback through the DP table. Test on two short DNA sequences.

Hint: after filling the DP table, walk back from `dp[m][n]` to `dp[0][0]` and record which operation was used at each step.

```python
def edit_distance_with_ops(s1, s2):
    """
    Edit distance with traceback returning the list of operations.
    Operations: ('Match', c), ('Substitute', c1, c2), ('Delete', c), ('Insert', c)
    """
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])

    # Traceback
    ops = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and s1[i-1] == s2[j-1]:
            ops.append(('Match', s1[i-1]))
            i -= 1
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + 1:
            ops.append(('Substitute', s1[i-1], s2[j-1]))
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + 1:
            ops.append(('Delete', s1[i-1]))
            i -= 1
        else:
            ops.append(('Insert', s2[j-1]))
            j -= 1

    ops.reverse()
    return dp[m][n], ops


# YOUR CODE HERE
# Test on two short DNA sequences
distance, operations = edit_distance_with_ops("ATCG", "TACG")
print(f"Edit distance: {distance}")
print("Operations:")
for op in operations:
    print(f"  {op}")
```

---

## Summary

✅ Tabulation builds solutions bottom-up iteratively  
✅ Avoids recursion stack overflow  
✅ Enables space optimization (rolling arrays)  
✅ Edit distance and LCS are foundations of sequence alignment  
✅ LCS similarity score: a quick conservation metric without full alignment  
✅ Edit distance traceback recovers the exact sequence of operations

**Next:** [03_knapsack.ipynb](03_knapsack.ipynb) - Resource optimization with knapsack DP
