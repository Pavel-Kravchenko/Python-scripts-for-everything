# Dynamic Programming

This module covers dynamic programming techniques for solving optimization problems by breaking them into overlapping subproblems.

## Core Concepts

### What is Dynamic Programming?
Dynamic Programming (DP) is an algorithmic technique that solves complex problems by:
1. Breaking them into simpler **overlapping subproblems**
2. Storing solutions to avoid **redundant computation**
3. Building up solutions to larger problems from smaller ones

### Two Approaches

| Approach | Description | Direction | When to Use |
|----------|-------------|-----------|-------------|
| **Memoization** (Top-Down) | Recursive with caching | Start from problem, recurse down | Natural recursive structure |
| **Tabulation** (Bottom-Up) | Iterative with table | Start from base cases, build up | Clear subproblem ordering |

## Topics Covered

### 1. Foundational Problems
- **Fibonacci Sequence** - Classic DP introduction
- **Climbing Stairs** - Counting paths
- **Coin Change** - Minimum coins needed

### 2. Sequence Problems
- **Longest Common Subsequence (LCS)** - Two sequences
- **Longest Increasing Subsequence (LIS)** - Single sequence
- **Edit Distance** - String transformation

### 3. Knapsack Problems
- **0/1 Knapsack** - Each item used once
- **Unbounded Knapsack** - Items can be reused
- **Subset Sum** - Can we make a target sum?

### 4. Advanced Problems
- **Matrix Chain Multiplication** - Optimal parenthesization
- **Longest Palindromic Subsequence** - String palindromes
- **Partition Equal Subset Sum** - Divide into equal halves

## Complexity Reference

| Problem | Time Complexity | Space Complexity | State |
|---------|-----------------|------------------|-------|
| Fibonacci | O(n) | O(1) | 1D |
| Climbing Stairs | O(n) | O(1) | 1D |
| Coin Change | O(n × amount) | O(amount) | 1D |
| LCS | O(m × n) | O(m × n) or O(min(m,n)) | 2D |
| LIS | O(n²) or O(n log n) | O(n) | 1D |
| Edit Distance | O(m × n) | O(m × n) or O(n) | 2D |
| 0/1 Knapsack | O(n × W) | O(W) | 2D → 1D |
| Matrix Chain | O(n³) | O(n²) | 2D |

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_intro_memoization.ipynb](01_intro_memoization.ipynb) | Fibonacci, recursion vs memoization |
| [02_tabulation.ipynb](02_tabulation.ipynb) | Bottom-up approach, LCS, edit distance |
| [03_knapsack.ipynb](03_knapsack.ipynb) | 0/1 knapsack and subset sum |
| [04_sequence_alignment.ipynb](04_sequence_alignment.ipynb) | Needleman-Wunsch and Smith-Waterman |

## Interactive Visualizations

- [DP Table Visualization](../interactive/dp/table-fill.html) - Watch DP tables fill step by step
- [Knapsack Solver](../interactive/dp/knapsack.html) - Interactive 0/1 knapsack
- [LCS/LIS Visualizer](../interactive/dp/sequences.html) - Sequence alignment visualization

## Problem-Solving Framework

### Step 1: Define the State
What information do we need to represent a subproblem?
- `dp[i]` = answer for first i elements
- `dp[i][j]` = answer for elements i to j

### Step 2: Identify Base Cases
What are the trivial solutions?
- `dp[0] = 0` (empty input)
- `dp[1] = 1` (single element)

### Step 3: Write the Recurrence
How do we build larger solutions from smaller ones?
- `dp[i] = dp[i-1] + dp[i-2]` (Fibonacci)
- `dp[i][j] = max(dp[i-1][j], dp[i-1][j-w[i]] + v[i])` (Knapsack)

### Step 4: Determine Order of Computation
Which subproblems must be solved first?
- Usually left-to-right, bottom-to-top

### Step 5: Extract the Answer
Where is the final answer?
- Often `dp[n]` or `dp[n][target]`

## Common Patterns

### Pattern 1: Linear DP
```
dp[i] depends on dp[i-1], dp[i-2], ...
Examples: Fibonacci, Climbing Stairs, House Robber
```

### Pattern 2: Two Sequence DP
```
dp[i][j] = f(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
Examples: LCS, Edit Distance, Regex Matching
```

### Pattern 3: Interval DP
```
dp[i][j] = f(dp[i][k], dp[k+1][j]) for k in [i, j)
Examples: Matrix Chain, Burst Balloons, Palindrome Partitioning
```

### Pattern 4: Knapsack DP
```
dp[i][w] = max(include item i, exclude item i)
Examples: 0/1 Knapsack, Subset Sum, Coin Change
```

## Tips for Interviews

1. **Identify DP problems**: Optimization, counting, or decision problems with overlapping subproblems
2. **Start with recursion**: Write the brute-force recursive solution first
3. **Add memoization**: Cache results to avoid recomputation
4. **Convert to tabulation**: If space optimization is needed
5. **Trace back**: To reconstruct the actual solution (not just the value)
