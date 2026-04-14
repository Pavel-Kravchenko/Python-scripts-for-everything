---
name: algo-knapsack
description: "Knapsack DP variants — 0/1, unbounded, subset sum with traceback and space optimization"
tool_type: python
primary_tool: Python
---

# Knapsack Problems

## 0/1 Knapsack — O(nW) time, O(nW) space

Each item used at most once. `dp[i][w]` = max value using first i items with capacity w.

```python
def knapsack_01(weights, values, capacity):
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for w in range(capacity + 1):
            dp[i][w] = dp[i-1][w]
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i][w], dp[i-1][w - weights[i-1]] + values[i-1])
    return dp[n][capacity]
```

## 0/1 Knapsack with Item Traceback

```python
def knapsack_01_with_items(weights, values, names, capacity):
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for w in range(capacity + 1):
            dp[i][w] = dp[i-1][w]
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i][w], dp[i-1][w - weights[i-1]] + values[i-1])
    selected = []
    w = capacity
    for i in range(n, 0, -1):
        if dp[i][w] != dp[i-1][w]:
            selected.append(names[i-1])
            w -= weights[i-1]
    return dp[n][capacity], selected[::-1]
```

## Space-Optimized 0/1 — O(W) space

Traverse capacity **backwards** to avoid using updated values from the same row.

```python
def knapsack_01_space_optimized(weights, values, capacity):
    dp = [0] * (capacity + 1)
    for i in range(len(weights)):
        for w in range(capacity, weights[i] - 1, -1):
            dp[w] = max(dp[w], dp[w - weights[i]] + values[i])
    return dp[capacity]
```

## Unbounded Knapsack

Items can be reused. Iterate capacity **forward** so current item's updates are available.

```python
def knapsack_unbounded(weights, values, capacity):
    dp = [0] * (capacity + 1)
    for w in range(1, capacity + 1):
        for i in range(len(weights)):
            if weights[i] <= w:
                dp[w] = max(dp[w], dp[w - weights[i]] + values[i])
    return dp[capacity]
```

## Subset Sum

Special case: can we select items to exactly reach a target?

```python
def subset_sum(nums, target):
    dp = [False] * (target + 1)
    dp[0] = True
    for num in nums:
        for t in range(target, num - 1, -1):
            dp[t] = dp[t] or dp[t - num]
    return dp[target]
```

## Bitmask DP (Minimum Cost Set Cover)

```python
def min_cost_coverage(costs, coverage, required):
    INF = float('inf')
    dp = [INF] * (required + 1)
    dp[0] = 0
    for mask in range(required + 1):
        if dp[mask] == INF:
            continue
        for i, cov in enumerate(coverage):
            new_mask = mask | cov
            if new_mask <= required:
                dp[new_mask] = min(dp[new_mask], dp[mask] + costs[i])
    return dp[required] if dp[required] != INF else -1
```

## Pitfalls

- 0/1 vs unbounded: the direction of the inner loop (backward vs forward) is the critical difference
- Space-optimized 0/1 loses traceback capability — use the full 2D table if you need selected items
- Subset sum is a decision problem (boolean DP), not an optimization
