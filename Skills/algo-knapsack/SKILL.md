---
name: algo-knapsack
description: "1. Understand 0/1 and unbounded knapsack problems 2. Apply knapsack DP to resource optimization 3. Solve subset sum problems"
tool_type: python
source_notebook: "Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/03_knapsack.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# 🎒 Knapsack Problems

*Source: Course notebook `Tier_4_Algorithms_and_Data_Structures/10_Dynamic_Programming/03_knapsack.ipynb`*


## Learning Objectives

1. Understand 0/1 and unbounded knapsack problems
2. Apply knapsack DP to resource optimization
3. Solve subset sum problems

---

## Biological Applications

- **Experiment design:** Select genes to study within budget constraints
- **Drug combination:** Maximize efficacy within toxicity limits
- **Feature selection:** Select markers within computational budget

## The Knapsack Problem

**Formal definition:** Given n items, each with weight w_i and value v_i, and a knapsack of capacity W, choose a subset of items that maximizes total value without exceeding total weight W.

**Worked example** (capacity W = 5):

| Item | Weight | Value |
|------|--------|-------|
| A    | 2      | 3     |
| B    | 3      | 4     |
| C    | 4      | 5     |

- Take A + B: weight 2+3=5, value 3+4=7
- Take A + C: weight 2+4=6 — exceeds capacity
- Take B only: weight 3, value 4
- Take A only: weight 2, value 3

Optimal: A + B with total value **7**.

The DP table `dp[i][w]` stores the best value using the first i items with capacity w.

```python
def knapsack_verbose(weights, values, names, capacity):
    """0/1 Knapsack with step-by-step table printout."""
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for w in range(capacity + 1):
            dp[i][w] = dp[i-1][w]
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i][w], dp[i-1][w - weights[i-1]] + values[i-1])

    # Print table with labels
    header = "     " + "".join(f"{w:>4}" for w in range(capacity + 1))
    print(header)
    print("     " + "----" * (capacity + 1))
    for i in range(n + 1):
        label = f"{'':>4}" if i == 0 else f"{names[i-1]:>4}"
        row = "".join(f"{dp[i][w]:>4}" for w in range(capacity + 1))
        print(f"{label} |{row}")
    return dp[n][capacity]

# Small example
items = ['A', 'B', 'C']
weights = [2, 3, 4]
values = [3, 4, 5]
result = knapsack_verbose(weights, values, items, 5)
print(f"\nOptimal value: {result}")
```python

```python
def knapsack_01(weights, values, capacity):
    """
    0/1 Knapsack: Each item can be used at most once.
    Returns maximum value achievable.
    """
    n = len(weights)
    # dp[i][w] = max value using items 0..i-1 with capacity w
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        for w in range(capacity + 1):
            # Don't take item i-1
            dp[i][w] = dp[i-1][w]
            # Take item i-1 if it fits
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i][w], dp[i-1][w - weights[i-1]] + values[i-1])
    
    return dp[n][capacity]

# Example: Gene panel selection
# Each gene has a cost (sequencing depth) and value (clinical relevance)
genes = ['BRCA1', 'TP53', 'EGFR', 'KRAS', 'MYC', 'PTEN']
costs = [50, 30, 40, 35, 45, 25]   # Sequencing cost units
relevance = [90, 85, 70, 75, 60, 65]  # Clinical relevance score
budget = 100

max_value = knapsack_01(costs, relevance, budget)
print(f"Maximum clinical relevance within budget: {max_value}")
```python

```python
def knapsack_01_with_items(weights, values, names, capacity):
    """0/1 Knapsack with item traceback."""
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        for w in range(capacity + 1):
            dp[i][w] = dp[i-1][w]
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i][w], dp[i-1][w - weights[i-1]] + values[i-1])
    
    # Traceback
    selected = []
    w = capacity
    for i in range(n, 0, -1):
        if dp[i][w] != dp[i-1][w]:
            selected.append(names[i-1])
            w -= weights[i-1]
    
    return dp[n][capacity], selected[::-1]

value, selected = knapsack_01_with_items(costs, relevance, genes, budget)
print(f"\nOptimal gene panel: {selected}")
print(f"Total relevance: {value}")
print(f"Total cost: {sum(costs[genes.index(g)] for g in selected)}")
```python

```python
def knapsack_01_space_optimized(weights, values, capacity):
    """Space-optimized 0/1 Knapsack using O(capacity) space."""
    n = len(weights)
    dp = [0] * (capacity + 1)
    
    for i in range(n):
        # Traverse backwards to avoid using updated values
        for w in range(capacity, weights[i] - 1, -1):
            dp[w] = max(dp[w], dp[w - weights[i]] + values[i])
    
    return dp[capacity]

print(f"\nSpace-optimized result: {knapsack_01_space_optimized(costs, relevance, budget)}")
```python

## Unbounded Knapsack

In the unbounded variant, each item can be selected any number of times.

**Bioinformatics analogy:** Selecting PCR primer pairs where you can reuse the same primer design. Each primer pair has a synthesis cost (weight) and an expected amplification efficiency (value). You want to maximize coverage within a reagent budget, and the same primer pair can be ordered multiple times for different reactions.

**Key difference from 0/1 knapsack:** When filling the DP table, iterate forward (not backward) through capacities so that updates from the current item are available for larger capacities — this allows an item to be picked more than once.

```python
def knapsack_unbounded(weights, values, capacity):
    """Unbounded knapsack: items can be reused."""
    dp = [0] * (capacity + 1)
    for w in range(1, capacity + 1):
        for i in range(len(weights)):
            if weights[i] <= w:
                dp[w] = max(dp[w], dp[w - weights[i]] + values[i])
    return dp[capacity]

# PCR primer pairs: synthesis cost and amplification efficiency
primer_costs = [2, 3, 4]      # Synthesis cost units
primer_efficiency = [3, 4, 5] # Amplification efficiency score
reagent_budget = 10

result = knapsack_unbounded(primer_costs, primer_efficiency, reagent_budget)
print(f"Max efficiency with budget {reagent_budget}: {result}")

# Compare with 0/1 knapsack on the same input
print(f"0/1 knapsack result: {knapsack_01_space_optimized(primer_costs, primer_efficiency, reagent_budget)}")
```python

## Subset Sum

Special case: Can we select items to exactly reach a target sum?

```python
def subset_sum(nums, target):
    """Check if subset sums to target."""
    dp = [False] * (target + 1)
    dp[0] = True
    
    for num in nums:
        for t in range(target, num - 1, -1):
            dp[t] = dp[t] or dp[t - num]
    
    return dp[target]

# Can we select genes costing exactly 80 units?
print(f"Can select genes costing exactly 80? {subset_sum(costs, 80)}")
print(f"Can select genes costing exactly 77? {subset_sum(costs, 77)}")
```python

## Exercises

### Exercise 1 (*)

**Experiment Budget**

Given a list of experiments with costs and expected data quality scores, find the best set of experiments within a budget of 500 units. Print the selected experiments and their total score.

```python
experiments = [
    'RNA-seq replicates',
    'Whole-genome sequencing',
    'ChIP-seq H3K27ac',
    'ATAC-seq',
    'Proteomics panel',
    'Single-cell RNA-seq',
    'Hi-C chromatin',
]
exp_costs = [120, 300, 80, 90, 200, 250, 180]
exp_quality = [75, 95, 60, 65, 80, 98, 70]
budget = 500


def knapsack_with_traceback(weights, values, names, capacity):
    """0/1 Knapsack returning max value and selected item names."""
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


# YOUR CODE HERE
# Call knapsack_with_traceback with the data above and print results
best_score, chosen = knapsack_with_traceback(exp_costs, exp_quality, experiments, budget)
print(f"Best data quality score: {best_score}")
print(f"Selected experiments: {chosen}")
print(f"Total cost: {sum(exp_costs[experiments.index(e)] for e in chosen)}")
```python

### Exercise 2 (**)

**Minimum Cost Coverage**

Given a set of sequencing library prep methods, each covering a subset of genomic regions and having a cost, find the minimum-cost subset of methods that covers all required regions.

This is a variant of set cover. For small inputs, enumerate all subsets using bitmask DP.

Hint: use bitmask DP where `dp[mask]` = minimum cost to cover the regions in `mask`.

```python
# Library prep methods: each covers a subset of regions (encoded as a bitmask)
# Regions: 0=Exome, 1=Promoters, 2=Enhancers, 3=Introns
methods = ['Standard WGS', 'ATAC-seq', 'ChIP-seq', 'RNA-seq', 'Targeted panel']
method_costs = [300, 90, 80, 120, 150]
# Each method covers certain regions (bitmask over 4 regions)
method_coverage = [
    0b1111,  # Standard WGS: covers all
    0b0110,  # ATAC-seq: Promoters + Enhancers
    0b0100,  # ChIP-seq: Enhancers
    0b1001,  # RNA-seq: Exome + Introns
    0b1011,  # Targeted panel: Exome + Promoters + Introns
]
required = 0b1111  # Must cover all 4 regions


def min_cost_coverage(costs, coverage, required):
    """
    Bitmask DP for minimum cost set cover.
    dp[mask] = minimum cost to cover exactly the regions in mask.
    """
    # YOUR CODE HERE
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


min_cost = min_cost_coverage(method_costs, method_coverage, required)
print(f"Minimum cost to cover all regions: {min_cost}")
```python

---

## Summary

✅ 0/1 Knapsack: Each item used once, O(nW) time and space  
✅ Unbounded Knapsack: Items reusable, iterate forward over capacities  
✅ Space optimization: O(W) using a single rolling array  
✅ Subset sum: Special case of knapsack for exact target matching  
✅ Applications: Gene panel selection, experiment design, coverage planning

**Next:** [04_sequence_alignment.ipynb](04_sequence_alignment.ipynb) - Global and local sequence alignment

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
