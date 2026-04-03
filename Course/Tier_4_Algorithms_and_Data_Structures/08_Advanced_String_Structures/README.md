# 08. Advanced String Structures

**Tier 4: Algorithms & Data Structures**

Index structures that make genome-scale search practical. Tries enable prefix queries in O(m) time regardless of index size; Aho-Corasick extends this to simultaneous multi-pattern search; suffix arrays and suffix trees index every substring of a text, forming the basis of the BWT/FM-index used by BWA and Bowtie2.

## Topics Covered

- Trie construction and insertion: O(m) per operation
- Prefix search and autocomplete on a trie
- Compressed tries (Patricia/radix trees)
- Aho-Corasick automaton: failure links, dictionary links, multi-pattern O(n+m+z) search
- Suffix array construction (naive O(n² log n) and SA-IS O(n))
- LCP (Longest Common Prefix) array and its applications
- Suffix tree construction: Ukkonen's O(n) online algorithm
- Generalized suffix trees and their use in longest common substring

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [01_tries.ipynb](01_tries.ipynb) | Trie construction, prefix search, autocomplete |
| 2 | [02_aho_corasick.ipynb](02_aho_corasick.ipynb) | Multi-pattern automaton with failure and dictionary links |
| 3 | [03_suffix_arrays.ipynb](03_suffix_arrays.ipynb) | Suffix array construction, LCP array, applications |
| 4 | [04_suffix_trees.ipynb](04_suffix_trees.ipynb) | Ukkonen's algorithm, generalized suffix trees |

## Complexity Reference

| Structure | Build Time | Build Space | Pattern Query | Notes |
|-----------|-----------|-------------|---------------|-------|
| Trie | O(n \|Σ\|) | O(n \|Σ\|) | O(m) | n = total chars inserted |
| Aho-Corasick | O(n \|Σ\|) | O(n \|Σ\|) | O(text + matches) | n = sum of pattern lengths |
| Suffix Array | O(n log n) | O(n) | O(m log n) | with LCP |
| Suffix Tree | O(n) | O(n) | O(m) | Ukkonen's algorithm |

m = query pattern length, n = text length, |Σ| = alphabet size

## Bioinformatics Connections

- Genome indexing for short-read alignment (BWA, Bowtie2) is built on the Burrows-Wheeler Transform, which is computed from the suffix array. See [NGS Fundamentals](../../Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb).
- de Bruijn graphs for genome assembly use a trie-like structure over k-mers. See [Comparative Genomics](../../Tier_2_Core_Bioinformatics/12_Comparative_Genomics/01_comparative_genomics.ipynb).
- Aho-Corasick multi-pattern scanning is used for simultaneous motif search across a genome, avoiding repeated single-pattern passes.
- Suffix trees support finding the longest common substring between two sequences, useful for dot-plot style comparisons.

## Prerequisites

- [07. String Algorithms](../07_String_Algorithms/README.md) -- DFA construction and KMP are conceptual prerequisites for Aho-Corasick
- [05. Tree Structures](../05_Tree_Structures/README.md) -- tree traversal algorithms

---

[<< 07. String Algorithms](../07_String_Algorithms/README.md) | [Tier 4 Overview](../README.md) | [09. Graph Algorithms >>](../09_Graph_Algorithms/README.md)
