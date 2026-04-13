# Tier 4: Algorithms & Data Structures

The computer science foundation behind bioinformatics tools. 10 modules with 30 Jupyter notebooks, 30 interactive HTML5 visualizations, and animated GIF demonstrations.

**Study this tier alongside Tier 2 and 3** to understand _why_ the algorithms work, not just _how_ to call them.

---

## Modules

| # | Module | Notebooks | Topics |
|---|--------|-----------|--------|
| 00 | [Skills Check](00_Skills_Check/) | 1 | Self-assessment for algorithm prerequisites |
| 01 | [Complexity Analysis](01_Complexity_Analysis/) | 2 | Big-O, Big-Omega, Big-Theta, recurrences |
| 02 | [Sorting Algorithms](02_Sorting_Algorithms/) | 2 | Bubble, Selection, Insertion, Merge, Quick, Counting, Radix |
| 03 | [Searching Algorithms](03_Searching_Algorithms/) | 1 | Linear, Binary search, two-pointer |
| 04 | [Linear Data Structures](04_Linear_Data_Structures/) | 3 | Linked lists, stacks, queues, dynamic arrays |
| 05 | [Tree Structures](05_Tree_Structures/) | 3 | BST, AVL trees, Red-Black trees |
| 06 | [Hash-Based Structures](06_Hash_Based_Structures/) | 1 | Hash tables, collision resolution, Bloom filters |
| 07 | [String Algorithms](07_String_Algorithms/) | 4 | Naive, KMP, Rabin-Karp, DFA matching |
| 08 | [Advanced String Structures](08_Advanced_String_Structures/) | 4 | Tries, Aho-Corasick, suffix arrays, suffix trees |
| 09 | [Graph Algorithms](09_Graph_Algorithms/) | 5 | BFS, DFS, Dijkstra, Bellman-Ford, MST, topological sort |
| 10 | [Dynamic Programming](10_Dynamic_Programming/) | 4 | Memoization, tabulation, LCS, edit distance, knapsack |

### All Notebooks

- [00 Skills Check](00_Skills_Check/00_skills_check.ipynb)
- [01a Complexity Analysis](01_Complexity_Analysis/01_complexity_analysis.ipynb) | [01b Basic Algorithms](01_Complexity_Analysis/02_basic_algorithms.ipynb)
- [02a Comparison Sorts](02_Sorting_Algorithms/01_comparison_sorts.ipynb) | [02b Linear Sorts](02_Sorting_Algorithms/02_linear_sorts.ipynb)
- [03 Linear & Binary Search](03_Searching_Algorithms/01_linear_binary_search.ipynb)
- [04a Linked Lists](04_Linear_Data_Structures/01_linked_lists.ipynb) | [04b Stacks & Queues](04_Linear_Data_Structures/02_stacks_queues.ipynb) | [04c Dynamic Arrays](04_Linear_Data_Structures/03_dynamic_arrays.ipynb)
- [05a Binary Search Trees](05_Tree_Structures/01_binary_search_trees.ipynb) | [05b AVL Trees](05_Tree_Structures/02_avl_trees.ipynb) | [05c Red-Black Trees](05_Tree_Structures/03_red_black_trees.ipynb)
- [06 Hash Tables & Bloom Filters](06_Hash_Based_Structures/01_hash_tables_bloom.ipynb)
- [07a Naive Matching](07_String_Algorithms/01_naive_pattern_matching.ipynb) | [07b KMP](07_String_Algorithms/02_kmp_algorithm.ipynb) | [07c Rabin-Karp](07_String_Algorithms/03_rabin_karp.ipynb) | [07d DFA Matching](07_String_Algorithms/04_dfa_matching.ipynb)
- [08a Tries](08_Advanced_String_Structures/01_tries.ipynb) | [08b Aho-Corasick](08_Advanced_String_Structures/02_aho_corasick.ipynb) | [08c Suffix Arrays](08_Advanced_String_Structures/03_suffix_arrays.ipynb) | [08d Suffix Trees](08_Advanced_String_Structures/04_suffix_trees.ipynb)
- [09a Graph Representations](09_Graph_Algorithms/01_graph_representations.ipynb) | [09b BFS & DFS](09_Graph_Algorithms/02_bfs_dfs.ipynb) | [09c Dijkstra](09_Graph_Algorithms/03_dijkstra.ipynb) | [09d MST](09_Graph_Algorithms/04_mst_kruskal_prim.ipynb) | [09e Topological Sort](09_Graph_Algorithms/05_topological_sort.ipynb)
- [10a Memoization](10_Dynamic_Programming/01_intro_memoization.ipynb) | [10b Tabulation](10_Dynamic_Programming/02_tabulation.ipynb) | [10c Knapsack](10_Dynamic_Programming/03_knapsack.ipynb) | [10d Sequence Alignment](10_Dynamic_Programming/04_sequence_alignment.ipynb)

---

## Interactive Visualizations

Open `interactive/index.html` in a browser for the full hub. Includes:

- **Sorting**: bubble, merge, quicksort race, linear sorts
- **Trees**: BST operations, AVL rotations, heap operations
- **Data structures**: linked lists, stacks/queues, hash tables
- **Strings**: pattern matching step-through, trie insertion
- **Graphs**: BFS/DFS traversal, Dijkstra shortest path, MST (Kruskal/Prim)
- **DP**: table fill visualization
- **Exercises**: 6 practice sets + 4 graded assignments (400 points total)

---

## Bioinformatics Connections

Every module here directly underpins tools used in the bioinformatics tiers:

| Algorithm Module | Bioinformatics Application | Course Notebook |
|---|---|---|
| **Dynamic Programming** (10) | Needleman-Wunsch & Smith-Waterman sequence alignment | [Pairwise Alignment](../Tier_2_Core_Bioinformatics/03_Pairwise_Sequence_Alignment/01_pairwise_sequence_alignment.ipynb) |
| **String Matching** (07) | BLAST seed-and-extend heuristic | [BLAST Searching](../Tier_2_Core_Bioinformatics/04_BLAST_Searching/01_blast_searching.ipynb) |
| **Suffix Trees/Arrays** (08) | Genome indexing (BWT, FM-index) for read alignment | [NGS Fundamentals](../Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb) |
| **Tries** (08) | k-mer counting, de Bruijn graphs for assembly | [Comparative Genomics](../Tier_2_Core_Bioinformatics/12_Comparative_Genomics/01_comparative_genomics.ipynb) |
| **Hash Tables & Bloom Filters** (06) | k-mer counting (Jellyfish), read deduplication | [NGS Fundamentals](../Tier_3_Applied_Bioinformatics/01_NGS_Fundamentals/01_ngs_fundamentals.ipynb) |
| **Trees** (05) | Phylogenetic tree construction and traversal | [Phylogenetics](../Tier_2_Core_Bioinformatics/06_Phylogenetics/01_phylogenetics.ipynb) |
| **Graphs** (09) | Metabolic pathways, gene regulatory networks | [GO and Pathways](../Tier_2_Core_Bioinformatics/11_Gene_Ontology_and_Pathways/01_gene_ontology_and_pathways.ipynb) |
| **Sorting** (02) | BAM coordinate sorting, variant prioritization | [Variant Calling](../Tier_3_Applied_Bioinformatics/02_Variant_Calling_and_SNP_Analysis/01_variant_calling_and_snp_analysis.ipynb) |
| **Complexity Analysis** (01) | Choosing efficient tools for large genomes | All Tier 2--3 notebooks |

---

## Quick Reference

See the [Data Structure Complexity Reference](assets/data_structure_complexity.md) for a complete cheat sheet of time/space complexity across all data structures and algorithms covered in this tier.

---

## Kodomo Reference Implementations

Alternative implementations from the Kodomo Bioinformatics Program coursework are available in [Solutions/Kodomo_Algorithm_Implementations/](../../Solutions/Kodomo_Algorithm_Implementations/). These show student approaches to the same algorithms covered in the course notebooks.

---

## Animated Demonstrations

<p align="center">
  <img src="assets/gifs/quicksort.gif" width="320" alt="QuickSort"/>
  <img src="assets/gifs/binary-search.gif" width="320" alt="Binary Search"/>
</p>
<p align="center">
  <img src="assets/gifs/merge-sort.gif" width="320" alt="Merge Sort"/>
  <img src="assets/gifs/counting-sort.gif" width="320" alt="Counting Sort"/>
</p>
<p align="center">
  <img src="assets/gifs/avl-rotation.gif" width="320" alt="AVL Rotation"/>
  <img src="assets/gifs/bst-insert.gif" width="320" alt="BST Insert"/>
</p>
