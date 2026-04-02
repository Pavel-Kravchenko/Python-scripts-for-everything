# 🧬 Module 6: An Introduction to Applied Bioinformatics (IAB)

This module contains the excellent **"An Introduction to Applied Bioinformatics"** course 
by Greg Caporaso and collaborators (Caporaso Lab, Northern Arizona University).

📚 **Original source**: [caporaso-lab/An-Introduction-to-Applied-Bioinformatics](https://github.com/caporaso-lab/An-Introduction-to-Applied-Bioinformatics)

---

## 📖 Course Contents

### 01_Getting_Started
- Reading An Introduction to Applied Bioinformatics
- Who should read IAB?
- Using the IPython Notebook
- Recommended reading list

### 02_Fundamentals ⭐ Core Algorithms
| Notebook | Topic | Key Concepts |
|----------|-------|--------------|
| `01_pairwise_alignment.ipynb` | Pairwise Sequence Alignment | Needleman-Wunsch, Smith-Waterman, scoring matrices |
| `02_homology_search.ipynb` | Sequence Homology Searching | BLAST concepts, k-mer heuristics, p-values |
| `03_multiple_alignment.ipynb` | Multiple Sequence Alignment | Progressive alignment, guide trees |
| `04_phylogenetics.ipynb` | Phylogenetic Reconstruction | UPGMA, neighbor-joining, bootstrap |
| `05_sequence_clustering.ipynb` | Sequence Mapping & Clustering | OTUs, centroid clustering |

### 03_Applications
- **Studying Microbial Diversity**
  - Alpha diversity (observed species, PD)
  - Beta diversity (Bray-Curtis, UniFrac)
  - Ordination (PCoA, PCA)
  - Working with feature tables

### 04_Exercises
- Local sequence alignment exercises
- Multiple sequence alignment exercises

### 05_Glossary
- Key bioinformatics terminology

---

## 🔑 Key Algorithms Covered

```
┌─────────────────────────────────────────────────────────────┐
│              SEQUENCE ALIGNMENT ALGORITHMS                  │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  GLOBAL (Needleman-Wunsch)     LOCAL (Smith-Waterman)      │
│  ┌─────────────────────┐       ┌─────────────────────┐     │
│  │ Aligns entire       │       │ Finds best local    │     │
│  │ sequences end-to-end│       │ similarity regions  │     │
│  │                     │       │                     │     │
│  │ ACGT---ACG          │       │    ACGTACG          │     │
│  │ ||||   |||          │       │       |||           │     │
│  │ ACGTTTTACG          │       │    ---ACG           │     │
│  └─────────────────────┘       └─────────────────────┘     │
│                                                             │
│  Dynamic Programming: O(n×m) time and space                │
└─────────────────────────────────────────────────────────────┘
```

```
┌─────────────────────────────────────────────────────────────┐
│              PHYLOGENETIC METHODS                           │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Distance-Based              Character-Based                │
│  ├── UPGMA                   ├── Maximum Parsimony          │
│  └── Neighbor-Joining        ├── Maximum Likelihood         │
│                              └── Bayesian Inference         │
│                                                             │
│         ┌──── Species A                                     │
│     ┌───┤                                                   │
│     │   └──── Species B                                     │
│  ───┤                                                       │
│     │   ┌──── Species C                                     │
│     └───┤                                                   │
│         └──── Species D                                     │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

```
┌─────────────────────────────────────────────────────────────┐
│              DIVERSITY METRICS                              │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  ALPHA DIVERSITY (within-sample)                            │
│  ├── Observed Species: Count unique taxa                    │
│  ├── Shannon Index: H' = -Σ(pi × ln(pi))                   │
│  └── Phylogenetic Diversity: Sum of branch lengths         │
│                                                             │
│  BETA DIVERSITY (between-sample)                            │
│  ├── Bray-Curtis: Abundance-weighted distance               │
│  ├── Jaccard: Presence/absence distance                     │
│  └── UniFrac: Phylogeny-weighted distance                   │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

---

## 🚀 Prerequisites

Before starting this module, you should be comfortable with:
- Python basics (Modules 1-2 of this course)
- BioPython basics (Module 5)
- Basic understanding of DNA/RNA/protein sequences

## 📦 Dependencies

```python
pip install scikit-bio pandas numpy matplotlib scipy
```

---

## 📝 Suggested Study Path

```
Week 1: 01_Getting_Started + 02_Fundamentals/01_pairwise_alignment
Week 2: 02_Fundamentals/02_homology_search + 03_multiple_alignment
Week 3: 02_Fundamentals/04_phylogenetics
Week 4: 02_Fundamentals/05_sequence_clustering
Week 5: 03_Applications (Microbial Diversity)
Week 6: 04_Exercises (hands-on practice)
```

---

## 🙏 Attribution

This content is from the IAB project:
- **Authors**: Greg Caporaso and contributors
- **License**: See original repository
- **Website**: http://readiab.org
