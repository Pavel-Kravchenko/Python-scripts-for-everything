# 🧬 Python for Bioinformatics - Course Overview

## Welcome!

This comprehensive course teaches Python programming with a focus on bioinformatics applications. 
Whether you're a biology student learning to code or a programmer entering the life sciences,
this course will give you practical skills for computational biology.

**📓 Start Here:** Open `Course/00_MASTER_INDEX.ipynb` for the interactive course navigator!

---

## 🏗️ Course Architecture

The course follows a **tiered progression** from foundational skills to advanced applications:

```
┌─────────────────────────────────────────────────────────────┐
│                    LEARNING PATHWAY                          │
├─────────────────────────────────────────────────────────────┤
│  TIER 0: Foundations     → Linux, Git, Bash, Statistics     │
│          ↓                                                   │
│  TIER 1: Python          → Variables → OOP → Data Analysis  │
│          ↓                                                   │
│  TIER 2: Core Bio        → BioPython, Alignment, Structure  │
│          ↓                                                   │
│  TIER 3: Applied Bio     → NGS, RNA-seq, ML, Clinical       │
│          ↓                                                   │
│  TIER 4: Algorithms      → Sorting, Trees, Graphs, DP       │
└─────────────────────────────────────────────────────────────┘
```

---

## 📚 Tier Overview

### Tier 0: Computational Foundations (7 modules)
*Prerequisites - skip if you're comfortable with command line*

| Module | Topic | Key Skills |
|--------|-------|------------|
| 0.0 | Skills Check | Test your readiness |
| 0.1 | Linux Basics | Navigation, file operations |
| 0.2 | Git Version Control | Commits, branches, collaboration |
| 0.3 | Bash Scripting | Automation, pipelines |
| 0.4 | File Encodings | UTF-8, ASCII, binary |
| 0.5 | R Basics | Statistical computing |
| 0.6 | Biostatistics | Hypothesis testing |

### Tier 1: Python for Bioinformatics (18 modules)
*Core Python programming with biological examples*

| Part | Modules | Topics |
|------|---------|--------|
| A: Fundamentals | 1-7 | Variables, control flow, functions, files |
| B: Data Structures | 8-11 | Lists, dicts, iterators, generators |
| C: Advanced Python | 12-15 | Regex, OOP, decorators, error handling |
| D: Data Science | 16-18 | NumPy, Pandas, visualization |

### Tier 2: Core Bioinformatics (12 modules)
*Essential bioinformatics algorithms and tools*

| Module | Topic | Key Skills |
|--------|-------|------------|
| 2.1 | Biological Databases | NCBI, UniProt, Ensembl |
| 2.2 | BioPython Essentials | Seq, SeqIO, Entrez |
| 2.3-5 | Sequence Alignment | Pairwise, BLAST, MSA |
| 2.6 | Phylogenetics | Trees, evolution |
| 2.7-8 | Molecular Structure | Proteins, nucleic acids |
| 2.9-12 | Advanced Analysis | Motifs, GO, genomics |

### Tier 3: Applied Bioinformatics (11 modules)
*Real-world applications and advanced topics*

| Module | Topic | Applications |
|--------|-------|--------------|
| 3.1-2 | NGS & Variants | FASTQ, VCF, SNP calling |
| 3.3 | RNA-seq | Differential expression |
| 3.4 | Metagenomics | Microbial diversity |
| 3.5-6 | Statistics | Regulatory analysis, FDR |
| 3.7 | Machine Learning | Classification, clustering |
| 3.8 | Capstone Project | Integration |
| 3.9-11 | Advanced | Docking, deep learning, clinical |

### Tier 4: Algorithms & Data Structures (10 modules)
*Computer science fundamentals for efficient bioinformatics*

| Module | Topic | Key Skills |
|--------|-------|------------|
| 4.1 | Complexity Analysis | Big-O notation, algorithm analysis |
| 4.2 | Sorting Algorithms | QuickSort, MergeSort, Counting Sort |
| 4.3 | Searching Algorithms | Binary search, interpolation |
| 4.4 | Linear Data Structures | Linked lists, stacks, queues |
| 4.5 | Tree Structures | BST, AVL trees, Red-Black trees |
| 4.6 | Hash-Based Structures | Hash tables, Bloom filters |
| 4.7 | String Algorithms | KMP, Rabin-Karp pattern matching |
| 4.8 | Advanced Strings | Tries, suffix arrays, suffix trees |
| 4.9 | Graph Algorithms | BFS, DFS, shortest paths |
| 4.10 | Dynamic Programming | Memoization, sequence alignment |

---

## 🎯 Learning Objectives

By the end of this course, you will be able to:

1. **Write Python code** for bioinformatics tasks
2. **Parse biological file formats** (FASTA, GenBank, PDB, CSV)
3. **Analyze DNA/RNA/protein sequences** programmatically
4. **Visualize biological data** with matplotlib and seaborn
5. **Use BioPython** for complex bioinformatics workflows
6. **Build reusable tools** for your research
7. **Analyze algorithm complexity** using Big-O notation
8. **Implement efficient data structures** for biological data

---

## 🔬 Practical Applications

### What You'll Build:

- GC content calculator with sliding window
- DNA to protein translator
- FASTA file parser and analyzer
- Restriction site finder
- Open Reading Frame (ORF) detector
- Gene expression heatmaps
- Genome GC landscape visualizer
- **Sequence alignment tools** with BLOSUM scoring
- **BLAST result analyzer** with homology assessment
- **Molecular visualization scripts** (Jmol/PyMol)
- **DNA structure models** (A/B/Z forms)

---

## 📚 How to Use This Course

### Learning Paths

**🆕 Complete Beginner:**
```
Tier 0 → Tier 1 (all) → Tier 2 → Tier 3 → Tier 4
Estimated time: 16-20 weeks
```

**🐍 Know Programming, New to Python:**
```
Tier 1 (Part A quick review) → Tier 1 (B,C,D) → Tier 2 → Tier 3
Estimated time: 8-12 weeks
```

**🧬 Python User, New to Bioinformatics:**
```
Take Tier 1 Skills Check → Tier 2 → Tier 3
Estimated time: 6-8 weeks
```

**🔬 Experienced Bioinformatician:**
```
Take Tier 2 Skills Check → Tier 3 (selected topics)
Estimated time: 2-4 weeks
```

### For Self-Study:

1. **Take the Skills Check** - At the start of each tier to test if you can skip
2. **Run all code cells** - Don't just read, execute and experiment
3. **Complete assignments** - Practice problems in `Assignments/` folder
4. **Build your own tools** - Apply concepts to your own data

### For Instructors:

- Each module includes lecture notebooks with explanations
- Assignments have starter code and full solutions in `Solutions/`
- Data files provided in `Data/` for all exercises
- Difficulty levels marked (⭐ to ⭐⭐⭐)

---

## 🛠️ Setup Checklist

- [ ] Linux/macOS terminal or WSL on Windows
- [ ] Git installed and configured
- [ ] Python 3.8+ installed
- [ ] Jupyter Notebook working
- [ ] BioPython installed (`pip install biopython`)
- [ ] Pandas installed (`pip install pandas`)
- [ ] Matplotlib installed (`pip install matplotlib`)
- [ ] Course repository cloned

### Optional Tools (for Module 5):
- [ ] Jmol (molecular visualization): https://jmol.sourceforge.net
- [ ] PyMol (publication figures): https://pymol.org
- [ ] 3DNA (nucleic acid modeling): http://x3dna.org
- [ ] py3Dmol (`pip install py3Dmol`) for Jupyter visualization

---

## 📖 Quick Reference

### Common DNA Operations

```python
# Count nucleotides
seq.count('G')

# GC content
(seq.count('G') + seq.count('C')) / len(seq) * 100

# Reverse complement
seq.translate(str.maketrans('ATGC', 'TACG'))[::-1]

# Extract codons
[seq[i:i+3] for i in range(0, len(seq)-2, 3)]
```

### BioPython Essentials

```python
from Bio.Seq import Seq
from Bio import SeqIO

# Create sequence
dna = Seq("ATGCGATC")

# Transcribe & translate
rna = dna.transcribe()
protein = dna.translate()

# Parse FASTA
for record in SeqIO.parse("file.fasta", "fasta"):
    print(record.id, len(record))
```

---

## 🎓 Prerequisites

- **Required**: Basic biology knowledge (DNA, RNA, proteins)
- **Helpful**: Some programming experience
- **Not Required**: Advanced math or statistics

---

## 📞 Getting Help

1. Check the course notebooks for explanations
2. Review the solutions for guidance
3. Consult the BioPython documentation
4. Search Stack Overflow for common issues

---

## 🚀 Let's Begin!

Start with: `Course/00_Computational_Foundations/01_Linux_Basics/01_linux_fundamentals.ipynb`

Or if you're already comfortable with command line and Git, jump to:
`Course/01_Python_Basics/01_Variables_DataTypes/01_variables_datatypes.ipynb`

Good luck and happy coding! 🧬
