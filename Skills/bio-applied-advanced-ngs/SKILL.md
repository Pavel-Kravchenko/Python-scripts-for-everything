---
name: bio-applied-advanced-ngs
description: "This notebook builds directly on **3.01 NGS Fundamentals** (sequencing platforms, FASTQ, basic QC, SAM/BAM, BWA/HISAT2 alignment). Here we go deeper: assembling genomes *de novo* when no reference exi"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/17_Genome_Assembly_and_Advanced_NGS/02_advanced_ngs.ipynb"
---

# Genome Assembly and Advanced NGS

*Source: Course notebook `Tier_3_Applied_Bioinformatics/17_Genome_Assembly_and_Advanced_NGS/02_advanced_ngs.ipynb`*

# Genome Assembly and Advanced NGS

## Tier 3 - Applied Bioinformatics

This notebook builds directly on **3.01 NGS Fundamentals** (sequencing platforms, FASTQ, basic QC, SAM/BAM, BWA/HISAT2 alignment). Here we go deeper: assembling genomes *de novo* when no reference exists, understanding the algorithms that make alignment and assembly work, assessing assembly quality, and working with long-read technologies.

### Learning Objectives
- Understand when and why *de novo* assembly is used instead of reference mapping
- Implement a toy Overlap-Layout-Consensus (OLC) assembler
- Implement a toy de Bruijn graph assembler with k-mer decomposition
- Calculate N50, L50, and NG50 assembly quality metrics from scratch
- Understand BUSCO, QUAST, and coverage-based quality assessment
- Explain scaffolding: paired-end, mate-pair, Hi-C
- Implement Burrows-Wheeler Transform and understand FM-index backward search
- Understand k-mer frequency analysis for genome size estimation (GenomeScope)
- Compare Oxford Nanopore and PacBio HiFi long-read technologies
- Design hybrid assembly strategies

---

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict, Counter
import random
import math
import itertools

random.seed(42)
np.random.seed(42)

print("Libraries loaded.")
```

---

## Part 1: Genome Assembly Concepts

### 1.1 Reference Mapping vs. De Novo Assembly

After sequencing, you have millions of short reads (100–300 bp for Illumina). Two fundamentally different analysis routes exist:

| Approach | When to use | Key requirement |
|----------|------------|----------------|
| **Reference mapping** | Resequencing a known species | Reference genome |
| **De novo assembly** | New organism, no reference | Sufficient coverage (~50–100×) |

**De novo assembly** is needed for:
- Newly discovered organisms (environmental metagenomics, rare species)
- Highly diverged strains where mapping fails
- Structural variant discovery (inversions, translocations)
- Telomere-to-telomere "complete" reference genomes

**The shotgun assembly problem** (illustrated by the lecture's newspaper analogy):
- Many copies of the genome are shredded into fragments 200–200,000 bp long
- Each fragment is sequenced from one or both ends (100–1000 bp reads)
- No position information is retained — reads arrive in random order
- **Goal**: reconstruct the original sequence from overlapping fragments

The challenge grows exponentially with genome size and repeat content. For a bacterial genome (~5 Mb) there may be 2^200 possible reconstructions; for a mammalian genome the number is astronomical.

### 1.2 The Overlap-Layout-Consensus (OLC) Approach

OLC is the classical approach, originally used for Sanger reads and now the basis of long-read assemblers (Canu, Flye):

1. **Overlap**: Find all pairs of reads that share a suffix-prefix overlap
2. **Layout**: Build an overlap graph; reads are nodes, edges weighted by overlap length. Find a Hamiltonian path (visits each node once)
3. **Consensus**: Derive a consensus sequence by majority-vote over aligned reads

**Why it works for long reads**: overlaps are long and unique enough to resolve most repeats. **Why it fails for short reads**: O(n²) pairwise comparisons are computationally prohibitive for billions of short reads.

```python
def find_overlaps(reads: list[str], min_overlap: int = 4) -> list[tuple]:
    """Find all suffix-prefix overlaps between reads."""
    overlaps = []
    for i, r1 in enumerate(reads):
        for j, r2 in enumerate(reads):
            if i == j:
                continue
            # Find longest suffix of r1 that is a prefix of r2
            max_overlap = min(len(r1), len(r2))
            for k in range(max_overlap, min_overlap - 1, -1):
                if r1[-k:] == r2[:k]:
                    overlaps.append((i, j, k))
                    break
    return overlaps


def greedy_assemble_olc(reads: list[str], min_overlap: int = 4) -> str:
    """Greedy OLC assembler: repeatedly merge the pair with the longest overlap."""
    contigs = list(reads)
    while True:
        overlaps = find_overlaps(contigs, min_overlap)
        if not overlaps:
            break
        # Take the best overlap
        best = max(overlaps, key=lambda x: x[2])
        i, j, ov = best
        merged = contigs[i] + contigs[j][ov:]
        # Remove merged reads and add new contig
        new_contigs = [c for idx, c in enumerate(contigs) if idx not in (i, j)]
        new_contigs.append(merged)
        contigs = new_contigs
    return contigs


# Demonstration with the lecture example
reads = [
    "AGCTACAGTATGCT",
    "TACAGTATGCTTAT",
    "GTATGCTTATCTGA",
    "TGCTTATCTGATAC",
    "TGATACCTTAGCCA",
]

print("Reads:")
for r in reads:
    print(f"  {r}")

overlaps = find_overlaps(reads, min_overlap=5)
print(f"\nOverlaps found (min_overlap=5):")
print(f"{'Read i':>8} {'Read j':>8} {'Overlap':>10} {'Seq':>20}")
for i, j, ov in sorted(overlaps, key=lambda x: -x[2]):
    print(f"{i:>8} {j:>8} {ov:>10}  {reads[i][-ov:]}")

result = greedy_assemble_olc(reads, min_overlap=5)
print(f"\nAssembly result: {result[0]}")
print(f"Expected (lecture):  AGCTACAGTATGCTATCTGATACCTTAGCCA")
```

### 1.3 De Bruijn Graph Assembly for Short Reads

For Illumina reads (billions of 150 bp reads), OLC is too slow. The de Bruijn graph approach solves this by reducing the comparison problem:

**Key idea**: instead of comparing whole reads, decompose each read into overlapping k-mers.

For a read of length L with k-mer size k, we extract **L − k + 1** k-mers, each shifting by one position.

**De Bruijn graph construction**:
- **Nodes**: (k−1)-mers (the "left" and "right" ends of each k-mer)
- **Edges**: k-mers (an edge from the left (k−1)-mer to the right (k−1)-mer)
- Assembly = finding an **Eulerian path** (traverses each *edge* once)

This reduces Hamiltonian path (NP-hard) to Eulerian path (solvable in O(E) by Hierholzer's algorithm).

**Why k matters** (from lecture slides):
- Small k (e.g., k=2): many nodes with high in/out-degree → ambiguous, can't resolve repeats
- Large k (e.g., k=10): very specific edges, but harder to connect at repeat boundaries; mismatches break the graph
- Practical choice: k=31–127 for Illumina; SPAdes uses multiple k values simultaneously

**What complicates real graphs** (from lecture):
- Sequencing errors → spurious k-mers → "tips" and erroneous edges
- Diploid/polyploid organisms → heterozygous SNPs create "bubbles" in the graph
- Repeats → multiple paths through the same node ("complex nodes")

```python
def build_debruijn_graph(reads: list[str], k: int) -> dict[str, list[str]]:
    """Build a de Bruijn graph from reads.
    
    Nodes are (k-1)-mers. Edges are k-mers.
    Returns adjacency list: {left_kmer: [right_kmer, ...]}
    """
    graph = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            left = kmer[:-1]   # (k-1)-mer prefix
            right = kmer[1:]   # (k-1)-mer suffix
            graph[left].append(right)
    return dict(graph)


def get_kmers(sequence: str, k: int) -> list[str]:
    """Extract all k-mers from a sequence."""
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def eulerian_path(graph: dict[str, list[str]]) -> list[str] | None:
    """Find an Eulerian path in a directed graph using Hierholzer's algorithm."""
    # Make a mutable copy
    adj = {node: list(edges) for node, edges in graph.items()}
    # Count in-degree and out-degree
    out_deg = {node: len(edges) for node, edges in adj.items()}
    in_deg = defaultdict(int)
    for node, edges in adj.items():
        for dest in edges:
            in_deg[dest] += 1
    
    all_nodes = set(out_deg) | set(in_deg)
    
    # Find start node: out_deg - in_deg == 1 (or any node with edges if Eulerian circuit)
    start = None
    for node in all_nodes:
        diff = out_deg.get(node, 0) - in_deg.get(node, 0)
        if diff == 1:
            start = node
            break
    if start is None:
        start = next(iter(out_deg))
    
    # Hierholzer
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if adj.get(v):
            stack.append(adj[v].pop())
        else:
            path.append(stack.pop())
    
    path.reverse()
    return path


def path_to_sequence(path: list[str]) -> str:
    """Reconstruct sequence from a de Bruijn path of (k-1)-mers."""
    if not path:
        return ""
    return path[0] + "".join(node[-1] for node in path[1:])


# Demonstrate k-mer decomposition (lecture slides 22-35)
read1 = "AGCTACAGTATGC"
read2 = "TATGCTTATCTGA"

for k in [5, 10]:
    kmers1 = get_kmers(read1, k)
    kmers2 = get_kmers(read2, k)
    shared = set(kmers1) & set(kmers2)
    print(f"k={k}: {len(kmers1)} k-mers from read1, {len(kmers2)} from read2, {len(shared)} shared")
    if shared:
        print(f"  Shared (junction): {shared}")
print()

# Build graph for k=5 (from lecture)
reads_toy = ["AGCTACAGTATGC", "TATGCTTATCTGA"]
g = build_debruijn_graph(reads_toy, k=5)
print("De Bruijn graph (k=5), edges (k-mer = node_left -> node_right):")
for src in sorted(g.keys()):
    for dst in g[src]:
        kmer = src + dst[-1]
        print(f"  {src} -> {dst}  [{kmer}]")
```

```python
# Fix the adjacency copy (use items() properly) and assemble a small sequence
def assemble_debruijn(reads: list[str], k: int) -> list[str]:
    """Assemble reads using de Bruijn graph + Eulerian path.
    Returns list of contigs (one per connected component with valid path).
    """
    graph = build_debruijn_graph(reads, k)
    
    # Count in/out degrees
    out_deg = defaultdict(int)
    in_deg = defaultdict(int)
    for src, dsts in graph.items():
        for dst in dsts:
            out_deg[src] += 1
            in_deg[dst] += 1

    all_nodes = set(out_deg) | set(in_deg)
    print(f"Graph: {len(all_nodes)} nodes, {sum(out_deg.values())} edges")

    # Show degree imbalance (indicates start/end of path)
    for node in sorted(all_nodes):
        o = out_deg.get(node, 0)
        i = in_deg.get(node, 0)
        if o != i:
            label = "START" if o > i else "END"
            print(f"  {node}: out={o} in={i} [{label}]")

    # Greedy path: follow unambiguous chains (simplified)
    adj = defaultdict(list)
    for src, dsts in graph.items():
        adj[src].extend(dsts)

    visited_edges = set()
    contigs = []

    # Find start nodes (out > in, or any unvisited)
    starts = [n for n in all_nodes if out_deg.get(n, 0) > in_deg.get(n, 0)]
    if not starts:
        starts = list(out_deg.keys())

    for start in starts:
        if not adj[start]:
            continue
        path = [start]
        current = start
        while adj[current]:
            nxt = adj[current].pop(0)
            edge_key = (current, nxt, len(path))
            path.append(nxt)
            current = nxt
        if len(path) > 1:
            contigs.append(path_to_sequence(path))

    return contigs


# Test with a known sequence fragmented into short reads
target = "AGCTACAGTATGCTTATCTGATAC"
read_len = 10
step = 3
test_reads = [target[i:i+read_len] for i in range(0, len(target) - read_len + 1, step)]

print(f"Target:    {target} (len={len(target)})")
print(f"Reads ({len(test_reads)} reads, length {read_len}, step {step}):")
for r in test_reads:
    print(f"  {r}")
print()

contigs = assemble_debruijn(test_reads, k=7)
print(f"\nAssembled contigs:")
for c in contigs:
    match = "MATCH" if c == target else ("PARTIAL" if c in target else "MISMATCH")
    print(f"  {c} [{match}]")
```

### 1.4 Real-World Assemblers

| Assembler | Read Type | Algorithm | Best For |
|-----------|-----------|-----------|----------|
| **SPAdes** | Illumina (short) | Multi-k de Bruijn | Bacteria, small eukaryotes, metagenomes |
| **Velvet** | Illumina (short) | De Bruijn | Early short-read assembler |
| **Canu** | PacBio CLR / ONT | OLC + correction | Large genomes, high repeat content |
| **Flye** | PacBio CLR / ONT / HiFi | Repeat graph | Fast, handles high error rates |
| **hifiasm** | PacBio HiFi | Graph-based | Haplotype-resolved assembly, state of the art |
| **verkko** | HiFi + ONT ultra-long | Graph-based | T2T-quality assemblies |

**SPAdes** improvements over naive de Bruijn:
1. Runs multiple k values (e.g., 21, 33, 55, 77) and merges the graphs
2. Uses paired-end information to resolve the graph
3. Error correction via k-mer count thresholds (low-count k-mers are errors)
4. Performs graph simplification: tip clipping, bubble popping, low-coverage edge removal

```bash
# SPAdes basic usage
spades.py -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o spades_output/ -t 8

# For metagenomes
spades.py --meta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o metaspades_output/

# Flye for long reads
flye --nano-raw reads.fastq.gz --genome-size 5m --out-dir flye_output/ --threads 8

# hifiasm for PacBio HiFi
hifiasm -o assembly -t 16 hifi_reads.fastq.gz
```
