---
name: bio-applied-genome-assembly
description: "Genome assembly algorithms (OLC and de Bruijn graph), k-mer selection, assembler comparison table, and SPAdes/Flye/hifiasm usage."
tool_type: python
primary_tool: Python
---

# Genome Assembly

## Pitfalls

- **k-mer size is critical:** Too small (k=2) → ambiguous graph, can't resolve repeats. Too large → breaks at sequencing errors, poor connectivity. Practical: k=31–127 for Illumina; SPAdes uses multiple k simultaneously.
- **Repeat content breaks assemblies:** Repeats longer than the read length create irresolvable branches in the de Bruijn graph. Long reads (PacBio/ONT) are essential for repeat-rich genomes.
- **Short reads cannot use OLC:** O(n²) pairwise overlaps are computationally prohibitive for billions of Illumina reads — use de Bruijn graph instead.
- **N50 is not a measure of correctness:** High N50 can result from chimeric contigs. Always run BUSCO or CheckM after assembly.
- **Coverage matters:** De novo assembly requires ~50–100x coverage. Lower coverage creates gaps and breaks contigs at low-coverage regions.

## Algorithm Comparison

| Algorithm | Reads | Key idea | Complexity |
|---|---|---|---|
| OLC (Overlap-Layout-Consensus) | Long reads (PacBio, ONT, Sanger) | Hamiltonian path in overlap graph | O(n²) overlaps |
| De Bruijn graph | Short reads (Illumina) | Eulerian path in k-mer graph | O(n·L) |

### OLC — Python skeleton
```python
from collections import defaultdict

def find_overlaps(reads: list[str], min_overlap: int = 4) -> list[tuple]:
    overlaps = []
    for i, r1 in enumerate(reads):
        for j, r2 in enumerate(reads):
            if i == j: continue
            for k in range(min(len(r1), len(r2)), min_overlap - 1, -1):
                if r1[-k:] == r2[:k]:
                    overlaps.append((i, j, k))
                    break
    return overlaps
```

### De Bruijn graph — Python skeleton
```python
def build_debruijn_graph(reads: list[str], k: int) -> dict[str, list[str]]:
    """Nodes are (k-1)-mers; edges are k-mers."""
    graph = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            graph[kmer[:-1]].append(kmer[1:])
    return dict(graph)

def path_to_sequence(path: list[str]) -> str:
    return path[0] + "".join(node[-1] for node in path[1:])
```

**de Bruijn graph complications in practice:**
- Sequencing errors → spurious k-mers → "tip" branches (dead ends)
- Diploid SNPs → "bubbles" (two paths with same entry/exit)
- Repeats → complex nodes with multiple in/out edges

## Assembler Reference

| Assembler | Read type | Algorithm | Best for |
|---|---|---|---|
| **SPAdes** | Illumina | Multi-k de Bruijn | Bacteria, small eukaryotes, metagenomes |
| **Canu** | PacBio CLR / ONT | OLC + correction | Large genomes, high repeat content |
| **Flye** | PacBio CLR / ONT / HiFi | Repeat graph | Fast, tolerates high error rates |
| **hifiasm** | PacBio HiFi | Graph-based | Haplotype-resolved assembly |
| **verkko** | HiFi + ONT ultra-long | Graph-based | Telomere-to-telomere assemblies |

## CLI Usage

```bash
# SPAdes — Illumina paired-end
spades.py -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o spades_out/ -t 8

# SPAdes — metagenome
spades.py --meta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o metaspades_out/

# Flye — long reads
flye --nano-raw reads.fastq.gz --genome-size 5m --out-dir flye_out/ --threads 8

# hifiasm — PacBio HiFi
hifiasm -o assembly -t 16 hifi_reads.fastq.gz

# Assembly QC
quast.py contigs.fasta -r reference.fasta -o quast_out/
busco -i contigs.fasta -l bacteria_odb10 -o busco_out/ -m genome
```

## SPAdes Graph Simplification Steps

1. Runs multiple k values (e.g., 21, 33, 55, 77) and merges graphs
2. Error correction: k-mers below coverage threshold are discarded
3. Tip clipping: dead-end branches shorter than threshold
4. Bubble popping: merges alternative paths caused by SNPs
5. Paired-end scaffolding: links contigs using insert size information
