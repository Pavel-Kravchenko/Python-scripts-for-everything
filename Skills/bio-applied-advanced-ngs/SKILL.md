---
name: bio-applied-advanced-ngs
description: "Advanced NGS: genome assembly algorithms (OLC, de Bruijn graph), k-mer theory, assembler selection table, and SPAdes/Flye/hifiasm CLI usage."
tool_type: python
primary_tool: Python
---

# Genome Assembly and Advanced NGS

## Pitfalls

- **k-mer size is critical:** Too small → ambiguous graph, unresolvable repeats. Too large → breaks at sequencing errors. Practical: k=31–127 for Illumina; SPAdes uses multiple k simultaneously.
- **OLC is O(n²) — not for short reads:** Pairwise overlap computation is prohibitive for billions of Illumina reads. Use de Bruijn graph for short reads.
- **Repeat content breaks assemblies:** Repeats longer than read length create irresolvable branches. Long reads (ONT/HiFi) are essential for repeat-rich genomes.
- **N50 is not correctness:** High N50 can result from chimeric contigs. Always run BUSCO or CheckM.
- **Coverage floors:** De novo assembly typically requires ≥50x coverage; below this, contigs break at low-coverage regions.

## Assembly Approach Selection

| Approach | Use when | Key requirement |
|---|---|---|
| Reference mapping | Known species resequencing | Reference genome |
| De novo assembly | New organism, structural variants, T2T | ≥50–100x coverage |

## Algorithm Comparison

| Algorithm | Reads | Key idea |
|---|---|---|
| OLC (Overlap-Layout-Consensus) | Long reads (Sanger, PacBio, ONT) | Hamiltonian path through overlap graph; O(n²) |
| De Bruijn graph | Short reads (Illumina) | Eulerian path through k-mer graph; O(n·L) |

### De Bruijn graph construction
```python
from collections import defaultdict

def build_debruijn_graph(reads: list[str], k: int) -> dict[str, list[str]]:
    """Nodes: (k-1)-mers. Edges: k-mers."""
    graph = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            graph[kmer[:-1]].append(kmer[1:])
    return dict(graph)

def path_to_sequence(path: list[str]) -> str:
    return path[0] + "".join(node[-1] for node in path[1:])
```

**Graph complications in practice:**
- Sequencing errors → spurious k-mers → "tip" dead ends
- Diploid SNPs → "bubbles" (two paths with same entry/exit node)
- Repeats → complex nodes with multiple in/out edges

## Assembler Reference

| Assembler | Read type | Algorithm | Best for |
|---|---|---|---|
| **SPAdes** | Illumina | Multi-k de Bruijn | Bacteria, small eukaryotes, metagenomes |
| **Velvet** | Illumina | De Bruijn | Legacy; superseded by SPAdes |
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

## SPAdes Improvements Over Naive De Bruijn

1. Runs multiple k values (e.g., 21, 33, 55, 77) and merges graphs
2. Error correction: k-mers below count threshold are discarded
3. Tip clipping: removes dead-end branches shorter than threshold
4. Bubble popping: merges heterozygous SNP paths
5. Paired-end scaffolding: links contigs using insert size distribution
