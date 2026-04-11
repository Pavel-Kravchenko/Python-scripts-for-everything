---
name: genome-assembly
description: De novo genome assembly — N50/L50/BUSCO metrics, de Bruijn graphs, OLC approach, assembler selection (SPAdes/Flye/Hifiasm), assembly graph concepts, contamination checking
---

# Genome Assembly

## When to Use
- De novo genome assembly or evaluating existing assemblies (N50, L50, BUSCO, assembly graph)
- Understanding de Bruijn graph vs overlap-layout-consensus (OLC) approaches
- Choosing an assembler for short reads, long reads, or hybrid assemblies
- Assessing completeness and contiguity of a new genome
- Detecting contamination or heterozygosity in an assembly

---

## Quick Reference

### Assembly Quality Metrics
| Metric | Definition | Good value (typical eukaryote) |
|--------|-----------|-------------------------------|
| **N50** | Length L such that 50% of assembled bases are in contigs ≥ L | > 1 Mb (chromosome-scale) |
| **L50** | Number of contigs whose combined length = N50 | Lower is better |
| **NG50** | N50 normalised to estimated genome size | Comparable across species |
| **BUSCO** | Fraction of conserved single-copy orthologs present | > 95% complete |
| **QV score** | Phred-like per-base accuracy of assembly | > 40 (Merqury) |
| **k-mer completeness** | Fraction of k-mers from reads present in assembly | > 95% |

### Assembler Comparison
| Assembler | Read type | Algorithm | Best for |
|-----------|-----------|-----------|----------|
| **SPAdes** | Illumina | Multi-k de Bruijn | Bacteria, small eukaryotes, metagenomes |
| **Velvet** | Illumina | De Bruijn | Legacy short-read assembler |
| **Canu** | PacBio CLR / ONT | OLC + correction | Large genomes, high repeat content |
| **Flye** | PacBio / ONT | Repeat graph | Long-read, incl. highly repetitive genomes |
| **Hifiasm** | PacBio HiFi | String graph | Haplotype-resolved assemblies |
| **MaSuRCA** | Hybrid | Super-reads | Short + long reads combined |
| **Verkko** | HiFi + ONT | Multiplex de Bruijn | Telomere-to-telomere assemblies |

### Short-Read vs Long-Read Assembly
| Property | Short-read (Illumina) | Long-read (PacBio/ONT) |
|----------|----------------------|----------------------|
| Read length | 150–300 bp | 10 kb – 100 kb |
| Error rate | ~0.1% | 0.1% (HiFi) – 5–15% (CLR/ONT) |
| Repeat resolution | Poor (reads < repeat length) | Excellent |
| Phasing | Requires linked-reads | Direct for HiFi |
| Cost/Gb | Low | Higher |

---

## Key Patterns

### N50 Calculation
```python
def n50(contig_lengths: list[int],
        genome_size: int | None = None) -> tuple[int, int]:
    """Return (N50, L50). If genome_size given, computes NG50."""
    lengths = sorted(contig_lengths, reverse=True)
    target = (genome_size or sum(lengths)) / 2
    cumsum = 0
    for i, length in enumerate(lengths, 1):
        cumsum += length
        if cumsum >= target:
            return length, i   # (N50/NG50, L50)
    return lengths[-1], len(lengths)
```

### De Bruijn Graph Construction
```python
from collections import defaultdict

def build_debruijn_graph(reads: list[str], k: int) -> dict[str, list[str]]:
    """Nodes = (k-1)-mers; edges = k-mers. Returns adjacency list."""
    graph: dict[str, list[str]] = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            graph[kmer[:-1]].append(kmer[1:])
    return dict(graph)
```

### Euler Path Assembly (toy)
```python
def eulerian_path_contigs(graph: dict[str, list[str]]) -> list[str]:
    """Assemble contigs by following non-branching paths (unitigs)."""
    from copy import deepcopy
    g = deepcopy(graph)
    visited_edges: set[tuple] = set()
    contigs = []

    def extend(start):
        path = [start]
        while True:
            node = path[-1]
            neighbors = [v for v in g.get(node, [])
                         if (node, v, len(path)) not in visited_edges]
            if not neighbors:
                break
            nxt = neighbors[0]
            visited_edges.add((node, nxt, len(path)))
            path.append(nxt)
        return path

    all_nodes = set(graph) | {v for vs in graph.values() for v in vs}
    for node in all_nodes:
        path = extend(node)
        if len(path) > 1:
            contigs.append(path[0] + ''.join(p[-1] for p in path[1:]))
    return contigs
```

---

## Code Templates

### Assembly Statistics Report
```python
import numpy as np

def assembly_stats(contig_lengths: list[int],
                   genome_size: int | None = None) -> dict:
    lengths = sorted(contig_lengths, reverse=True)
    total = sum(lengths)
    n50_val, l50_val = n50(lengths, genome_size)
    return {
        'n_contigs':     len(lengths),
        'total_bases':   total,
        'largest_contig': lengths[0],
        'N50':           n50_val,
        'L50':           l50_val,
        'mean_length':   int(np.mean(lengths)),
        'median_length': int(np.median(lengths)),
    }
```

### SPAdes Assembly (shell)
```bash
# Short-read assembly
spades.py -1 R1.fastq.gz -2 R2.fastq.gz -o spades_out/ --threads 16

# Hybrid assembly (Illumina + ONT)
spades.py -1 R1.fastq.gz -2 R2.fastq.gz \
          --nanopore long_reads.fastq.gz \
          -o hybrid_out/ --threads 16
```

### Flye Long-Read Assembly (shell)
```bash
# PacBio HiFi
flye --pacbio-hifi reads.fastq.gz --out-dir flye_out/ --threads 16

# ONT (high error rate)
flye --nano-raw reads.fastq.gz --out-dir flye_out/ \
     --genome-size 3g --threads 16
```

### BUSCO Completeness Check (shell)
```bash
busco -i assembly.fa -l vertebrata_odb10 \
      -o busco_out -m genome --cpu 8

# Summarise
cat busco_out/short_summary.specific.vertebrata_odb10.busco_out.txt
```

### Contamination Screen with Blobtools (shell)
```bash
# Step 1: align reads back to assembly for coverage
minimap2 -ax sr assembly.fa R1.fastq R2.fastq | samtools sort -o reads.bam
samtools index reads.bam

# Step 2: BLAST assembly against nt (or use Kraken2)
blastn -query assembly.fa -db nt -out blast.out \
       -outfmt "6 qseqid staxids bitscore" \
       -max_target_seqs 10 -evalue 1e-25 -num_threads 8

# Step 3: Blobtools
blobtools create -i assembly.fa -b reads.bam -t blast.out -o blob
blobtools view -i blob.blobDB.json -r species -o blob_view
```

---

## Common Pitfalls

- **k-mer size selection**: too small → repetitive regions collapse; too large → sparse coverage. For Illumina 150 bp reads, k = 51–77 is typical; use assembler auto-tuning (SPAdes) when unsure.
- **Coverage requirement**: minimum ~30× for reliable de Bruijn graph; < 10× means many k-mers appear only once and are discarded as errors.
- **N50 alone is insufficient**: pair N50 with BUSCO; a high N50 with < 90% BUSCO completeness suggests misassembly or highly fragmented biology.
- **Heterozygosity**: diploid organisms produce bubble structures in the graph. Use `--haplotype` modes (Hifiasm) or purge haplotigs (purge_dups) before final assembly.
- **Contamination**: always run BlobTools or Kraken2 before submission; contaminant contigs inflate assembly size and confuse annotation.
- **Long-read error models**: Flye handles raw ONT reads; Hifiasm requires high-accuracy HiFi reads (> 99% accuracy). Do not mix parameters.
- **Circular genomes** (bacteria, organelles): use `--plasmid` mode in SPAdes or check for circular paths in the assembly graph (Bandage).

---

## Related Skills
- `ngs-variant-calling` — NGS fundamentals, short-read alignment, SAM/BAM (prerequisite)
- `long-read-sequencing` — ONT/PacBio basecalling, Minimap2, Flye/Hifiasm, SV calling
- `metagenomics-shotgun` — MEGAHIT metagenome assembly, MetaBAT2 binning, CheckM MAGs
- `structural-bioinformatics` — downstream structural analysis of assembled proteins
