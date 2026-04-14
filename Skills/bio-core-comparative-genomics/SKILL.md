---
name: bio-core-comparative-genomics
description: Dot plots, synteny analysis, genomic rearrangement detection, ortholog/paralog distinction, pan-genome concepts, and tool selection for pairwise genome alignment.
tool_type: python
primary_tool: NumPy
---

## When to Use
- Visualising sequence similarity and structural rearrangements (dot plots)
- Detecting conserved gene order between species (synteny blocks)
- Identifying inversions, translocations, duplications from genome alignments
- Choosing a whole-genome alignment tool

## Key Concepts

**Ortholog vs. paralog**: Orthologs diverged via speciation (same function across species). Paralogs diverged via duplication (may gain new functions). Reliable ortholog inference requires sequence similarity + syntenic context.

**Dot plot diagonals encode structure**:
- Main diagonal (bottom-left → top-right) = conserved linear segment
- Anti-diagonal = inverted segment (reverse complement match)
- Off-main-diagonal block = transposition
- Repeated parallel diagonals = tandem repeats

**Core genome vs. pan-genome**: Core = genes in all strains. Pan = core + accessory (some strains) + unique (single strain). Bacterial pan-genomes are "open" (grow with more sequencing); eukaryotic are "closed".

**Identity thresholds are heuristic**: >90% ANI = same species; >70% identity usable with LASTZ; <30% = twilight zone (synteny required). ANI < 95% across strains = different species.

## Tool Selection

| Scenario | Tool |
|----------|------|
| Closely related genomes (>90% ANI), fast | MUMmer/nucmer |
| Moderate divergence (>70% identity) | LASTZ |
| Any pairwise comparison, long reads | minimap2 |
| Multiple sequence alignment | MAFFT, MUSCLE |
| Synteny visualisation | MCScanX, SynMap, JCVI |

## Dot Plot Implementation

```python
import numpy as np

def filtered_dotplot(seq1, seq2, window=5, stringency=0.6):
    """Sliding-window dot plot. Returns 2D match-score matrix."""
    threshold = int(window * stringency)
    matrix = np.zeros((len(seq2), len(seq1)))
    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            matches = sum(seq1[j+k] == seq2[i+k] for k in range(window))
            if matches >= threshold:
                matrix[i + window//2, j + window//2] = matches / window
    return matrix

def dna_dotplot(seq1, seq2, window=11, stringency=0.7):
    """Dot plot with reverse-complement detection. Returns (fwd, rev) matrices."""
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    rc   = lambda s: ''.join(comp.get(c,c) for c in reversed(s))
    threshold = int(window * stringency)
    fwd = np.zeros((len(seq2), len(seq1)))
    rev = np.zeros((len(seq2), len(seq1)))
    seq2_rc = rc(seq2)
    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            fm = sum(seq1[j+k] == seq2[i+k] for k in range(window))
            if fm >= threshold:
                fwd[i + window//2, j + window//2] = fm / window
            rm = sum(seq1[j+k] == seq2_rc[i+k] for k in range(window))
            if rm >= threshold:
                rev[len(seq2)-1-i - window//2, j + window//2] = rm / window
    return fwd, rev
```

## Synteny Analysis

```python
def find_ortholog_pairs(genes_a, genes_b):
    """Match genes by name (proxy for ortholog identification)."""
    name_to_b = {g['name']: g for g in genes_b}
    return [(ga, name_to_b[ga['name']]) for ga in genes_a if ga['name'] in name_to_b]

def detect_synteny_blocks(pairs):
    """Consecutive gene pairs in conserved order and strand relationship."""
    blocks, current = [], [pairs[0]]
    for i in range(1, len(pairs)):
        prev_a, prev_b = pairs[i-1]
        curr_a, curr_b = pairs[i]
        same_order     = curr_b['start'] > prev_b['start']
        same_strand_rel = ((curr_a['strand'] == curr_b['strand']) ==
                           (prev_a['strand'] == prev_b['strand']))
        if same_order and same_strand_rel:
            current.append(pairs[i])
        else:
            blocks.append(current)
            current = [pairs[i]]
    blocks.append(current)
    return blocks
```

Synteny block interpretation:
- Collinear block: genes in same order and relative strand — conserved region
- Inverted block: same genes, reversed order and flipped strands — inversion event
- Orphan pairs: genes not in any block — transposition or loss event

## Rearrangement Detection from Dot Plots

```python
# Build reference with known rearrangement
def random_dna(n):
    return ''.join(np.random.choice(list('ACGT'), n))

seg_A, seg_B, seg_C, seg_D = [random_dna(150) for _ in range(4)]
rc = lambda s: ''.join({'A':'T','T':'A','C':'G','G':'C'}.get(c,c) for c in reversed(s))

reference  = seg_A + seg_B + seg_C + seg_D        # A-B-C-D
rearranged = seg_A + rc(seg_C) + seg_B + seg_D    # A-C'-B-D (inversion of C, shift of B)

fwd, rev = dna_dotplot(reference, rearranged, window=11, stringency=0.7)
# fwd: blue diagonal for A (top-left) and D (bottom-right), off-diagonal for B
# rev: red anti-diagonal for inverted C
```

## Pitfalls

- **Dot plot window size trade-off**: small windows (1–3) show all matches including noise; large windows (>15) miss short conserved regions. Start with window=11, stringency=0.7 for genomic DNA.
- **Self-comparison diagonal is uninformative**: mask the main diagonal ± window/2 to reveal internal repeats.
- **Ortholog thresholds are heuristic**: >30% for bacterial orthologs, >70% for strain-level, but synteny + phylogeny are required for confident calls.
- **ANI vs percent identity**: ANI averages over all shared genomic regions; BLAST percent identity is local. Use ANI (fastANI/PyANI) for species-level comparisons.
- **Inversions appear as anti-diagonals only when reverse complement is checked**: a dot plot using only forward strand matches will show a gap for inverted regions, not a diagonal.
- **Pan-genome openness**: never extrapolate a bacterial core genome from <10 strains — the core shrinks significantly as more diverse strains are added.
- **Coordinate systems**: MUMmer uses 1-based; BED tools use 0-based. Always check which coordinate system a tool outputs before downstream merging.
