---
name: bio-applied-assembly-sv
description: Long-Read Assembly  Structural Variants with NumPy
tool_type: python
primary_tool: NumPy
---

# Long-Read Assembly  Structural Variants

- [Flye documentation](https://github.com/fenderglass/Flye)
- [Hifiasm documentation](https://github.com/chhylp123/hifiasm)
- [Sniffles2 documentation](https://github.com/fritzsedlazeck/Sniffles)

## De Novo Assembly with Flye

Flye uses a **repeat graph** assembly algorithm. Unlike overlap-layout-consensus (OLC) assemblers, Flye constructs an assembly graph where repeated sequences are collapsed into shared nodes. This handles the complex repeat structure of eukaryotic genomes better than traditional OLC approaches, because repeats that would cause an OLC graph to collapse into a tangled hairball are instead represented as traversed-multiple-times edges. The underlying graph is iterative: Flye first builds an approximate repeat graph from a disjointig set, then resolves it by checking which paths through repeat nodes are supported by bridging reads.

Input modes control how Flye models the error profile of your reads. Use `--nano-hq` for current ONT R10.4.1 reads (Q20+), `--nano-raw` for older R9.4.1 data (Q < 15), `--nano-corr` for ONT reads that have already been error-corrected by another tool, `--pacbio-hifi` for PacBio CCS/HiFi reads (Q20+), and `--pacbio-raw` for older PacBio CLR reads. Choosing the wrong mode significantly degrades assembly quality because Flye tunes its overlap finding and error tolerance per mode.

The `--genome-size` flag is required (e.g., `5m` for 5 Mb bacterial, `3g` for 3 Gb human). Flye uses this estimate to set expected coverage depth and to calibrate repeat detection thresholds. It does not need to be exact — values within ±50% are generally fine. The Flye output directory contains: `assembly.fasta` (final contigs), `assembly_graph.gfa` (assembly graph for visualization in Bandage), and `assembly_info.txt` (per-contig length, coverage, circularity, and multiplicity statistics).

| Organism | Input | Expected contigs | Typical N50 |
|---|---|---|---|
| Bacterial (5 Mb) | ONT R10.4.1, 100× | 1–5 contigs | ~1–5 Mb (near-complete) |
| Yeast (12 Mb) | ONT R10.4.1, 50× | 20–50 contigs | ~500 kb – 1 Mb |
| Human (3 Gb) | ONT R10.4.1, 30× | 1,000–3,000 contigs | ~40–80 Mb |
| Human (3 Gb) | PacBio HiFi, 30× | 500–1,500 contigs | ~80–150 Mb |
| Metagenome | ONT R10.4.1, mixed | depends on diversity | 100 kb – 5 Mb (per species) |


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate a bacterial assembly result (5 Mb genome, ONT R10.4.1)
# Flye typically produces near-complete assemblies for bacteria
n_contigs = rng.poisson(3)  # typically 1-5 contigs
n_contigs = max(n_contigs, 1)

genome_size = 5_000_000
# Simulate contig lengths that sum to approximately genome size
raw_lengths = rng.exponential(genome_size / n_contigs, n_contigs)
contig_lengths = (raw_lengths / raw_lengths.sum() * genome_size).astype(int)
# Fix rounding
contig_lengths[-1] = genome_size - contig_lengths[:-1].sum()
contig_lengths = np.sort(contig_lengths)[::-1]

# Generate coverage (depth  100x)
mean_cov = 100.0
coverages = rng.normal(mean_cov, mean_cov * 0.05, n_contigs)

# Build assembly_info.txt style table
df_info = pd.DataFrame({
    "seq_name": [f"contig_{i+1:03d}" for i in range(n_contigs)],
    "length": contig_lengths,
    "cov.": coverages.round(1),
    "circ.": ["Y" if l > 1_000_000 else "N" for l in contig_lengths],
    "repeat": ["N"] * n_contigs,
    "mult.": [1] * n_contigs,
    "alt_group": ["*"] * n_contigs,
    "graph_path": [f"edge_{i+1}" for i in range(n_contigs)],
})

print("=== Flye assembly_info.txt (simulated) ===")
print(df_info.to_string(index=False))
print(f"\nTotal assembly size: {contig_lengths.sum():,} bp")
print(f"Number of contigs: {n_contigs}")
print(f"Largest contig:    {contig_lengths[0]:,} bp ({contig_lengths[0]/genome_size*100:.1f}% of genome)")
print(f"Mean coverage:     {coverages.mean():.1f}\u00d7")

# N50 calculation
def compute_n50(lengths):
    sorted_l = np.sort(lengths)[::-1]
    cumsum = np.cumsum(sorted_l)
    idx = np.searchsorted(cumsum, cumsum[-1] / 2)
    return sorted_l[idx]

n50 = compute_n50(contig_lengths)
print(f"Assembly N50:      {n50:,} bp")
```python

## HiFi Assembly with Hifiasm

Hifiasm is the leading assembler for PacBio HiFi reads and can also assemble Hi-C+HiFi or HiFi+ONT ultra-long combinations. It uses a **haplotype-resolved** assembly graph that produces **phased** assemblies by default. Rather than collapsing heterozygous loci into a consensus, Hifiasm preserves both haplotypes as separate paths through the graph, enabling downstream allele-specific analyses that are impossible with collapsed assemblers.

HiFi reads are ~15–20 kb with Q20+ accuracy (>99%). This combination allows Hifiasm to construct assembly graphs with far fewer gaps than ONT-only assemblies. The output format is **GFA (Graph Fragment Assembly)** files, not FASTA directly — the primary assembly contigs are in `*.bp.p_ctg.gfa`. Convert to FASTA with a simple awk one-liner: `awk '/^S/{print ">"$2"\n"$3}'`. The key output files are: `sample.asm.bp.p_ctg.gfa` (primary contigs, one representative haplotype), `sample.asm.bp.hap1.p_ctg.gfa` (haplotype 1), and `sample.asm.bp.hap2.p_ctg.gfa` (haplotype 2).

**Trio binning** is available when parental Illumina reads are provided. Hifiasm uses yak to build k-mer frequency databases from each parent, then assigns each HiFi read to a parental haplotype before assembly. This produces two fully phased haplotype assemblies (paternal and maternal) — critical for allele-specific gene expression, imprinting analysis, and compound heterozygous variant phasing in clinical genetics.

Compared to Flye, Hifiasm typically produces larger contigs (higher N50) and better phasing from HiFi data. Flye remains more versatile — it accepts many read types including ONT, CLR, and mixed inputs — and is the preferred choice for metagenomics and samples where only ONT data is available. For human genomes with HiFi data, Hifiasm is generally the first choice.

| Feature | Flye | Hifiasm |
|---|---|---|
| Read types | ONT, HiFi, CLR, mixed | HiFi primary; HiFi+ONT UL |
| Phasing | Collapsed (one haplotype) | Haplotype-resolved by default |
| Metagenomics | Yes (`--meta` flag) | No |
| Typical human N50 (HiFi 30×) | ~80–120 Mb | ~120–180 Mb |
| Trio binning | No | Yes (with yak) |
| Hi-C scaffolding | Via 3rd-party | Native support |


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(10)

# Simulate contig length distributions for different assemblers (human genome, chr1-scale)
def sim_contigs(n_contigs, mean_len, cv=0.8, genome_size=3e9):
    """Simulate contig lengths."""
    lengths = rng.lognormal(np.log(mean_len), cv, n_contigs)
    lengths = (lengths / lengths.sum() * genome_size).astype(int)
    return np.sort(lengths)[::-1]

assemblers = {
    "Flye (ONT R10.4.1)":        sim_contigs(1200, 2_000_000, 1.0),
    "Hifiasm (PacBio HiFi)":     sim_contigs(400,  6_000_000, 0.9),
    "Hifiasm (HiFi + ONT UL)":   sim_contigs(150, 15_000_000, 0.8),
}

def n50(lengths):
    s = np.sort(lengths)[::-1]
    cs = np.cumsum(s)
    return s[np.searchsorted(cs, cs[-1]/2)]

print(f"{'Assembler':<28} {'Contigs':>8} {'Total (Gb)':>10} {'N50 (Mb)':>10} {'Largest (Mb)':>13}")
print("-" * 72)
for name, lens in assemblers.items():
    print(f"{name:<28} {len(lens):>8,} {lens.sum()/1e9:>10.2f} {n50(lens)/1e6:>10.1f} {lens[0]/1e6:>13.1f}")

# Plot NG50 plots (cumulative contig length vs sorted contig length)
fig, ax = plt.subplots(figsize=(9, 5))
colors = ["steelblue", "darkorange", "seagreen"]
for (name, lens), color in zip(assemblers.items(), colors):
    sorted_l = np.sort(lens)[::-1]
    cumulative = np.cumsum(sorted_l) / 3e9 * 100  # % of genome size
    ax.plot(cumulative, sorted_l / 1e6, label=name, color=color, linewidth=2)

ax.axhline(1, color="gray", linestyle=":", alpha=0.5, label="1 Mb reference line")
ax.set_xlabel("Cumulative genome coverage (%)")
ax.set_ylabel("Contig length (Mb)")
ax.set_title("Assembly contig length distributions (simulated, human genome scale)")
ax.legend()
ax.set_xlim(0, 100)
ax.set_ylim(0)
plt.tight_layout()
plt.show()
```python

## Assembly Polishing with Medaka

Polishing corrects remaining systematic errors in the draft assembly. For ONT assemblies, Flye's draft typically has ~0.1–1% error rate, dominated by homopolymer indels — insertions or deletions within runs of the same base (e.g., AAAAAAA → AAAAAA). These are characteristic of nanopore sequencing because the electrical signal changes only when a new base enters the sensing region, making homopolymer runs inherently ambiguous. Short-read polishing tools like Pilon can also address these errors if Illumina data is available.

Medaka is ONT's neural-network polisher. It maps the original reads back to the draft assembly using Minimap2, then uses a recurrent neural network (RNN) — specifically a gated recurrent unit (GRU) architecture — to call consensus corrections from the per-position pileup. Critically, Medaka requires that the model matches the flow cell chemistry used during sequencing. The model naming convention encodes chemistry: `r1041_e82_400bps_hac_v4.2.0` = pore R10.4.1, chemistry E8.2, 400 bps translocation speed, HAC basecalling model, version 4.2.0. Using the wrong model degrades polishing performance because the error signatures differ between chemistries.

Medaka is specifically designed for ONT reads. For PacBio HiFi assemblies, polishing is usually unnecessary because HiFi reads are already Q20+ and the assembly consensus is already very accurate. For ONT assemblies with <30× coverage, Medaka may not help much — insufficient depth means the RNN cannot distinguish true variants from sequencing errors. One round of polishing is usually sufficient for high-depth (>50×) ONT data; additional rounds give diminishing returns and can introduce systematic biases.

| Error type | Flye draft | After Medaka | Improvement |
|---|---|---|---|
| Substitutions | ~0.03% | ~0.005% | ~6× |
| Insertions | ~0.15% | ~0.02% | ~7× |
| Deletions | ~0.10% | ~0.015% | ~7× |
| Homopolymer errors | ~0.5–1% | ~0.1–0.2% | ~5× |
| Overall QV | ~QV30 | ~QV40–45 | significant |

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
