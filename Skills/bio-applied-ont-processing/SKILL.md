---
name: bio-applied-ont-processing
description: "**Tier 3 — Applied Bioinformatics | Module 25 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/25_Long_Read_Sequencing/01_ont_processing.ipynb"
---

# ONT Data Processing

*Source: Course notebook `Tier_3_Applied_Bioinformatics/25_Long_Read_Sequencing/01_ont_processing.ipynb`*

# ONT Data Processing

**Tier 3 — Applied Bioinformatics | Module 25 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 17 (Genome Assembly)*

---

**By the end of this notebook you will be able to:**
1. Explain ONT sequencing chemistry and FAST5/POD5 raw signal format
2. Run basecalling with Dorado and assess read quality with NanoStat/NanoPlot
3. Align long reads to a reference with Minimap2
4. Detect base modifications (5mC methylation) from ONT signal
5. Compare long-read vs short-read alignment statistics

**Key resources:**
- [Oxford Nanopore community tutorials](https://community.nanoporetech.com/)
- [Awesome Nanopore training](https://github.com/GenomicsAotearoa/Genomics-Aotearoa-Nanopore-training)
- [Minimap2 documentation](https://github.com/lh3/minimap2)

## 1. ONT Technology Overview

Oxford Nanopore Technology (ONT) sequences DNA by threading a single strand through a protein nanopore embedded in an electrically resistive membrane. A constant voltage drives ionic current through the pore; as each k-mer of bases occupies the constriction zone, it causes a characteristic disruption to that current. The resulting time-series of picoampere values — often called a "squiggle" — encodes the sequence. The R9.4.1 chemistry uses a 5-mer sensing region, while the newer R10.4.1 dual-reader pore uses a 9-mer sensing window, dramatically improving homopolymer resolution and raw-read accuracy.

**FAST5 format** is an HDF5-based container that stores raw signal traces (picoampere arrays), channel metadata, and — if basecalling was run — the basecalled sequence and quality scores. Originally each FAST5 file held data from a single channel/read. Multi-read FAST5 was later introduced but I/O performance remained a bottleneck. FAST5 is now considered legacy.

**POD5 format** is Oxford Nanopore's modern, purpose-built binary format that replaces FAST5. It uses Apache Arrow for columnar on-disk storage, enabling much faster random-access I/O and smaller file sizes. Dorado natively reads POD5. ONT provides `pod5 convert fast5` to migrate legacy data. One POD5 file typically contains reads from an entire batch (not per-channel), simplifying data management.

| Feature | ONT R9.4.1 | ONT R10.4.1 | PacBio HiFi (Revio) |
|---|---|---|---|
| Read length | 50 bp – 4 Mb | 50 bp – 4 Mb | 10–25 kb (CCS) |
| Modal read length | ~8–12 kb | ~10–20 kb | ~15–18 kb |
| Raw accuracy | ~95% | ~97–99% | ~99.9% (CCS) |
| Consensus accuracy | 99.5%+ | 99.9%+ | 99.9%+ |
| Throughput per flow cell | ~30–50 Gb | ~50–120 Gb | ~90 Gb (Revio) |
| Native 5mC detection | Yes (retrained model) | Yes (dual-base calling) | No (requires bisulfite or 6mA) |
| Instrument | MinION / PromethION | MinION R10 / P2 / PromethION | Revio |
| Typical N50 (human WGS) | 15–25 kb | 20–40 kb | 15–20 kb |

## 2. Basecalling with Dorado

**Dorado** is ONT's current GPU-accelerated basecaller, replacing the older Guppy. It uses a convolutional + recurrent neural network trained to translate raw squiggle signals into nucleotide sequences and per-base quality scores. Dorado is open-source and supports both NVIDIA CUDA and Apple Silicon (Metal) backends, making it usable on workstations as well as HPC clusters.

Dorado offers three model accuracy tiers for each chemistry. The **`fast`** model runs at maximum throughput (~95–97% accuracy) and is suitable for rapid screening or real-time adaptive sampling. The **`hac`** (high-accuracy) model provides a good balance of speed and quality (~97–99%) and is the standard choice for most genomic applications. The **`sup`** (super-accuracy) model achieves ~99%+ accuracy but is 5–10× slower than `hac`, making it best suited for post-run reprocessing of critical samples on a GPU cluster.

Output from Dorado is an **unaligned BAM** by default — this is the recommended format because it preserves all per-read metadata (move tables, basecalling model version, sequencing summary statistics) as BAM tags. FASTQ can be obtained by piping through `samtools fastq`. For modification-aware basecalling, append the modification code to the model name (e.g. `hac,5mCG_5hmCG`); probabilities are stored in the MM and ML BAM tags per the SAM specification. **Duplex mode** pairs the template and complement strands from the same DNA molecule, achieving ~Q30 (99.9%) accuracy on the duplex fraction of reads.

```bash
# Quick reference: Dorado model tiers
# fast  → ~95-97% accuracy, highest throughput
# hac   → ~97-99% accuracy, standard choice
# sup   → ~99%+  accuracy, 5-10x slower than hac
```

```python
# Dorado basecalling — simplex HAC model
# !dorado basecaller hac pod5_data/ > calls.bam
# Explanation: 'hac' = high-accuracy model; output is unaligned BAM

# Convert to FASTQ if needed
# !samtools fastq calls.bam | gzip > calls.fastq.gz

# Duplex basecalling (pairs complementary strands, ~Q30)
# !dorado duplex hac pod5_data/ > calls_duplex.bam

# Methylation-aware: add modification tag to model name
# !dorado basecaller hac,5mCG_5hmCG pod5_data/ > calls_modcall.bam

# Check basecalling summary stats
# !dorado summary calls.bam | head -5

# List available models (requires internet or local model cache)
# !dorado download --list
```

## 3. Read Quality Assessment

**NanoStat** produces a concise text summary of key run metrics: total reads, total bases sequenced, mean and median read length, N50 read length (the length at which 50% of all sequenced bases are in reads of that length or longer), mean and median Phred quality score, and the percentage of reads passing Q10, Q15, and Q20 thresholds. This is the first QC step after basecalling — a healthy R10.4.1 run typically shows N50 > 15 kb and median quality > Q15.

**NanoPlot** generates a suite of interactive HTML visualizations including a read-length histogram, a quality-score distribution, a yield-over-time plot, and a length-vs-quality scatter ("dot") plot. The dot plot is especially useful for identifying distinct populations (e.g. short adapter-only fragments, ultra-long reads) and for comparing runs. NanoPlot can accept FASTQ, BAM, or a sequencing summary CSV as input.

**NanoFilt** is a streaming FASTQ filter. Typical parameters for whole-genome sequencing assembly are `-q 8 -l 3000` (retain reads ≥ Q8 and ≥ 3 kb). For variant calling where accuracy matters more than depth, `-q 10 -l 1000` removes low-quality and very short reads while retaining the majority of bases. Quality thresholds on the Phred scale: **Q10** = 90% accuracy per base, **Q15** = 96.8%, **Q20** = 99%.

| Quality | Per-base accuracy | Typical ONT chemistry |
|---|---|---|
| Q10 | 90.0% | R9.4.1 minimum |
| Q12 | 93.7% | R9.4.1 median |
| Q15 | 96.8% | R10.4.1 median |
| Q20 | 99.0% | R10.4.1 sup / duplex |
| Q30 | 99.9% | ONT duplex / PacBio CCS |

```python
# Run NanoStat on raw FASTQ
# !NanoStat --fastq calls.fastq.gz --outdir nanostat_out/ --threads 4

# NanoPlot (produces interactive HTML)
# !NanoPlot --fastq calls.fastq.gz --outdir nanoplot_out/ --plots dot --N50

# Filter reads: min quality Q10, min length 1000 bp
# !NanoFilt -q 10 -l 1000 calls.fastq.gz | gzip > calls_filtered.fastq.gz

# Check number of reads before/after filtering
# !echo "Before: $(zcat calls.fastq.gz | wc -l | awk '{print $1/4}') reads"
# !echo "After:  $(zcat calls_filtered.fastq.gz | wc -l | awk '{print $1/4}') reads"
```

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

# Simulate ONT R10.4.1 read length distribution (log-normal) and quality
n_reads = 5000
# Read lengths: log-normal with mean ~15 kb
lengths = rng.lognormal(mean=9.6, sigma=0.9, size=n_reads).astype(int)
lengths = np.clip(lengths, 200, 500_000)

# Quality scores: approximately normal around Q16 for R10.4.1
mean_q = rng.normal(16.5, 2.5, size=n_reads)
mean_q = np.clip(mean_q, 5, 30)

df_reads = pd.DataFrame({"length_bp": lengths, "mean_q": mean_q})

# --- compute NanoStat-style summary ---
sorted_len = np.sort(lengths)[::-1]
cumsum = np.cumsum(sorted_len)
n50_val = sorted_len[np.searchsorted(cumsum, cumsum[-1] / 2)]

print("=== Simulated NanoStat Summary ===")
print(f"Total reads:          {n_reads:>10,}")
print(f"Total bases (Gb):     {lengths.sum()/1e9:>10.2f}")
print(f"Mean read length:     {lengths.mean():>10,.0f} bp")
print(f"Median read length:   {np.median(lengths):>10,.0f} bp")
print(f"N50 read length:      {n50_val:>10,} bp")
print(f"Mean read quality:    {mean_q.mean():>10.1f}")
print(f"Reads > Q10:          {(mean_q >= 10).mean()*100:>9.1f}%")
print(f"Reads > Q15:          {(mean_q >= 15).mean()*100:>9.1f}%")
print(f"Reads > Q20:          {(mean_q >= 20).mean()*100:>9.1f}%")

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# Panel 1: read length histogram (log scale)
axes[0].hist(lengths / 1000, bins=60, color="steelblue", edgecolor="white", linewidth=0.3)
axes[0].axvline(n50_val / 1000, color="crimson", linestyle="--", label=f"N50 = {n50_val/1000:.1f} kb")
axes[0].set_xlabel("Read length (kb)")
axes[0].set_ylabel("Count")
axes[0].set_title("Read length distribution")
axes[0].set_xscale("log")
axes[0].legend()

# Panel 2: quality distribution
axes[1].hist(mean_q, bins=40, color="darkorange", edgecolor="white", linewidth=0.3)
axes[1].axvline(10, color="red", linestyle="--", label="Q10 threshold")
axes[1].axvline(15, color="green", linestyle="--", label="Q15 threshold")
axes[1].set_xlabel("Mean read quality (Phred)")
axes[1].set_ylabel("Count")
axes[1].set_title("Quality score distribution")
axes[1].legend()

# Panel 3: length vs quality scatter
sc = axes[2].scatter(lengths / 1000, mean_q, alpha=0.05, s=4, c=mean_q, cmap="RdYlGn")
axes[2].set_xlabel("Read length (kb)")
axes[2].set_ylabel("Mean quality")
axes[2].set_title("Length vs. quality")
axes[2].set_xscale("log")
plt.colorbar(sc, ax=axes[2], label="Q score")

plt.suptitle("ONT R10.4.1 — simulated read QC (5,000 reads)", fontweight="bold")
plt.tight_layout()
plt.show()
```

## 4. Alignment with Minimap2

**Minimap2** is the standard aligner for long reads. Rather than aligning every position like short-read aligners (BWA, Bowtie2), Minimap2 uses **minimizer-based seeding**: it extracts a sparse set of representative k-mers (minimizers) from both query and reference, finds collinear chains of matching minimizers, then performs banded dynamic programming only in the seed-chain regions. This makes alignment of 10–100 kb reads orders of magnitude faster than classical approaches.

**Alignment presets** (`-ax` flag) configure Minimap2's parameters for different data types. `map-ont` is optimised for Oxford Nanopore genomic reads (higher mismatch tolerance, no splice-aware scoring). `map-hifi` tunes for PacBio HiFi's higher accuracy and read length distribution. `splice` enables intron-aware alignment for RNA/cDNA reads, recognising the N CIGAR operation for intron skipping. `ava-ont` is used for all-vs-all overlap computation during de novo assembly. The preset can make a 5–10× difference in sensitivity and specificity.

Long reads typically achieve **95–99% alignment rate** for genomic DNA sequencing when the sample is the same species as the reference. The supplementary alignment rate (reads split across two loci — a signature of structural variants) is ~1–5%, compared with <0.1% for short reads. After alignment, always sort with `samtools sort` and index with `samtools index` for efficient random access. Use `samtools flagstat` for a quick summary of mapping statistics.

```
Key samtools flagstat fields for long-read BAMs:
  - total reads            → all primary + supplementary + secondary
  - mapped (%)             → primary mapped reads
  - supplementary          → chimeric/split reads (SV signatures)
  - paired-in-sequencing   → always 0 for long reads (single-end)
```

```python
# Align ONT reads (map-ont preset)
# !minimap2 -ax map-ont -t 8 hg38.fa calls_filtered.fastq.gz \
#     | samtools sort -o ont_aligned.bam -@ 8
# !samtools index ont_aligned.bam

# Check alignment statistics
# !samtools flagstat ont_aligned.bam

# For PacBio HiFi reads use map-hifi preset:
# !minimap2 -ax map-hifi -t 8 hg38.fa hifi_reads.fastq.gz \
#     | samtools sort -o hifi_aligned.bam
```
