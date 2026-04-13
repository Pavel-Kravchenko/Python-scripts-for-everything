---
name: bio-applied-chipseq-pipeline
description: "**Tier 3 — Applied Bioinformatics | Module 24 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/24_ChIP_seq_Epigenomics/01_chipseq_pipeline.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# ChIP-seq Processing Pipeline

*Source: Course notebook `Tier_3_Applied_Bioinformatics/24_ChIP_seq_Epigenomics/01_chipseq_pipeline.ipynb`*

# ChIP-seq Processing Pipeline

**Tier 3 — Applied Bioinformatics | Module 24 · Notebook 1**

*Prerequisites: Module 01 (NGS Fundamentals), Module 23 (TF Footprinting & ATAC-seq)*

---

**By the end of this notebook you will be able to:**
1. Describe the ChIP-seq experimental workflow and the role of Input/IgG controls
2. Run read trimming, alignment, and duplicate removal on ChIP-seq data
3. Call peaks with MACS3 for transcription factors and histone marks
4. Compute FRiP score and cross-correlation quality metrics
5. Visualize signal enrichment with deepTools heatmaps and profile plots

**Key resources:**
- [ENCODE ChIP-seq pipeline](https://www.encodeproject.org/chip-seq/)
- [Harvard HBC ChIP-seq training](https://hbctraining.github.io/Intro-ChIP-seq/)
- [MACS3 documentation](https://macs3-project.github.io/MACS/)
- [deepTools documentation](https://deeptools.readthedocs.io/)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Optional dependencies — caught gracefully
try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False

try:
    import matplotlib_venn as venn_mod
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11
print('Setup complete.')
```

## 1. ChIP-seq Experimental Design

Chromatin immunoprecipitation followed by sequencing (ChIP-seq) maps the genome-wide binding sites of transcription factors (TFs) or the locations of histone modifications. The core workflow involves crosslinking proteins to DNA, sonicating chromatin into short fragments, immunoprecipitating with an antibody specific to the protein or histone mark of interest, reversing crosslinks, and sequencing the recovered DNA. A matching **Input** (or IgG control) sample undergoes the same protocol minus the immunoprecipitation step, capturing background chromatin accessibility bias without any enrichment. The ratio of ChIP signal to Input signal forms the basis of all downstream peak calling.

| Experiment Type | Peak Shape | Recommended Reads | Example Marks |
|---|---|---|---|
| Transcription Factor (TF) | Narrow (< 500 bp) | 20–40 M | CTCF, GATA1, p53 |
| Active histone mark | Narrow | 30–50 M | H3K4me3, H3K9ac |
| Broad histone mark | Broad (> 5 kb) | 40–80 M | H3K27me3, H3K9me3 |
| Enhancer mark | Narrow/Mixed | 30–50 M | H3K27ac, H3K4me1 |

Antibody quality is the single most important variable in ChIP-seq. Only **ChIP-grade** antibodies (tested by the supplier for ChIP, ideally with a knockout validation) should be used; polyclonal antibodies validated for western blot alone are insufficient. For sequencing strategy, **paired-end** reads are strongly preferred over single-end because they allow accurate fragment-size estimation and more reliable removal of PCR duplicates (duplicate pairs must share both 5′ coordinates). ENCODE requires **≥ 2 biological replicates** for all ChIP-seq experiments to enable assessment of reproducibility via IDR (Irreproducibility Discovery Rate). Pooled pseudoreplicate analysis is also recommended to benchmark self-consistency.

## 2. Read QC and Trimming

**FastQC** generates per-read quality reports for each FASTQ file. The metrics most relevant to ChIP-seq are: (1) per-base quality scores — a drop below Q20 in the last 10 bp of 3′ ends is common and tolerable; (2) adapter content — high adapter content (> 20%) indicates short fragments and must be trimmed before alignment; (3) sequence duplication level — elevated duplication at the FastQC stage reflects library complexity, though most duplicates are removed by Picard later; (4) GC content distribution — an unexpected bimodal GC peak suggests contamination or strong PCR bias. **MultiQC** aggregates FastQC reports across all samples into a single interactive HTML report, making it easy to spot outliers.

**Trim Galore** wraps Cutadapt and automatically detects the most common Illumina adapter sequences (TruSeq, Nextera). For paired-end data use `--paired`; it will trim both mates and discard read pairs where either mate is shorter than the minimum length after trimming. The `--fastqc` flag reruns FastQC after trimming so both pre- and post-trim reports are available immediately. **fastp** is a faster C++-based alternative that handles adapter detection, quality filtering, and per-sample HTML reports in a single tool call; it is increasingly preferred for high-throughput pipelines.

Warning signs to act on before proceeding: adapter content > 20% means fragment sizes are very short (re-evaluate sonication); per-base quality dropping below Q20 at 3′ end beyond 10 bp may benefit from aggressive trimming (`--quality 25`); an unexpected secondary peak in GC content distribution often signals PCR amplification of a small subset of fragments (library complexity issue). If > 30% of read pairs are discarded after trimming, the library likely needs to be re-prepared.

```python
# Note: requires bioinformatics tools installed (conda environment)

# --- Step 1: Quality control with FastQC ---
# Run FastQC on raw FASTQ files (paired-end)
!fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc/raw/ -t 4

# Aggregate reports with MultiQC
!multiqc qc/raw/ -o qc/raw/multiqc/

# --- Step 2: Adapter trimming with Trim Galore ---
# Automatically detects adapters; --fastqc re-runs QC after trimming
!trim_galore \
    --paired \
    --fastqc \
    --cores 4 \
    --quality 20 \
    --length 20 \
    --output_dir trimmed/ \
    sample_R1.fastq.gz sample_R2.fastq.gz

# Alternative: fastp (faster, similar features)
!fastp \
    -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o trimmed/sample_R1_trimmed.fastq.gz \
    -O trimmed/sample_R2_trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --thread 4 \
    -h trimmed/fastp_report.html
```

## 3. Alignment with Bowtie2

**Bowtie2** is the standard short-read aligner for ChIP-seq. It builds an FM-index from the reference genome, enabling fast gapped alignment of short reads (typically 50–150 bp). For human experiments use the **hg38** reference assembly; for mouse use mm10 or mm39. After alignment, reads mapping to ENCODE blacklisted regions — repetitive elements, centromeres, and telomeres that produce artifactually high signal in nearly all ChIP-seq experiments — must be removed using `bedtools intersect -v`. The hg38 blacklist is available from the ENCODE portal.

A good ChIP-seq experiment should achieve **> 80% overall alignment rate** for trimmed paired-end reads against hg38. `samtools flagstat` reports total reads, mapped, properly paired, and duplicate counts at a glance. If alignment rate falls below 70%, first check the trimming log (adapter content), then verify the correct genome version was used. Low mapping can also reflect cross-species contamination or heavily degraded DNA.

Post-alignment filtering is critical for reducing noise: remove unmapped reads (`-F 4`), secondary alignments (`-F 256`), and low-confidence mappings with MAPQ < 30 (`-q 30`). Keep only properly paired reads (`-f 2`) to ensure both mates in a pair mapped concordantly. Combining these flags: `samtools view -F 4 -F 256 -q 30 -f 2 -b`. For paired-end data, `--no-mixed` and `--no-discordant` flags in Bowtie2 itself prevent mixed/discordant alignments from entering the BAM.

```python
# Note: requires bioinformatics tools installed (conda environment)

# --- Download or build Bowtie2 index (hg38) ---
# If not already indexed:
# !bowtie2-build hg38.fa hg38_index/hg38

# --- Align paired-end trimmed reads ---
!bowtie2 \
    -x /data/indices/hg38 \
    -1 trimmed/sample_R1_val_1.fq.gz \
    -2 trimmed/sample_R2_val_2.fq.gz \
    -p 8 \
    --no-mixed \
    --no-discordant \
    2> logs/sample_bowtie2.log | \
    samtools view -bS -F 4 -F 256 -q 30 -f 2 | \
    samtools sort -@ 4 -o aligned/sample.bam

# --- Index BAM ---
!samtools index aligned/sample.bam

# --- QC: Alignment statistics ---
!samtools flagstat aligned/sample.bam | tee logs/sample_flagstat.txt

# --- Remove blacklisted regions (ENCODE) ---
!bedtools intersect \
    -a aligned/sample.bam \
    -b hg38-blacklist.v2.bed \
    -v > aligned/sample_filtered.bam
!samtools index aligned/sample_filtered.bam
```

## 4. Duplicate Removal with Picard

Duplicate reads in ChIP-seq arise from two distinct sources. **PCR duplicates** are identical fragments amplified multiple times during library preparation; they reflect low library complexity and disproportionately inflate peak signal at highly enriched loci, creating false peaks. **Optical duplicates** are reads from physically proximate clusters on the flow cell that are mis-called as separate molecules; they are a technical artifact of the sequencer and are more prevalent on newer patterned flow cells (NovaSeq). In high-quality ChIP-seq experiments a total duplication rate above **50%** is a warning sign of library complexity problems; for TF ChIP the rate should ideally be < 30%.

**Picard MarkDuplicates** identifies duplicates by comparing the 5′ mapping coordinates of read pairs. Setting `REMOVE_DUPLICATES=true` physically removes duplicate reads from the output BAM rather than just flagging them. The `OPTICAL_DUPLICATE_PIXEL_DISTANCE` parameter should be set appropriately for the sequencer: **100** for HiSeq (non-patterned flow cell) and **2500** for NovaSeq (patterned flow cell). The output metrics file contains the `PERCENT_DUPLICATION` value and the estimated library size.

ENCODE uses two complementary complexity metrics computed from the duplication metrics file. The **PCR Bottleneck Coefficient (PBC1)** is the ratio of genomic positions with exactly one unique read to positions with at least one read; values > 0.7 are acceptable, > 0.9 are ideal. The **Non-Redundant Fraction (NRF)** is the number of distinct alignments divided by total alignments; NRF > 0.8 is acceptable, > 0.9 is ideal. After deduplication, re-run `samtools flagstat` to compare pre- and post-dedup read counts.

```python
# Note: requires bioinformatics tools installed (conda environment)

# --- Mark and remove PCR/optical duplicates with Picard ---
!picard MarkDuplicates \
    I=aligned/sample_filtered.bam \
    O=dedup/sample_dedup.bam \
    M=dedup/sample_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=./ \
    2> logs/picard_dedup.log

!samtools index dedup/sample_dedup.bam

# --- Compare flagstat before and after deduplication ---
!echo '=== BEFORE deduplication ===' && samtools flagstat aligned/sample_filtered.bam
!echo '=== AFTER deduplication ===' && samtools flagstat dedup/sample_dedup.bam
```

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Simulate duplication metrics (as would be output by Picard)
# In practice, parse dedup/sample_dup_metrics.txt
np.random.seed(42)

metrics_data = {
    'LIBRARY': ['ChIP_sample', 'Input_sample'],
    'UNPAIRED_READS_EXAMINED': [0, 0],
    'READ_PAIRS_EXAMINED': [35_420_000, 28_100_000],
    'UNMAPPED_READS': [1_200_000, 850_000],
    'UNPAIRED_READ_DUPLICATES': [0, 0],
    'READ_PAIR_DUPLICATES': [8_500_000, 3_200_000],
    'READ_PAIR_OPTICAL_DUPLICATES': [420_000, 180_000],
    'PERCENT_DUPLICATION': [0.240, 0.114],
    'ESTIMATED_LIBRARY_SIZE': [4_800_000, 11_200_000],
}
df_metrics = pd.DataFrame(metrics_data)

df_metrics['PCR_DUPS'] = (
    df_metrics['READ_PAIR_DUPLICATES'] - df_metrics['READ_PAIR_OPTICAL_DUPLICATES']
)
df_metrics['NRF'] = 1 - df_metrics['PERCENT_DUPLICATION']  # simplified approximation

print('Duplication Metrics Summary')
print('=' * 50)
cols = ['LIBRARY', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION',
        'NRF', 'ESTIMATED_LIBRARY_SIZE']
print(df_metrics[cols].to_string(index=False))
print('\nENCODE thresholds: NRF > 0.8 (acceptable), > 0.9 (ideal)')

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
libraries = df_metrics['LIBRARY']
dup_pct = df_metrics['PERCENT_DUPLICATION'] * 100

axes[0].bar(libraries, dup_pct, color=['steelblue', 'coral'])
axes[0].axhline(20, color='red', linestyle='--', label='20% threshold')
axes[0].set_ylabel('Duplication Rate (%)')
axes[0].set_title('PCR Duplication Rate')
axes[0].legend()

# Stacked: optical vs PCR duplicates
optical = df_metrics['READ_PAIR_OPTICAL_DUPLICATES'] / 1e6
pcr = df_metrics['PCR_DUPS'] / 1e6
x = np.arange(len(libraries))
axes[1].bar(x, optical, label='Optical duplicates', color='tomato')
axes[1].bar(x, pcr, bottom=optical, label='PCR duplicates', color='steelblue')
axes[1].set_xticks(x)
axes[1].set_xticklabels(libraries)
axes[1].set_ylabel('Read pairs (millions)')
axes[1].set_title('Duplicate Breakdown')
axes[1].legend()

plt.tight_layout()
plt.show()
```
