---
name: bio-applied-assembly-binning
description: "Metagenomic assembly with MEGAHIT, contig binning with MetaBAT2, and MAG quality assessment with CheckM. Includes binning signals, multi-sample strategy, and MIMAG quality tiers."
tool_type: python
primary_tool: NumPy
---

# Metagenomic Assembly, Binning, and MAGs

## References
- [MEGAHIT](https://github.com/voutcn/megahit)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat)
- [CheckM](https://github.com/Ecogenomics/CheckM)

## Pitfalls

- **Discard short contigs before binning:** Contigs < 1500 bp carry insufficient tetranucleotide signal for reliable binning. Use `--minContig 1500` in MetaBAT2.
- **Multi-sample binning dramatically improves quality:** Mapping all samples to the same assembly and providing all BAMs to `jgi_summarize_bam_contig_depths` adds differential abundance as a binning signal — do this whenever you have ≥2 related samples.
- **Strain heterogeneity ≠ contamination:** High CheckM contamination with high strain heterogeneity (≥90% AA identity between duplicates) indicates strain mixing, not foreign contamination.
- **Run all three binners, then use DAS_Tool:** MetaBAT2 alone leaves quality on the table. DAS_Tool dereplicated bins from MetaBAT2 + CONCOCT + MaxBin2 consistently outperforms any single binner.
- **MEGAHIT vs metaSPAdes trade-off:** MEGAHIT is RAM-efficient and fast; metaSPAdes is more accurate for low-coverage organisms but requires 100–400 GB RAM for complex communities.
- **N50 expectations:** Typical gut metagenome (5 Gbp data) → 100k–500k contigs, N50 of 1–10 kb.

## Workflow Overview

```
reads → MEGAHIT assembly → filter ≥500 bp contigs
      → map reads back → jgi_summarize_bam_contig_depths
      → MetaBAT2/CONCOCT/MaxBin2 binning
      → DAS_Tool dereplication
      → CheckM quality assessment
      → filter HQ MAGs (≥90% complete, <5% contamination)
```

## MEGAHIT Assembly

```bash
megahit \
    -1 decontam_1.fastq.gz \
    -2 decontam_2.fastq.gz \
    -o megahit_assembly/ \
    --min-contig-len 500 \
    -t 16 \
    -m 0.5              # use at most 50% of available RAM

# Check stats
seqkit stats megahit_assembly/final.contigs.fa
```

**Key parameter:** `--k-list 21,29,39,59,79,99,119,141` (default). Multi-k handles extreme coverage variation (1x rare organism vs 1000x abundant organism in the same community).

## Coverage Estimation

```bash
# Index assembly and map reads back
bowtie2-build megahit_assembly/final.contigs.fa contigs_index

bowtie2 -x contigs_index -1 decontam_1.fastq.gz -2 decontam_2.fastq.gz \
    -p 16 --no-unal | samtools sort -@ 8 -o contigs_mapped.bam
samtools index contigs_mapped.bam

# Generate depth file (required by MetaBAT2)
jgi_summarize_bam_contig_depths \
    --outputDepth depths.txt \
    contigs_mapped.bam
# Output columns: contigName, contigLen, totalAvgDepth, sampleDepth, sampleDepthVar
```

## Contig Binning

```bash
# MetaBAT2
metabat2 \
    -i megahit_assembly/final.contigs.fa \
    -a depths.txt \
    -o bins/bin \
    --minContig 1500 \
    --minClsSize 100000 \   # minimum 100 kb bin size
    -t 8 \
    --saveCls

# Run all three binners then combine with DAS_Tool
das_tool \
    -i metabat2_bins,concoct_bins,maxbin2_bins \
    -l MetaBAT2,CONCOCT,MaxBin2 \
    -c megahit_assembly/final.contigs.fa \
    -o dastool_output/ \
    -t 8 \
    --write_bins
```

**Binning signals used by MetaBAT2:**
1. Tetranucleotide frequency (TNF) — 256-dimensional composition vector, genome-specific
2. Coverage depth across samples — organism-specific abundance

## MAG Quality Assessment (CheckM)

```bash
checkm lineage_wf bins/ checkm_out/ -x fa -t 16 --pplacer_threads 4
checkm qa checkm_out/lineage.ms checkm_out/ -o 2 > checkm_summary.txt
```

### MIMAG Quality Tiers (Bowers et al. 2017, *Nature Biotechnology*)

| Tier | Completeness | Contamination |
|---|---|---|
| High quality | ≥ 90% | < 5% |
| Medium quality | ≥ 50% | < 10% |
| Low quality | < 50% | — |

CheckM uses ~100 conserved single-copy marker genes expected exactly once in any complete bacterial genome. Marker gene set is selected based on phylogenetic placement of each bin.

## Binner Comparison

| Binner | Signal used | Strengths |
|---|---|---|
| MetaBAT2 | TNF + coverage | Fast, widely used |
| CONCOCT | TNF + coverage (PCA) | Good for low-coverage data |
| MaxBin2 | Marker gene abundance + coverage | Marker-guided |
| DAS_Tool | Combines all above | Best final MAG set |
