---
name: bio-applied-variant-surveillance
description: SARS-CoV-2 lineage classification (Pango), Freyja wastewater deconvolution, spike protein mutation tracking, and surveillance pipeline tools
tool_type: python
primary_tool: Matplotlib
---

# Variant Surveillance and Wastewater Epidemiology

## Key Tools

| Tool | Purpose | Output |
|------|---------|--------|
| Pangolin | Lineage assignment from consensus | CSV with lineage + confidence |
| Nextclade | Clade + QC + amino acid changes | TSV with quality flags |
| Freyja | Wastewater lineage deconvolution | Lineage fractions per sample |
| GISAID / outbreak.info | Lineage frequency tracking | Web API + download |

## Pango Nomenclature

- Root: A, B (early Wuhan)
- Sub-lineages: B.1 → B.1.1 → B.1.1.7 (Alpha)
- After 3 sub-levels: alias issued (BA, BQ, XBB) or recombinant prefix X
- WHO categories: VOC (Concern) > VOI (Interest) > VUM (Monitoring)

```bash
# Lineage assignment
pangolin sequences.fasta --outfile lineages.csv --threads 4

# Clade + quality flags
nextclade run --input-fasta sequences.fasta \
    --input-dataset sars-cov-2 --output-tsv nextclade.tsv
```

## Spike Protein Key Regions

| Region | Positions | Mutation effects |
|--------|-----------|-----------------|
| NTD | 13–305 | Antibody binding |
| RBD | 319–541 | ACE2 contact; immune escape (E484K, K417N, F486V); affinity (N501Y) |
| Furin cleavage | 681–685 | P681H/R increases fitness |

## Freyja Wastewater Pipeline

```bash
# 1. Call variants from wastewater BAM
freyja variants wastewater.bam --variants variants.tsv --depths depths.tsv

# 2. Deconvolve lineages (uses UShER barcodes)
freyja demix variants.tsv depths.tsv --output lineages.csv

# 3. Aggregate time-series
freyja aggregate --inputdir ./samples/ --output aggregated.tsv

# 4. Plot
freyja plot aggregated.tsv --output lineage_plot.pdf
```

**Freyja model**: linear mixture model solving for lineage fractions given observed variant frequencies and a barcode matrix of lineage-defining mutations.

## Wastewater Epidemiology Facts

- Detects ~1 infected person per 100,000 catchment population
- Provides 3–7 day **lead time** before clinical case surge
- Population-level signal unaffected by clinical testing access gaps
- Signal = viral RNA copies/L (log-scale); smooth before correlation analysis

## Python Patterns

```python
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter1d

# Smooth wastewater signal before correlating with clinical cases
ww_smooth = gaussian_filter1d(viral_load, sigma=3)

# Lead-time correlation (WW leads cases by ~5 days/weeks)
r, p = pearsonr(ww_smooth[:-lag], cases[lag:])

# Lineage frequency normalization (stacked area plot)
freq_matrix = np.array(list(lineage_freqs.values()))
freq_norm = freq_matrix / freq_matrix.sum(axis=0) * 100

# Dominant lineage per time point
dominant_idx = np.argmax(freq_norm, axis=0)
```

## Pitfalls

- **Amplicon dropout**: primer mismatches in divergent lineages → missing mutations → incorrect lineage calls; update primer schemes with each major variant wave
- **Freyja barcode staleness**: lineage barcodes must be updated regularly (`freyja update`); stale barcodes misclassify new sublineages as "Other"
- **Wastewater normalization**: viral load varies with flow rate and precipitation — normalize to pepper mild mottle virus (PMMoV) or crAssphage as fecal indicator
- **Consensus vs quasispecies**: Pangolin assigns one lineage to a consensus; co-infections or recombinants require manual inspection of per-site variant frequencies
- **Coordinate systems**: BED is 0-based; VCF/GFF are 1-based — off-by-one errors when mixing
- **Multiple testing**: apply Benjamini-Hochberg FDR when testing associations across many positions/lineages
