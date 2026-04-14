---
name: bio-applied-isoform-analysis
description: Isoform analysis with long reads — Minimap2 splice alignment, bambu isoform discovery, DRIMSeq differential isoform usage.
tool_type: python
primary_tool: NumPy
---

# Isoform Analysis with Long Reads

- [FLAMES](https://github.com/LuyiTian/FLAMES)
- [IsoSeq3 (PacBio)](https://github.com/PacificBiosciences/IsoSeq)
- [bambu](https://github.com/GoekeLab/bambu)

## Technology Comparison

| Feature | Short-read (Illumina) | ONT cDNA/direct-RNA | PacBio IsoSeq (HiFi) |
|---|---|---|---|
| Read length | 75–300 bp | 1–20 kb | 1–30 kb (CCS) |
| Per-read error | ~0.1% | ~1–3% (R10.4.1) | ~0.1% (HiFi) |
| Isoform resolution | Inference required | Direct | Direct |
| Modification detection | No | Yes (direct-RNA only) | No |

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors.
- **PCR-cDNA amplification bias**: PCR amplification distorts isoform frequency ratios — avoid when input is sufficient.
- **NDR threshold (bambu)**: Default NDR=1 accepts all novel isoforms. Use NDR=0.1 for strict filtering (90% confidence a transcript is genuine). Too lenient → many false positives.
- **Isoform quantification is compositional**: Counts per gene sum to a total. Use Dirichlet-multinomial models (DRIMSeq), not simple t-tests, for differential isoform usage.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.

## Isoform Categories (bambu)

| Category | Meaning |
|---|---|
| `annotated` | Exact match to reference transcript |
| `novel_in_catalog` | New combination of known exons |
| `novel_splice_site` | New 5' or 3' splice donor/acceptor |
| `novel_exon` | Entirely new exon |
| `intergenic` | In unannotated region — usually filter out |

## Alignment

```bash
# ONT cDNA
minimap2 -ax splice --secondary=no -C5 --cs \
    hg38.fa cdna_reads.fastq.gz \
    | samtools sort -o cdna_aligned.bam -@ 8
samtools index cdna_aligned.bam

# PacBio HiFi / IsoSeq (splice:hq for high-accuracy reads)
minimap2 -ax splice:hq --secondary=no -C5 \
    hg38.fa isoseq_reads.fastq.gz \
    | samtools sort -o isoseq_aligned.bam -@ 8
samtools index isoseq_aligned.bam
```

**Key Minimap2 flags:**
- `-ax splice`: long-read spliced alignment preset
- `--secondary=no`: one alignment per read
- `-C5`: extra cost for non-canonical splice sites (GT-AG = 0, others penalized)
- `--cs`: output alignment difference string

## bambu Isoform Discovery (R)

```r
library(bambu)

annotations <- prepareAnnotations('gencode.v44.annotation.gtf')

se <- bambu(
  reads       = c('ctrl_rep1.bam', 'ctrl_rep2.bam', 'treat_rep1.bam', 'treat_rep2.bam'),
  annotations = annotations,
  genome      = 'hg38.fa',
  NDR         = 0.1   # novel discovery rate; lower = stricter
)

writeBambuOutput(se, path = 'bambu_output/')
# counts_transcript.txt: transcript-level counts (rows=transcripts, cols=samples)
# counts_gene.txt:       gene-level counts
# extended_annotations.gtf: all transcripts including novel ones
```

## Differential Isoform Usage — DRIMSeq (R)

DRIMSeq models isoform counts per gene as a **Dirichlet-multinomial distribution** (counts are compositional: they sum to the gene total).

```r
library(DRIMSeq)

counts_tx <- read.table('bambu_output/counts_transcript.txt', header=TRUE)
sample_info <- data.frame(
  sample_id = colnames(counts_tx)[-c(1,2)],
  condition = c('Control', 'Control', 'Treatment', 'Treatment')
)

d <- dmDSdata(counts = counts_tx, samples = sample_info)

# Filter: require ≥2 isoforms and adequate expression
d <- dmFilter(d,
  min_samps_gene_expr    = 2,
  min_samps_feature_expr = 2,
  min_gene_expr          = 10,
  min_feature_expr       = 5
)

d <- dmPrecision(d)
d <- dmFit(d)
d <- dmTest(d, coef = 'conditionTreatment')

res_gene <- results(d, level = 'gene')    # omnibus test
res_tx   <- results(d, level = 'feature') # per-transcript
sig      <- res_gene[res_gene$adj_pvalue < 0.05, ]
```

## QC Metrics to Check

- Mapping rate (expect >90% for good library)
- Multi-exon read fraction
- NFR (<150 bp) fragment proportion
- Isoform category breakdown (high `intergenic` fraction = quality concern)
- D50: top N clones covering 50% of reads (lower = more diversity in expression)
