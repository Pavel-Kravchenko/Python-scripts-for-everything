---
name: bio-applied-taxonomic-profiling
description: "Taxonomic profiling of shotgun metagenomes: host decontamination, Kraken2 classification, Bracken abundance re-estimation"
tool_type: python
primary_tool: NumPy
---

# Taxonomic Profiling of Shotgun Metagenomes

## 16S vs Shotgun

| Feature | 16S amplicon | Shotgun |
|---|---|---|
| Resolution | Genus level | Species/strain level |
| Functional info | No | Yes |
| Host contamination | N/A | Major issue |
| Cost/sample | ~$30-80 | ~$200-500 |

## Host Read Removal

Tool: Bowtie2 with `--un-conc-gz` to output non-host read pairs.

Typical host contamination: gut 1-10%, skin 10-30%, BAL 40-80%, blood 60-95%.

## Kraken2 Classification

Assigns reads to LCA of all genomes sharing k-mer sequences.

```bash
kraken2 --db kraken2_db/ --paired --gzip-compressed --threads 16 \
    --confidence 0.1 --minimum-hit-groups 3 \
    --report kraken2_report.txt --output kraken2_output.txt \
    decontam_1.fastq.gz decontam_2.fastq.gz
```

**Key parameters:**
- `--confidence 0.1`: min fraction of k-mers for classification (reduces false assignments)
- `--minimum-hit-groups 3`: min minimizer groups before assignment
- Database: standard (~55 GB RAM), PlusPF (~100 GB, adds protozoa+fungi)

**Report columns:** % reads in clade, reads at exact taxon, reads in clade, rank code (S/G/F/O/C/P/K), taxonomy ID, scientific name.

## Bracken Abundance Re-estimation

Redistributes reads assigned to higher LCA nodes back to species level.

```bash
bracken -d kraken2_db/ -i kraken2_report.txt \
    -o bracken_species.txt -w bracken_species_report.txt \
    -r 150 -l S -t 10
```

| Parameter | Meaning |
|---|---|
| `-r 150` | Read length (match your data) |
| `-l S` | Level: S=species, G=genus, F=family |
| `-t 10` | Min reads threshold at species level |

Use `fraction_total_reads` column for downstream diversity/differential abundance analysis.

## Pitfalls

- **Database choice matters**: organisms not in the reference DB will be missed or misclassified
- **LCA inflation without Bracken**: reads sharing k-mers across species get pushed to genus/family level
- **Host contamination**: if not removed, host reads inflate counts and cause false taxonomic assignments
