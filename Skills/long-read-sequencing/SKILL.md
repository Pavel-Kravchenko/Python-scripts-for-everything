---
name: long-read-sequencing
description: Oxford Nanopore and PacBio long-read sequencing — basecalling (Dorado), QC (NanoStat), alignment (Minimap2), assembly (Flye, Hifiasm), SV calling (Sniffles2), methylation, isoform analysis
primary_tool: Pandas
---

# Long-Read Sequencing

## When to Use

Use this skill when:
- Processing Oxford Nanopore (ONT) or PacBio HiFi reads
- Performing de novo genome or metagenome assembly
- Detecting structural variants (deletions, insertions, inversions, duplications)
- Analyzing CpG methylation from ONT signal
- Performing full-length transcript isoform analysis

## Quick Reference

| Task | Tool | Key Command |
|------|------|-------------|
| ONT basecalling | Dorado | `dorado basecaller hac pod5/` |
| Read QC | NanoStat | `NanoStat --fastq reads.fq.gz` |
| Read filtering | NanoFilt | `NanoFilt -q 10 -l 1000` |
| Alignment | Minimap2 | `minimap2 -ax map-ont ref.fa reads.fq.gz` |
| HiFi alignment | Minimap2 | `minimap2 -ax map-hifi ref.fa reads.fq.gz` |
| cDNA alignment | Minimap2 | `minimap2 -ax splice reads.fq.gz` |
| De novo assembly (ONT) | Flye | `flye --nano-hq reads.fq.gz --genome-size 5m` |
| De novo assembly (HiFi) | Hifiasm | `hifiasm -o out.asm -t 16 reads.fq.gz` |
| Assembly polishing | Medaka | `medaka_consensus -i reads.fq.gz -d assembly.fa` |
| Assembly QC | QUAST | `quast.py assembly.fa -r ref.fa` |
| SV calling | Sniffles2 | `sniffles --input aligned.bam --vcf svs.vcf` |
| Methylation | modkit | `modkit pileup aligned.bam meth.bed --cpg` |
| Isoform analysis | bambu (R) | `bambu(reads='aln.bam', annotations=gtf)` |

## Key Patterns

**Pattern 1: ONT processing pipeline**
```bash
# Basecalling with Dorado (high accuracy model)
dorado basecaller hac pod5_data/ | samtools fastq > reads.fastq.gz

# QC
NanoStat --fastq reads.fastq.gz --outdir nanostat/
NanoFilt -q 10 -l 1000 reads.fastq.gz > reads_filtered.fastq.gz

# Alignment
minimap2 -ax map-ont hg38.fa reads_filtered.fastq.gz | \
    samtools sort -o aligned.bam && samtools index aligned.bam
```

**Pattern 2: De novo assembly**
```bash
# Flye for ONT genome assembly
flye --nano-hq reads_filtered.fastq.gz --genome-size 3g \
    --out-dir flye_out/ --threads 16

# Hifiasm for PacBio HiFi
hifiasm -o sample.asm -t 16 hifi_reads.fastq.gz
awk '/^S/{print ">"$2"\n"$3}' sample.asm.bp.p_ctg.gfa > assembly.fasta

# Polish ONT assembly with Medaka
medaka_consensus -i reads_filtered.fastq.gz -d flye_out/assembly.fasta \
    -o medaka/ -t 8 -m r1041_e82_400bps_hac_v4.2.0
```

**Pattern 3: Structural variant calling**
```bash
sniffles --input aligned.bam --vcf svs.vcf \
    --reference hg38.fa --threads 8 --minsupport 5
```

**Pattern 4: Methylation detection**
```bash
# Basecall with methylation model
dorado basecaller hac,5mCG_5hmCG pod5/ > calls_mod.bam

# Extract CpG methylation
modkit pileup calls_mod.bam methylation.bed \
    --ref hg38.fa --cpg --threads 8
```

**Pattern 5: Isoform analysis (R)**
```r
library(bambu)
se <- bambu(reads='aligned.bam',
            annotations=gencode_gtf,
            genome=hg38_fa)
writeBambuOutput(se, path='bambu_output/')
# serowRanges contains isoform coordinates
# assay(se, 'CPM') contains isoform-level expression
```

## Technology Comparison

| Property | ONT R9/R10 | PacBio HiFi |
|----------|-----------|-------------|
| Read length | Typically 5–50 kb | Typically 15–25 kb |
| Raw accuracy | 97–99% (R10) | >99.9% |
| Throughput | High (P2 Solo: 80 Gb) | Moderate (Sequel II: 160 Gb) |
| Methylation | Direct (native DNA) | 5mC with Kinetics |
| Cost per Gb | Low | Higher |

## Pitfalls

- **Basecall model selection** — match model to flow cell (R9.4 vs R10.4) and kit chemistry
- **Assembly genome size** — always provide `--genome-size` to Flye for ploidy-aware assembly
- **Coverage for assembly** — aim for ≥50× for Flye; ≥30× for Hifiasm HiFi
- **SV minimum support** — default `--minsupport 5` for Sniffles2; lower for low-coverage data
- **Medaka model** — use the correct model matching your basecaller version and flow cell

## Code Templates

### NanoStat QC Summary Parser
```python
import subprocess
import re

def parse_nanostat(fastq_path):
    """Run NanoStat and return key metrics as dict."""
    out = subprocess.check_output(
        ['NanoStat', '--fastq', fastq_path, '-t', '4'],
        stderr=subprocess.DEVNULL).decode()
    metrics = {}
    for line in out.splitlines():
        m = re.match(r'(.+?):\s+([\d.,]+)', line.strip())
        if m:
            key = m.group(1).strip().lower().replace(' ', '_')
            metrics[key] = float(m.group(2).replace(',', ''))
    return metrics

stats = parse_nanostat('reads.fastq.gz')
print(f"N50: {stats.get('read_length_n50', 'N/A')} bp")
print(f"Mean Q: {stats.get('mean_read_quality', 'N/A')}")
```

### Parse Sniffles2 VCF for SVs
```python
import pandas as pd

def parse_sv_vcf(vcf_path):
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, sv_id, ref, alt = fields[:5]
            info = dict(
                kv.split('=', 1) if '=' in kv else (kv, True)
                for kv in fields[7].split(';')
            )
            records.append({
                'chrom': chrom, 'pos': int(pos),
                'svtype': info.get('SVTYPE', ''),
                'svlen': abs(int(info.get('SVLEN', 0))),
                'support': int(info.get('SUPPORT', 0)),
                'af': float(info.get('AF', 0)),
            })
    return pd.DataFrame(records)

svs = parse_sv_vcf('sniffles_svs.vcf')
deletions = svs[svs['svtype'] == 'DEL']
large_dels = deletions[deletions['svlen'] >= 1000]
print(f"Large deletions (≥1 kb): {len(large_dels)}")
```

### Assembly N50 Calculator
```python
from Bio import SeqIO

def assembly_stats(fasta_path):
    lengths = sorted([len(r.seq) for r in SeqIO.parse(fasta_path, 'fasta')],
                     reverse=True)
    total = sum(lengths)
    cumsum = 0
    n50 = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total * 0.5:
            n50 = l
            break
    return {
        'num_contigs': len(lengths),
        'total_length': total,
        'largest': lengths[0],
        'n50': n50,
    }

stats = assembly_stats('assembly.fasta')
print(f"N50 = {stats['n50']:,} bp | Total = {stats['total_length']:,} bp")
```

## Related Skills
- `ngs-variant-calling` — short-read alignment, GATK variant calling, VCF format
- `genome-assembly` / `proteomics` — assembly algorithms, quality metrics, annotation
- `sequence-alignment` — pairwise alignment, BLAST, Minimap2 PAF format
- `structural-bioinformatics` — interpreting structural variants in protein context
