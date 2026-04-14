---
name: foundations-linux-fundamentals
description: "Linux command-line essentials for bioinformatics: navigation, file ops, pipes, grep, file format inspection, and vim survival guide."
tool_type: python
primary_tool: Python
---

# Linux Fundamentals for Bioinformatics

## Pitfalls

- **Spaces around `=` in Bash:** `var = value` is wrong; `var=value` is right. The space makes Bash interpret `var` as a command.
- **`rm` is permanent:** No trash bin. Always double-check before `rm -rf`. Never run `rm -rf /` or `rm -rf ~`.
- **Pipes discard stderr:** `cmd1 | cmd2` only passes stdout. Add `2>&1` to also pipe error messages.
- **Wildcards expand before the command runs:** `rm *.fastq` — the shell does the glob expansion, not `rm`. If no files match, you get an error (not a no-op).
- **Relative vs absolute paths:** Scripts run from cron or other directories need absolute paths.
- **Compressed files:** Use `zcat`/`zgrep`/`zless` for `.gz` files; never decompress just to inspect.

## Navigation

| Symbol | Meaning |
|---|---|
| `/` | Root of filesystem |
| `~` | Home directory |
| `.` | Current directory |
| `..` | Parent directory |
| `-` | Previous directory (`cd -`) |

```bash
pwd           # where am I
ls -lah       # long + hidden + human sizes — daily driver
cd -          # back to previous directory
```

## File Operations

```bash
# Create project structure
mkdir -p RNA_Seq/{00_raw,01_qc,02_trimmed,03_aligned,04_counts,scripts,logs}

# Copy, move, delete
cp -r src_dir/ dest_dir/        # -r required for directories
mv old.txt new.txt              # rename
rm -rf directory/               # permanent — no undo

# Inspect without decompressing
zcat sample.fastq.gz | head -8  # first 2 reads (4 lines each)
zgrep "PASS" variants.vcf.gz    # grep compressed file
```

## Viewing Files

```bash
head -n 4 reads.fastq    # first FASTQ read (4 lines per read)
tail -f pipeline.log     # follow log in real time
wc -l file.txt           # count lines
wc -c file.txt           # count bytes
less large.bam           # interactive scroll (q=quit, /=search)
```

## Bioinformatics File Formats

| Format | Extension | Content |
|---|---|---|
| FASTA | `.fa`, `.fasta` | Sequences (genome, protein) |
| FASTQ | `.fq`, `.fastq.gz` | Reads + quality scores (4 lines/read) |
| SAM/BAM | `.sam`, `.bam` | Alignments (BAM = binary SAM) |
| BED | `.bed` | Genomic intervals (0-based half-open) |
| VCF | `.vcf` | Variant calls (1-based) |
| GFF/GTF | `.gff`, `.gtf` | Gene annotations (1-based) |

## Pipes and Redirection

```bash
# Operators
>    # redirect stdout (overwrite)
>>   # redirect stdout (append)
2>   # redirect stderr
2>&1 # merge stderr into stdout
|    # pipe stdout to next command's stdin

# Common bioinformatics pipelines
grep -c "^>" proteins.fasta            # count FASTA sequences
zcat sample.fastq.gz | wc -l | awk '{print $1/4, "reads"}'  # count reads
grep -v "^#" variants.vcf | cut -f1 | sort | uniq -c        # variants per chr
cut -f1 regions.bed | sort | uniq -c | sort -rn              # regions per chr
```

## grep Quick Reference

```bash
grep -c "^>" genome.fa          # count sequences in FASTA
grep -v "^#" variants.vcf       # strip VCF header
grep "BRCA1" gencode.gtf        # find gene in annotation
grep -B 1 "GAATTC" seqs.fasta   # show header above EcoRI site
zgrep "PASS" variants.vcf.gz    # grep inside .gz
grep -r "error" logs/           # recursive search
```

## Vim Survival Guide

```
vim filename    open file
i               insert mode (type text)
Esc             back to normal mode
:w              save
:q              quit
:wq             save and quit
:q!             quit without saving
/pattern        search forward
n / N           next / previous match
dd              delete line
u               undo
:%s/old/new/g   replace all occurrences
```
