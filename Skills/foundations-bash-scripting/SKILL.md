---
name: foundations-bash-scripting
description: "Bash scripting essentials for bioinformatics: variables, conditionals, loops, and pipeline patterns."
tool_type: python
primary_tool: Python
---

# Bash Scripting for Bioinformatics

## Script Template

```bash
#!/bin/bash
set -euo pipefail

# Configuration
THREADS=8
INPUT_DIR="${1:?Usage: $0 <input_dir> <output_dir>}"
OUTPUT_DIR="${2:?Usage: $0 <input_dir> <output_dir>}"

log() { echo "[$(date '+%H:%M:%S')] $1" >&2; }

mkdir -p "$OUTPUT_DIR"
log "Starting..."
```

## Variables

```bash
sample="BRCA_001"
threads=8
out="${results_dir}/${sample}.bam"   # ${} safer when appending

current_date=$(date +%Y-%m-%d)       # command substitution
num_cpus=$(nproc)

# Special variables
$0   # script name
$1   # first argument
$@   # all arguments (array)
$#   # argument count
$?   # exit code of last command
$$   # current PID
```

## Conditionals

```bash
# File tests (most common in bioinformatics)
[[ -f "$file" ]]   # is a regular file
[[ -d "$dir" ]]    # is a directory
[[ -s "$file" ]]   # file exists and is non-empty
[[ -z "$var" ]]    # string is empty
[[ -n "$var" ]]    # string is non-empty

# Check tool availability
command -v samtools &>/dev/null || { echo "samtools not found"; exit 1; }

# Input validation pattern
if [[ -z "$input_file" ]]; then
    echo "ERROR: No input file" >&2; exit 1
fi
if [[ ! -f "$input_file" ]]; then
    echo "ERROR: Not found: $input_file" >&2; exit 1
fi
```

Numeric: `-eq -ne -lt -le -gt -ge`. String: `== != -z -n`.
Always use `[[ ]]` (double brackets), not `[ ]`.

## Loops

```bash
# Process all FASTQ files
for fastq in "${INPUT_DIR}"/*.fastq.gz; do
    [[ -f "$fastq" ]] || { echo "No .fastq.gz files found"; exit 1; }
    sample=$(basename "$fastq" .fastq.gz)
    fastqc -t "$THREADS" -o "$OUTPUT_DIR" "$fastq"
done

# Read sample sheet (tab-separated)
while IFS=$'\t' read -r sample_id condition; do
    echo "Processing $sample_id ($condition)"
done < sample_sheet.tsv

# C-style with index
for ((i=1; i<=22; i++)); do
    echo "chr${i}"
done
```

## Case Statement

```bash
case "$file" in
    *.fastq.gz | *.fq.gz) fastqc "$file" ;;
    *.bam)                 samtools flagstat "$file" ;;
    *.vcf | *.vcf.gz)      bcftools stats "$file" ;;
    *.fasta | *.fa)        grep -c "^>" "$file" ;;
    *) echo "Unknown: $file"; exit 1 ;;
esac
```

## Comparison Table

| Goal | Command |
|------|---------|
| Default value | `${var:-default}` |
| Require arg | `${1:?Usage: ...}` |
| Strip extension | `$(basename "$f" .gz)` |
| Dir of file | `$(dirname "$f")` |
| Count lines | `wc -l < "$file"` |
| Redirect stderr | `cmd 2>/dev/null` |

## Pitfalls

- **`set -euo pipefail` is non-negotiable**: without it, a failed command silently continues; `-u` catches unset variables; `pipefail` catches failures inside pipes.
- **No spaces around `=`**: `var=value` is correct; `var = value` runs a command named `var`.
- **Always double-quote variables**: `"$var"` not `$var`; a filename with a space breaks unquoted expansions into multiple words.
- **`$()` not backticks**: backticks cannot be nested and are harder to read.
- **`local` in functions**: undeclared variables leak into global scope and may clobber outer variables with the same name.
- **Glob matching empty directories**: `for f in dir/*.gz` — when no files match, the loop body runs once with the literal glob string; check `[[ -f "$f" ]]` at the top of the loop.
