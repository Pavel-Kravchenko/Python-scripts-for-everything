---
name: linux-git-bash
description: Linux shell commands, git workflows, bash scripting, and file encoding handling for bioinformatics data processing.
tool_type: python
primary_tool: Python
---

# Linux, Git & Bash for Bioinformatics

## When to Use
- Processing FASTA/FASTQ/BAM/VCF files from the command line
- Setting up reproducible analysis projects with git
- Writing batch pipeline scripts that loop over samples

## Shell Quick Reference

### Navigation & File Operations
| Command | Purpose |
|---------|---------|
| `ls -lah` | Long listing with human-readable sizes |
| `mkdir -p project/{data/{raw,processed},results,scripts}` | Nested dirs via brace expansion |
| `find . -name "*.fastq.gz" -size +1G` | Find large FASTQ files |
| `find . -mtime -7 -name "*.bam"` | BAM files modified in last 7 days |
| `chmod +x script.sh` | Make executable |

### Text Processing
| Command | Example / key flags |
|---------|---------------------|
| `grep -c "^>" seqs.fasta` | Count FASTA sequences |
| `grep -v "^#" variants.vcf` | Strip VCF headers |
| `grep -B 1 "GAATTC" seqs.fasta` | Match + preceding header |
| `cut -f1,3 genes.bed` | Tab-delimited columns 1 and 3 |
| `sort -k1,1 -k2,2n` | Multi-key sort (chr then pos) |
| `awk 'NR%4==2'` | FASTQ sequence lines only |
| `sed 's/chr//'` | Strip chr prefix |
| `sed -i 's/old/new/g' file` | In-place global substitution |

### samtools (never use cat/grep on BAM)
```bash
samtools view aligned.bam | head -5           # View as SAM
samtools view -c -F 4 aligned.bam             # Count aligned reads
samtools index aligned.bam                    # Required before random access
samtools view aligned.bam chr17:7571720-7590868  # Region extract
samtools flagstat aligned.bam                 # Alignment statistics
samtools sort -o sorted.bam unsorted.bam      # Sort by coordinate
```

### Git Commands
| Command | Purpose |
|---------|---------|
| `git log --oneline --graph` | Visual history |
| `git diff --staged` | Staged vs last commit |
| `git restore file` | Discard working-dir changes |
| `git restore --staged file` | Unstage |
| `git reset --soft HEAD~1` | Undo commit, keep staged |
| `git revert <hash>` | Safe undo (new commit) |
| `git stash` / `git stash pop` | Shelve uncommitted changes |
| `git log -S "alpha"` | Find commits that changed a string |
| `git tag -a v1.0 -m "msg"` + `git push --tags` | Annotated release tag |

## Key Shell Patterns

```bash
# Count bio-file records
grep -c "^>" proteins.fasta
zcat sample.fastq.gz | wc -l | awk '{print $1/4}'
grep -v "^#" variants.vcf | wc -l

# VCF chromosome distribution
grep -v "^#" variants.vcf | cut -f1 | sort | uniq -c | sort -rn

# FASTQ → FASTA
sed -n '1~4s/^@/>/p;2~4p' reads.fastq > reads.fasta

# Average read length
awk 'NR%4==2 {sum+=length($0); count++} END {print sum/count}' reads.fastq

# Extract gene names from GTF
awk -F'\t' '$3=="gene"' gencode.gtf \
  | grep -o 'gene_name "[^"]*"' \
  | sed 's/gene_name "//;s/"//' | sort -u

# BED feature lengths
awk -F'\t' '{print $0 "\t" $3-$2}' regions.bed

# Parallel FastQC
find data/ -name "*.fastq.gz" | xargs -P 4 -I {} fastqc {} -o results/qc/

# Paired-end R2 from R1 path
r2="${r1/_R1/_R2}"; sample=$(basename "$r1" _R1.fastq.gz)
```

## Code Templates

### Robust bash script skeleton
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="${1:-}"
OUTPUT_DIR="${2:-results}"
LOGFILE="${OUTPUT_DIR}/pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"; }
cleanup() { rm -f /tmp/pipeline_*.tmp 2>/dev/null || true; }
trap cleanup EXIT

[[ -z "$INPUT_DIR" ]] && { echo "Usage: $0 <input_dir> [output_dir]"; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: Not a directory: $INPUT_DIR"; exit 1; }
command -v samtools &>/dev/null || { log "ERROR: samtools not installed"; exit 1; }

mkdir -p "$OUTPUT_DIR"
log "Starting pipeline. Input: $INPUT_DIR"
```

### Batch FASTQ processing
```bash
#!/bin/bash
set -euo pipefail
INPUT_DIR="${1:-.}"; OUTPUT_DIR="${2:-qc_results}"
mkdir -p "$OUTPUT_DIR"
count=0
for fastq in "${INPUT_DIR}"/*.fastq.gz; do
    [[ -f "$fastq" ]] || { echo "No .fastq.gz files found"; exit 1; }
    sample=$(basename "$fastq" .fastq.gz)
    echo "[$(( ++count ))] Processing: $sample"
    fastqc "$fastq" -o "$OUTPUT_DIR" -t 4
done
echo "Done. Processed $count files."
```

### Sample sheet generator (paired-end)
```bash
#!/bin/bash
set -euo pipefail
input_dir="${1:-.}"; output_file="${2:-sample_sheet.tsv}"
echo -e "sample_id\tR1_path\tR2_path" > "$output_file"
for r1 in "${input_dir}"/*_R1.fastq.gz; do
    [[ -f "$r1" ]] || { echo "No *_R1.fastq.gz files"; exit 1; }
    r2="${r1/_R1/_R2}"; sample=$(basename "$r1" _R1.fastq.gz)
    [[ -f "$r2" ]] || { echo "WARNING: Missing R2 for $sample"; continue; }
    echo -e "${sample}\t${r1}\t${r2}" >> "$output_file"
done
```

### Bioinformatics .gitignore
```gitignore
*.fastq *.fastq.gz *.fq.gz
*.bam *.bam.bai *.sam *.cram
*.vcf *.vcf.gz *.bcf *.sra
*.fa *.fasta *.fa.fai *.dict
data/raw/ results/ *.log *.tmp
__pycache__/ *.pyc .ipynb_checkpoints/
.Rhistory .RData .DS_Store .vscode/ .idea/
```

## File Encoding Handling

| Scenario | Solution |
|----------|---------|
| Read any text file | `open(f, encoding='utf-8')` |
| Windows file with BOM | `open(f, encoding='utf-8-sig')` |
| Windows line endings | `text.replace('\r\n', '\n')` |
| Garbled characters | Re-decode: `raw.decode('latin-1')` |
| Unknown encoding | `chardet.detect(raw_bytes)` then try UTF-8 → Latin-1 |
| Binary formats (BAM, gzip) | Always `'rb'` mode |

FASTQ Phred+33: `phred = ord(char) - 33`, `P_error = 10 ** (-phred / 10)`. Valid range: ASCII 33 (`!`) to 126 (`~`).

```python
import unicodedata

def sanitize_sequence(seq: str, valid_chars: str = 'ATGCNatgcn') -> str:
    cleaned, removed = [], []
    for char in seq:
        if char in valid_chars:
            cleaned.append(char)
        elif char not in ('\n', '\r', ' ', '\t'):
            removed.append(f"'{char}' ({unicodedata.name(char, f'U+{ord(char):04X}')})")
    if removed:
        print(f"WARNING: Removed {removed}")
    return ''.join(cleaned)
```

## Pitfalls

- **`set -euo pipefail` omitted**: silent failures cascade — pipelines produce garbage without error messages.
- **Unquoted variables**: `ls $file` breaks on spaces; always use `"$file"`.
- **`git add .` in large repos**: accidentally stages `.bam`/`.fastq.gz`; use `git add <specific files>`.
- **Spaces around `=` in bash**: `var = "value"` is a syntax error; `var="value"` is correct.
- **`cat large.bam` or `grep pattern file.bam`**: BAM is binary — use `samtools view` instead.
- **`for f in *.fastq.gz` with no matches**: `$f` becomes the literal string `*.fastq.gz`; guard with `[[ -f "$f" ]]`.
- **`cleanup` trap failing**: use `|| true` so cleanup errors don't trigger `set -e` exit inside the trap.
- **Committing large data files**: GitHub rejects files >100 MB; configure `.gitignore` before first commit.
- **`git reset --hard`**: permanently destroys uncommitted work; prefer `git restore` or `git reset --soft`.
- **Windows `\r\n` line endings in FASTA**: `\r` at end of sequence corrupts parsers; run `dos2unix`.
- **Non-breaking space U+00A0 in sequences**: looks like space, breaks parsers when copy-pasted from PDF/Word.
