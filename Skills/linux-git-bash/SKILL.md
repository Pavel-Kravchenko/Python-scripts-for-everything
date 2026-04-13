---
name: linux-git-bash
description: Linux shell commands, git workflows, bash scripting, and file encodings for bioinformatics data processing
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Linux, Git & Bash for Bioinformatics

## When to Use
- Processing FASTA/FASTQ/BAM/VCF files from the command line (counting, filtering, format conversion)
- Setting up reproducible analysis projects with git version control and proper `.gitignore`
- Writing batch pipeline scripts that loop over samples with robust error handling

## Quick Reference

### Navigation & File Operations
| Command | Purpose |
|---------|---------|
| `ls -lah` | Long listing with human-readable sizes |
| `mkdir -p project/{data/{raw,processed},results,scripts,logs}` | Nested dirs via brace expansion |
| `cp -r src/ dst/` | Recursive copy |
| `rm -r dir/` | Recursive delete (no trash bin) |
| `chmod +x script.sh` | Make executable |
| `find . -name "*.fastq.gz" -size +1G` | Find large FASTQ files |
| `find . -mtime -7 -name "*.bam"` | BAM files modified in last 7 days |

### Viewing Files
| Command | Purpose |
|---------|---------|
| `less file` | Interactive pager (q=quit, /=search) |
| `head -n 4 reads.fastq` | First read (4 lines per read) |
| `tail -f pipeline.log` | Follow log in real time |
| `wc -l file` | Line count |
| `zcat file.fastq.gz \| head -8` | View compressed FASTQ |

### Redirection & Pipes
| Operator | Meaning |
|----------|---------|
| `>` | Redirect stdout (overwrite) |
| `>>` | Redirect stdout (append) |
| `2>` | Redirect stderr |
| `2>&1` | Merge stderr into stdout |
| `\|` | Pipe stdout to next command's stdin |

### Text Processing
| Command | Key flags / example |
|---------|---------------------|
| `grep -c "^>" seqs.fasta` | Count FASTA sequences |
| `grep -v "^#" variants.vcf` | Strip VCF headers |
| `grep -B 1 "GAATTC" seqs.fasta` | Match + preceding header |
| `cut -f1,3 genes.bed` | Tab-delimited columns 1 and 3 |
| `cut -d',' -f2 data.csv` | CSV column 2 |
| `sort -k2,2n file` | Numeric sort on column 2 |
| `sort -k1,1 -k2,2n` | Multi-key sort (chr then pos) |
| `uniq -c` | Count consecutive duplicates |
| `tr -d '\n'` | Delete newlines (join sequences) |
| `sed 's/chr//'` | Strip chr prefix in-place |
| `sed -i 's/old/new/g' file` | In-place global substitution |
| `sed '/^#/d'` | Delete comment lines |
| `awk -F'\t' '{print $1,$3-$2}' bed` | Feature lengths |
| `awk 'NR%4==2'` | FASTQ sequence lines only |

### Compression
| Command | Purpose |
|---------|---------|
| `gzip -k file.fastq` | Compress, keep original |
| `gunzip file.fastq.gz` | Decompress |
| `zcat file.fastq.gz \| wc -l` | Count lines without decompressing |
| `tar -czvf archive.tar.gz dir/` | Create compressed archive |
| `tar -xzvf archive.tar.gz` | Extract |

### samtools (BAM files — never use cat/grep on BAM)
```bash
samtools view aligned.bam | head -5          # View as SAM
samtools view -H aligned.bam                 # Header only
samtools view -c -F 4 aligned.bam            # Count aligned reads
samtools index aligned.bam                   # Required before random access
samtools view aligned.bam chr17:7571720-7590868  # Region extract (TP53)
samtools flagstat aligned.bam                # Alignment statistics
samtools sort -o sorted.bam unsorted.bam     # Sort by coordinate
```

### Git Commands
| Command | Purpose |
|---------|---------|
| `git init` | New repo |
| `git clone URL` | Clone remote |
| `git status` | Working tree status |
| `git add file` | Stage specific file |
| `git commit -m "msg"` | Commit staged changes |
| `git log --oneline --graph` | Visual history |
| `git diff --staged` | Staged vs last commit |
| `git diff HEAD~1` | Changes since previous commit |
| `git show HEAD` | Latest commit + diff |
| `git restore file` | Discard working-dir changes |
| `git restore --staged file` | Unstage |
| `git reset --soft HEAD~1` | Undo commit, keep changes staged |
| `git revert <hash>` | Safe undo (new commit) |
| `git stash` / `git stash pop` | Shelve uncommitted changes |
| `git tag -a v1.0 -m "msg"` | Annotated tag |
| `git push --tags` | Push tags to remote |
| `git log -S "alpha"` | Find commits that changed a string |

### Branch Workflow
```bash
git checkout -b feature-batch-correction   # Create + switch
git merge feature-batch-correction         # Merge into current
git branch -d feature-batch-correction     # Delete merged branch
git push -u origin feature-name            # Push new branch to remote
```

---

## Key Patterns

**1. Count bio-file records**
```bash
grep -c "^>" proteins.fasta                         # FASTA sequences
zcat sample.fastq.gz | wc -l | awk '{print $1/4}'  # FASTQ reads
grep -v "^#" variants.vcf | wc -l                   # VCF variants
```

**2. VCF chromosome distribution**
```bash
grep -v "^#" variants.vcf | cut -f1 | sort | uniq -c | sort -rn
```

**3. FASTQ → FASTA conversion**
```bash
sed -n '1~4s/^@/>/p;2~4p' reads.fastq > reads.fasta
```

**4. Average read length from FASTQ**
```bash
awk 'NR%4==2 {sum+=length($0); count++} END {print sum/count}' reads.fastq
```

**5. Extract gene names from GTF**
```bash
awk -F'\t' '$3=="gene"' gencode.gtf \
  | grep -o 'gene_name "[^"]*"' \
  | sed 's/gene_name "//;s/"//' | sort -u
```

**6. GC content of a FASTA genome**
```bash
grep -v "^>" genome.fa | tr -d '\n' \
  | awk '{gc=gsub(/[GC]/,"",$0); at=gsub(/[AT]/,"",$0); print gc/(gc+at)*100"%"}'
```

**7. BED feature length with awk**
```bash
awk -F'\t' '{print $0 "\t" $3-$2}' regions.bed
awk -F'\t' '{len=$3-$2; sum+=len; n++} END {print sum/n}' regions.bed  # average
```

**8. Filter high-quality VCF variants**
```bash
grep -v "^#" sample.vcf | awk -F'\t' '$6 > 30' | sort -k6,6rn
```

**9. Parallel FastQC with xargs**
```bash
find data/ -name "*.fastq.gz" | xargs -P 4 -I {} fastqc {} -o results/qc/
```

**10. Paired-end R1 → R2 derivation**
```bash
r2="${r1/_R1/_R2}"
sample=$(basename "$r1" _R1.fastq.gz)
```

---

## Code Templates

### Bioinformatics .gitignore
```gitignore
# Large sequencing data
*.fastq *.fastq.gz *.fq.gz
*.bam *.bam.bai *.sam *.cram
*.vcf *.vcf.gz *.bcf
*.sra
data/raw/

# Reference genomes
*.fa *.fasta *.fa.fai *.dict

# Outputs (regenerated by scripts)
results/
*.log *.tmp *.out

# Python
__pycache__/ *.pyc
.ipynb_checkpoints/

# R
.Rhistory .RData

# OS / IDE
.DS_Store Thumbs.db
.vscode/ .idea/
```

### Robust bash script skeleton
```bash
#!/bin/bash
# ============================================================
# Script: pipeline.sh
# Usage:  ./pipeline.sh <input_dir> [output_dir]
# ============================================================
set -euo pipefail

INPUT_DIR="${1:-}"
OUTPUT_DIR="${2:-results}"
THREADS=8
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

### Batch FASTQ processing loop
```bash
#!/bin/bash
set -euo pipefail
INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-qc_results}"
mkdir -p "$OUTPUT_DIR"

count=0
for fastq in "${INPUT_DIR}"/*.fastq.gz; do
    [[ -f "$fastq" ]] || { echo "No .fastq.gz files in $INPUT_DIR"; exit 1; }
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
input_dir="${1:-.}"
output_file="${2:-sample_sheet.tsv}"
echo -e "sample_id\tR1_path\tR2_path" > "$output_file"

for r1 in "${input_dir}"/*_R1.fastq.gz; do
    [[ -f "$r1" ]] || { echo "No *_R1.fastq.gz in $input_dir"; exit 1; }
    r2="${r1/_R1/_R2}"
    sample=$(basename "$r1" _R1.fastq.gz)
    [[ -f "$r2" ]] || { echo "WARNING: Missing R2 for $sample"; continue; }
    echo -e "${sample}\t${r1}\t${r2}" >> "$output_file"
done
```

### FASTA report
```bash
#!/bin/bash
set -euo pipefail
file="${1:?Usage: $0 <fasta_file>}"
[[ -f "$file" ]]  || { echo "ERROR: Not found: $file"; exit 1; }
[[ -s "$file" ]]  || { echo "ERROR: Empty file: $file"; exit 1; }

num_seqs=$(grep -c "^>" "$file")
total_nt=$(grep -v "^>" "$file" | tr -d '\n' | wc -c | tr -d ' ')
echo "Sequences:  $num_seqs"
echo "Total nt:   $total_nt"
echo "Avg length: $(( total_nt / num_seqs ))"
```

### Case statement: route by file type
```bash
case "$file" in
    *.fastq.gz | *.fq.gz) fastqc "$file" ;;
    *.bam)                 samtools index "$file"; samtools flagstat "$file" ;;
    *.vcf | *.vcf.gz)      bcftools stats "$file" ;;
    *.fasta | *.fa)        grep -c "^>" "$file" ;;
    *)                     echo "Unknown format: $file"; exit 1 ;;
esac
```

### Git workflow for paper submission
```bash
# Daily work cycle
git pull
# ... edit scripts ...
git status
git add scripts/deseq2_analysis.R
git commit -m "Tighten significance threshold from 0.05 to 0.01"
git push

# Tag paper submission
git tag -a v1.0 -m "Submitted to Nature Genetics 2024-01-15"
git push --tags

# Feature branch for experimental analysis
git checkout -b experiment-limma-vs-deseq2
# ... work ...
git checkout main
git merge experiment-limma-vs-deseq2
git branch -d experiment-limma-vs-deseq2
```

---

## File Encoding Handling

### Core rules
| Scenario | Solution |
|----------|---------|
| Read any text file | `open(f, encoding='utf-8')` |
| Windows file with BOM | `open(f, encoding='utf-8-sig')` |
| Windows line endings | `text.replace('\r\n', '\n')` or `open(f, newline='')` |
| Garbled characters (mojibake) | Re-decode: `raw.decode('latin-1')` |
| Unknown encoding | `chardet.detect(raw_bytes)` then try UTF-8 → Latin-1 |
| Inspect raw bytes | `open(f, 'rb')` |
| Binary formats (BAM, gzip) | Always `'rb'` mode |

### FASTQ Phred+33 quality encoding
- Quality char → Phred score: `phred = ord(char) - 33`
- Phred score → error probability: `P = 10 ** (-phred / 10)`
- Valid range: ASCII 33 (`!`) to 126 (`~`) → Phred 0–93

### Invisible character pitfalls (copy-paste from PDF/Word)
- Non-breaking space U+00A0 — looks like space, breaks sequence parsers
- Zero-width space U+200B — completely invisible
- Soft hyphen U+00AD — may silently corrupt sequence data

```python
import unicodedata

def sanitize_sequence(seq, valid_chars='ATGCNatgcn'):
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

### Fix FASTA bytes (BOM + line endings + non-ASCII sequences)
```python
def fix_fasta(raw_bytes):
    text = raw_bytes.decode('utf-8-sig')                   # strips BOM
    text = text.replace('\r\n', '\n').replace('\r', '\n')  # normalize endings
    lines = []
    for line in text.split('\n'):
        if line.startswith('>'):
            lines.append(line)                             # keep header intact
        else:
            lines.append(''.join(c for c in line if ord(c) < 128))
    return '\n'.join(lines)
```

### Shell: detect and fix line endings
```bash
file sample.txt           # reports "CRLF line terminators" for Windows files
dos2unix sample.txt       # convert \r\n -> \n in-place
sed -i 's/\r//' file.txt  # alternative without dos2unix
```

---

## Common Pitfalls

- **`set -euo pipefail` omitted** — silent failures cascade; pipelines produce garbage without error
- **Unquoted variables** — `ls $file` breaks on spaces; always use `"$file"`
- **`git add .` in large repos** — accidentally stages `.bam`/`.fastq.gz`; use `git add <specific files>`
- **Spaces around `=` in bash** — `var = "value"` is a syntax error; use `var="value"`
- **`cat large.bam`** — BAM is binary; use `samtools view` instead
- **`grep pattern file.bam`** — same; grep cannot parse binary BAM
- **Forgetting `2>/dev/null` or `|| true` in cleanup traps** — cleanup itself can cause non-zero exit inside `set -e`
- **Committing large data files** — GitHub rejects files >100 MB; use `.gitignore` before first commit
- **`git reset --hard`** — permanently destroys uncommitted work; prefer `git restore` or `git reset --soft`
- **Windows `\r\n` line endings in FASTA** — `\r` appears at end of sequence, corrupts parsers; run `dos2unix`
- **Latin-1 fallback silently succeeds** — always logs a WARNING; UTF-8 is the correct target
- **`for f in *.fastq.gz`** — if no files match, `$f` is the literal string `*.fastq.gz`; guard with `[[ -f "$f" ]]`

---

## Related Skills
- `python-core-bio` — File I/O patterns complement shell processing; Python for complex parsing
- `ngs-variant-calling` — Shell pipelines for NGS data: BWA → samtools → GATK
