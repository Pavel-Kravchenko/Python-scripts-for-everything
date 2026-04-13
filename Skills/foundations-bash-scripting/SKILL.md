---
name: foundations-bash-scripting
description: "By the end of this module, you will be able to: - Write and execute Bash scripts to automate repetitive tasks - Use variables, conditionals, loops, and functions - Process multiple files in batch (FAS"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/03_Bash_Scripting/01_bash_scripting.ipynb"
---

# Module 0.3: Bash Scripting for Bioinformatics

*Source: Course notebook `Tier_0_Computational_Foundations/03_Bash_Scripting/01_bash_scripting.ipynb`*

# Module 0.3: Bash Scripting for Bioinformatics

---

### Learning Objectives

By the end of this module, you will be able to:
- Write and execute Bash scripts to automate repetitive tasks
- Use variables, conditionals, loops, and functions
- Process multiple files in batch (FASTQ, BAM, VCF)
- Build robust scripts with error handling and input validation
- Apply text processing tools (awk, sed, cut) inside scripts
- Create reusable analysis pipeline scripts

**Prerequisites:** Linux command line (Module 0.1)

**Estimated time:** 3-4 hours

---

## How to use this notebook
1. Read each script before running it — the comments explain what each section does and why.
2. Scripts are written to `/tmp/` so they do not pollute your project. After running a cell, try editing the script and re-running.
3. Pay particular attention to the error-handling sections — silent failures in pipelines are the most dangerous bugs.
4. The final exercise builds a complete pipeline script. Work through it even if it looks intimidating.

## Common stumbling points

- **`set -euo pipefail` is non-negotiable**: Without it, a script will blissfully continue after a command fails, producing garbage results with no warning.
- **No spaces around `=`**: `var=value` is correct; `var = value` makes Bash run a command named `var`.
- **Always double-quote variables**: `"$var"` not `$var`. If the value contains a space (a common filename problem), unquoted variables split into multiple words, breaking the command.
- **`$()` vs backticks**: Use `$(command)` for command substitution, not `` `command` ``. Backticks are harder to read and cannot be nested.
- **`local` in functions**: Always declare function variables with `local`. Without it, they leak into global scope and can overwrite variables with the same name elsewhere.

---

## 1. Your First Bash Script

A Bash script is just a text file containing commands that the shell executes in order.

```python
%%bash
# Create a simple script
cat > /tmp/hello_bioinfo.sh << 'EOF'
#!/bin/bash
# My first bioinformatics script
#
# The first line (#!/bin/bash) is called the "shebang".
# It tells the system which interpreter to use.

echo "Welcome to Bioinformatics!"
echo "Today is $(date '+%Y-%m-%d %H:%M')"
echo "User: $USER"
echo "Working directory: $(pwd)"
echo "System: $(uname -s)"
EOF

# Make it executable
chmod +x /tmp/hello_bioinfo.sh

# Run it
/tmp/hello_bioinfo.sh
```

### Script structure conventions

```bash
#!/bin/bash
# ============================================================
# Script: pipeline_step1.sh
# Description: Run FastQC on all raw FASTQ files
# Author: Your Name
# Date: 2024-01-15
# Usage: ./pipeline_step1.sh <input_dir> <output_dir>
# ============================================================

set -euo pipefail    # Strict error handling (explained later)

# --- Configuration ---
THREADS=8
MIN_QUALITY=20

# --- Functions ---
log() { echo "[$(date '+%H:%M:%S')] $1"; }

# --- Main logic ---
log "Starting analysis..."
# ... your commands ...
log "Done!"
```

Two ways to run a script:
```bash
chmod +x script.sh    # Make executable (once)
./script.sh           # Run directly

bash script.sh        # Run with bash (does not need chmod)
```

---

## 2. Variables

Variables store values that you reference throughout your script. This avoids
hard-coding paths and parameters, making scripts reusable.

```python
%%bash
# === Variable assignment ===
# CRITICAL: NO SPACES around the = sign!
# RIGHT:  name="value"
# WRONG:  name = "value"   (bash interprets this as running a command called 'name')

sample_name="BRCA_patient_001"
num_threads=8
genome_file="/references/GRCh38.fa"
output_dir="results/aligned"

# === Using variables ===
echo "Processing: $sample_name"
echo "Threads: ${num_threads}"           # ${} is safer when appending text
echo "Output: ${output_dir}/${sample_name}.bam"

# === Command substitution -- capture command output in a variable ===
current_date=$(date +%Y-%m-%d)
num_cpus=$(nproc)
echo "Date: $current_date"
echo "Available CPUs: $num_cpus"
```

```python
%%bash
# === Special variables ===
# These are set automatically by Bash

cat > /tmp/show_args.sh << 'EOF'
#!/bin/bash
echo "Script name:         $0"
echo "First argument:      $1"
echo "Second argument:     $2"
echo "All arguments:       $@"
echo "Number of arguments: $#"
echo "Exit code of last command: $?"
echo "Process ID:          $$"
EOF
chmod +x /tmp/show_args.sh

/tmp/show_args.sh sample_001.fastq.gz /data/output
```

### Quoting rules (this trips up everyone)

```bash
name="World"

echo "Hello $name"     # Double quotes: variables ARE expanded  -> Hello World
echo 'Hello $name'     # Single quotes: variables NOT expanded  -> Hello $name
echo Hello $name       # No quotes: works, but breaks on spaces in values!
```

**Golden rule: Always double-quote your variables.** `"$variable"` not `$variable`.

Why? If `file="my data.bam"`, then `ls $file` becomes `ls my data.bam` (two arguments),
while `ls "$file"` correctly becomes `ls "my data.bam"` (one argument).

### Bioinformatics example: Configurable pipeline

Using variables to make a pipeline configurable is the first step toward
reproducibility. Anyone can see and change the parameters at the top of the script.

```python
%%bash
cat > /tmp/pipeline_config.sh << 'EOF'
#!/bin/bash
# === RNA-Seq Pipeline Configuration ===

# Directories
PROJECT_DIR="/home/user/rna_seq_project"
RAW_DATA="${PROJECT_DIR}/00_raw_data"
QC_DIR="${PROJECT_DIR}/01_qc"
TRIMMED="${PROJECT_DIR}/02_trimmed"
ALIGNED="${PROJECT_DIR}/03_aligned"
COUNTS="${PROJECT_DIR}/04_counts"

# Reference files
GENOME="/references/GRCh38.fa"
ANNOTATION="/references/gencode.v38.gtf"
STAR_INDEX="/references/STAR_index"

# Parameters
THREADS=16
MIN_READ_LENGTH=35
QUALITY_THRESHOLD=20
ADAPTER="AGATCGGAAGAGC"

echo "Pipeline configuration:"
echo "  Project:     ${PROJECT_DIR}"
echo "  Genome:      ${GENOME}"
echo "  Threads:     ${THREADS}"
echo "  Min length:  ${MIN_READ_LENGTH}"
echo "  Min quality: ${QUALITY_THRESHOLD}"
EOF
bash /tmp/pipeline_config.sh
```

---

## 3. Conditionals (if/elif/else)

Conditionals let your script make decisions: does this file exist? Did the last
command succeed? Are there enough arguments?

```python
%%bash
# Basic syntax
number=42

if [[ $number -gt 100 ]]; then
    echo "Greater than 100"
elif [[ $number -gt 10 ]]; then
    echo "Between 11 and 100"
else
    echo "10 or less"
fi
```

### Comparison operators

```
NUMERIC                       STRING
  -eq   equals                  ==    equals
  -ne   not equals              !=    not equals
  -lt   less than               -z    string is empty
  -le   less or equal           -n    string is not empty
  -gt   greater than
  -ge   greater or equal

FILE TESTS                    LOGICAL
  -e file    exists             &&    AND
  -f file    is regular file    ||    OR
  -d file    is directory       !     NOT
  -s file    size > 0
  -r file    is readable
  -w file    is writable
  -x file    is executable
```

**Important:** Use `[[ ]]` (double brackets) in Bash, not `[ ]`. Double brackets
are safer with spaces and special characters.

```python
%%bash
# File existence checks -- you will use these constantly

test_file="/tmp/hello_bioinfo.sh"

if [[ -f "$test_file" ]]; then
    echo "File exists: $test_file"
else
    echo "File NOT found: $test_file"
fi

# Check if a directory exists, create it if not
output_dir="/tmp/test_results"
if [[ ! -d "$output_dir" ]]; then
    echo "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
fi

# Check if a command/tool is available
if command -v python3 &> /dev/null; then
    echo "Python3 is installed: $(python3 --version)"
else
    echo "Python3 is NOT installed"
fi
```

### Bioinformatics example: Input validation script

```python
%%bash
cat > /tmp/validate_input.sh << 'SCRIPT'
#!/bin/bash
# Validate input files before running a pipeline

input_file=$1

# Check if argument was provided
if [[ -z "$input_file" ]]; then
    echo "ERROR: No input file specified"
    echo "Usage: $0 <fastq_file>"
    exit 1
fi

# Check if file exists
if [[ ! -f "$input_file" ]]; then
    echo "ERROR: File not found: $input_file"
    exit 1
fi

# Check if file is not empty
if [[ ! -s "$input_file" ]]; then
    echo "ERROR: File is empty: $input_file"
    exit 1
fi

# Check file extension
if [[ "$input_file" != *.fastq && "$input_file" != *.fastq.gz && "$input_file" != *.fq.gz ]]; then
    echo "WARNING: Non-standard extension for FASTQ file: $input_file"
fi

echo "PASS: Input validation succeeded for $input_file"
SCRIPT
chmod +x /tmp/validate_input.sh

# Test with a valid file
echo "@read1" > /tmp/test.fastq
/tmp/validate_input.sh /tmp/test.fastq

# Test with missing file
/tmp/validate_input.sh /tmp/nonexistent.fastq

# Test with no arguments
/tmp/validate_input.sh
```

---

## 4. Loops

Loops are the reason you write Bash scripts. They let you process hundreds of files
with the same commands.

### for loop -- iterate over a list

```python
%%bash
# Iterate over explicit list
echo "=== Explicit list ==="
for sample in sample_01 sample_02 sample_03; do
    echo "Processing: $sample"
done

# Iterate over files matching a pattern
echo ""
echo "=== Files in /tmp matching *.sh ==="
for script in /tmp/*.sh; do
    echo "Found script: $(basename "$script")"
done

# Iterate over a range of numbers
echo ""
echo "=== Number range ==="
for i in {1..5}; do
    echo "  Chromosome $i"
done

# C-style for loop (useful when you need the index)
echo ""
echo "=== C-style loop ==="
for ((i=0; i<3; i++)); do
    echo "  Index: $i"
done
```

### Bioinformatics example: Batch processing FASTQ files

```python
%%bash
cat > /tmp/batch_qc.sh << 'SCRIPT'
#!/bin/bash
set -euo pipefail

INPUT_DIR="${1:-.}"       # First argument, default to current dir
OUTPUT_DIR="${2:-qc_results}"
THREADS=4

mkdir -p "$OUTPUT_DIR"

echo "Running QC on all FASTQ files in: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "---"

count=0
for fastq in "${INPUT_DIR}"/*.fastq.gz; do
    # Check if the glob actually matched any files
    if [[ ! -f "$fastq" ]]; then
        echo "No .fastq.gz files found in $INPUT_DIR"
        exit 1
    fi

    sample=$(basename "$fastq" .fastq.gz)
    echo "[$((++count))] Processing: $sample"

    # In a real pipeline, you would run:
    # fastqc -t $THREADS -o "$OUTPUT_DIR" "$fastq"
    echo "    fastqc -t $THREADS -o $OUTPUT_DIR $fastq"
done

echo "---"
echo "Processed $count files"
SCRIPT

# Create test files
mkdir -p /tmp/test_fastq
for i in 1 2 3 4 5; do
    echo "@read" | gzip > /tmp/test_fastq/sample_${i}.fastq.gz
done

bash /tmp/batch_qc.sh /tmp/test_fastq /tmp/qc_output
```

### while loop -- condition-based iteration

```python
%%bash
# Basic while loop
count=1
while [[ $count -le 5 ]]; do
    echo "Iteration $count"
    ((count++))
done

echo ""

# Read a file line by line -- extremely common pattern
cat > /tmp/sample_list.txt << 'EOF'
sample_001	control
sample_002	control
sample_003	treatment
sample_004	treatment
EOF

echo "=== Reading sample list ==="
while IFS=$'\t' read -r sample_id condition; do
    echo "Sample: $sample_id  |  Condition: $condition"
done < /tmp/sample_list.txt
```

### Loop control: break and continue

```python
%%bash
# continue -- skip the rest of the current iteration
# break    -- exit the loop entirely

echo "=== Processing files (skip empty, stop if error) ==="

# Create some test files
echo "data" > /tmp/good_file.txt
touch /tmp/empty_file.txt              # Empty file
echo "more data" > /tmp/another_good.txt

for f in /tmp/good_file.txt /tmp/empty_file.txt /tmp/another_good.txt; do
    # Skip empty files
    if [[ ! -s "$f" ]]; then
        echo "  SKIP (empty): $(basename $f)"
        continue
    fi

    echo "  Processing: $(basename $f)"
done
```

---

## 5. Case Statements

When you need to handle multiple conditions based on a single value, `case` is
cleaner than a long chain of if/elif.

```python
%%bash
cat > /tmp/process_file.sh << 'SCRIPT'
#!/bin/bash
# Route bioinformatics files to the appropriate processing tool

file=$1

if [[ -z "$file" ]]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

case "$file" in
    *.fastq.gz | *.fq.gz)
        echo "FASTQ file -> Running FastQC"
        # fastqc "$file"
        ;;
    *.bam)
        echo "BAM file -> Generating index and stats"
        # samtools index "$file"
        # samtools flagstat "$file"
        ;;
    *.vcf | *.vcf.gz)
        echo "VCF file -> Running bcftools stats"
        # bcftools stats "$file"
        ;;
    *.fasta | *.fa)
        echo "FASTA file -> Counting sequences"
        # grep -c "^>" "$file"
        ;;
    *.bed)
        echo "BED file -> Sorting"
        # sort -k1,1 -k2,2n "$file"
        ;;
    *)
        echo "Unknown file type: $file"
        exit 1
        ;;
esac
SCRIPT
chmod +x /tmp/process_file.sh

/tmp/process_file.sh reads.fastq.gz
/tmp/process_file.sh aligned.bam
/tmp/process_file.sh variants.vcf
/tmp/process_file.sh genome.fa
```
