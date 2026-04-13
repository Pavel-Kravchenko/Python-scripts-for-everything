---
name: foundations-linux-fundamentals
description: "By the end of this module, you will be able to: - Navigate the Linux file system with confidence - Create, move, copy, and delete files and directories - Use pipes and redirections to build data-proce"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/01_Linux_Basics/01_linux_fundamentals.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Module 0.1: Linux Fundamentals for Bioinformatics

*Source: Course notebook `Tier_0_Computational_Foundations/01_Linux_Basics/01_linux_fundamentals.ipynb`*

# Module 0.1: Linux Fundamentals for Bioinformatics

---

### Learning Objectives

By the end of this module, you will be able to:
- Navigate the Linux file system with confidence
- Create, move, copy, and delete files and directories
- Use pipes and redirections to build data-processing pipelines
- Search files with `grep`, `find`, and regular expressions
- Process structured text with `awk` and `sed`
- Work with compressed files (gzip, tar) -- the standard in genomics
- Manage processes and connect to remote servers
- Apply all of the above to real bioinformatics file formats (FASTA, FASTQ, BED, VCF, BAM)

**Prerequisites:** None. This module starts from scratch.

**Estimated time:** 3-4 hours

---

## How to use this notebook
1. Run cells top-to-bottom — later cells depend on earlier ones (e.g., files created in one cell are used in the next).
2. For every `%%bash` cell: read the command first, understand what it does, then run it.
3. **Experiment**: modify commands and re-run. Breaking things is the best way to learn Linux.
4. If a cell fails, check that previous cells were run and that the working directory is what you expect.

## Common stumbling points

- **Spaces around `=` in Bash**: `var = value` is wrong; `var=value` is right. The space makes Bash interpret `var` as a command name.
- **Relative vs. absolute paths**: `cat data/file.txt` works only from the right directory; `cat /home/user/project/data/file.txt` works from anywhere.
- **`rm` is permanent**: Linux has no trash bin. Always double-check before running `rm -rf`.
- **Pipes discard stderr**: `command1 | command2` only passes stdout. Add `2>&1` to also pipe error messages.
- **Wildcards expand before the command runs**: `rm *.fastq` becomes `rm sample1.fastq sample2.fastq ...` — the shell does the expansion, not `rm`.

---

## 1. Getting Help

Before we start learning commands, the most important skill: **how to look things up.**

Every standard Linux command comes with a manual page (`man`) and usually a `--help` flag.

```python
%%bash
# The man command opens the manual page for any command
# (Press 'q' to quit, '/' to search within the manual)
# man ls

# Quick help -- most commands support --help
ls --help 2>&1 | head -20

# whatis gives a one-line description
whatis grep
whatis awk
```

---

## 2. Navigation: Where Am I, What Is Here?

The Linux file system is a tree rooted at `/`. Everything -- devices, processes, your
home directory -- is a file or directory somewhere in this tree.

```
                     / (root)
                     |
     +---------------+----------------+
     |               |                |
    home            bin              usr
     |                                |
     +-- username                    local/bin
          |
          +-- projects/
          |     +-- rna_seq_analysis/
          |     +-- genome_assembly/
          +-- data/
                +-- raw_reads/
                +-- references/
```

### Key path concepts

| Symbol | Meaning | Example |
|--------|---------|---------|
| `/` | Root directory (top of the tree) | `cd /` |
| `~` | Your home directory | `cd ~` (same as `cd /home/username`) |
| `.` | Current directory | `./run_script.sh` |
| `..` | Parent directory (one level up) | `cd ..` |
| `-` | Previous directory (where you just were) | `cd -` |

```python
%%bash
# pwd -- Print Working Directory: where am I right now?
pwd

# ls -- List directory contents
ls              # Basic listing
ls -l           # Long format: permissions, owner, size, date
ls -a           # Include hidden files (starting with .)
ls -lah         # Long + all + human-readable sizes -- the combo you will use most
ls -1           # One file per line -- useful for piping
```

```python
%%bash
# cd -- Change Directory
cd ~              # Go to home directory
cd ..             # Go up one level
cd ../..          # Go up two levels
cd /usr/local/bin # Absolute path (starts with /)
cd projects       # Relative path (from current directory)
cd -              # Go back to the previous directory

# Tip: Press Tab to autocomplete paths -- saves enormous time
```

---

## 3. File and Directory Operations

### Creating

```python
%%bash
# mkdir -- Make Directory
mkdir my_project                    # Create a single directory
mkdir -p project/{data/{raw,processed},results,scripts,logs}
#   -p creates parent directories as needed AND does not error if they exist
#   Brace expansion {a,b} creates multiple siblings at once

# touch -- Create an empty file (or update timestamp of existing file)
touch README.md
touch data/sample_{01,02,03}.txt    # Creates three files via brace expansion
```

### Copying, Moving, Renaming, Deleting

```python
%%bash
# cp -- Copy
cp source.txt destination.txt       # Copy a file
cp -r source_dir/ dest_dir/         # Copy a directory recursively (-r is required)

# mv -- Move or Rename (same command for both)
mv old_name.txt new_name.txt        # Rename a file
mv file.txt target_directory/       # Move a file into a directory

# rm -- Remove (WARNING: there is no trash bin -- deletion is permanent)
rm file.txt                         # Delete a file
rm -r directory/                    # Delete a directory and everything inside
rm -rf directory/                   # Force-delete without confirmation
#   NEVER run: rm -rf /   or   rm -rf ~   -- you will destroy your system/data
```

### Bioinformatics context: Project directory conventions

Most bioinformatics projects follow a numbered-folder convention to make the
analysis order obvious:

```bash
mkdir -p RNA_Seq_Project/{00_raw_data,01_qc,02_trimmed,03_aligned,04_counts,05_analysis,scripts,logs}
```

This mirrors the actual pipeline steps: raw reads go in `00_raw_data`, quality
control reports in `01_qc`, and so on. Anyone opening the project folder immediately
sees the workflow.

---

## 4. Viewing File Contents

You will spend a huge amount of time inspecting files -- checking whether a FASTQ file
downloaded correctly, previewing a VCF, or reviewing pipeline logs.

```python
%%bash
# cat -- Print entire file to terminal (fine for small files)
cat small_file.txt

# less -- Interactive viewer with scrolling (press q to quit, / to search)
# less large_file.txt

# head / tail -- View the beginning or end of a file
head -n 10 file.txt         # First 10 lines
tail -n 10 file.txt         # Last 10 lines
head -n 4 reads.fastq       # First read in a FASTQ file (4 lines per read)

# tail -f -- Follow a file in real time (indispensable for monitoring logs)
# tail -f alignment.log     # Watch a log file as it grows

# wc -- Word/line/character count
wc -l file.txt              # Count lines
wc -c file.txt              # Count bytes
```

### Bioinformatics file formats you will encounter

| Format | Extension | What it contains | Typical size |
|--------|-----------|------------------|--------------|
| FASTA | `.fa`, `.fasta` | Sequences (genome, proteins) | MB to GB |
| FASTQ | `.fq`, `.fastq` | Sequencing reads + quality scores | GB (compressed) |
| SAM/BAM | `.sam`, `.bam` | Read alignments (BAM = binary compressed SAM) | GB |
| BED | `.bed` | Genomic intervals (chr, start, end) | KB to MB |
| VCF | `.vcf` | Variant calls | MB |
| GFF/GTF | `.gff`, `.gtf` | Gene annotations | MB |

**FASTA format example:**
```
>seq1 Homo sapiens BRCA1 mRNA
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAAT
GCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGG
>seq2 Mus musculus TP53 mRNA
ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCCCTC
```

**FASTQ format example (4 lines per read):**
```
@SEQ_ID_001
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### Text Editors in the Terminal

| Editor | Difficulty | Best for |
|--------|------------|----------|
| `nano` | Beginner | Quick edits. Ctrl+O to save, Ctrl+X to exit. |
| `vim`  | Steep learning curve | Power editing, scripting, servers. |
| `emacs`| Steep learning curve | Full IDE-like experience. |

**Vim survival guide** (you WILL encounter vim on servers):

```
vim filename          # Open file
i                     # Enter INSERT mode (now you can type)
Esc                   # Return to NORMAL mode
:w                    # Save (write)
:q                    # Quit
:wq                   # Save and quit
:q!                   # Quit without saving
/pattern              # Search forward for pattern
n / N                 # Next / previous search match
dd                    # Delete current line
u                     # Undo
Ctrl+R                # Redo
gg / G                # Go to beginning / end of file
:21                   # Go to line 21
:%s/old/new/g         # Replace all occurrences of 'old' with 'new'
```

---

## 5. Pipes and Redirections

This is the Unix superpower. Each command does one thing well; pipes (`|`) chain them
together into powerful data-processing pipelines.

```
Data flow:

  command1 | command2 | command3 > output.txt
      |         |          |          |
      v         v          v          v
  [stdout] -> [stdin] -> [stdin] -> [file]

Redirection operators:
  >    Redirect stdout to file (overwrite)
  >>   Redirect stdout to file (append)
  2>   Redirect stderr to file
  2>>  Redirect stderr to file (append)
  2>&1 Redirect stderr to same place as stdout
  <    Use file as stdin
  |    Pipe stdout of left command to stdin of right command
```

```python
%%bash
# Basic redirection
echo "Sample_001" > sample_list.txt      # Create/overwrite file
echo "Sample_002" >> sample_list.txt     # Append to file

# Separate stdout and stderr
# command > output.log 2> error.log

# Combine stdout and stderr into one file
# command > all_output.log 2>&1

# Pipe chains -- the bread and butter of command-line bioinformatics
# Count unique genes in a list:
# cat genes.txt | sort | uniq | wc -l

# Count reads in a FASTQ file (4 lines per read):
# cat sample.fastq | wc -l | awk '{print $1/4}'

# Extract sequence IDs from FASTA, clean them up, save to file:
# grep "^>" sequences.fasta | cut -d' ' -f1 | sed 's/>//' > sequence_ids.txt
```

### Bioinformatics pipeline examples

```bash
# Count variants per chromosome from a VCF file
grep -v "^#" variants.vcf | cut -f1 | sort | uniq -c | sort -rn

# Find unique chromosomes in a BED file, sorted by frequency
cut -f1 regions.bed | sort | uniq -c | sort -rn

# Estimate GC content of a genome
grep -v "^>" genome.fa | tr -d '\n' | \
  awk '{gc=gsub(/[GC]/,"",$0); at=gsub(/[AT]/,"",$0); print gc/(gc+at)*100"%"}'

# Count reads in a compressed FASTQ
zcat sample.fastq.gz | wc -l | awk '{print $1/4, "reads"}'
```

---

## 6. Searching: grep and find

### grep -- Search Within Files

`grep` is your go-to tool for finding patterns in text. The name stands for
"Global Regular Expression Print."

```python
%%bash
# Basic grep usage
# grep "pattern" filename

# Useful flags:
# grep -c "pattern" file      # Count matching lines (not matches)
# grep -i "pattern" file      # Case-insensitive
# grep -v "pattern" file      # Invert: show lines that do NOT match
# grep -r "pattern" dir/      # Recursive search through all files in directory
# grep -n "pattern" file      # Show line numbers
# grep -w "gene" file         # Match whole word only (not "genes" or "genebank")
# grep -l "pattern" *.txt     # List filenames that contain the pattern
# grep -A 3 "pattern" file    # Show 3 lines After each match
# grep -B 2 "pattern" file    # Show 2 lines Before each match

# Search compressed files without decompressing
# zgrep "PASS" variants.vcf.gz
```

#### grep in bioinformatics

```bash
# Count sequences in a FASTA file
grep -c "^>" proteins.fasta

# Extract header lines from FASTA
grep "^>" proteins.fasta

# Find all lines mentioning BRCA1 in an annotation file
grep "BRCA1" gencode.gtf

# Count TATA boxes in promoter sequences
grep -c "TATAAA" promoters.fa

# Remove comment lines from a VCF to get only data rows
grep -v "^#" variants.vcf

# Find which samples contain a specific variant
grep -l "rs12345" *.vcf

# Search for a motif and show the full FASTA header above it
grep -B 1 "GAATTC" sequences.fasta    # EcoRI recognition site
```
