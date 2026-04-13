---
name: foundations-git-version-control
description: "By the end of this module, you will be able to: - Explain why version control is essential for reproducible research - Initialize a Git repository and track changes - Navigate project history (log, di"
tool_type: python
source_notebook: "Tier_0_Computational_Foundations/02_Git_Version_Control/01_git_version_control.ipynb"
---

# Module 0.2: Git Version Control for Scientists

*Source: Course notebook `Tier_0_Computational_Foundations/02_Git_Version_Control/01_git_version_control.ipynb`*

# Module 0.2: Git Version Control for Scientists

---

### Learning Objectives

By the end of this module, you will be able to:
- Explain why version control is essential for reproducible research
- Initialize a Git repository and track changes
- Navigate project history (log, diff, show)
- Use branches to develop features and experiments in isolation
- Collaborate with others using GitHub (push, pull, pull requests)
- Undo mistakes safely
- Write a proper `.gitignore` for bioinformatics projects

**Prerequisites:** Basic command-line familiarity (Module 0.1)

**Estimated time:** 2-3 hours

---

## How to use this notebook
1. Run cells top-to-bottom — each section builds a Git repository that later sections continue to use.
2. The `%%bash` cells run real git commands. Read the git output carefully — it tells you exactly what happened.
3. After each commit, run `git log --oneline` to see the growing history.
4. To experiment with undoing things, create a separate directory so you can make mistakes without worrying.

## Common stumbling points

- **Staging vs. committing**: `git add` does not save your work permanently — it just marks files for the *next* commit. Files in the working directory that are not staged will not be included in the commit.
- **Commit messages are forever**: Write messages that explain *why* you made the change, not just *what* changed. "Fix bug" is useless six months later; "Fix off-by-one error in exon boundary parsing" is useful.
- **Never commit large data files**: Git is not a data store. A single BAM file committed to a repo can make it permanently bloated and impossible to clone quickly.
- **`git reset --hard` is destructive**: Unlike most Git operations, `--hard` discards your changes with no undo. Use `--soft` to undo commits while keeping your edits.
- **`HEAD~1` means "one commit before the current one"**: `git diff HEAD~1` shows what changed in the last commit.

---

## 1. Setup and Configuration

First-time Git setup. You only need to do this once per computer.

```python
%%bash
# Tell Git who you are (used in commit metadata)
git config --global user.name "Your Name"
git config --global user.email "your.email@university.edu"

# Set the default branch name to 'main' (modern convention)
git config --global init.defaultBranch main

# Set your preferred text editor for commit messages
git config --global core.editor "nano"    # or "vim", or "code --wait" for VS Code

# Verify your configuration
git config --list
```

---

## 2. Core Concepts: The Three States

Every file in a Git repo exists in one of three states. Understanding this model
is the key to understanding everything else in Git.

```
+----------------+     git add     +----------------+    git commit    +----------------+
|    WORKING     | -------------> |    STAGING      | -------------> |   REPOSITORY   |
|   DIRECTORY    |                |   AREA (index)  |                |   (commits)    |
+----------------+                +----------------+                +----------------+

  Your files as         Files marked           Permanent snapshots
  you edit them         "ready to commit"      of your project
```

**Why a staging area?** It lets you commit a subset of your changes. For example,
you fixed a bug AND started a new feature in the same session. You can stage and
commit the bug fix first, then commit the new feature separately. This keeps your
history clean and each commit focused.

---

## 3. Creating a Repository

Two ways to start: initialize a new one, or clone an existing one.

```python
%%bash
# Option 1: Start a new repository from scratch
mkdir my_rna_seq_analysis
cd my_rna_seq_analysis
git init
# This creates a hidden .git/ directory that stores all version history

# Option 2: Clone an existing repository from GitHub
# git clone https://github.com/username/repository.git
# git clone https://github.com/username/repository.git my_local_name

# Check the status of your repo
git status
```

### Setting Up a Bioinformatics Project with Git

Let's create a real project structure and initialize it properly.

```python
%%bash
# Create project structure
mkdir -p /tmp/gene_expression_study/{data,scripts,results,docs}
cd /tmp/gene_expression_study
git init

# Create a README
cat > README.md << 'EOF'
# Gene Expression Study

Analysis of differential gene expression in breast cancer subtypes.

## Directory Structure
- `data/` - Raw and processed data (not tracked by Git)
- `scripts/` - Analysis scripts and pipelines
- `results/` - Generated results (not tracked by Git)
- `docs/` - Documentation and notes
EOF

echo "Project created and initialized!"
git status
```

---

## 4. The .gitignore File

In bioinformatics, most files should NOT be in Git. Sequence data, alignment files,
and results are too large (Git has a ~100 MB file limit on GitHub) and can be
regenerated from scripts + raw data.

**Rule of thumb: Track code and configuration. Do NOT track data and outputs.**

```python
%%bash
cd /tmp/gene_expression_study

# Create a comprehensive .gitignore for bioinformatics
cat > .gitignore << 'GITIGNORE'
# ==========================================
# Bioinformatics .gitignore
# ==========================================

# --- Large data files ---
*.fastq
*.fastq.gz
*.fq.gz
*.bam
*.bam.bai
*.sam
*.cram
*.bcf
*.sra
data/raw/

# --- Reference genomes ---
*.fa
*.fasta
*.fa.fai
*.dict

# --- Results and outputs (regenerated by scripts) ---
results/
*.log
*.tmp
*.out

# --- Python ---
__pycache__/
*.pyc
.ipynb_checkpoints/
*.egg-info/

# --- R ---
.Rhistory
.RData

# --- OS files ---
.DS_Store
Thumbs.db

# --- IDE files ---
.vscode/
.idea/
GITIGNORE

echo ".gitignore created"
cat .gitignore
```

**What TO track in Git:**

| Track | Do NOT Track |
|-------|-------------|
| Analysis scripts (`.py`, `.R`, `.sh`) | Raw data (`.fastq.gz`, `.bam`) |
| Pipeline definitions (Snakemake, Nextflow) | Reference genomes |
| Configuration files | Generated results |
| Documentation (README, methods) | Log files |
| Environment specs (`environment.yml`) | OS/IDE metadata |
| Small metadata files (sample sheets) | Anything >50 MB |

---

## 5. The Basic Workflow: add, commit, push

This is the workflow you will use dozens of times per day.

```python
%%bash
cd /tmp/gene_expression_study

# Step 1: Check what has changed
git status

# Step 2: Stage files for commit
git add README.md               # Stage a specific file
git add .gitignore              # Stage another file
# git add .                     # Stage ALL changes (use with caution)
# git add scripts/              # Stage all files in a directory

# Step 3: Commit with a descriptive message
git commit -m "Initialize project with README and .gitignore"

# Verify
git log --oneline
```

```python
%%bash
cd /tmp/gene_expression_study

# Let's add an analysis script and commit it
cat > scripts/deseq2_analysis.R << 'EOF'
# DESeq2 Differential Expression Analysis
# Author: Your Name
# Date: 2024-01-15

library(DESeq2)

# Load count matrix
counts <- read.csv("data/counts_matrix.csv", row.names=1)
metadata <- read.csv("data/sample_metadata.csv")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

# Run differential expression
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)

# Save results
write.csv(as.data.frame(results), "results/deseq2_results.csv")
EOF

git add scripts/deseq2_analysis.R
git commit -m "Add DESeq2 differential expression analysis script"
git log --oneline
```

### Writing Good Commit Messages

Commit messages are for your future self and your collaborators. They should explain
the *why*, not just the *what*.

```
BAD commit messages:             GOOD commit messages:

  "fix"                            "Fix off-by-one error in codon extraction"
  "update"                         "Add support for GFF3 annotation format"
  "stuff"                          "Optimize GC content calc (2x faster)"
  "final version"                  "Complete DESeq2 analysis with batch correction"
  "asdf"                           "Filter low-quality variants (QUAL < 20)"

Template:  <verb> <what> [context/reason]

Common verbs: Add, Fix, Update, Remove, Refactor, Optimize, Implement
```

For longer messages, use the editor instead of `-m`:

```bash
git commit    # Opens your editor for a multi-line message
```

Multi-line message convention:
```
Short summary (50 chars or less)

Detailed explanation of what changed and why. Wrap at 72 characters.
Reference relevant issues: Fixes #42
```

---

## 6. Viewing History and Differences

Git's real power is in navigating history -- comparing versions, finding when a bug
was introduced, or recovering old code.

```python
%%bash
cd /tmp/gene_expression_study

# View commit history
git log                          # Full history (press q to exit)
git log --oneline                # Compact: one line per commit
git log --oneline -5             # Last 5 commits
git log --graph --oneline        # Visual branch history
git log --stat                   # Show which files changed per commit

# View changes
git diff                         # Unstaged changes (working dir vs staging)
git diff --staged                # Staged changes (staging vs last commit)
git diff HEAD~1                  # Changes since the previous commit
# git diff abc123 def456         # Changes between two specific commits

# Show details of a specific commit
git show HEAD                    # Latest commit with full diff
# git show abc123                # Specific commit by hash
```

### Practical example: Tracking pipeline changes

```python
%%bash
cd /tmp/gene_expression_study

# Modify the analysis script -- change the significance threshold
sed -i '' 's/alpha = 0.05/alpha = 0.01/' scripts/deseq2_analysis.R 2>/dev/null ||
sed -i 's/alpha = 0.05/alpha = 0.01/' scripts/deseq2_analysis.R

# See what changed
echo "=== What changed? ==="
git diff scripts/deseq2_analysis.R

# Commit the change
git add scripts/deseq2_analysis.R
git commit -m "Tighten significance threshold from 0.05 to 0.01"

echo ""
echo "=== Commit history ==="
git log --oneline
```

---

## 7. Undoing Mistakes

Everyone makes mistakes. Git gives you multiple safety nets, from gentle to aggressive.

```python
%%bash
# === DISCARD CHANGES in working directory ===
# (File is modified but not yet staged)
# git restore filename.py               # Modern syntax (Git 2.23+)
# git checkout -- filename.py           # Older syntax

# === UNSTAGE a file ===
# (File is staged but not yet committed)
# git restore --staged filename.py      # Modern syntax
# git reset HEAD filename.py            # Older syntax

# === UNDO the last commit (keep changes) ===
# git reset --soft HEAD~1               # Changes remain staged
# git reset HEAD~1                      # Changes become unstaged

# === UNDO the last commit (DISCARD changes) ===
# git reset --hard HEAD~1               # DANGEROUS: changes are lost!

# === SAFELY undo a specific commit ===
# (Creates a NEW commit that reverses the old one -- safe for shared repos)
# git revert abc123

echo 'Summary of undo operations:'
echo '  restore file      = discard local changes'
echo '  restore --staged  = unstage'
echo '  reset --soft      = undo commit, keep changes staged'
echo '  reset --hard      = undo commit, DELETE changes (dangerous)'
echo '  revert            = new commit that undoes an old one (safe)'
```

### When to use each undo method

| Situation | Command | Safe? |
|-----------|---------|-------|
| Accidentally edited a file, want to restore it | `git restore file` | Yes (loses edits) |
| Staged wrong file, want to unstage | `git restore --staged file` | Yes |
| Just committed, forgot a file | `git reset --soft HEAD~1` | Yes |
| Just committed, want to discard everything | `git reset --hard HEAD~1` | DANGEROUS |
| Need to undo an old commit in shared repo | `git revert <hash>` | Yes |

---

## 8. Branches

Branches let you work on different things simultaneously without interfering with
each other. This is essential when:

- You want to try a different normalization method without breaking the working pipeline
- Multiple lab members are working on different analyses
- You need to fix a bug while a new feature is half-done

```
main -------*--------*---------------------*---*---> (stable code)
             \                             /   /
              *---*---*---*--- feature ---/   /      (new analysis method)
                   \                        /
                    *---*--- bugfix -------/         (quick fix)
```
