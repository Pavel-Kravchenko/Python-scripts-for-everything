---
name: foundations-git-version-control
description: Git for bioinformatics — setup, staging workflow, .gitignore, history navigation, undo operations, and branches.
tool_type: bash
primary_tool: git
---

## When to Use
- Version-controlling analysis scripts, pipelines, and configs
- Collaborating on code across lab members
- Tracking which script version produced which result
- Safely experimenting with alternative methods (branches)

## Setup (one-time)

```bash
git config --global user.name "Your Name"
git config --global user.email "your@email.com"
git config --global init.defaultBranch main
git config --global core.editor "nano"   # or "vim", "code --wait"
```

## Core Workflow

```bash
git status                    # what's changed
git diff                      # unstaged changes
git diff --staged             # staged vs last commit
git add scripts/deseq2.R      # stage specific file (preferred over git add .)
git commit -m "Fix off-by-one in exon boundary parsing"
git log --oneline -10         # compact history
git log --graph --oneline     # visual branch history
```

## Bioinformatics .gitignore

```gitignore
# Large data — never track
*.fastq *.fastq.gz *.fq.gz *.bam *.bam.bai *.sam *.cram *.bcf *.sra
data/raw/
*.fa *.fasta *.fa.fai *.dict

# Generated outputs
results/ *.log *.tmp *.out

# Python
__pycache__/ *.pyc .ipynb_checkpoints/ *.egg-info/

# R
.Rhistory .RData

# OS / IDE
.DS_Store Thumbs.db .vscode/ .idea/
```

**Track**: scripts, pipeline definitions, configs, README, `environment.yml`, small sample sheets.
**Never track**: raw data, reference genomes, generated results, anything >50 MB.

## Undo Operations

| Situation | Command | Destructive? |
|-----------|---------|-------------|
| Discard file edits (unstaged) | `git restore file.py` | Loses edits |
| Unstage a file | `git restore --staged file.py` | No |
| Undo last commit, keep changes staged | `git reset --soft HEAD~1` | No |
| Undo last commit, keep changes unstaged | `git reset HEAD~1` | No |
| Undo last commit, discard changes | `git reset --hard HEAD~1` | **Yes** |
| Undo old commit in shared repo | `git revert <hash>` | No (new commit) |

## Branches

```bash
git branch feature/normalize-rpkm    # create
git checkout feature/normalize-rpkm  # switch (or: git switch feature/...)
git checkout -b feature/normalize-rpkm   # create + switch in one step

git merge feature/normalize-rpkm     # merge into current branch
git branch -d feature/normalize-rpkm # delete after merge
```

Use branches when:
- Testing a different normalization without breaking the working pipeline
- Multiple lab members work on different analyses simultaneously
- Fixing a bug while a new feature is half-done

## Commit Message Style

```
# Format: <verb> <what> [context]
# Verbs: Add, Fix, Update, Remove, Refactor, Optimize

Fix off-by-one error in exon boundary parsing
Add DESeq2 analysis with batch correction (LRT test)
Update STAR alignment to use 2-pass mode
Remove deprecated RPKM normalization function
```

Bad: `"fix"`, `"update"`, `"stuff"`, `"final version"`.

Multi-line (for significant changes):
```
Short summary (50 chars max)

Detailed explanation of what changed and why.
Wrap at 72 characters. Reference issues: Fixes #42.
```

## Pitfalls

- **Staging vs. committing**: `git add` marks files for the *next* commit — it does not save your work. Unstaged files are excluded from the commit.
- **Never commit large data files**: a single BAM committed to a repo permanently bloats it and makes `git clone` slow. Use `.gitignore` and data management tools (DVC, lfs) instead.
- **`git reset --hard` is irreversible**: unlike most Git operations, it discards working-directory changes with no undo. Use `--soft` unless you are certain.
- **`git revert` is safe for shared repos**: it creates a new commit that undoes the old one, preserving history. `reset --hard` on pushed commits rewrites history and breaks collaborators' clones.
- **`.gitignore` must be committed**: it has no effect until tracked. `git add .gitignore && git commit` before adding data files.
- **Commit messages are forever**: "Fix bug" is useless six months later. Write messages that explain *why*, not just *what*.
