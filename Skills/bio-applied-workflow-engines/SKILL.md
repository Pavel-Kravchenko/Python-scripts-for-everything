---
name: bio-applied-workflow-engines
description: "1. Explain the DAG-based execution model shared by all modern workflow engines 2. Write **production-quality Snakefiles** with config files, wildcards, params, logs, benchmarks, conda environments, an"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/12_Modern_Workflows/02_workflow_engines.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Workflow Engines: Snakemake and nf-core / Nextflow

*Source: Course notebook `Tier_3_Applied_Bioinformatics/12_Modern_Workflows/02_workflow_engines.ipynb`*


**Tier 3 – Applied Bioinformatics**

---

## Learning Objectives

By the end of this notebook you will be able to:

1. Explain the DAG-based execution model shared by all modern workflow engines
2. Write **production-quality Snakefiles** with config files, wildcards, params, logs, benchmarks, conda environments, and containers
3. Run, profile, and scale Snakemake on a **cluster / cloud** (SLURM, AWS, Google Cloud)
4. Understand the **Nextflow DSL2** dataflow model and write basic processes
5. Use **nf-core** community pipelines, explore their parameter space, and customise them with profiles and config files
6. Choose between Snakemake and Nextflow for a given project

---

## Contents

1. [Why Workflow Engines?](#1.-Why-Workflow-Engines?)
2. [Snakemake Fundamentals](#2.-Snakemake-Fundamentals)
3. [Snakemake — Advanced Features](#3.-Snakemake-Advanced-Features)
4. [Snakemake — Cluster and Cloud](#4.-Snakemake-Cluster-and-Cloud)
5. [Nextflow / DSL2 Fundamentals](#5.-Nextflow-/-DSL2-Fundamentals)
6. [nf-core Community Pipelines](#6.-nf-core-Community-Pipelines)
7. [Snakemake vs Nextflow Cheat-Sheet](#7.-Snakemake-vs-Nextflow)
8. [Exercises](#8.-Exercises)
9. [Summary](#9.-Summary)

---

## 1. Why Workflow Engines?

### 1.1 The shell-script problem

A typical genomics project evolves into a tangle of shell scripts:

```bash
#!/bin/bash
fastqc raw/*.fastq.gz -o qc/
trimmomatic PE raw/s1_R1.fastq.gz raw/s1_R2.fastq.gz     trimmed/s1_R1.fq.gz trimmed/s1_R1_unpaired.fq.gz     trimmed/s1_R2.fq.gz trimmed/s1_R2_unpaired.fq.gz     ILLUMINACLIP:adapters.fa:2:30:10 MINLEN:36
bwa mem -t 8 ref/genome.fa trimmed/s1_R1.fq.gz trimmed/s1_R2.fq.gz |     samtools sort -o aligned/s1.bam
...
```python

Problems:
- **No resume**: reruns everything even if step 2 failed halfway through step 5
- **No parallelism**: samples processed sequentially
- **No resource control**: cannot request 8 CPUs for BWA and 1 for sorting
- **Not reproducible**: hard-coded paths, tool versions not pinned
- **No provenance**: no record of which command produced which file

### 1.2 What workflow engines add

| Feature | Shell script | Snakemake | Nextflow |
|---------|-------------|-----------|----------|
| DAG-based dependency tracking | ✗ | ✓ | ✓ |
| Automatic resume / re-run | ✗ | ✓ | ✓ |
| Parallelism | manual | automatic | automatic |
| Resource management | ✗ | ✓ | ✓ |
| Conda / containers | manual | built-in | built-in |
| Cluster / cloud | manual | ✓ | ✓ |
| Provenance / audit trail | ✗ | partial | ✓ (full) |
| Community pipelines | ✗ | catalog | nf-core |

```python
# ── Illustrate a minimal pipeline DAG (text representation) ──────────────
pipeline_dag = {
    "raw FASTQ": ["FastQC", "Trimmomatic"],
    "FastQC":    ["MultiQC"],
    "Trimmomatic": ["BWA MEM"],
    "BWA MEM":   ["samtools sort"],
    "samtools sort": ["samtools index", "GATK HaplotypeCaller"],
    "samtools index": [],
    "GATK HaplotypeCaller": ["VQSR", "VEP annotation"],
    "VQSR": ["VEP annotation"],
    "VEP annotation": ["final VCF"],
    "MultiQC": ["final QC report"],
}

def topo_sort(graph):
    """Kahn's topological sort for a DAG dict."""
    in_degree = {n: 0 for n in graph}
    for node, neighbors in graph.items():
        for nb in neighbors:
            if nb not in in_degree:
                in_degree[nb] = 0
            in_degree[nb] += 1
    queue = [n for n, d in in_degree.items() if d == 0]
    order = []
    while queue:
        n = queue.pop(0)
        order.append(n)
        for nb in graph.get(n, []):
            in_degree[nb] -= 1
            if in_degree[nb] == 0:
                queue.append(nb)
    return order

order = topo_sort(pipeline_dag)
print("Execution order (topological sort):")
for i, step in enumerate(order, 1):
    print(f"  {i:>2}. {step}")
```python

---

## 2. Snakemake Fundamentals

### 2.1 Installation

```bash
# Recommended: conda / mamba
conda create -n snake -c conda-forge -c bioconda snakemake
conda activate snake

# pip (no conda dependencies)
pip install snakemake
```python

### 2.2 Core concepts

| Concept | Description |
|---------|-------------|
| **Rule** | Unit of work: defines input, output, and the shell/script command |
| **Wildcard** | `{sample}` – Snakemake infers which rule to run by matching filenames |
| **DAG** | Directed Acyclic Graph built by backward-chaining from the target |
| **Snakefile** | The recipe file (Python + rule syntax) |
| **`all` rule** | The default target rule — lists all final outputs |

### 2.3 Minimal Snakefile anatomy

```python
# Snakefile

SAMPLES = ["sample1", "sample2", "sample3"]

rule all:
    input:
        expand("results/vcf/{sample}.vcf.gz", sample=SAMPLES),
        "results/qc/multiqc_report.html"

rule fastqc:
    input:  "data/{sample}_R1.fastq.gz"
    output: "results/qc/{sample}_fastqc.html"
    shell:  "fastqc {input} --outdir results/qc/"

rule trim:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        r1 = "trimmed/{sample}_R1.fq.gz",
        r2 = "trimmed/{sample}_R2.fq.gz"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule bwa_mem:
    input:
        r1 = "trimmed/{sample}_R1.fq.gz",
        r2 = "trimmed/{sample}_R2.fq.gz",
        ref = "ref/hg38.fa"
    output: "aligned/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} "
        "| samtools sort -o {output}"
```python

```python
# ── Simulate Snakemake's rule-matching logic ─────────────────────────────
import re

def match_rule_output(output_pattern, target_file):
    """
    Check if a target filename matches a rule's output pattern,
    and return the wildcard dict if it does.
    Handles {wildcard} patterns like Snakemake.
    """
    # Convert Snakemake pattern to regex: {name} → (?P<name>.+)
    regex = re.sub(r'\{(\w+)\}', r'(?P<\1>[^/]+)', re.escape(output_pattern))
    regex = regex.replace('\/', '/')  # unescape slashes
    m = re.fullmatch(regex, target_file)
    if m:
        return m.groupdict()
    return None

# Define some rules (simplified)
rules = [
    ("fastqc",  "results/qc/{sample}_fastqc.html"),
    ("trim",    "trimmed/{sample}_R1.fq.gz"),
    ("bwa_mem", "aligned/{sample}.bam"),
    ("call",    "results/vcf/{sample}.vcf.gz"),
]

# Test targets
targets = [
    "results/qc/sample1_fastqc.html",
    "aligned/sample2.bam",
    "results/vcf/sample3.vcf.gz",
    "some/unknown/file.txt",
]

print(f"{'Target':<40}  {'Matching rule':<12}  Wildcards")
print('-' * 75)
for target in targets:
    matched = None
    for rule_name, pattern in rules:
        wc = match_rule_output(pattern, target)
        if wc is not None:
            matched = (rule_name, wc)
            break
    if matched:
        print(f"{target:<40}  {matched[0]:<12}  {matched[1]}")
    else:
        print(f"{target:<40}  (no match)")
```python

---

## 3. Snakemake — Advanced Features

### 3.1 Config files

Separate parameters from logic using a YAML config file:

```yaml
# config/config.yaml
samples:
  - sample1
  - sample2
  - sample3

reference: ref/hg38.fa
adapters:  ref/adapters.fa

trimming:
  min_length: 36
  quality:    20

calling:
  min_base_quality: 20
  min_mapping_quality: 30
```python

Load in Snakefile:

```python
configfile: "config/config.yaml"

SAMPLES = config["samples"]
REF     = config["reference"]
```python

### 3.2 params, log, benchmark

```python
rule bwa_mem:
    input:
        r1  = "trimmed/{sample}_R1.fq.gz",
        r2  = "trimmed/{sample}_R2.fq.gz",
        ref = config["reference"]
    output:
        bam = "aligned/{sample}.bam"
    params:
        rg  = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/bwa_mem/{sample}.tsv"
    threads: 8
    resources:
        mem_mb = 16000
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.ref} "
        "{input.r1} {input.r2} | samtools sort -o {output.bam}) "
        "2> {log}"
```python

**log** files capture stderr/stdout for debugging without breaking Snakemake's
output-tracking.  
**benchmark** writes a TSV with wall-clock time, CPU time, and peak RSS per rule
invocation — invaluable for profiling.

### 3.3 Conda environments per rule

```python
rule fastqc:
    input:  "data/{sample}_R1.fastq.gz"
    output: "results/qc/{sample}_fastqc.html"
    conda:  "envs/qc.yaml"     # ← per-rule conda env
    shell:  "fastqc {input} --outdir results/qc/"
```python

```yaml
# envs/qc.yaml
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - fastqc=0.12.1
  - multiqc=1.19
```python

Run with `--use-conda` to activate:
```bash
snakemake --use-conda --cores 8
```python

### 3.4 Singularity / Apptainer containers

```python
rule gatk_haplotypecaller:
    input:  bam="aligned/{sample}.bam", ref=config["reference"]
    output: gvcf="gvcf/{sample}.g.vcf.gz"
    container:
        "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} "
        "-O {output.gvcf} -ERC GVCF"
```python

Run with:
```bash
snakemake --use-singularity --cores 8
```python

### 3.5 Checkpoints (dynamic output)

When you don't know the number of output files in advance (e.g., genome
assembly splitting), use checkpoints:

```python
checkpoint split_by_chromosome:
    input: "assembly.fa"
    output: directory("chromosomes/")
    shell: "csplit assembly.fa ... "

def aggregate_chromosomes(wildcards):
    checkpoint_output = checkpoints.split_by_chromosome.get(**wildcards).output[0]
    chroms = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.fa")).chrom
    return expand("annotated/{chrom}.gff3", chrom=chroms)

rule all:
    input: aggregate_chromosomes
```python

```python
# ── Build a Snakemake-style workflow graph from rules ────────────────────
from dataclasses import dataclass, field
from typing import Dict, List, Optional

@dataclass
class Rule:
    name: str
    inputs: List[str]
    outputs: List[str]
    threads: int = 1
    mem_mb: int = 1000
    has_log: bool = False
    has_benchmark: bool = False
    conda_env: Optional[str] = None
    container: Optional[str] = None

def build_dependency_graph(rules: List[Rule]) -> Dict[str, List[str]]:
    """
    Build a rule → [upstream rules] dependency graph by matching
    output patterns to input patterns (simplified exact-match version).
    """
    # Map each output pattern to its rule
    output_to_rule = {}
    for rule in rules:
        for out in rule.outputs:
            output_to_rule[out] = rule.name
    # For each rule, find which rules produce its inputs
    graph = {rule.name: [] for rule in rules}
    for rule in rules:
        for inp in rule.inputs:
            producer = output_to_rule.get(inp)
            if producer and producer != rule.name:
                graph[rule.name].append(producer)
    return graph

# ── Example variant-calling pipeline ─────────────────────────────────────
variant_pipeline = [
    Rule("fastqc",       inputs=["data/{sample}_R1.fastq.gz"],
                         outputs=["results/qc/{sample}_fastqc.html"],
                         conda_env="envs/qc.yaml"),
    Rule("trim",         inputs=["data/{sample}_R1.fastq.gz",
                                 "data/{sample}_R2.fastq.gz"],
                         outputs=["trimmed/{sample}_R1.fq.gz",
                                  "trimmed/{sample}_R2.fq.gz"],
                         conda_env="envs/fastp.yaml"),
    Rule("bwa_mem",      inputs=["trimmed/{sample}_R1.fq.gz",
                                 "trimmed/{sample}_R2.fq.gz"],
                         outputs=["aligned/{sample}.bam"],
                         threads=8, mem_mb=16000,
                         has_log=True, has_benchmark=True),
    Rule("mark_dups",    inputs=["aligned/{sample}.bam"],
                         outputs=["dedup/{sample}.bam"],
                         has_log=True, container="docker://broadinstitute/gatk:4.5.0.0"),
    Rule("haplotype",    inputs=["dedup/{sample}.bam"],
                         outputs=["gvcf/{sample}.g.vcf.gz"],
                         threads=4, has_log=True,
                         container="docker://broadinstitute/gatk:4.5.0.0"),
]

print(f"Pipeline: {len(variant_pipeline)} rules")
print()
print(f"{'Rule':<15}  {'Threads':>7}  {'RAM':>6}  {'Conda':>5}  {'Container':>9}  {'Log':>4}  {'Bench':>5}")
print('-' * 70)
for r in variant_pipeline:
    print(f"{r.name:<15}  {r.threads:>7}  {r.mem_mb:>5}M  "          f"{'✓' if r.conda_env else '':>5}  "          f"{'✓' if r.container else '':>9}  "          f"{'✓' if r.has_log else '':>4}  "          f"{'✓' if r.has_benchmark else '':>5}")
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
