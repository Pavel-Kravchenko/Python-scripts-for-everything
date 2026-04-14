---
name: bioinformatics-workflows-cicd
description: Snakemake, Nextflow DSL2, GitHub Actions CI, and pytest patterns for bioinformatics pipelines.
tool_type: python
primary_tool: snakemake
---

## When to Use
- Scaling analyses across samples on HPC (Snakemake, Nextflow)
- Adding CI/CD to bioinformatics repos (GitHub Actions + pytest)
- Containerising pipeline processes (Docker/Singularity via Nextflow)

## Snakemake

### Core Structure
```python
SAMPLES = config['samples']   # from config.yaml

rule all:
    input:
        "results/counts_matrix.csv",
        expand("qc/{sample}_fastqc.html", sample=SAMPLES)

rule align:
    input:
        fastq = "data/{sample}.fastq.gz",
        index = "reference/star_index"
    output: bam = "aligned/{sample}.bam"
    threads: 8
    resources: mem_mb=32000
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq} "
        "--readFilesCommand zcat --runThreadN {threads} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix aligned/{wildcards.sample}_ && "
        "mv aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}"
```

### CLI Quick Reference
```bash
snakemake -n                                          # dry run
snakemake --cores 8 --use-conda
snakemake --dag | dot -Tpng > dag.png
snakemake --cluster "sbatch -c {threads} --mem {resources.mem_mb}" --jobs 100
snakemake --forcerun align --cores 8                  # force re-run one rule
```

## Nextflow DSL2

```groovy
nextflow.enable.dsl=2
params.reads  = "data/*_{1,2}.fastq.gz"
params.outdir = "results"

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc"
    cpus 2; memory '4 GB'
    container 'biocontainers/fastqc:v0.11.9'
    input:  tuple val(sample_id), path(reads)
    output: path "*_fastqc.{html,zip}"
    script: "fastqc ${reads} --threads ${task.cpus}"
}

workflow {
    Channel.fromFilePairs(params.reads, checkIfExists: true) | FASTQC
}
```

```bash
nextflow run main.nf -with-docker
nextflow run main.nf -resume              # checkpoint restart
nextflow run nf-core/rnaseq -r 3.12.0 \
    --input samplesheet.csv --genome GRCh38 -profile docker
```

## GitHub Actions CI

```yaml
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v5
      with: { python-version: "3.11" }
    - run: pip install pytest pytest-cov && pip install -r requirements.txt
    - run: pytest --cov=src --cov-report=xml
```

## pytest Patterns for Bioinformatics

```python
# conftest.py
import pytest
from pathlib import Path

@pytest.fixture
def sample_fasta_file(tmp_path):
    content = ">seq1\nATGCATGC\n>seq2\nGCGCGCGC\n"
    p = tmp_path / "test.fasta"
    p.write_text(content)
    return p

# test_bio_utils.py
@pytest.mark.parametrize("seq,expected", [
    ("GGGG", 100.0), ("AAAA", 0.0), ("ATGC", 50.0),
])
def test_gc_parametrized(seq, expected):
    assert gc_content(seq) == expected

def test_reverse_complement_palindrome():
    assert reverse_complement("GAATTC") == "GAATTC"   # EcoRI site

def test_motif_overlapping():
    assert find_motif("AAAA", "AA") == [1, 2, 3]      # 1-based, overlapping
```

## Pitfalls

- **Snakemake wildcards**: Rule `all` must list final targets explicitly; Snakemake resolves the DAG backward from outputs.
- **Snakemake `config`**: Always pass `--configfile config.yaml`; `config['samples']` raises `KeyError` if the file is missing.
- **Nextflow `-resume`**: Requires the same working directory and unchanged process signatures; always use it after partial failures — rerunning without it reruns everything.
- **Nextflow DSL2 channels**: `fromFilePairs` expects exactly two files matching the glob; odd counts silently drop the sample.
- **Container versions**: Pin container tags in `publishDir` processes (`v0.11.9` not `latest`) — `latest` can break reproducibility silently.
- **pytest coordinates**: Test 0-based vs 1-based boundaries explicitly — off-by-one is the most common bioinformatics bug in tests.
- **CI secrets**: Never hard-code API keys; use `${{ secrets.MY_KEY }}` in GitHub Actions env vars.

## Related Skills
- `clinical-modeling-workflows` — ACMG, Scanpy, docking (split from this skill)
- `bio-workflow-management-snakemake-workflows` — advanced Snakemake patterns
- `bio-workflow-management-nextflow-pipelines` — advanced Nextflow/nf-core
