---
name: bio-applied-workflow-engines
description: Snakemake and Nextflow/nf-core workflow patterns for genomics pipelines — rules, wildcards, config, cluster execution, and comparison table.
tool_type: python
primary_tool: Python
---

# Workflow Engines: Snakemake and Nextflow

## Why Workflow Engines?

| Feature | Shell script | Snakemake | Nextflow |
|---------|-------------|-----------|----------|
| DAG-based dependency tracking | No | Yes | Yes |
| Automatic resume on failure | No | Yes | Yes |
| Parallelism | Manual | Automatic | Automatic |
| Resource management (CPU/RAM) | No | Yes | Yes |
| Conda / containers | Manual | Built-in | Built-in |
| Cluster / cloud | Manual | Yes | Yes |
| Provenance / audit trail | No | Partial | Full |
| Community pipelines | No | Catalog | nf-core |

## Snakemake

### Minimal Snakefile
```python
SAMPLES = ["sample1", "sample2", "sample3"]

rule all:
    input:
        expand("results/vcf/{sample}.vcf.gz", sample=SAMPLES),
        "results/qc/multiqc_report.html"

rule fastqc:
    input:  "data/{sample}_R1.fastq.gz"
    output: "results/qc/{sample}_fastqc.html"
    conda:  "envs/qc.yaml"
    shell:  "fastqc {input} --outdir results/qc/"

rule bwa_mem:
    input:
        r1  = "trimmed/{sample}_R1.fq.gz",
        r2  = "trimmed/{sample}_R2.fq.gz",
        ref = config["reference"]
    output: "aligned/{sample}.bam"
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    log:       "logs/bwa_mem/{sample}.log"
    benchmark: "benchmarks/bwa_mem/{sample}.tsv"
    threads: 8
    resources: mem_mb = 16000
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.ref} "
        "{input.r1} {input.r2} | samtools sort -o {output}) 2> {log}"
```

### Config file pattern
```yaml
# config/config.yaml
samples: [sample1, sample2, sample3]
reference: ref/hg38.fa
trimming:
  min_length: 36
  quality: 20
```
```python
# Snakefile top
configfile: "config/config.yaml"
SAMPLES = config["samples"]
```

### Conda environment per rule
```yaml
# envs/qc.yaml
channels: [bioconda, conda-forge]
dependencies:
  - fastqc=0.12.1
  - multiqc=1.19
```
```bash
snakemake --use-conda --cores 8
```

### Container per rule
```python
rule gatk_haplotypecaller:
    container: "docker://broadinstitute/gatk:4.5.0.0"
    shell: "gatk HaplotypeCaller ..."
```
```bash
snakemake --use-singularity --cores 8
```

### Checkpoint (dynamic output count)
```python
checkpoint split_by_chromosome:
    input:  "assembly.fa"
    output: directory("chromosomes/")
    shell:  "csplit assembly.fa ..."

def aggregate_chromosomes(wildcards):
    out = checkpoints.split_by_chromosome.get(**wildcards).output[0]
    chroms = glob_wildcards(os.path.join(out, "{chrom}.fa")).chrom
    return expand("annotated/{chrom}.gff3", chrom=chroms)

rule all:
    input: aggregate_chromosomes
```

### Cluster / cloud execution
```bash
# SLURM profile
snakemake --profile slurm --jobs 100

# Google Life Sciences
snakemake --google-lifesciences --default-remote-prefix mybucket/results
```

## Nextflow / DSL2

### Process + workflow skeleton
```nextflow
nextflow.enable.dsl=2

process FASTQC {
    conda 'bioconda::fastqc=0.12.1'
    input:  path reads
    output: path "*.html"
    script:
    """
    fastqc ${reads} --outdir .
    """
}

process BWA_MEM {
    cpus   8
    memory '16 GB'
    input:
        tuple val(sample_id), path(r1), path(r2)
        path ref
    output: tuple val(sample_id), path("${sample_id}.bam")
    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} \
        | samtools sort -o ${sample_id}.bam
    """
}

workflow {
    reads_ch = Channel.fromFilePairs("data/*_R{1,2}.fastq.gz")
    FASTQC(reads_ch.map { id, files -> files }.flatten())
    BWA_MEM(reads_ch, file(params.reference))
}
```

### nf-core pipelines
```bash
# Install and run curated pipelines
pip install nf-core
nf-core launch nf-core/rnaseq   # interactive config
nf-core launch nf-core/sarek    # somatic/germline variant calling
nf-core launch nf-core/chipseq

# Run with test profile
nextflow run nf-core/rnaseq -profile test,docker
```

## Pitfalls

- **Snakemake runs from the Snakefile directory**: use relative paths or `workflow.basedir` for portability.
- **Wildcards are greedy by default**: `{sample}` matches slashes; use `{sample,[^/]+}` to restrict.
- **`rule all` must list all final outputs**: Snakemake builds the DAG backward from this target.
- **Missing `log:` directive**: without it, failed jobs produce no traceable stderr; always capture logs.
- **`benchmark:` is per-rule not per-run**: benchmark files are overwritten on re-runs; include `{sample}` in the path.
- **Nextflow channel consumed once**: a channel is exhausted after use; use `.into{}` or `.multiMap{}` to branch.
- **nf-core pipelines require `--input` samplesheet**: format varies by pipeline; check the docs for the exact TSV/CSV schema.
- **Coordinate systems**: BED is 0-based half-open; VCF/GFF are 1-based inclusive — mixing them causes off-by-one errors in custom rules.
