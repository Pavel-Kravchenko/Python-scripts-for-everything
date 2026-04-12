---
name: single-cell-scanpy
description: Single-cell analysis workflows with Scanpy/AnnData.
---

## When to Use

Use this atomic skill for focused work on **single-cell-scanpy** without bundling unrelated topics.

## Quick Reference

This skill was split from `clinical-modeling-workflows.md` to keep topics independent and self-contained.

## Core Patterns

Use the parent material below as the source reference, then keep implementations specific to this topic.

## Source Reference (from merged skill)

---
name: clinical-modeling-workflows
description: Molecular docking, clinical variant interpretation (ACMG), single-cell analysis (Scanpy), workflow engines (Snakemake/Nextflow), and CI/CD for bioinformatics
---

# Clinical Genomics, Modeling & Modern Workflows

## When to Use
- Classifying germline/somatic variants for clinical reporting (ACMG/AMP rules)
- Virtual screening or pose prediction for small-molecule drug candidates (AutoDock Vina)
- Building homology models or running MD simulation setup (GROMACS)
- Single-cell RNA-seq clustering and cell-type annotation (Scanpy/AnnData)
- Scaling bioinformatics pipelines across samples or HPC (Snakemake, Nextflow)
- Adding automated tests and GitHub Actions CI to bioinformatics code

## Quick Reference

### ACMG/AMP 5-Tier Classification
`Pathogenic (P) > Likely Pathogenic (LP) > VUS > Likely Benign (LB) > Benign (B)`

**Pathogenic criteria (accumulate points):**
| Strength | Codes | Key examples |
|----------|-------|--------------|
| Very Strong | PVS1 | Null variant (nonsense/frameshift/splice) in LoF-disease gene |
| Strong | PS1-PS4 | Same AA change as established P (PS1); de novo confirmed (PS2); functional assay damaging (PS3); elevated prevalence in cases (PS4) |
| Moderate | PM1-PM6 | Mutational hotspot (PM1); absent from gnomAD (PM2); protein length change (PM4); assumed de novo (PM6) |
| Supporting | PP1-PP5 | Co-segregation (PP1); computational evidence PP3; patient phenotype fits gene (PP4) |

**Benign criteria:**
| Strength | Codes | Key examples |
|----------|-------|--------------|
| Stand-alone | BA1 | MAF > 5% in any population |
| Strong | BS1-BS4 | Freq higher than expected for disorder (BS1); healthy adult (BS2); benign functional study (BS3) |
| Supporting | BP1-BP7 | In silico benign (BP4); synonymous, no splice effect (BP7) |

**Classification combining rules (Richards 2015, Table 5):**
- **Pathogenic**: PVS1 + (≥1 PS OR ≥2 PM OR 1 PM+1 PP OR ≥2 PP); or ≥2 PS; or 1 PS + (≥3 PM OR 2 PM+≥2 PP)
- **Likely Pathogenic**: PVS1+1 PM; or 1 PS+1-2 PM; or 1 PS+≥2 PP; or ≥3 PM; or 2 PM+≥2 PP; or 1 PM+≥4 PP
- **Benign**: BA1 alone; or ≥2 BS; or 1 BS+1 BP; or ≥2 BP → Likely Benign

### In Silico Predictors
| Tool | Score range | Damaging threshold |
|------|-------------|-------------------|
| SIFT | 0–1 (lower = damaging) | < 0.05 |
| PolyPhen-2 | 0–1 (higher = damaging) | > 0.908 probably; > 0.446 possibly |
| CADD (Phred) | higher = more deleterious | ≥ 20 (top 1%); ≥ 25 (top 0.3%) |
| REVEL | 0–1 | > 0.5 possible; > 0.75 likely pathogenic |
| AlphaMissense | 0–1 | > 0.564 likely pathogenic; < 0.34 likely benign |
| SpliceAI | delta 0–1 | > 0.2 suggestive; > 0.5 high confidence splice impact |

PP3 applies when ≥4/5 tools call damaging; BP4 applies when ≤1/5 call damaging.

### gnomAD Constraint Metrics
- **pLI > 0.9**: gene intolerant to heterozygous LoF (haploinsufficient)
- **LOEUF < 0.35**: strong LoF constraint
- **Missense o/e < 0.5**: missense-constrained gene

### ClinVar Review Stars
`****` Expert panel | `***` Expert review | `**` Multi-submitter no conflict | `*` Single submitter | `(0)` No criteria

### Molecular Docking (AutoDock Vina)
```bash
# Prepare receptor (add H, remove water, assign charges)
prepare_receptor4.py -r protein.pdb -o receptor.pdbqt
# Prepare ligand (3D coords + H)
prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt
# Run docking (config.txt defines center_x/y/z, size_x/y/z)
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
     --config config.txt --out output.pdbqt
# Scores in kcal/mol; more negative = stronger predicted binding
```
Scoring functions: force-field (AutoDock), empirical (Vina/ChemPLP), knowledge-based (DrugScore).

### Force Field Terms
`E_total = E_bond + E_angle + E_dihedral + E_electrostatic + E_vdW`
- Bond: harmonic `½ k_b(r−r₀)²`
- vdW: Lennard-Jones `4ε[(σ/r)¹²−(σ/r)⁶]`; equilibrium at r = 2^(1/6)σ
- Common FFs: AMBER (proteins/NA), CHARMM (proteins/lipids), OPLS-AA (organic molecules)

### GROMACS MD Workflow
```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p -ff amber99sb-ildn
gmx editconf -f protein.gro -o box.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
# Add ions → energy minimization (steepest descent) → NVT equil → NPT equil → production MD
```
Analysis: RMSD (`gmx rms`), RMSF (`gmx rmsf`), Rg (`gmx gyrate`), H-bonds (`gmx hbond`)

### Homology Modeling Quality Thresholds
- > 50% identity: high-quality model; 30–50%: reasonable backbone; < 30%: "twilight zone"
- Ramachandran: > 90% favored, < 0.5% outliers
- pLDDT (AlphaFold): > 90 high confidence; 70–90 moderate; < 50 likely disordered

## Key Patterns

### Scanpy Single-Cell Workflow
```
AnnData (cells × genes) → QC → Normalize → HVG → Scale → PCA → Neighbors → UMAP → Leiden → Markers
```

### Snakemake Core Structure
```python
# Snakefile
SAMPLES = ["s1", "s2", "s3"]

rule all:
    input: expand("results/{sample}.txt", sample=SAMPLES)

rule process:
    input:  "data/{sample}.fastq.gz"
    output: "results/{sample}.txt"
    threads: 8
    resources: mem_mb=16000
    shell: "mytool --threads {threads} {input} > {output}"
```

### GitHub Actions CI Template
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

## Code Templates

### ACMG Classifier
```python
CRITERIA_STRENGTH = {
    'PVS1': 'very_strong',
    'PS1': 'strong', 'PS2': 'strong', 'PS3': 'strong', 'PS4': 'strong',
    'PM1': 'moderate', 'PM2': 'moderate', 'PM4': 'moderate', 'PM6': 'moderate',
    'PP1': 'supporting', 'PP3': 'supporting', 'PP4': 'supporting',
    'BA1': 'stand_alone',
    'BS1': 'strong', 'BS2': 'strong', 'BS3': 'strong',
    'BP4': 'supporting', 'BP7': 'supporting',
}

def classify_variant(criteria: list[str]) -> str:
    pvs = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'very_strong')
    ps  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'strong'
              and c.startswith('PS'))
    pm  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'moderate')
    pp  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'supporting'
              and c.startswith('P'))
    ba  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'stand_alone')
    bs  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'strong'
              and c.startswith('BS'))
    bp  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'supporting'
              and c.startswith('BP'))

    if ba >= 1 or bs >= 2: return 'Benign'
    if (bs >= 1 and bp >= 1) or bp >= 2: return 'Likely Benign'
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
        return 'Pathogenic'
    if ps >= 2: return 'Pathogenic'
    if pvs >= 1 and pm >= 1: return 'Likely Pathogenic'
    if ps >= 1 and pm >= 1: return 'Likely Pathogenic'
    if pm >= 3 or (pm >= 2 and pp >= 2): return 'Likely Pathogenic'
    if pvs >= 1: return 'Likely Pathogenic'
    return 'VUS'
```

### ClinVar API Query
```python
import requests

def query_clinvar(term: str, retmax: int = 5) -> list[dict]:
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    ids = requests.get(f'{base}/esearch.fcgi',
                       params={'db': 'clinvar', 'term': term,
                               'retmax': retmax, 'retmode': 'json'}).json()
    id_list = ids['esearchresult']['idlist']
    if not id_list:
        return []
    summaries = requests.get(f'{base}/esummary.fcgi',
                             params={'db': 'clinvar',
                                     'id': ','.join(id_list),
                                     'retmode': 'json'}).json()
    result = summaries.get('result', {})
    return [{'uid': uid,
             'title': result[uid].get('title'),
             'sig': result[uid].get('clinical_significance', {}).get('description')}
            for uid in id_list if uid in result]
```

### Scanpy Standard Pipeline
```python
import scanpy as sc

adata = sc.datasets.pbmc3k()                      # or sc.read_10x_mtx(...)

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20]

# Normalize + HVG selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

# Dimensionality reduction
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)

# Clustering + markers
sc.tl.leiden(adata, resolution=0.5)               # adjust resolution: 0.2 coarse–2.0 fine
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=5)

# Cell-type annotation
cluster_annotations = {'0': 'CD4+ T cells', '1': 'CD14+ Monocytes', ...}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
```

### Snakemake RNA-seq Pipeline (key rules)
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
```bash
snakemake -n                          # dry run
snakemake --cores 8 --use-conda
snakemake --dag | dot -Tpng > dag.png
snakemake --cluster "sbatch -c {threads} --mem {resources.mem_mb}" --jobs 100
snakemake --forcerun align --cores 8  # force re-run one rule
```

### Nextflow Process (DSL2)
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
nextflow run main.nf -resume             # checkpoint restart
nextflow run nf-core/rnaseq -r 3.12.0 \
    --input samplesheet.csv --genome GRCh38 -profile docker
```

### pytest for Bioinformatics
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
import pytest
from bio_utils import gc_content, reverse_complement, find_motif

class TestGCContent:
    def test_balanced(self):       assert gc_content("ATGC") == 50.0
    def test_empty(self):          assert gc_content("") == 0.0
    def test_lowercase(self):      assert gc_content("atgc") == 50.0

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

### Capstone End-to-End Analysis Template
```python
# Pipeline: Unknown sequences → QC → BLAST → MSA → Phylogeny → Structure → Domains → GO
from Bio import SeqIO, Entrez, Phylo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

Entrez.email = "your@email.com"

# 1. QC: remove Ns, check GC, filter short seqs
clean = [r for r in records if len(r.seq) > 100 and r.seq.count('N') < 10]

# 2. BLAST identification
result = NCBIWWW.qblast("blastn", "nt", clean[0].seq)
records_blast = list(NCBIXML.parse(result))

# 3. MSA → phylogeny
cmd = ClustalOmegaCommandline(infile="seqs.fasta", outfile="aligned.fasta", auto=True)
cmd()
aln = AlignIO.read("aligned.fasta", "fasta")
calc = DistanceCalculator('identity')
dm = calc.get_distance(aln)
tree = DistanceTreeConstructor().nj(dm)
Phylo.draw(tree)
```

## Common Pitfalls

- **ACMG**: BA1 alone is stand-alone benign — no additional pathogenic criteria override it. VUS is the default; LP requires >90% certainty.
- **Variant predictors**: Concordance across ≥4 tools before applying PP3/BP4; no single tool is authoritative.
- **gnomAD filtering**: Use population-specific popmax, not global AF; recessive disorders tolerate higher carrier frequency (e.g., sickle cell *HBB* afr MAF ~6%).
- **Force fields**: Never mix parameters from different FFs. Always energy-minimize before MD.
- **Docking scores**: kcal/mol estimates rank poses for one ligand reliably; cross-ligand ranking is unreliable — validate top hits experimentally.
- **Scanpy**: Save `adata.raw = adata` *before* HVG subsetting — `rank_genes_groups` needs unscaled counts. Scale after `raw` assignment.
- **Leiden resolution**: Default 1.0 over-splits PBMC; 0.3–0.6 typically gives interpretable major types.
- **Snakemake wildcards**: Rule `all` must specify final targets explicitly; Snakemake works backward from outputs.
- **Nextflow `-resume`**: Requires same working directory and unchanged process signatures; always use it after partial failures.
- **pytest coordinates**: Test 0-based vs 1-based explicitly — off-by-one is the most common bioinformatics bug.
- **Capstone workflow**: QC and BLAST before alignment; never align raw uncleaned sequences — gap inflation ruins phylogeny.

---

## Related Skills
- `phylogenetics-evolution` — MSA methods, distance/ML tree construction, bootstrap support
- `python-advanced-sql` — SQL queries for clinical variant databases and annotation lookups
- `biopython-databases` — NCBI Entrez, BioPython SeqRecord, PDB structure parsing
- `biostatistics-r` — statistical tests for differential expression; PRS regression models
- `numpy-pandas-wrangling` — AnnData/DataFrame operations, variant annotation table joins


## Related Skills

- `single-cell-scanpy` (this file)
- `clinical-modeling-workflows` (legacy merged skill)
