---
name: clinical-modeling-workflows
description: ACMG/AMP variant classification, in silico predictors, AutoDock Vina docking, GROMACS MD setup, and Scanpy single-cell analysis.
tool_type: python
primary_tool: scanpy
---

## When to Use
- Classifying germline/somatic variants for clinical reporting (ACMG/AMP)
- Virtual screening / pose prediction (AutoDock Vina)
- MD simulation setup (GROMACS)
- Single-cell RNA-seq clustering and cell-type annotation (Scanpy/AnnData)

## ACMG/AMP Classification

`Pathogenic (P) > Likely Pathogenic (LP) > VUS > Likely Benign (LB) > Benign (B)`

### Criteria Strengths
| Strength | Pathogenic codes | Key triggers |
|----------|-----------------|-------------|
| Very Strong | PVS1 | Null variant (nonsense/frameshift/splice) in LoF-disease gene |
| Strong | PS1–PS4 | Same AA change as established P; de novo confirmed; damaging functional assay |
| Moderate | PM1–PM6 | Mutational hotspot; absent gnomAD; length change; assumed de novo |
| Supporting | PP1–PP5 | Co-segregation; computational PP3; phenotype fits gene |
| Stand-alone benign | BA1 | MAF > 5% any population |
| Strong benign | BS1–BS4 | Freq > expected; healthy adult; benign functional study |
| Supporting benign | BP1–BP7 | In silico benign BP4; synonymous no splice effect BP7 |

### Classification Rules (Richards 2015 Table 5)
- **Pathogenic**: PVS1 + (≥1 PS OR ≥2 PM OR PM+PP OR ≥2 PP); or ≥2 PS
- **Likely Pathogenic**: PVS1+PM; or PS+1-2 PM; or PS+≥2 PP; or ≥3 PM; or 2 PM+≥2 PP
- **Benign**: BA1 alone; or ≥2 BS; or BS+BP → Likely Benign

### In Silico Predictors
| Tool | Range | Damaging threshold |
|------|-------|--------------------|
| SIFT | 0–1 (lower=damaging) | < 0.05 |
| PolyPhen-2 | 0–1 (higher=damaging) | > 0.908 probably; > 0.446 possibly |
| CADD (Phred) | higher=worse | ≥ 20 (top 1%); ≥ 25 (top 0.3%) |
| REVEL | 0–1 | > 0.75 likely pathogenic |
| AlphaMissense | 0–1 | > 0.564 likely P; < 0.34 likely B |
| SpliceAI | 0–1 | > 0.2 suggestive; > 0.5 high confidence |

PP3: ≥4/5 tools damaging. BP4: ≤1/5 damaging.

### gnomAD Constraint
- **pLI > 0.9**: haploinsufficient (LoF intolerant)
- **LOEUF < 0.35**: strong LoF constraint
- **Missense o/e < 0.5**: missense-constrained

ClinVar stars: `****` Expert panel | `**` Multi-submitter no conflict | `*` Single submitter

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
    ps  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'strong' and c.startswith('PS'))
    pm  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'moderate')
    pp  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'supporting' and c.startswith('P'))
    ba  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'stand_alone')
    bs  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'strong' and c.startswith('BS'))
    bp  = sum(1 for c in criteria if CRITERIA_STRENGTH.get(c) == 'supporting' and c.startswith('BP'))

    if ba >= 1 or bs >= 2: return 'Benign'
    if (bs >= 1 and bp >= 1) or bp >= 2: return 'Likely Benign'
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2): return 'Pathogenic'
    if ps >= 2: return 'Pathogenic'
    if pvs >= 1 and pm >= 1: return 'Likely Pathogenic'
    if ps >= 1 and pm >= 1: return 'Likely Pathogenic'
    if pm >= 3 or (pm >= 2 and pp >= 2): return 'Likely Pathogenic'
    if pvs >= 1: return 'Likely Pathogenic'
    return 'VUS'
```

## Molecular Docking (AutoDock Vina)

```bash
prepare_receptor4.py -r protein.pdb -o receptor.pdbqt
prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
     --config config.txt --out output.pdbqt
# config.txt: center_x/y/z, size_x/y/z of search box
# Scores in kcal/mol — more negative = stronger predicted binding
```

Force field terms: `E_total = E_bond + E_angle + E_dihedral + E_electrostatic + E_vdW`
Common FFs: AMBER (proteins/NA), CHARMM (proteins/lipids), OPLS-AA (organic molecules)

## GROMACS MD Workflow

```bash
gmx pdb2gmx -f protein.pdb -o protein.gro -water tip3p -ff amber99sb-ildn
gmx editconf -f protein.gro -o box.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
# Then: add ions → energy minimization → NVT equil → NPT equil → production MD
```
Analysis: `gmx rms` (RMSD), `gmx rmsf` (RMSF), `gmx gyrate` (Rg), `gmx hbond`

Homology model quality: >50% identity = reliable; 30–50% = reasonable; <30% = twilight zone.
pLDDT (AlphaFold): >90 high confidence; 70–90 moderate; <50 likely disordered.

## Scanpy Single-Cell Pipeline

```python
import scanpy as sc

adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')  # or sc.datasets.pbmc3k()

# QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20]

# Normalize → HVG → raw snapshot → scale
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata                              # save before subsetting
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# Dimensionality reduction + clustering
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)            # 0.2 coarse – 2.0 fine

# Marker genes + annotation
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
cluster_annotations = {'0': 'CD4+ T cells', '1': 'CD14+ Monocytes'}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
```

## Pitfalls

- **ACMG**: BA1 alone is stand-alone benign — no pathogenic criteria override it. VUS is the default when evidence is insufficient.
- **Variant predictors**: require ≥4/5 tools concordant before applying PP3/BP4; no single tool is authoritative.
- **gnomAD popmax**: use population-specific MAF, not global AF; recessive disorders tolerate higher carrier frequency.
- **Force fields**: never mix parameters from different FFs. Always energy-minimize before MD.
- **Docking scores**: reliable for ranking poses of one ligand; unreliable for cross-ligand ranking without rescoring.
- **Scanpy `adata.raw`**: must be set *before* HVG subsetting — `rank_genes_groups` uses the raw unscaled counts.
- **Leiden resolution**: default 1.0 over-splits PBMC data; 0.3–0.6 typically gives interpretable major cell types.

## Related Skills
- `bioinformatics-workflows-cicd` — Snakemake, Nextflow, GitHub Actions CI
- `bio-applied-clinical-genomics` — VCF annotation pipelines, ClinVar lookup
- `bio-applied-single-cell-scanpy` — advanced Scanpy/AnnData patterns
