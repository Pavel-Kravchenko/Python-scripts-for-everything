---
name: bio-applied-cell-type-annotation
description: "scRNA-seq cell type annotation — manual marker scoring, SingleR reference-based, and CellTypist automated classification"
tool_type: python
primary_tool: NumPy
---

# scRNA-seq: Cell Type Annotation

- [SingleR documentation](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
- [CellTypist](https://www.celltypist.org/)
- [CellMarker database](http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
- [PanglaoDB](https://panglaodb.se/)

## PBMC Canonical Markers

| Cell Type | Marker Genes | Notes |
|-----------|-------------|-------|
| T cells | CD3D, CD3E, CD3G | Pan-T marker |
| CD4+ T cells | CD4, IL7R | Helper T |
| CD8+ T cells | CD8A, CD8B | Cytotoxic T |
| B cells | CD79A, CD79B, MS4A1 (CD20) | MS4A1 = rituximab target |
| NK cells | GNLY, NKG7, KLRD1 | No CD3 distinguishes from T |
| Monocytes (classical) | LYZ, S100A8, S100A9, CST3 | High LYZ = phagocytic |
| Monocytes (non-classical) | FCGR3A (CD16), MS4A7 | Patrolling |
| Dendritic cells | FCER1A, CLEC10A | Low count, high HLA-DR |
| Platelets | PPBP, PF4 | Often contaminants |

## Manual Annotation Workflow

1. Plot known markers on UMAP
2. Check violin plots per cluster
3. Assign if >= 2 concordant markers strongly expressed
4. Ambiguous clusters: check for doublet markers (two cell type programs)

```python
marker_dict = {
    'T_cell':   ['CD3D', 'CD3E', 'CD3G'],
    'B_cell':   ['CD79A', 'MS4A1'],
    'Monocyte': ['LYZ', 'S100A8'],
    'NK_cell':  ['GNLY', 'NKG7'],
    'Dendritic':['FCER1A', 'CLEC10A'],
}

X = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
gene_names = adata.var_names.tolist()
clusters = adata.obs['leiden'].values
cluster_ids = sorted(adata.obs['leiden'].unique())

# Score each cluster for each cell type
scores = {}
for ct, markers in marker_dict.items():
    valid_markers = [m for m in markers if m in gene_names]
    if valid_markers:
        marker_idx = [gene_names.index(m) for m in valid_markers]
        scores[ct] = {cid: X[clusters == cid][:, marker_idx].mean() for cid in cluster_ids}

score_df = pd.DataFrame(scores, index=cluster_ids)
cluster_annotation = score_df.idxmax(axis=1).to_dict()
adata.obs['cell_type_annotated'] = adata.obs['leiden'].map(cluster_annotation)
```

## SingleR (R, Reference-Based)

```r
library(SingleR)
library(celldex)
ref <- MonacoImmuneData()  # 29 immune subtypes, best for PBMC
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
plotScoreHeatmap(pred)
table(pred$pruned.labels)  # NA = low confidence — inspect these cells
```

Available references: `HumanPrimaryCellAtlasData()` (37 types), `MonacoImmuneData()` (29 immune), `BlueprintEncodeData()` (hematopoietic).

`pred$pruned.labels` NA cells may be novel types, transitional states, or poor-quality cells.

## CellTypist (Python)

```python
import celltypist
from celltypist import models

models.download_models(force_update=False)
# Models: Immune_All_Low.pkl (36 subtypes), Immune_All_High.pkl (9 broad), Pan_Fetal_Human.pkl
model = models.Model.load(model='Immune_All_Low.pkl')

# IMPORTANT: expects log1p-normalized data with target_sum=1e4
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata = predictions.to_adata()
```

`majority_voting=True`: pools per-cluster predictions, reduces noise from individual cell misclassification.

**When CellTypist fails**: non-immune cells need tissue-specific models; non-human species need retraining; if all clusters predict same type, check normalization (`target_sum=1e4`).

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
