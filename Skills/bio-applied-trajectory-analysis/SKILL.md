---
name: bio-applied-trajectory-analysis
description: "scRNA-seq trajectory analysis: pseudotime (DPT), PAGA graph abstraction, and RNA velocity (scVelo). Decision guide, key parameters, and pitfalls."
tool_type: python
primary_tool: NumPy
---

# scRNA-seq: Trajectory Analysis and RNA Velocity

## Method Overview

| Method | Tool | Input | Output |
|--------|------|-------|--------|
| Diffusion pseudotime (DPT) | scanpy `sc.tl.dpt` | k-NN graph + root cell | Pseudotime per cell |
| PAGA | scanpy `sc.tl.paga` | Cluster labels + k-NN graph | Cluster connectivity graph |
| RNA velocity (steady-state) | scVelo `mode='stochastic'` | Spliced + unspliced counts | Velocity arrows |
| RNA velocity (dynamical) | scVelo `mode='dynamical'` | Spliced + unspliced counts | Velocity + kinetics |

## When Pseudotime is Valid / Invalid

| Valid | Invalid |
|-------|---------|
| Continuous developmental transitions sampled | Only start + end states, no intermediates |
| Trajectory visible as smooth continuum in UMAP | Discrete blobs in UMAP |
| Differentiation is the dominant variation source | Batch effects / cell cycle dominate PCs |
| ≥50 cells in intermediate stages | Too few intermediate cells |

## DPT Implementation Pattern

```python
import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import dijkstra

def compute_dpt(X_pca, root_idx, n_neighbors=20):
    """Diffusion pseudotime via Gaussian-kernel affinity + Dijkstra."""
    knn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    knn.fit(X_pca)
    distances, indices = knn.kneighbors(X_pca)

    sigma = np.median(distances[:, 1:])
    affinities = np.exp(-distances**2 / (2 * sigma**2))

    n = X_pca.shape[0]
    A = sp.lil_matrix((n, n))
    for i in range(n):
        for j_local, j in enumerate(indices[i]):
            A[i, j] = affinities[i, j_local]
    A = sp.csr_matrix(A)
    A = (A + A.T) / 2  # symmetrize

    # Convert affinities → distances for shortest-path
    with np.errstate(divide='ignore'):
        dist_graph = sp.csr_matrix(1.0 / A.toarray())
    dist_graph = sp.csr_matrix(np.where(A.toarray() > 0, dist_graph.toarray(), 0))

    graph_dist = dijkstra(dist_graph, directed=False, indices=root_idx)
    graph_dist[np.isinf(graph_dist)] = graph_dist[~np.isinf(graph_dist)].max()
    dpt = (graph_dist - graph_dist.min()) / (graph_dist.max() - graph_dist.min())
    return dpt

# Root selection
# Option 1: use a known marker
root_idx = int(np.argmax(adata.X[:, adata.var_names.get_loc('SOX2')]))
# Option 2: most extreme on diffusion component 1
root_idx = int(np.argmin(adata.obsm['X_diffmap'][:, 0]))
```

## scanpy DPT Workflow

```python
import scanpy as sc

sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.diffmap(adata)
adata.uns['iroot'] = root_idx
sc.tl.dpt(adata)
# result: adata.obs['dpt_pseudotime']
```

## PAGA

```python
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, threshold=0.03, layout='fr')

# PAGA-initialized UMAP (more reproducible, preserves global topology)
sc.tl.umap(adata, init_pos='paga')
```

**Use PAGA when:** >20 clusters, trajectory topology is unclear, or individual cell UMAP is too cluttered to interpret.

## RNA Velocity (scVelo)

```python
import scvelo as scv

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Dynamical model: slower but captures transient states
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
```

**Data requirements:** spliced + unspliced count matrices; ≥50,000 reads/cell. Obtain via STARsolo (`--soloFeatures Gene Velocyto`), Velocyto CLI, or alevin-fry.

## Pitfalls

- **Root cell choice dominates results**: wrong root flips pseudotime direction. Validate using known markers or experimental time points.
- **Pseudotime ≠ real time**: ordering reflects transcriptional similarity, not clock time. Do not compare pseudotime values across datasets.
- **RNA velocity on fully differentiated cells is noise**: velocity requires cells in transition. Stable mature types show circular/short arrows — this is correct, not a bug.
- **UMAP distances are not biological distances**: cells far apart in UMAP may not be transcriptionally distant. Use PC space for quantitative comparisons.
- **Batch effects corrupt trajectories**: correct batch effects before trajectory inference; a batch effect can look like a developmental branch.
- **Coordinate systems**: BED = 0-based half-open; VCF/GFF = 1-based inclusive
