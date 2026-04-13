---
name: bio-applied-trajectory-analysis
description: "**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 4**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/04_trajectory_analysis.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: anndata 0.10+, matplotlib 3.8+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+, scipy 1.12+, scvelo 0.3+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# scRNA-seq: Trajectory Analysis and RNA Velocity

*Source: Course notebook `Tier_3_Applied_Bioinformatics/30_Single_Cell_RNA_seq/04_trajectory_analysis.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 30 · Notebook 4**

*Prerequisites: Notebook 3 (Cell Type Annotation)*

---

**By the end of this notebook you will be able to:**
1. Explain pseudotime as a proxy for developmental progression
2. Build diffusion maps and compute pseudotime with DPT (diffusion pseudotime)
3. Run Monocle 3 trajectory inference for branched differentiation trees
4. Compute spliced/unspliced ratios and run scVelo for RNA velocity
5. Interpret velocity arrows in UMAP space and identify driver genes



**Key resources:**
- [scVelo documentation](https://scvelo.readthedocs.io/)
- [Monocle 3 documentation](https://cole-trapnell-lab.github.io/monocle3/)
- [PAGA (Wolf et al. 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1663-x)

## Why trajectory analysis matters

Clustering gives you a static snapshot: these are the cell types present. Trajectory analysis asks: how did cells transition between states? In development, stem cells differentiate into mature cell types through continuous intermediate states. Pseudotime recapitulates these transitions by ordering cells along a one-dimensional (or branched) continuum that correlates with developmental progression. This enables identification of driver genes, regulatory switches at branch points, and the kinetics of gene expression changes.

```python
# Setup: create a differentiation-like dataset (stem -> two branches)
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(42)
n_cells = 400
n_genes = 500

# Simulate three differentiation stages: stem -> progenitor -> {typeA, typeB}
n_stem, n_prog, n_typeA, n_typeB = 80, 120, 100, 100
stages = (['Stem'] * n_stem + ['Progenitor'] * n_prog + 
          ['Type_A'] * n_typeA + ['Type_B'] * n_typeB)
pseudotime_truth = np.concatenate([
    rng.uniform(0, 0.3, n_stem),
    rng.uniform(0.3, 0.6, n_prog),
    rng.uniform(0.6, 1.0, n_typeA),
    rng.uniform(0.6, 1.0, n_typeB)
])

counts = rng.negative_binomial(2, 0.9, (n_cells, n_genes)).astype(np.float32)

# Gene programs along pseudotime
for i, pt in enumerate(pseudotime_truth):
    # Early genes (expressed in stem cells, decrease)
    counts[i, :50] += (rng.poisson(20 * (1 - pt), 50)).astype(np.float32)
    # Late genes (expressed in mature cells, increase)
    counts[i, 50:100] += (rng.poisson(15 * pt, 50)).astype(np.float32)
    # Branch A specific
    if stages[i] == 'Type_A':
        counts[i, 100:130] += rng.poisson(25, 30).astype(np.float32)
    # Branch B specific
    if stages[i] == 'Type_B':
        counts[i, 130:160] += rng.poisson(25, 30).astype(np.float32)

# Normalize
cell_totals = counts.sum(axis=1, keepdims=True)
X_norm = np.log1p(counts / cell_totals * 1e4)

adata = ad.AnnData(
    X=X_norm,
    obs=pd.DataFrame({'stage': stages, 'pseudotime_truth': pseudotime_truth},
                     index=[f'CELL_{i:04d}' for i in range(n_cells)]),
    var=pd.DataFrame(index=[f'GENE{i:04d}' for i in range(n_genes)])
)

# PCA + UMAP
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
X_scaled = StandardScaler().fit_transform(X_norm)
X_scaled = np.clip(X_scaled, -10, 10)
X_pca = PCA(n_components=20, random_state=42).fit_transform(X_scaled)
adata.obsm['X_pca'] = X_pca

try:
    from umap import UMAP
    X_umap = UMAP(n_neighbors=20, min_dist=0.2, random_state=42).fit_transform(X_pca[:, :15])
    adata.obsm['X_umap'] = X_umap
except ImportError:
    adata.obsm['X_umap'] = X_pca[:, :2]

print(f"Trajectory dataset: {adata.shape}")
print(f"Stages: {pd.Series(stages).value_counts().to_dict()}")

# Visualize true structure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors_map = {'Stem':'#1b7837','Progenitor':'#fee08b','Type_A':'#d73027','Type_B':'#4575b4'}
coords = adata.obsm['X_umap']

for st in ['Stem','Progenitor','Type_A','Type_B']:
    mask = np.array(stages) == st
    axes[0].scatter(coords[mask,0], coords[mask,1], c=colors_map[st], label=st, s=12, alpha=0.8)
axes[0].set_title('True differentiation stages')
axes[0].legend()

sc = axes[1].scatter(coords[:,0], coords[:,1], c=pseudotime_truth, cmap='viridis', s=10)
axes[1].set_title('True pseudotime')
plt.colorbar(sc, ax=axes[1], label='Pseudotime')
for ax in axes:
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.savefig('trajectory_data.png', dpi=100, bbox_inches='tight')
plt.show()
```python

## 1. Pseudotime Concepts

Pseudotime is not real time — it is a computational ordering of cells along a trajectory that maximally aligns with the dominant direction of transcriptional change.

### Key assumptions
1. **Continuous transitions**: cells at intermediate states must be sampled. If you only have stem cells and mature cells with no intermediates, pseudotime ordering is unreliable.
2. **The dominant source of variation is developmental**: if batch effects or cell cycle dominate variation, pseudotime will capture those instead of differentiation.
3. **Root selection**: you must specify where the trajectory starts (the "root" cell). This is the only prior knowledge required and significantly affects interpretation. Root cells are typically the least differentiated (lowest expression of maturation markers, highest stem cell marker expression).

### When pseudotime is valid
- Datasets explicitly designed to capture a developmental process (embryonic time courses, ex vivo differentiation assays)
- Clear trajectory structure visible in UMAP (cells form a smooth continuum, not discrete blobs)

### When pseudotime is unreliable
- Mature, terminally differentiated cell types (macrophages, neurons) — no transitions exist
- Too few cells sampled along the trajectory (< 50 cells in intermediate stages)
- Strong batch effects that weren't corrected before trajectory analysis

## 2. Diffusion Maps and DPT (Diffusion Pseudotime)

### Diffusion maps
Diffusion maps (Coifman & Lafon 2006) model the data as a diffusion process on the cell-cell similarity graph. Each diffusion component is an eigenvector of the diffusion operator, capturing progressively finer structure in the manifold. Unlike PCA, diffusion maps respect the manifold geometry.

In scanpy: `sc.tl.diffmap(adata)` stores components in `adata.obsm['X_diffmap']`.

### DPT algorithm (Haghverdi et al. 2016)
DPT assigns pseudotime by measuring the diffusion distance from a root cell to every other cell:
1. Compute the diffusion distance matrix using diffusion map components
2. Select a root cell (usually the most extreme cell on diffusion component 1)
3. Pseudotime of each cell = diffusion distance from root

The root is specified via `adata.uns['iroot']` — the integer index of the root cell, not its name.

### Selecting the root cell
```python
# Option 1: Use a known marker to find the most stem-like cell
stem_marker_idx = adata.var_names.tolist().index('SOX2')
stem_scores = adata.X[:, stem_marker_idx]
adata.uns['iroot'] = int(np.argmax(stem_scores))

# Option 2: Use the cell at the extreme end of diffusion component 1
adata.uns['iroot'] = int(np.argmin(adata.obsm['X_diffmap'][:, 0]))
```python

```python
# Implement diffusion pseudotime
import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

X_pca = adata.obsm['X_pca'][:, :15]
coords = adata.obsm['X_umap']
stages = adata.obs['stage'].values
pseudotime_truth = adata.obs['pseudotime_truth'].values

# Step 1: Build k-NN graph
knn = NearestNeighbors(n_neighbors=20, metric='euclidean')
knn.fit(X_pca)
distances, indices = knn.kneighbors(X_pca)

# Step 2: Gaussian kernel to convert distances to affinities
sigma = np.median(distances[:, 1:])  # bandwidth = median nearest-neighbor distance
affinities = np.exp(-distances**2 / (2 * sigma**2))

# Build affinity matrix
n_cells = adata.n_obs
A = sp.lil_matrix((n_cells, n_cells))
for i in range(n_cells):
    for j_local, j in enumerate(indices[i]):
        A[i, j] = affinities[i, j_local]
A = sp.csr_matrix(A)
A = (A + A.T) / 2  # symmetrize

# Step 3: Normalize to Markov transition matrix
D = sp.diags(np.array(A.sum(axis=1)).flatten())
T = sp.linalg.spsolve(D, A)  # row-normalized = transition probabilities

# Step 4: Diffusion pseudotime from root cell
# Root = stem cell with highest early gene expression
early_genes_idx = list(range(50))
early_expr = adata.X[:, early_genes_idx].sum(axis=1)
root_idx = int(np.argmax(early_expr))
adata.uns['iroot'] = root_idx
print(f"Root cell: {adata.obs_names[root_idx]} (stage: {stages[root_idx]})")

# Compute pseudotime via shortest path distances in affinity graph
# (approximation of diffusion distance)
from scipy.sparse.csgraph import shortest_path, dijkstra

# Convert affinities to distances for Dijkstra
dist_matrix = sp.csr_matrix(1.0 / (A.toarray() + 1e-10))
np.fill_diagonal(dist_matrix.toarray(), 0)
graph_distances = dijkstra(sp.csr_matrix(dist_matrix), 
                            directed=False, indices=root_idx)
# Handle unreachable cells
graph_distances[np.isinf(graph_distances)] = graph_distances[~np.isinf(graph_distances)].max()

# Normalize to [0, 1]
dpt = (graph_distances - graph_distances.min()) / (graph_distances.max() - graph_distances.min())
adata.obs['dpt_pseudotime'] = dpt

# Evaluate: correlation with true pseudotime
from scipy.stats import spearmanr
corr, p = spearmanr(dpt, pseudotime_truth)
print(f"Spearman correlation with true pseudotime: {corr:.3f} (p={p:.2e})")

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(16, 4))

sc1 = axes[0].scatter(coords[:,0], coords[:,1], c=pseudotime_truth, cmap='viridis', s=10)
axes[0].scatter(coords[root_idx,0], coords[root_idx,1], c='red', s=100, marker='*', 
                zorder=5, label='Root cell')
axes[0].set_title('True pseudotime')
plt.colorbar(sc1, ax=axes[0])
axes[0].legend()

sc2 = axes[1].scatter(coords[:,0], coords[:,1], c=dpt, cmap='plasma', s=10)
axes[1].set_title(f'DPT pseudotime (r={corr:.2f})')
plt.colorbar(sc2, ax=axes[1])

axes[2].scatter(pseudotime_truth, dpt, s=5, alpha=0.4)
axes[2].set_xlabel('True pseudotime')
axes[2].set_ylabel('DPT pseudotime')
axes[2].set_title(f'Correlation: Spearman r={corr:.3f}')

for ax in axes[:2]:
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.savefig('dpt_pseudotime.png', dpi=100, bbox_inches='tight')
plt.show()
```python

## 3. PAGA Graph Abstraction

PAGA (Wolf et al. 2019, *Genome Biology*) provides a coarse-grained view of the data topology: instead of individual cells, it shows connectivity between clusters. This is especially useful for complex datasets with many cell types where visualizing individual cell-level trajectories is cluttered.

### How PAGA works
1. For each pair of clusters, count the number of inter-cluster edges in the k-NN graph
2. Compute expected inter-cluster edges under a null model (random cluster assignment with same sizes)
3. PAGA score = (observed - expected) / expected
4. High PAGA score → clusters are more connected than expected → likely a real transition

### Using PAGA to initialize UMAP
PAGA positions clusters in a meaningful layout, which can then be used as initial coordinates for UMAP (PAGA-initialized UMAP). This makes UMAP more reproducible and preserves global topology better than random initialization.

```python
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, threshold=0.03, layout='fr')  # Fruchterman-Reingold layout

# PAGA-initialized UMAP
sc.tl.umap(adata, init_pos='paga')
```python

### When to use PAGA
- Complex datasets: 20+ clusters where point cloud UMAP is uninterpretable
- Trajectory topology questions: does differentiation follow a tree? A DAG? Cycles?
- Multi-sample datasets: PAGA edges reflect shared transitions across all cells/samples

## 4. RNA Velocity with scVelo

### The RNA velocity concept (La Manno et al. 2018)
Every gene has two observable states: **unspliced pre-mRNA** (still in the nucleus, recently transcribed) and **spliced mRNA** (mature, cytoplasmic). 

- If a gene is being actively upregulated: unspliced > expected ratio → velocity points toward higher expression
- If a gene is being downregulated: unspliced < expected ratio → velocity points toward lower expression

By modeling this for thousands of genes simultaneously, you can infer the direction of transcriptional change for each cell — pointing it toward its future state.

### The two models in scVelo
1. **Steady-state / ratio model** (`mode='stochastic'` in scVelo): assumes transcription and splicing are in steady state. Fast, but less accurate for transient dynamics.
2. **Dynamical model** (`mode='dynamical'`): explicitly models the kinetics (transcription, splicing, and degradation rates). Slower (EM algorithm) but captures transient states accurately.

### Data requirements
- Requires BOTH spliced and unspliced count matrices
- Obtained from: STARsolo (with `--soloFeatures Gene Velocyto`), Velocyto CLI, or alevin-fry
- Rule of thumb: need ≥ 50,000 reads/cell for reliable velocity estimates

### Key scVelo workflow
```python
import scvelo as scv
scv.pp.filter_and_normalize(adata)      # filter genes, normalize
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)  # compute first/second moments
scv.tl.recover_dynamics(adata)          # fit kinetic model (dynamical mode)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)            # transition probabilities from velocities
scv.pl.velocity_embedding_stream(adata, basis='umap')  # streamline plot
```python

### Interpreting velocity arrows
- Arrow direction: predicted future state
- Arrow length: speed of transcriptional change
- Circular arrows: cells in a stable attractor state (no net change)
- **Caveat**: velocity requires sufficient cells in transition states. If your dataset only has mature cell types, velocity will be noisy and uninterpretable.

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
