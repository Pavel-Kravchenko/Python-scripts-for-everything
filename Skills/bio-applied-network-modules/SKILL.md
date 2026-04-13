---
name: bio-applied-network-modules
description: "**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 2**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/28_Network_Biology/02_network_modules.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, networkx 3.2+, numpy 1.26+, pandas 2.1+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Network Modules, Enrichment & Visualization

*Source: Course notebook `Tier_3_Applied_Bioinformatics/28_Network_Biology/02_network_modules.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 2**

*Prerequisites: Notebook 1 (PPI Network Construction)*

---

**By the end of this notebook you will be able to:**
1. Apply Louvain community detection to find functional network modules
2. Run GO enrichment on each detected module
3. Export networks to Cytoscape format for interactive exploration
4. Build protein co-expression networks from RNA-seq data (WGCNA concept)
5. Integrate network modules with differential expression results



**Key resources:**
- [Cytoscape documentation](https://cytoscape.org/documentation_users.html)
- [py4cytoscape documentation](https://py4cytoscape.readthedocs.io/)
- [NetworkX community detection](https://networkx.org/documentation/stable/reference/algorithms/community.html)

## 1. Community Detection in Biological Networks

Biological networks are organized into **modules** (communities) — groups of densely interconnected nodes that are sparsely connected to the rest of the network. These modules often correspond to biological functions or complexes.

### Why modules matter
- **Functional coherence**: nodes in the same module tend to share biological function
- **Robustness**: modular organization prevents cascading failures
- **Drug targeting**: disrupt specific modules to affect specific pathways

### Community detection algorithms

| Algorithm | Key idea | Complexity | Notes |
|---|---|---|---|
| **Louvain** | Greedy modularity maximization + refinement | O(n log n) | Fast, widely used |
| **Leiden** | Fixed partition bug in Louvain | O(n log n) | Preferred over Louvain |
| **Girvan-Newman** | Edge betweenness removal | O(m²n) | Slow, interpretable |
| **Spectral clustering** | Eigenvectors of Laplacian | O(n³) | Principled, no resolution limit |
| **Infomap** | Random walk compression | O(m) | Good for directed networks |

### Modularity score Q
Q measures quality of a partition: fraction of edges within modules minus expected fraction under random placement.
- Q > 0.3: meaningful community structure
- Q > 0.5: strong modularity
- Q = 0: no better than random
- Q = 1: perfect (no between-module edges)

```python
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from collections import defaultdict

np.random.seed(42)

# ----- Generate a modular PPI network -----
# 5 biological modules: Cell cycle, DNA repair, Signaling, Metabolism, Immune
module_names = ['Cell_cycle', 'DNA_repair', 'MAPK_signaling', 'Metabolism', 'Immune']
module_sizes = [12, 10, 11, 9, 10]

module_gene_sets = {
    'Cell_cycle':    ['CDK1','CDK2','CDK4','CDK6','CCND1','CCNE1','RB1','E2F1','E2F3','CDKN1A','CDKN2A','CDC25A'],
    'DNA_repair':    ['BRCA1','BRCA2','ATM','ATR','CHEK1','CHEK2','RAD51','FANCD2','PALB2','MLH1'],
    'MAPK_signaling':['KRAS','BRAF','MAP2K1','MAPK1','MAPK3','EGFR','ERBB2','GRB2','SOS1','RAF1','HRAS','PIK3CA'],
    'Metabolism':    ['LDHA','PKM','HK2','G6PD','FASN','ACLY','HMGCR','SLC2A1'],
    'Immune':        ['CD8A','CD274','PDCD1','CTLA4','LAG3','TIM3','TIGIT','FOXP3','IL2','IFNG'],
}
all_genes = [g for gs in module_gene_sets.values() for g in gs]
n_total = len(all_genes)

# Build adjacency matrix with modular structure
# High intra-module connectivity, low inter-module connectivity
G = nx.Graph()
G.add_nodes_from(all_genes)

# Add node module labels
for mod, genes in module_gene_sets.items():
    for g in genes:
        G.nodes[g]['module_true'] = mod

# Intra-module edges (dense)
for mod, genes in module_gene_sets.items():
    n = len(genes)
    for i in range(n):
        for j in range(i+1, n):
            if np.random.random() < 0.55:  # 55% intra-module edge probability
                score = np.random.uniform(500, 999)
                G.add_edge(genes[i], genes[j], weight=score)

# Inter-module edges (sparse, known cross-talk)
cross_talk = [
    ('Cell_cycle', 'DNA_repair', 0.20),  # cell cycle-DNA damage checkpoint
    ('MAPK_signaling', 'Cell_cycle', 0.15),
    ('MAPK_signaling', 'Metabolism', 0.10),
    ('DNA_repair', 'Immune', 0.08),
    ('Immune', 'MAPK_signaling', 0.06),
]
for mod1, mod2, prob in cross_talk:
    for g1 in module_gene_sets[mod1]:
        for g2 in module_gene_sets[mod2]:
            if np.random.random() < prob:
                G.add_edge(g1, g2, weight=np.random.uniform(150, 500))

print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print(f"Modules: {module_names}")
print(f"Density: {nx.density(G):.3f}")
print(f"Connected: {nx.is_connected(G)}")
```python

## 2. Louvain Community Detection

The **Louvain algorithm** proceeds in two phases:
1. **Phase 1** (modularity gain): Move each node to the neighbor community that maximizes ΔQ
2. **Phase 2** (community contraction): Aggregate each community into a supernode; repeat

This is repeated until ΔQ = 0 (local maximum of modularity).

### Leiden algorithm (improved Louvain)
Louvain can produce **disconnected communities** (a bug). The Leiden algorithm guarantees:
- Well-connected communities
- Faster convergence
- Better modularity scores

```python
# Louvain (python-louvain package)
import community as community_louvain
partition = community_louvain.best_partition(G, random_state=42)

# Leiden (leidenalg package)
import leidenalg
import igraph as ig
ig_G = ig.Graph.from_networkx(G)
partition = leidenalg.find_partition(ig_G, leidenalg.ModularityVertexPartition, seed=42)
```python

```python
# ----- Community detection: Louvain algorithm -----
# Note: python-louvain package (community module) implements Louvain
# Here we implement a simplified version based on greedy modularity

try:
    import community as community_louvain
    partition_louvain = community_louvain.best_partition(G, random_state=42)
    method = 'Louvain (python-louvain)'
except ImportError:
    # Fallback: greedy modularity maximization from networkx
    communities = nx.community.greedy_modularity_communities(G)
    partition_louvain = {}
    for i, comm in enumerate(communities):
        for node in comm:
            partition_louvain[node] = i
    method = 'Greedy modularity (networkx)'

n_communities = len(set(partition_louvain.values()))
modularity = nx.community.modularity(G, 
    [{n for n, c in partition_louvain.items() if c == i} for i in range(n_communities)])

print(f"Community detection method: {method}")
print(f"Number of communities found: {n_communities}")
print(f"Modularity score Q = {modularity:.4f} (>0.3 is meaningful)")

# Compare detected communities with true modules
print("\nCommunity composition (detected vs true):")
df_comp = pd.DataFrame({'Gene': list(G.nodes()),
                         'Detected': [partition_louvain[n] for n in G.nodes()],
                         'True': [G.nodes[n]['module_true'] for n in G.nodes()]})
print(pd.crosstab(df_comp['True'], df_comp['Detected']))
```python

## 3. Visualizing Network Modules

Good module visualization requires:
1. **Layout that separates modules** (e.g., spring layout with high k)
2. **Color nodes by module membership**
3. **Show inter-module edges** as thin, different color
4. **Export to Cytoscape** for publication-quality figures

### Comparing detected vs. true modules
The **Normalized Mutual Information (NMI)** and **Adjusted Rand Index (ARI)** measure similarity between detected and true partitions:
- NMI = 1: perfect match
- NMI = 0: random match

```python
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
nmi = normalized_mutual_info_score(true_labels, detected_labels)
ari = adjusted_rand_score(true_labels, detected_labels)
```python

```python
# ----- Visualize communities -----

fig, axes = plt.subplots(1, 3, figsize=(18, 7))
fig.suptitle('Network Community Detection', fontsize=13, fontweight='bold')

pos = nx.spring_layout(G, seed=42, k=3.0)

# Panel 1: True modules (ground truth)
mod_color_map = cm.Set1(np.linspace(0, 1, len(module_names)))
mod_to_color = {mod: mod_color_map[i] for i, mod in enumerate(module_names)}
true_colors = [mod_to_color[G.nodes[n]['module_true']] for n in G.nodes()]
nx.draw_networkx(G, pos=pos, ax=axes[0], node_color=true_colors, node_size=150,
                 edge_color='gray', alpha=0.7, width=0.5, font_size=5, arrows=False)
axes[0].set_title('True Module Structure')
handles = [mpatches.Patch(color=mod_to_color[m], label=m.replace('_', ' '))
           for m in module_names]
axes[0].legend(handles=handles, fontsize=7, loc='lower left')
axes[0].axis('off')

# Panel 2: Detected communities
n_detected = len(set(partition_louvain.values()))
det_colors_map = cm.tab10(np.linspace(0, 1, n_detected))
det_colors = [det_colors_map[partition_louvain[n]] for n in G.nodes()]
nx.draw_networkx(G, pos=pos, ax=axes[1], node_color=det_colors, node_size=150,
                 edge_color='gray', alpha=0.7, width=0.5, font_size=5, arrows=False)
axes[1].set_title(f'Detected Communities (Q={modularity:.3f})')
axes[1].axis('off')

# Panel 3: Module size distribution
comm_sizes = pd.Series(partition_louvain).value_counts().sort_index()
axes[2].bar(comm_sizes.index, comm_sizes.values, color=det_colors_map[:len(comm_sizes)])
axes[2].set_xlabel('Community ID')
axes[2].set_ylabel('Number of nodes')
axes[2].set_title('Community Size Distribution')

plt.tight_layout()
plt.show()

print(f"\nCommunity sizes: {dict(comm_sizes)}")
print(f"Modularity Q={modularity:.4f}")
print("Q > 0.3: meaningful community structure")
print("Q > 0.5: strong modularity")
```python

## 4. WGCNA: Weighted Co-expression Network Analysis

**WGCNA** constructs co-expression modules from RNA-seq data (not from PPI databases):

### WGCNA workflow
1. **Pearson correlation matrix** of gene expression
2. **Soft thresholding**: raise adjacency to power β to emphasize strong correlations (scale-free test)
3. **Topological Overlap Measure (TOM)**: connectivity-based adjacency (more robust than simple correlation)
4. **Hierarchical clustering** of 1-TOM distance matrix
5. **Dynamic tree cut**: cut dendrogram to obtain modules
6. **Module eigengene (ME)**: first PC of module expression = summary of module activity

### Key WGCNA outputs
- **Module-trait correlation**: which modules are associated with clinical outcomes?
- **Module membership (kME)**: correlation of each gene with its module eigengene
- **Hub gene**: highest intra-module connectivity (kIM) within a module

### WGCNA vs. PPI network modules
| Aspect | WGCNA | STRING communities |
|---|---|---|
| Data | Expression matrix | Protein interactions |
| Edge | Co-expression r | PPI confidence |
| Interpretation | Co-regulated genes | Physical/functional partners |
| Context-specific | Yes | No |

```python
# ----- WGCNA concepts: co-expression network modules -----

# WGCNA (Weighted Gene Co-expression Network Analysis) works on expression data
# Step 1: Build adjacency matrix from correlation matrix
# Step 2: Raise to soft threshold power beta (ensures scale-free topology)
# Step 3: TOM (Topological Overlap Measure) = connectivity-based adjacency
# Step 4: Hierarchical clustering + dynamic tree cut → modules
# Step 5: Module eigengene = first PC of module expression

np.random.seed(5)

# Simulate expression data for the network genes
n_samples = 50
n_genes = len(all_genes)
gene_idx = {g: i for i, g in enumerate(all_genes)}

# Expression matrix: genes with module structure
expr_matrix = np.random.randn(n_samples, n_genes)
# Add module signal
module_signals = {mod: np.random.randn(n_samples) for mod in module_names}
for mod, genes in module_gene_sets.items():
    for g in genes:
        expr_matrix[:, gene_idx[g]] += module_signals[mod] * 1.5

# Correlation matrix
from sklearn.preprocessing import StandardScaler
expr_scaled = StandardScaler().fit_transform(expr_matrix)
corr_matrix = np.corrcoef(expr_scaled.T)  # gene x gene

# Module eigengenes
module_eigengenes = {}
from sklearn.decomposition import PCA
for mod, genes in module_gene_sets.items():
    mod_indices = [gene_idx[g] for g in genes]
    mod_expr = expr_scaled[:, mod_indices]
    pca = PCA(n_components=1)
    eigengene = pca.fit_transform(mod_expr).flatten()
    module_eigengenes[mod] = eigengene

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('WGCNA Co-expression Network Concepts', fontsize=12, fontweight='bold')

# 1: Gene correlation heatmap (subset)
# Order genes by module
ordered_genes = [g for mod in module_names for g in module_gene_sets[mod]]
ordered_idx = [gene_idx[g] for g in ordered_genes]
corr_ordered = corr_matrix[np.ix_(ordered_idx, ordered_idx)]

im = axes[0].imshow(corr_ordered, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
axes[0].set_title('Gene-Gene Correlation Matrix\n(ordered by module)')
axes[0].set_xticks([]); axes[0].set_yticks([])
plt.colorbar(im, ax=axes[0], label='Pearson r')

# Add module boundaries
cumsum = np.cumsum([0] + module_sizes)
for boundary in cumsum[1:-1]:
    axes[0].axhline(boundary - 0.5, color='black', linewidth=1)
    axes[0].axvline(boundary - 0.5, color='black', linewidth=1)

# 2: Module eigengene expression across samples
sample_labels = ['Type1'] * 25 + ['Type2'] * 25
type_order = sorted(range(n_samples), key=lambda i: sample_labels[i])
for i, (mod, eigengene) in enumerate(module_eigengenes.items()):
    axes[1].plot(range(n_samples), eigengene[type_order] + i*3,
                 label=mod.replace('_', ' '), alpha=0.8)
axes[1].axvline(25, color='black', linestyle='--', label='Type1|Type2')
axes[1].set_xlabel('Samples (sorted by type)')
axes[1].set_ylabel('Module eigengene (offset for clarity)')
axes[1].set_title('Module Eigengene Expression')
axes[1].legend(fontsize=7)

# 3: Module-module eigengene correlation
me_df = pd.DataFrame(module_eigengenes)
me_corr = me_df.corr()
im2 = axes[2].imshow(me_corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
axes[2].set_xticks(range(len(module_names)))
axes[2].set_yticks(range(len(module_names)))
axes[2].set_xticklabels([m.replace('_', ' ') for m in module_names], rotation=45, ha='right', fontsize=7)
axes[2].set_yticklabels([m.replace('_', ' ') for m in module_names], fontsize=7)
axes[2].set_title('Module-Module Eigengene Correlation')
plt.colorbar(im2, ax=axes[2])

for i in range(len(module_names)):
    for j in range(len(module_names)):
        axes[2].text(j, i, f'{me_corr.iloc[i,j]:.2f}', ha='center', va='center', fontsize=7)

plt.tight_layout()
plt.show()
```python

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
