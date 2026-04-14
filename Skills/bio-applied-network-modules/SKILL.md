---
name: bio-applied-network-modules
description: Community detection in biological networks — Louvain/Leiden algorithms, modularity Q, WGCNA co-expression modules, and Cytoscape export
tool_type: python
primary_tool: Matplotlib
---

# Network Modules, Enrichment and Visualization

## Community Detection Algorithms

| Algorithm | Complexity | Notes |
|-----------|------------|-------|
| Louvain | O(n log n) | Greedy modularity; can produce disconnected communities (bug) |
| Leiden | O(n log n) | Fixes Louvain bug; preferred for publication |
| Girvan-Newman | O(m²n) | Edge betweenness removal; slow, interpretable |
| Spectral clustering | O(n³) | Principled, no resolution limit |
| Infomap | O(m) | Best for directed networks |

### Modularity Score Q

- Q > 0.3: meaningful community structure
- Q > 0.5: strong modularity
- Q = 0: no better than random
- Q = 1: perfect (no between-module edges)

## Community Detection Pattern

```python
import networkx as nx

# Louvain (python-louvain package)
import community as community_louvain
partition = community_louvain.best_partition(G, random_state=42)

# Leiden (leidenalg + igraph)
import leidenalg, igraph as ig
ig_G = ig.Graph.from_networkx(G)
leiden_part = leidenalg.find_partition(ig_G, leidenalg.ModularityVertexPartition, seed=42)

# Fallback: NetworkX greedy modularity
communities = nx.community.greedy_modularity_communities(G)
partition = {node: i for i, comm in enumerate(communities) for node in comm}

# Compute Q
n_comm = len(set(partition.values()))
Q = nx.community.modularity(G,
    [{n for n, c in partition.items() if c == i} for i in range(n_comm)])

# Evaluate against ground truth
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
nmi = normalized_mutual_info_score(true_labels, detected_labels)
ari = adjusted_rand_score(true_labels, detected_labels)
```

## Building a Modular PPI Network

```python
import numpy as np
import networkx as nx

# Intra-module edges (dense, ~55% probability)
for mod, genes in module_gene_sets.items():
    for i, g1 in enumerate(genes):
        for g2 in genes[i+1:]:
            if np.random.random() < 0.55:
                G.add_edge(g1, g2, weight=np.random.uniform(500, 999))

# Inter-module edges (sparse cross-talk)
for mod1, mod2, prob in cross_talk:
    for g1 in module_gene_sets[mod1]:
        for g2 in module_gene_sets[mod2]:
            if np.random.random() < prob:
                G.add_edge(g1, g2, weight=np.random.uniform(150, 500))
```

## WGCNA Co-expression Modules

WGCNA constructs modules from RNA-seq expression data, not PPI databases.

### Workflow
1. Pearson correlation matrix of gene expression
2. Soft threshold: raise adjacency to power β (scale-free topology test → pick β where R² > 0.85)
3. Topological Overlap Measure (TOM): connectivity-based adjacency (more robust than raw correlation)
4. Hierarchical clustering of (1 − TOM) distance
5. Dynamic tree cut → modules
6. Module eigengene (ME) = first PC of module expression

### Key Outputs
- **Module-trait correlation**: which modules associate with clinical outcomes?
- **Module membership kME**: correlation of each gene with its module eigengene
- **Hub gene**: highest intra-module connectivity within a module

```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

# Module eigengene
def module_eigengene(expr_matrix, gene_indices):
    """First PC of module expression. expr_matrix: (samples × genes)."""
    mod_expr = StandardScaler().fit_transform(expr_matrix[:, gene_indices])
    return PCA(n_components=1).fit_transform(mod_expr).flatten()

# Soft thresholding (adjacency = |corr|^beta)
corr = np.corrcoef(expr_matrix.T)
beta = 6  # pick from scale-free topology plot
adjacency = np.abs(corr) ** beta
TOM = adjacency  # simplified; real TOM adds shared neighbor connectivity
```

### WGCNA vs PPI Modules

| Aspect | WGCNA | STRING/PPI |
|--------|-------|------------|
| Data source | Expression matrix | Protein interactions |
| Edge meaning | Co-expression | Physical/functional |
| Context-specific | Yes | No |
| Interpretation | Co-regulated genes | Partners/complexes |

## Visualization

```python
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.cm as cm

pos = nx.spring_layout(G, seed=42, k=3.0)  # high k separates modules
colors = [cm.tab10(partition[n] / max(partition.values())) for n in G.nodes()]
nx.draw_networkx(G, pos=pos, node_color=colors, node_size=150,
                 edge_color='gray', width=0.5, font_size=5)
```

## Cytoscape Export

```python
# Export for Cytoscape (py4cytoscape or save as edgelist)
import py4cytoscape as p4c
p4c.create_network_from_networkx(G)          # requires Cytoscape running
# Or: save edgelist + node attributes as TSV
nx.write_edgelist(G, 'network.edgelist', data=['weight'])
```

## Pitfalls

- **Resolution limit**: Louvain/Leiden cannot detect modules smaller than √(2m) nodes; use hierarchical approaches for small sub-communities
- **Disconnected communities (Louvain bug)**: use Leiden instead; or check `nx.is_connected(G.subgraph(comm))` after partitioning
- **WGCNA soft threshold selection**: always verify scale-free topology (R² of log(k) vs log(P(k)) > 0.85) before using a beta value
- **Module eigengene sign is arbitrary**: first PC can flip sign between runs; use absolute correlation for trait associations
- **Dense networks hide modules**: pre-filter STRING interactions to confidence ≥ 700 before community detection; all-connected graphs give Q ≈ 0
- **Batch effects**: always check for batch confounding before interpreting module-trait correlations
- **Multiple testing**: apply Benjamini-Hochberg FDR for module-trait correlation p-values
