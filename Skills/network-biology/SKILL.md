---
name: network-biology
description: Biological network analysis — PPI networks from STRING, NetworkX graph metrics, Louvain community detection, GO enrichment per module, Cytoscape export, gene regulatory network inference with GENIE3
primary_tool: NetworkX
---

# Network Biology

## When to Use

Use this skill when:
- Building protein-protein interaction (PPI) networks from STRING/BioGRID
- Computing network topology metrics (degree, centrality, clustering)
- Detecting functional modules in biological networks
- Inferring gene regulatory networks from expression data
- Integrating network modules with differential expression results

## Quick Reference

| Task | Tool | Key Method |
|------|------|-----------|
| Retrieve PPI data | STRING REST API | `requests.post(string_api, data=params)` |
| Build graph | NetworkX | `nx.from_pandas_edgelist(df)` |
| Degree centrality | NetworkX | `dict(G.degree())` |
| Betweenness | NetworkX | `nx.betweenness_centrality(G)` |
| Community detection | python-louvain | `best_partition(G)` |
| GO enrichment | gseapy | `gp.enrichr(gene_list, gene_sets='GO_BP')` |
| Export to Cytoscape | NetworkX | `nx.write_graphml(G, 'net.graphml')` |
| GRN inference | GENIE3 (R) | `GENIE3(exprMatrix, regulators=tfs)` |

## Key Patterns

**Pattern 1: STRING PPI network**
```python
import requests
import pandas as pd
import networkx as nx

genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']
string_api = 'https://string-db.org/api/json/network'
params = {
    'identifiers': '%0d'.join(genes),
    'species': 9606,
    'required_score': 700  # high confidence
}
resp = requests.post(string_api, data=params)
interactions = pd.DataFrame(resp.json())

G = nx.from_pandas_edgelist(
    interactions,
    source='preferredName_A',
    target='preferredName_B',
    edge_attr='score'
)
```

**Pattern 2: Network metrics**
```python
degree = dict(G.degree())
betweenness = nx.betweenness_centrality(G)
closeness = nx.closeness_centrality(G)
clustering = nx.clustering(G)

metrics = pd.DataFrame({
    'degree': degree,
    'betweenness': betweenness,
    'closeness': closeness,
    'clustering': clustering
})
# Hub genes: top 10 by degree
hubs = metrics[metrics['degree'] >= metrics['degree'].quantile(0.9)]
```

**Pattern 3: Louvain community detection**
```python
from community import best_partition  # pip install python-louvain

partition = best_partition(G, random_state=42)
# Map nodes to communities
modules = {}
for node, mod_id in partition.items():
    modules.setdefault(mod_id, []).append(node)
print(f'{len(modules)} modules detected')
```

**Pattern 4: GO enrichment per module**
```python
import gseapy as gp

for mod_id, genes in modules.items():
    if len(genes) < 5:
        continue
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets='GO_Biological_Process_2021',
        organism='Human',
        outdir=None
    )
    top = enr.results.head(5)[['Term', 'Adjusted P-value']]
    print(f'Module {mod_id} ({len(genes)} genes): {top.to_string()}')
```

**Pattern 5: GENIE3 GRN inference (R)**
```r
library(GENIE3)

# exprMatrix: genes  samples, rows  genes
tfs <- c('TP53', 'MYC', 'E2F1', 'FOXM1')
weight_matrix <- GENIE3(exprMatrix=expr_matrix, regulators=tfs,
                         nTrees=1000, nCores=8)
link_list <- getLinkList(weight_matrix, reportMax=20000)
# link_list: regulatoryGene  targetGene  weight
```

## Biological Network Properties

| Property | Random Network | Scale-Free (Biological) |
|----------|----------------|------------------------|
| Degree distribution | Poisson | Power law (P(k) ~ k^-γ) |
| Hubs | Absent | Present (oncogenes, TFs) |
| Robustness | Uniform | Robust to random, vulnerable to hub removal |
| Clustering | Low | High (modular) |
| Path length | √N | log(N) |

## Pitfalls

- **Literature bias** — STRING PPI data is biased toward well-studied genes; hub genes may reflect research attention, not true biology
- **Directionality** — PPI networks are undirected (physical); GRNs are directed (regulatory); do not mix
- **Community detection randomness** — Louvain is stochastic; use `random_state` and run multiple times
- **Score threshold** — STRING combined_score > 700 = high confidence; lower thresholds inflate network with false positives
- **GRN validation** — GENIE3 weights are correlational; validate against known TF binding data (JASPAR, TRRUST, ENCODE ChIP-seq)

## Code Templates

### Degree Distribution Plot
```python
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def plot_degree_distribution(G, title='Degree Distribution'):
    degrees = [d for _, d in G.degree()]
    counts = Counter(degrees)
    k = sorted(counts.keys())
    pk = [counts[ki] / len(degrees) for ki in k]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.bar(k, pk, color='steelblue', alpha=0.7)
    ax1.set(xlabel='Degree k', ylabel='P(k)', title='Linear scale')

    ax2.loglog(k, pk, 'o', color='steelblue', alpha=0.7)
    ax2.set(xlabel='Degree k', ylabel='P(k)', title='Log-log (power law check)')

    plt.suptitle(title)
    plt.tight_layout()
    return degrees

degs = plot_degree_distribution(G, title=f'Network ({G.number_of_nodes()} nodes)')
```

### Identify Hub Genes and Bottlenecks
```python
import networkx as nx
import pandas as pd

def find_network_key_nodes(G, top_n=20):
    metrics = pd.DataFrame({
        'degree':      dict(G.degree()),
        'betweenness': nx.betweenness_centrality(G, normalized=True),
        'closeness':   nx.closeness_centrality(G),
        'pagerank':    nx.pagerank(G),
    })
    metrics['composite'] = (
        metrics['degree'] / metrics['degree'].max() +
        metrics['betweenness'] / metrics['betweenness'].max() +
        metrics['pagerank'] / metrics['pagerank'].max()
    ) / 3
    return metrics.sort_values('composite', ascending=False).head(top_n)

key_nodes = find_network_key_nodes(G)
print(key_nodes[['degree', 'betweenness', 'pagerank']].head(10))
```

### Overlap Network Modules with DEGs
```python
import pandas as pd

def module_deg_overlap(modules, degs, background_size=20000):
    """Fisher's exact test for overlap between modules and DEG list."""
    from scipy.stats import fisher_exact
    results = []
    deg_set = set(degs)
    for mod_id, genes in modules.items():
        module_set = set(genes)
        a = len(module_set & deg_set)          # overlap
        b = len(module_set) - a                # in module, not DEG
        c = len(deg_set) - a                   # DEG, not in module
        d = background_size - a - b - c        # neither
        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        results.append({'module': mod_id, 'size': len(module_set),
                        'overlap': a, 'odds_ratio': odds, 'pvalue': pval})
    df = pd.DataFrame(results).sort_values('pvalue')
    df['padj'] = df['pvalue'] * len(df)  # Bonferroni
    return df

overlap = module_deg_overlap(modules, degs=significant_genes)
print(overlap[overlap['padj'] < 0.05])
```

### Export Network to GraphML (Cytoscape)
```python
import networkx as nx

def export_for_cytoscape(G, metrics_df, output_path):
    """Add node attributes from metrics and save as GraphML."""
    for node, row in metrics_df.iterrows():
        if node in G:
            G.nodes[node]['degree'] = int(row['degree'])
            G.nodes[node]['betweenness'] = float(row['betweenness'])
            G.nodes[node]['pagerank'] = float(row['pagerank'])
    nx.write_graphml(G, output_path)
    print(f"Saved {G.number_of_nodes()} nodes, {G.number_of_edges()} edges to {output_path}")

export_for_cytoscape(G, key_nodes, 'network.graphml')
```

## Related Skills
- `multi-omics-integration` — integrating network modules with expression data
- `ml-deep-learning-bio` — graph neural networks for molecular/biological graphs
- `data-visualization-bio` — heatmaps, volcano plots, pathway enrichment plots
- `python-core-bio` — data structures, dictionaries, set operations
