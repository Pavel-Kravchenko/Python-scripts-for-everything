---
name: network-biology
description: Biological network analysis — PPI networks from STRING, NetworkX graph metrics, Louvain community detection, GO enrichment per module, Cytoscape export, gene regulatory network inference with GENIE3
---

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
# Hub genes: top 10% by degree
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

# exprMatrix: genes × samples, rows = genes
tfs <- c('TP53', 'MYC', 'E2F1', 'FOXM1')
weight_matrix <- GENIE3(exprMatrix=expr_matrix, regulators=tfs,
                         nTrees=1000, nCores=8)
link_list <- getLinkList(weight_matrix, reportMax=20000)
# link_list: regulatoryGene | targetGene | weight
```

## Biological Network Properties

| Property | Random Network | Scale-Free (Biological) |
|----------|----------------|------------------------|
| Degree distribution | Poisson | Power law (P(k) ~ k^-γ) |
| Hubs | Absent | Present (oncogenes, TFs) |
| Robustness | Uniform | Robust to random, vulnerable to hub removal |
| Clustering | Low | High (modular) |
| Path length | √N | log(N) |

## Common Pitfalls

- **Literature bias** — STRING PPI data is biased toward well-studied genes; hub genes may reflect research attention, not true biology
- **Directionality** — PPI networks are undirected (physical); GRNs are directed (regulatory); do not mix
- **Community detection randomness** — Louvain is stochastic; use `random_state` and run multiple times
- **Score threshold** — STRING combined_score > 700 = high confidence; lower thresholds inflate network with false positives
- **GRN validation** — GENIE3 weights are correlational; validate against known TF binding data (JASPAR, TRRUST, ENCODE ChIP-seq)
