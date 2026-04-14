---
name: bio-applied-ppi-networks
description: "PPI network construction and analysis with NetworkX and STRING DB: centrality metrics, hub/bottleneck classification, scale-free properties, and community detection."
tool_type: python
primary_tool: Python
---

# PPI Network Construction and Analysis

## References
- [STRING DB API](https://string-db.org/cgi/help?subpage=api)
- [NetworkX docs](https://networkx.org/documentation/stable/)

## Graph Metrics Quick Reference

| Metric | Definition | Biological use |
|---|---|---|
| **Degree** | Direct neighbor count | Hub genes; essential proteins |
| **Betweenness** | Fraction of shortest paths through node | Bottleneck drug targets |
| **Closeness** | Inverse avg shortest path to all nodes | Signal propagation speed |
| **Eigenvector** | Importance weighted by neighbor importance | Quality over quantity of connections |
| **Clustering coeff** | Fraction of neighbors that are connected | Module membership |

### Hub/bottleneck classification
```
             High betweenness    Low betweenness
High degree  Hub + bottleneck    Hub only
Low degree   Bottleneck only     Peripheral
```

## Scale-Free Networks

PPI networks follow power-law degree distribution P(k) ~ k^-γ:
- Most nodes have few connections; a few hubs (TP53, EGFR, MYC) have many
- Hubs are essential — removal causes network collapse
- Verify with log-log degree distribution plot (should be linear)

## STRING Database

**Score channels:** neighborhood, gene fusion, co-occurrence, coexpression, experimental, database, text mining.
**Combined score** = 1 − ∏(1 − individual_scores).

| Score range | Confidence |
|---|---|
| 0–150 | Low (noise) |
| 400–700 | Medium (common threshold) |
| 700–900 | High |
| 900–999 | Highest (multiple independent sources) |

### Querying STRING API
```python
import requests, pandas as pd

params = {
    "identifiers": "TP53%0dMDM2%0dBRCA1%0dEGFR",
    "species": 9606,          # human
    "required_score": 400,
    "caller_identity": "my_app"
}
response = requests.get("https://string-db.org/api/json/network", params=params)
df = pd.DataFrame(response.json())
# Columns: preferredName_A, preferredName_B, combined_score, score_experimental, ...
```

## Building the NetworkX Graph

```python
import networkx as nx

# From pandas edge list
G = nx.from_pandas_edgelist(df, source='preferredName_A', target='preferredName_B',
                             edge_attr='combined_score')

# Add node attributes
for gene, annotation in gene_annotations.items():
    if gene in G.nodes:
        G.nodes[gene]['function'] = annotation

# Always work on the largest connected component for centrality
G_main = G.subgraph(max(nx.connected_components(G), key=len)).copy()
print(f"Nodes: {G_main.number_of_nodes()}, Edges: {G_main.number_of_edges()}, "
      f"Density: {nx.density(G_main):.3f}")
```

## Centrality Analysis

```python
import pandas as pd

degree_cent     = nx.degree_centrality(G_main)
betweenness     = nx.betweenness_centrality(G_main, normalized=True)
closeness       = nx.closeness_centrality(G_main)
eigenvector     = nx.eigenvector_centrality(G_main, max_iter=500, tol=1e-6)

df_cent = pd.DataFrame({
    'Gene':        list(G_main.nodes()),
    'Degree':      [G_main.degree(g) for g in G_main.nodes()],
    'Betweenness': [betweenness[g] for g in G_main.nodes()],
    'Closeness':   [closeness[g] for g in G_main.nodes()],
    'Eigenvector': [eigenvector[g] for g in G_main.nodes()],
}).sort_values('Degree', ascending=False)

# Classify hub vs bottleneck
deg_thresh = df_cent['Degree'].median()
bet_thresh = df_cent['Betweenness'].median()
df_cent['role'] = df_cent.apply(lambda r:
    'hub+bottleneck' if r.Degree > deg_thresh and r.Betweenness > bet_thresh
    else 'hub' if r.Degree > deg_thresh
    else 'bottleneck' if r.Betweenness > bet_thresh
    else 'peripheral', axis=1)
```

## Pitfalls

- **Always filter by confidence score** before building the graph — unfiltered STRING data includes many low-confidence edges that inflate connectivity.
- **Use largest connected component** for centrality metrics; isolates and small components distort global statistics.
- **Directed vs undirected:** Use `nx.DiGraph` for regulatory networks, `nx.Graph` for PPI/co-expression.
- **`eigenvector_centrality` may not converge** on disconnected graphs — always call on `G_main`, not the full graph.
- **Betweenness is O(VE)** — expensive on large graphs; use `nx.betweenness_centrality(G, k=500)` for approximation.
