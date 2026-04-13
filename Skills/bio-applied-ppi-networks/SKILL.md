---
name: bio-applied-ppi-networks
description: "**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 1**"
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/28_Network_Biology/01_ppi_networks.ipynb"
primary_tool: Matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, networkx 3.2+, numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# PPI Network Construction & Analysis

*Source: Course notebook `Tier_3_Applied_Bioinformatics/28_Network_Biology/01_ppi_networks.ipynb`*

# PPI Network Construction & Analysis

**Tier 3 — Applied Bioinformatics | Module 28 · Notebook 1**

*Prerequisites: Module 05 (Promoter & Regulatory Analysis), Module 07 (Machine Learning)*

---

**By the end of this notebook you will be able to:**
1. Query the STRING database API to retrieve PPI data for a gene set
2. Build a NetworkX graph from PPI interactions
3. Calculate and interpret key network metrics: degree, betweenness centrality, clustering coefficient
4. Identify hub genes (top-degree nodes) and bottleneck nodes
5. Visualize the network with biologically meaningful layouts



**Key resources:**
- [STRING DB API documentation](https://string-db.org/cgi/help?sessionId=&subpage=api)
- [NetworkX documentation](https://networkx.org/documentation/stable/)
- [Network biology review (Barabási & Oltvai, 2004)](https://www.nature.com/articles/nrg1272)

## 1. Graph Theory Basics for Biology

Biological networks are modeled as **graphs** G = (V, E) where:
- **V**: vertices (nodes) = genes, proteins, metabolites, etc.
- **E**: edges = interactions (physical, functional, regulatory)

### Key graph metrics

| Metric | Definition | Biological relevance |
|---|---|---|
| **Degree** | Number of direct neighbors | Hub genes = high degree = essential |
| **Betweenness centrality** | Fraction of shortest paths passing through node | Bottleneck proteins; drug targets |
| **Closeness centrality** | Inverse avg shortest path to all nodes | Global network position |
| **Eigenvector centrality** | Node importance weighted by neighbor importance | Quality > quantity of connections |
| **Clustering coefficient** | Fraction of node's neighbors that are connected | Local network density; module membership |
| **Path length** | Steps between two nodes | Network diameter = max path |

### Scale-free networks in biology
PPI networks often follow a **power-law degree distribution** P(k) ~ k^-γ:
- Most nodes have few connections
- A few "hub" genes have many connections
- Hubs are essential (removal = network collapse)
- Example: TP53, EGFR, MYC are hubs in cancer networks

### Types of biological networks
| Network type | Nodes | Edges | Database |
|---|---|---|---|
| PPI | Proteins | Physical interaction | STRING, BioGRID, IntAct |
| Co-expression | Genes | Correlated expression | WGCNA output |
| Metabolic | Metabolites/reactions | Enzymatic reactions | KEGG, Reactome |
| Gene regulatory | TFs + targets | TF binding | TRRUST, ChEA |

```python
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import defaultdict

np.random.seed(42)

# ----- Construct a synthetic PPI network from cancer driver genes -----
# Simulate STRING-like interaction data

# 25 cancer-related gene names
genes = ['TP53', 'EGFR', 'KRAS', 'MYC', 'PTEN', 'PIK3CA', 'AKT1', 'BRCA1',
         'BRCA2', 'CDK4', 'CDK6', 'RB1', 'ATM', 'CHEK1', 'CHEK2', 'MDM2',
         'CDKN2A', 'MAP2K1', 'BRAF', 'ERBB2', 'VHL', 'HIF1A', 'RAD51',
         'FANCD2', 'PALB2']

n_genes = len(genes)
n_interactions = 80

# Create interaction matrix with biologically inspired structure
# Hub genes (TP53, EGFR, MYC, BRCA1) have higher connectivity
hub_genes = ['TP53', 'EGFR', 'MYC', 'BRCA1']
hub_idx = [genes.index(g) for g in hub_genes]

# Generate edges with confidence scores
edges = []
edge_set = set()
for _ in range(n_interactions * 5):
    if len(edges) >= n_interactions:
        break
    # Hubs more likely to form edges
    probs = np.ones(n_genes)
    probs[hub_idx] = 5
    probs /= probs.sum()
    g1, g2 = np.random.choice(n_genes, size=2, replace=False, p=probs)
    if (g1, g2) not in edge_set and (g2, g1) not in edge_set:
        score = np.random.uniform(150, 999)
        edges.append({'protein1': genes[g1], 'protein2': genes[g2], 'combined_score': score})
        edge_set.add((g1, g2))

df_edges = pd.DataFrame(edges)
print(f"Interactions loaded: {len(df_edges)}")
print(df_edges.head(5).to_string(index=False))
```

## 2. STRING Database and Network Construction

**STRING** (Search Tool for the Retrieval of Interacting Genes/Proteins) integrates:
- Experimental PPI data (yeast 2-hybrid, co-IP, MS-based)
- Computational predictions (co-expression, phylogenetic profiling)
- Text mining (co-mentioned in abstracts)
- Curated databases (KEGG, Reactome)

### Confidence score interpretation
- **0-150**: low confidence (chance interactions)
- **150-400**: medium-low confidence
- **400-700**: medium confidence (commonly used threshold)
- **700-900**: high confidence
- **900-999**: highest confidence (multiple independent sources)

### Accessing STRING
```python
import requests

string_api = "https://string-db.org/api"
genes = "TP53%0dMDM2%0dBRCA1%0dEGFR%0dKRAS"
params = {
    "identifiers": genes,
    "species": 9606,    # human
    "required_score": 400,
    "caller_identity": "bioinformatics_course"
}
response = requests.get(f"{string_api}/json/network", params=params)
interactions = pd.DataFrame(response.json())
```

```python
# ----- Simulating STRING database query -----
import json

# In practice, use the STRING REST API:
# import requests
# string_api_url = "https://string-db.org/api"
# params = {"identifiers": "%0d".join(gene_list),
#           "species": 9606,  # human
#           "caller_identity": "my_app"}
# response = requests.post(string_api_url + "/json/network", data=params)
# df_interactions = pd.DataFrame(response.json())

# Here we demonstrate the response format:
example_response = [
    {"stringId_A": "9606.ENSP00000269305", "stringId_B": "9606.ENSP00000232424",
     "preferredName_A": "TP53", "preferredName_B": "MDM2",
     "combined_score": 998, "score_coexpression": 0.32,
     "score_experimental": 0.87, "score_database": 0.90},
    {"stringId_A": "9606.ENSP00000269305", "stringId_B": "9606.ENSP00000261509",
     "preferredName_A": "TP53", "preferredName_B": "CHEK2",
     "combined_score": 942, "score_coexpression": 0.15,
     "score_experimental": 0.78, "score_database": 0.85},
    {"stringId_A": "9606.ENSP00000396521", "stringId_B": "9606.ENSP00000299299",
     "preferredName_A": "BRCA1", "preferredName_B": "BRCA2",
     "combined_score": 997, "score_coexpression": 0.55,
     "score_experimental": 0.92, "score_database": 0.95},
]

df_example = pd.DataFrame(example_response)
print("Example STRING API response:")
print(df_example[['preferredName_A', 'preferredName_B', 'combined_score',
                    'score_experimental', 'score_database']].to_string(index=False))

# STRING score channels
print("\nSTRING confidence score channels:")
channels = {
    'Neighborhood': 'Genes found in close genomic proximity',
    'Gene fusion': 'Evidence of gene fusion events',
    'Co-occurrence': 'Similar phylogenetic profiles',
    'Coexpression': 'Correlated mRNA expression',
    'Experimental': 'Protein-protein interaction experiments (BiP, 2-hybrid)',
    'Database': 'Curated databases (KEGG, Reactome)',
    'Text mining': 'Co-mentioned in publications',
}
for channel, desc in channels.items():
    print(f"  {channel:<15}: {desc}")
print("\nCombined score = 1 - prod(1 - individual_scores)")
print("Threshold recommendations: 400 (medium), 700 (high), 900 (highest)")
```

## 3. Building the NetworkX Graph

After obtaining interactions from STRING (or other databases), we construct a NetworkX graph:

### Graph types
- **Undirected** (`nx.Graph`): PPI, co-expression
- **Directed** (`nx.DiGraph`): regulatory networks, signaling cascades
- **Weighted** (`weight` edge attribute): interaction confidence scores

### Node and edge attributes
```python
G = nx.from_pandas_edgelist(df, source='gene1', target='gene2',
                             edge_attr=['combined_score'])
# Add node attributes
for gene in G.nodes():
    G.nodes[gene]['degree'] = G.degree(gene)
    G.nodes[gene]['type'] = annotations[gene]
```

### Connected components
Real PPI networks often have:
- One giant connected component (~80% of nodes)
- Several small isolated components (poorly studied genes)
- Isolates (singleton nodes with no interactions)

Always analyze the **largest connected component** for centrality metrics.

```python
# ----- Build NetworkX graph -----

# Threshold at combined score >= 400 (STRING medium confidence)
threshold = 400
df_filtered = df_edges[df_edges['combined_score'] >= threshold].copy()

# Node annotations
gene_annotations = {
    'TP53': 'tumor_suppressor', 'EGFR': 'receptor_kinase', 'KRAS': 'oncogene',
    'MYC': 'transcription_factor', 'PTEN': 'tumor_suppressor', 'PIK3CA': 'kinase',
    'AKT1': 'kinase', 'BRCA1': 'DNA_repair', 'BRCA2': 'DNA_repair',
    'CDK4': 'kinase', 'CDK6': 'kinase', 'RB1': 'tumor_suppressor',
    'ATM': 'DNA_repair', 'CHEK1': 'kinase', 'CHEK2': 'kinase',
    'MDM2': 'ubiquitin_ligase', 'CDKN2A': 'tumor_suppressor', 'MAP2K1': 'kinase',
    'BRAF': 'kinase', 'ERBB2': 'receptor_kinase', 'VHL': 'tumor_suppressor',
    'HIF1A': 'transcription_factor', 'RAD51': 'DNA_repair',
    'FANCD2': 'DNA_repair', 'PALB2': 'DNA_repair'
}

# Build graph
G = nx.from_pandas_edgelist(df_filtered, source='protein1', target='protein2',
                             edge_attr='combined_score')

# Add node attributes
for gene, annotation in gene_annotations.items():
    if gene in G.nodes:
        G.nodes[gene]['function'] = annotation

print(f"Graph summary:")
print(f"  Nodes: {G.number_of_nodes()}")
print(f"  Edges: {G.number_of_edges()}")
print(f"  Density: {nx.density(G):.3f}")
print(f"  Is connected: {nx.is_connected(G)}")
if not nx.is_connected(G):
    components = list(nx.connected_components(G))
    print(f"  Connected components: {len(components)}")
    print(f"  Largest component: {len(max(components, key=len))} nodes")

# Use largest connected component for analysis
G_main = G.subgraph(max(nx.connected_components(G), key=len)).copy()
print(f"\nMain component: {G_main.number_of_nodes()} nodes, {G_main.number_of_edges()} edges")
```

## 4. Network Centrality Measures

Centrality metrics reveal different aspects of node importance:

### Degree centrality
Simple count of direct neighbors. Hub identification. Formula: k_i / (N-1).

### Betweenness centrality
Fraction of shortest paths between all pairs that pass through node i. Identifies **bottleneck** proteins — critical connectors between network modules. High betweenness + low degree = "bottleneck but not hub".

### Closeness centrality
How quickly a node can reach all others via shortest paths. Genes with high closeness are well-positioned to propagate signals globally.

### Eigenvector centrality (PageRank-like)
A node is important if its neighbors are important. Captures quality of connections, not just quantity. Basis of Google PageRank.

### Hub and bottleneck classification
```
             High betweenness    Low betweenness
High degree  Hub + bottleneck    Hub only
Low degree   Bottleneck only     Peripheral
```

```python
# ----- Network centrality analysis -----

# Compute multiple centrality measures
degree_cent = dict(nx.degree_centrality(G_main))
betweenness_cent = nx.betweenness_centrality(G_main, weight=None, normalized=True)
closeness_cent = nx.closeness_centrality(G_main)
eigenvector_cent = nx.eigenvector_centrality(G_main, max_iter=500, tol=1e-6)
degree_raw = dict(G_main.degree())

# Compile into DataFrame
df_centrality = pd.DataFrame({
    'Gene': list(G_main.nodes()),
    'Degree': [degree_raw[g] for g in G_main.nodes()],
    'Degree_centrality': [degree_cent[g] for g in G_main.nodes()],
    'Betweenness': [betweenness_cent[g] for g in G_main.nodes()],
    'Closeness': [closeness_cent[g] for g in G_main.nodes()],
    'Eigenvector': [eigenvector_cent[g] for g in G_main.nodes()],
}).sort_values('Degree', ascending=False)

print("Top 10 hub genes by degree centrality:")
print(df_centrality.head(10)[['Gene','Degree','Betweenness','Closeness','Eigenvector']].to_string(index=False))

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('PPI Network Centrality Analysis', fontsize=13, fontweight='bold')

# Degree distribution
degrees = sorted([d for n, d in G_main.degree()], reverse=True)
axes[0, 0].hist(degrees, bins=15, color='steelblue', alpha=0.8)
axes[0, 0].set_xlabel('Degree'); axes[0, 0].set_ylabel('Count')
axes[0, 0].set_title('Degree Distribution')

# Power law check (log-log)
if len(degrees) > 5:
    unique_deg = sorted(set(degrees))
    deg_count = [degrees.count(d) for d in unique_deg]
    axes[0, 1].loglog(unique_deg, deg_count, 'o-', color='steelblue')
    axes[0, 1].set_xlabel('Degree (log)'); axes[0, 1].set_ylabel('Count (log)')
    axes[0, 1].set_title('Degree Distribution (log-log)\nScale-free test')

# Betweenness vs degree
axes[0, 2].scatter(df_centrality['Degree'], df_centrality['Betweenness'],
                    s=60, alpha=0.7, color='orange')
for _, row in df_centrality.head(5).iterrows():
    axes[0, 2].annotate(row['Gene'], (row['Degree'], row['Betweenness']),
                         fontsize=8, ha='left')
axes[0, 2].set_xlabel('Degree'); axes[0, 2].set_ylabel('Betweenness centrality')
axes[0, 2].set_title('Degree vs Betweenness\n(bottleneck genes have high betweenness, low degree)')

# Top 10 by each metric
metrics_top = ['Degree', 'Betweenness', 'Eigenvector']
for ax, metric in zip(axes[1, :], metrics_top):
    top10 = df_centrality.nlargest(10, metric)
    bars = ax.barh(top10['Gene'], top10[metric], color='steelblue', alpha=0.8)
    ax.set_xlabel(metric)
    ax.set_title(f'Top 10 by {metric}')
    ax.invert_yaxis()

plt.tight_layout()
plt.show()
```
