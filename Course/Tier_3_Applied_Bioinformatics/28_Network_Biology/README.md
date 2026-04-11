# Module 28: Network Biology

**Tier 3 — Applied Bioinformatics | Module 28**

Protein-protein interaction networks, gene regulatory networks, and pathway topology analysis using NetworkX, Cytoscape, and STRING DB.

## What You'll Learn
- Graph theory fundamentals applied to biological networks
- Protein-protein interaction (PPI) data sources: STRING, BioGRID, IntAct
- Network construction in NetworkX from PPI databases
- Key network metrics: degree distribution, clustering coefficient, betweenness centrality, hub genes
- Network motif detection (feed-forward loops, bi-fan motifs)
- Community detection / module finding (Louvain, Girvan-Newman)
- Gene regulatory network (GRN) inference from expression data (GENIE3, ARACNE concepts)
- Visualization with Cytoscape (via py4cytoscape) and NetworkX layouts
- Pathway network analysis and enrichment in the network context

## Prerequisites
- Module 05 (Promoter & Regulatory Analysis) — regulatory elements
- Module 07 (Machine Learning for Biology) — graph-based methods

## Planned Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [PPI Network Construction](01_ppi_networks.ipynb) | STRING API, NetworkX, hub detection, centrality metrics |
| 2 | [Network Modules & Enrichment](02_network_modules.ipynb) | Community detection, module GO enrichment, visualization |
| 3 | [Gene Regulatory Networks](03_gene_regulatory_networks.ipynb) | GRN inference concepts, GENIE3, network motifs |

## Key Tools
| Tool | Purpose |
|------|---------|
| NetworkX | Graph construction and analysis (Python) |
| py4cytoscape | Cytoscape automation from Python |
| Cytoscape | Interactive network visualization |
| STRING DB API | PPI data retrieval |
| igraph | Fast graph algorithms (Python/R) |
| python-louvain | Community detection |

## Resources
- [Cytoscape user documentation](https://cytoscape.org/documentation_users.html)
- [NetworkX documentation](https://networkx.org/documentation/stable/)
- [STRING DB tutorials](https://string-db.org/cgi/help)
- [py4cytoscape documentation](https://py4cytoscape.readthedocs.io/)
- [Network biology review (Barabási & Oltvai, 2004)](https://www.nature.com/articles/nrg1272)
- [GENIE3 documentation](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)

## Related Skill
`network-biology.md` *(planned)*

---

[← Previous: 3.27 Multi-Omics Integration](../27_Multi_Omics_Integration/) | [Course Overview](../../README.md) | [Next: 3.29 Cheminformatics & Drug Discovery →](../29_Cheminformatics_Drug_Discovery/)
