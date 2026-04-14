---
name: bio-core-pathways
description: Gene Ontology and pathway enrichment — GO structure, hypergeometric/Fisher test, ORA vs GSEA, Benjamini-Hochberg FDR, KEGG/Reactome API patterns.
tool_type: python
primary_tool: NumPy
---

# Gene Ontology and Pathway Analysis

## Pitfalls

- **KEGG is dynamic**: Pathways are manually curated and updated. A gene list from 2020 may classify differently in 2024. Always note the KEGG release date.
- **ORA vs GSEA**: ORA requires a binary gene list (arbitrary threshold). GSEA uses a continuous ranking (fold change, -log10 p-value) and is more powerful — no threshold needed. For GSEA, the "leading edge" subset is the core signal.
- **True path rule**: GO annotations must be propagated to all ancestor terms before enrichment testing. A gene annotated to "apoptosis" is implicitly also annotated to "cell death" and "biological process." Failing to propagate underestimates term coverage.
- **Background gene set matters**: Use the genes actually tested (expressed in your assay), not the whole genome. Using the wrong background inflates or deflates enrichment.
- **Bonferroni is overly conservative**: For correlated GO terms (a gene's annotation propagates to many ancestors), Benjamini-Hochberg FDR is preferred.
- **IEA evidence codes are low-quality**: Electronically inferred annotations (IEA) are unreviewed. Filter them out for high-confidence analysis; keep only IDA, IMP, IPI, IGI, IEP (experimental) or ISS, IBA (computational, moderate confidence).

## GO Structure

Three orthogonal ontologies, each a Directed Acyclic Graph (DAG):

| Ontology | Abbreviation | Meaning |
|---|---|---|
| Molecular Function | MF | Biochemical activity of the gene product |
| Biological Process | BP | Pathway or larger biological process |
| Cellular Component | CC | Where in the cell the product is active |

## Evidence Code Quality

| Tier | Codes | Quality |
|---|---|---|
| Experimental | EXP, IDA, IPI, IMP, IGI, IEP | High |
| Computational (curated) | ISS, ISO, ISA, ISM, IBA | Medium |
| Traceable/Non-traceable | TAS, NAS | Low |
| Electronic | IEA | Lowest — filter for high-confidence work |

## Hypergeometric Test (ORA)

```
Contingency table:
                  | In GO term | Not in term | Total
------------------+------------+-------------+------
In gene list      |     k      |    n - k    |   n
Not in gene list  |   K - k    |  N-K-n+k    | N - n
Total             |     K      |    N - K    |   N

N = background size (expressed genes)
K = genes annotated to this GO term in background
n = genes in your list
k = overlap
P(X >= k) = hypergeom.sf(k-1, N, K, n)
```

```python
from scipy import stats

def hypergeom_enrichment(k, n, K, N):
    """One-tailed p-value: probability of seeing k or more by chance."""
    return stats.hypergeom.sf(k - 1, N, K, n)

fold_enrichment = (k / n) / (K / N)
```

## GO Enrichment Implementation

```python
from scipy import stats
import numpy as np

def go_enrichment(gene_list, term_to_genes, background_size=20000):
    """
    gene_list: list of gene symbols
    term_to_genes: dict of GO_ID -> set of annotated genes (propagated)
    Returns list of dicts sorted by p-value with BH FDR.
    """
    gene_set = set(g.upper() for g in gene_list)
    n = len(gene_set)
    N = background_size

    results = []
    for term_id, term_genes in term_to_genes.items():
        k = len(gene_set & {g.upper() for g in term_genes})
        if k == 0:
            continue
        K = len(term_genes)
        p = stats.hypergeom.sf(k - 1, N, K, n)
        results.append({'go_id': term_id, 'k': k, 'K': K,
                        'fold_enrichment': (k/n) / (K/N), 'p_value': p})

    results.sort(key=lambda r: r['p_value'])

    # Benjamini-Hochberg FDR
    m = len(results)
    for i, r in enumerate(results):
        r['fdr'] = min(r['p_value'] * m / (i + 1), 1.0)
    # Enforce monotonicity
    running_min = 1.0
    for r in reversed(results):
        running_min = min(running_min, r['fdr'])
        r['fdr'] = running_min

    return results
```

## Multiple Testing Correction

```python
def bh_correction(p_values):
    """Benjamini-Hochberg FDR correction."""
    m = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * m
    for rank, (orig_idx, p) in enumerate(indexed, 1):
        adjusted[orig_idx] = min(p * m / rank, 1.0)
    # Enforce monotonicity
    running_min = 1.0
    for _, (orig_idx, _) in reversed(list(enumerate(indexed))):
        running_min = min(running_min, adjusted[orig_idx])
        adjusted[orig_idx] = running_min
    return adjusted

# Or use statsmodels:
from statsmodels.stats.multitest import multipletests
_, adj_pvals, _, _ = multipletests(p_values, method='fdr_bh')
```

## KEGG REST API

```python
import urllib.request

def kegg_get(endpoint):
    url = f"https://rest.kegg.jp/{endpoint}"
    with urllib.request.urlopen(url, timeout=10) as r:
        return r.read().decode('utf-8')

# List all human pathways
pathways = kegg_get("list/pathway/hsa")

# Get genes in a specific pathway
genes_in_pathway = kegg_get("link/hsa/hsa05210")  # colorectal cancer

# Find pathways containing a gene
pathway_for_gene = kegg_get("link/pathway/hsa:7157")  # TP53
```

## Tool Decision Table

| Scenario | Tool |
|---|---|
| Gene list, binary significant/not | ORA (hypergeometric / Fisher exact) |
| Ranked gene list (all tested genes) | GSEA (fgsea, clusterProfiler) |
| Human pathways with reactions | Reactome (ReactomePA) |
| Community-curated, open-access | WikiPathways |
| Organism-specific metabolism | KEGG |
| All GO terms + FDR | goatools, clusterProfiler, g:Profiler |
