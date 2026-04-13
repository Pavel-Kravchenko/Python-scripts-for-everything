---
name: bio-core-pathways
description: "- Understand the Gene Ontology (GO) system: three ontologies, DAG structure, GO terms - Perform GO enrichment analysis using the hypergeometric and Fisher's exact tests - Use GOATOOLS and the g:Profil"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/11_Gene_Ontology_and_Pathways/02_pathways.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Gene Ontology and Pathway Analysis

*Source: Course notebook `Tier_2_Core_Bioinformatics/11_Gene_Ontology_and_Pathways/02_pathways.ipynb`*

# Gene Ontology and Pathway Analysis

---

## Learning Objectives

By the end of this notebook you will be able to:

- Understand the Gene Ontology (GO) system: three ontologies, DAG structure, GO terms
- Perform GO enrichment analysis using the hypergeometric and Fisher's exact tests
- Use GOATOOLS and the g:Profiler concept for functional enrichment
- Navigate the KEGG database: structure, pathway maps, identifiers
- Access KEGG programmatically through its REST API
- Understand Gene Set Enrichment Analysis (GSEA)
- Predict transmembrane protein topology from sequence
- Build a complete workflow from a gene list to biological interpretation

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section

## Complicated moments explained

- **KEGG pathway membership is dynamic**: KEGG pathways are manually curated and updated regularly. A gene list from 2020 may classify differently in 2024 as new interactions are discovered. Always note the KEGG database release date in your analysis.
- **ORA vs. GSEA**: Over-Representation Analysis (ORA) requires a binary gene list. GSEA uses a continuous ranking (e.g., by fold change or -log10(p-value)) and tests whether members of a gene set are enriched at either end of the ranked list. GSEA is more powerful because it does not require an arbitrary threshold.
- **GSEA leading edge**: The 'leading edge' subset consists of the genes that contribute most to the enrichment score — the core of the signal. These are the biologically most relevant genes to focus on.
- **Kyte-Doolittle hydropathy threshold**: The standard threshold for predicting transmembrane helices is mean hydropathy > 1.6 over a 19-residue window. This threshold was calibrated empirically and works reasonably well but is not universal. Modern tools (TMHMM, DeepTMHMM, Phobius) use machine learning and are far more accurate.
- **The positive-inside rule**: Cytoplasmic loops in membrane proteins tend to have more positively charged residues (Arg, Lys) than extracellular loops. This rule (established by Gunnar von Heijne) is used to predict membrane protein topology and achieves ~85% accuracy.

## Environment check (run this first)

```python
# Environment check
import numpy as np
import matplotlib.pyplot as plt
import urllib.request
import urllib.error

print("Core imports ready.")
print()

# Quick KEGG connectivity test
try:
    url = "https://rest.kegg.jp/list/pathway/hsa"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=5) as response:
        lines = response.read().decode('utf-8').strip().split('\n')
    print(f"KEGG API: connected. {len(lines)} human pathways available.")
except Exception as e:
    print(f"KEGG API: not reachable ({e})")
    print("Offline fallback data is used in this notebook.")
print()
print("Pathway databases:")
print("  KEGG      - curated biochemical/signaling pathways, organism-specific")
print("  Reactome  - mechanistic reactions and complexes, human-focused")
print("  WikiPathways - community-curated, CC0 license")
print("  MSigDB    - gene set collections for GSEA (Hallmarks, GO, KEGG, etc.)")
```

```python
# Gene-to-GO annotations with evidence codes (simplified real-world data)

GENE_ANNOTATIONS = {
    'TP53':  [('GO:0006915', 'IDA'), ('GO:0007049', 'IMP'), ('GO:0097193', 'TAS'),
             ('GO:0003677', 'IDA'), ('GO:0005634', 'IDA')],
    'CDK2':  [('GO:0004672', 'IDA'), ('GO:0000278', 'IDA'), ('GO:0000082', 'IMP'),
             ('GO:0005634', 'IDA'), ('GO:0005737', 'IEA')],
    'BCL2':  [('GO:0006915', 'IDA'), ('GO:0097193', 'IMP'), ('GO:0005739', 'IDA')],
    'CASP3': [('GO:0006915', 'IDA'), ('GO:0003824', 'IDA'), ('GO:0005737', 'IDA')],
    'CCNB1': [('GO:0000278', 'IDA'), ('GO:0005634', 'IDA'), ('GO:0005737', 'IDA')],
    'CCND1': [('GO:0000082', 'IDA'), ('GO:0000278', 'TAS'), ('GO:0005634', 'IDA')],
    'E2F1':  [('GO:0007049', 'IMP'), ('GO:0003677', 'IDA'), ('GO:0006260', 'TAS'),
             ('GO:0005634', 'IDA')],
    'RB1':   [('GO:0007049', 'IDA'), ('GO:0000082', 'IMP'), ('GO:0005634', 'IDA')],
    'CDKN1A':[('GO:0007049', 'IDA'), ('GO:0000082', 'IMP'), ('GO:0005634', 'IDA'),
             ('GO:0005737', 'IEA')],
    'BAX':   [('GO:0006915', 'IDA'), ('GO:0097193', 'IDA'), ('GO:0005739', 'IDA')],
    'MDM2':  [('GO:0006915', 'IMP'), ('GO:0005634', 'IDA'), ('GO:0005737', 'IDA')],
    'CASP9': [('GO:0006915', 'IDA'), ('GO:0097193', 'IDA'), ('GO:0003824', 'IDA'),
             ('GO:0005737', 'IDA')],
}

EVIDENCE_QUALITY = {
    'EXP': 5, 'IDA': 5, 'IPI': 5, 'IMP': 5, 'IGI': 5, 'IEP': 4,
    'ISS': 3, 'ISO': 3, 'ISA': 3, 'ISM': 3, 'IBA': 3,
    'TAS': 2, 'NAS': 1,
    'IEA': 1,
}


def filter_annotations_by_quality(gene_annotations, min_quality=3):
    """Keep only annotations supported by evidence at or above min_quality."""
    filtered = {}
    for gene, annots in gene_annotations.items():
        kept = [(go_id, ec) for go_id, ec in annots if EVIDENCE_QUALITY.get(ec, 0) >= min_quality]
        if kept:
            filtered[gene] = kept
    return filtered


# Show annotations for TP53
print("TP53 annotations (all):")
for go_id, ec in GENE_ANNOTATIONS['TP53']:
    term = GO_TERMS.get(go_id, {})
    print(f"  {go_id} {term.get('name', '?'):45s} [{term.get('domain', '?')}]  evidence={ec} (quality {EVIDENCE_QUALITY[ec]})")

# Filter to experimental evidence only
high_quality = filter_annotations_by_quality(GENE_ANNOTATIONS, min_quality=4)
print(f"\nGenes with high-quality annotations: {len(high_quality)}")
for gene in sorted(high_quality):
    terms = [go_id for go_id, _ in high_quality[gene]]
    print(f"  {gene}: {len(terms)} annotations")
```

---

## 2. GO Enrichment Analysis

Given a list of interesting genes (e.g., differentially expressed genes from an RNA-seq experiment), we ask: **are any GO terms over-represented in this list compared to the background?**

### 2.1 The Hypergeometric Test

The statistical model is the **hypergeometric distribution** (equivalent to Fisher's exact test for a 2x2 contingency table).

```
+--------------------+--------------+--------------+--------+
|                    | In GO term   | Not in term  | Total  |
+--------------------+--------------+--------------+--------+
| In gene list       |     k        |    n - k     |   n    |
| Not in gene list   |    K - k     |  N-K-n+k     |  N - n |
+--------------------+--------------+--------------+--------+
| Total              |     K        |    N - K     |   N    |
+--------------------+--------------+--------------+--------+

  N = total genes in background (e.g., genome)
  K = genes annotated to this GO term in the background
  n = genes in your list
  k = genes in your list AND annotated to this GO term

  P(X >= k) = sum_{i=k}^{min(n,K)} C(K,i)*C(N-K,n-i) / C(N,n)
```

A small p-value means the overlap `k` is larger than expected by chance -- the term is **enriched**.

```python
from scipy import stats
import numpy as np


def propagate_annotations(gene_annotations, go_terms):
    """
    Apply the true path rule: propagate each annotation to all ancestor terms.
    Returns dict  gene -> set of GO IDs (direct + propagated).
    """
    propagated = {}
    for gene, annots in gene_annotations.items():
        all_terms = set()
        for go_id, _ec in annots:
            all_terms.add(go_id)
            all_terms |= get_ancestors(go_id, go_terms)
        propagated[gene] = all_terms
    return propagated


def go_enrichment(gene_list, gene_annotations, go_terms, background_size=20000):
    """
    Perform GO enrichment analysis with the hypergeometric test.

    Parameters
    ----------
    gene_list : list of str
        Genes of interest (e.g., differentially expressed genes).
    gene_annotations : dict
        Gene -> list of (GO_ID, evidence_code).
    go_terms : dict
        GO_ID -> {'name', 'domain', 'parents'}.
    background_size : int
        Total number of genes in the background genome.

    Returns
    -------
    list of dict, sorted by p-value.
    """
    # Step 1: propagate annotations
    propagated = propagate_annotations(gene_annotations, go_terms)

    # Step 2: build term -> gene set mapping
    term_to_genes = {}
    for gene, terms in propagated.items():
        for t in terms:
            term_to_genes.setdefault(t, set()).add(gene)

    gene_set = set(g.upper() for g in gene_list)
    n = len(gene_set)
    N = background_size

    results = []
    for term_id, term_genes in term_to_genes.items():
        # Scale K for the full background (our annotation DB is small)
        annotated_in_list = gene_set & {g.upper() for g in term_genes}
        k = len(annotated_in_list)
        if k == 0:
            continue

        # Estimate K proportionally (in real analysis, use the full annotation DB)
        K = max(k, int(len(term_genes) / len(propagated) * N)) if propagated else k

        # Hypergeometric survival function: P(X >= k)
        p_value = stats.hypergeom.sf(k - 1, N, K, n)

        expected = n * K / N
        fold_enrichment = k / expected if expected > 0 else float('inf')

        info = go_terms.get(term_id, {})
        results.append({
            'go_id': term_id,
            'name': info.get('name', 'unknown'),
            'domain': info.get('domain', '?'),
            'k': k,
            'K': K,
            'p_value': p_value,
            'fold_enrichment': fold_enrichment,
            'genes': sorted(annotated_in_list),
        })

    results.sort(key=lambda r: r['p_value'])

    # Step 3: Benjamini-Hochberg FDR correction
    m = len(results)
    for i, r in enumerate(results):
        r['fdr'] = min(r['p_value'] * m / (i + 1), 1.0)

    return results
```

```python
# Run GO enrichment on an example gene list
# Scenario: RNA-seq identified these genes as differentially expressed in a cancer study

de_genes = ['TP53', 'BAX', 'CASP3', 'CASP9', 'BCL2', 'CDK2', 'CCNB1', 'CCND1', 'RB1', 'E2F1']

enrichment_results = go_enrichment(de_genes, GENE_ANNOTATIONS, GO_TERMS, background_size=20000)

print("GO Enrichment Results")
print("=" * 85)
print(f"{'GO ID':<14} {'Domain':<4} {'Name':<42} {'k':>3} {'p-value':>10} {'FDR':>10}")
print("-" * 85)
for r in enrichment_results[:12]:
    print(f"{r['go_id']:<14} {r['domain']:<4} {r['name'][:40]:<42} {r['k']:>3} {r['p_value']:>10.2e} {r['fdr']:>10.2e}")
    print(f"{'':>14} Genes: {', '.join(r['genes'])}")
```

---

### 2.2 Multiple Testing Correction

When testing thousands of GO terms simultaneously, many will appear significant by chance alone. We must correct for **multiple testing**.

Common methods:

| Method | Description | Stringency |
|--------|-------------|------------|
| Bonferroni | Multiply p by number of tests | Very strict |
| Benjamini-Hochberg (BH) | Controls False Discovery Rate (FDR) | Moderate |
| Benjamini-Yekutieli | FDR for dependent tests | Moderate-strict |

In practice, **BH FDR < 0.05** is the standard threshold for GO enrichment.

```python
def multiple_testing_correction(p_values, method='bh'):
    """
    Apply multiple testing correction.

    Parameters
    ----------
    p_values : list of float
    method : str
        'bonferroni' or 'bh' (Benjamini-Hochberg)

    Returns
    -------
    list of float  (adjusted p-values)
    """
    m = len(p_values)
    if method == 'bonferroni':
        return [min(p * m, 1.0) for p in p_values]

    elif method == 'bh':
        # Sort p-values, apply BH, then restore original order
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        adjusted = [0.0] * m
        prev_adj = 0.0
        for rank, (orig_idx, p) in enumerate(indexed, 1):
            adj = min(p * m / rank, 1.0)
            adjusted[orig_idx] = adj
        # Enforce monotonicity (larger rank => larger adjusted p)
        running_min = 1.0
        for rank in range(m, 0, -1):
            orig_idx = indexed[rank - 1][0]
            running_min = min(running_min, adjusted[orig_idx])
            adjusted[orig_idx] = running_min
        return adjusted

    raise ValueError(f"Unknown method: {method}")


# Demonstrate on example p-values
raw_p = [0.001, 0.03, 0.0001, 0.5, 0.04, 0.0005]
bonf = multiple_testing_correction(raw_p, method='bonferroni')
bh = multiple_testing_correction(raw_p, method='bh')

print(f"{'Raw p':>10} {'Bonferroni':>12} {'BH FDR':>10}")
print("-" * 35)
for p, b, f in zip(raw_p, bonf, bh):
    sig_raw = '*' if p < 0.05 else ' '
    sig_bh = '*' if f < 0.05 else ' '
    print(f"{p:>10.4f}{sig_raw} {b:>10.4f}   {f:>10.4f}{sig_bh}")
print("\n* = significant at alpha=0.05")
```
