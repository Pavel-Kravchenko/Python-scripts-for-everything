---
name: bio-core-gene-ontology
description: "Gene Ontology structure, evidence codes, enrichment analysis with hypergeometric test and BH correction"
tool_type: python
primary_tool: NumPy
---

# Gene Ontology

## Key Concepts

- **True Path Rule**: Annotation to a specific term implies annotation to all ancestors. Enrichment analysis must propagate up the DAG or it undercounts genes in broad terms.
- **Background set**: Use all expressed genes in your experiment, not the whole genome — wrong background inflates/deflates enrichment scores.
- **Multiple testing**: Testing ~20,000 GO terms at α=0.05 → ~1,000 false positives. Always apply BH FDR; report adjusted p-values.
- **Gene ID mapping**: GO DBs use Entrez/UniProt IDs; your data may use symbols or Ensembl IDs. Mapping errors silently shrink your list — report mapping success rate.
- **ORA vs GSEA**: ORA requires binary gene list. GSEA uses a ranked list (fold change, p-value) and is more powerful. Use GSEA when you have quantitative rankings.

```python
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# GO ontology structure:
# GO:0008150 - biological_process (BP)
# GO:0003674 - molecular_function (MF)
# GO:0005575 - cellular_component (CC)
```

## Evidence Codes

```
EXP, IDA, IPI, IMP, IGI, IEP  — Experimental (gold standard, quality 4–5)
ISS, ISO, ISA, ISM, IBA        — Computational/sequence-based (quality 3)
TAS                            — Traceable Author Statement (quality 2)
NAS                            — Non-traceable Author Statement (quality 1)
IEA                            — Electronic Annotation (automated, lowest, quality 1)
```

```python
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
    'TAS': 2, 'NAS': 1, 'IEA': 1,
}

def filter_annotations_by_quality(gene_annotations, min_quality=3):
    return {
        gene: [(go_id, ec) for go_id, ec in annots if EVIDENCE_QUALITY.get(ec, 0) >= min_quality]
        for gene, annots in gene_annotations.items()
        if any(EVIDENCE_QUALITY.get(ec, 0) >= min_quality for _, ec in annots)
    }
```

## GO Enrichment Analysis: Hypergeometric Test

2×2 contingency table:

```
                  | In GO term | Not in term | Total
In gene list      |     k      |    n - k    |   n
Not in gene list  |   K - k    |  N-K-n+k    | N - n
Total             |     K      |    N - K    |   N

N = background size; K = genes annotated to term in background
n = genes in your list; k = overlap
P(X >= k) = hypergeometric survival function
```

```python
def propagate_annotations(gene_annotations, go_terms):
    """Apply True Path Rule: propagate each annotation to all ancestor terms."""
    propagated = {}
    for gene, annots in gene_annotations.items():
        all_terms = set()
        for go_id, _ec in annots:
            all_terms.add(go_id)
            all_terms |= get_ancestors(go_id, go_terms)
        propagated[gene] = all_terms
    return propagated


def go_enrichment(gene_list, gene_annotations, go_terms, background_size=20000):
    propagated = propagate_annotations(gene_annotations, go_terms)

    term_to_genes = {}
    for gene, terms in propagated.items():
        for t in terms:
            term_to_genes.setdefault(t, set()).add(gene)

    gene_set = set(g.upper() for g in gene_list)
    n, N = len(gene_set), background_size

    results = []
    for term_id, term_genes in term_to_genes.items():
        annotated_in_list = gene_set & {g.upper() for g in term_genes}
        k = len(annotated_in_list)
        if k == 0:
            continue
        K = max(k, int(len(term_genes) / len(propagated) * N)) if propagated else k
        p_value = stats.hypergeom.sf(k - 1, N, K, n)
        expected = n * K / N
        info = go_terms.get(term_id, {})
        results.append({
            'go_id': term_id,
            'name': info.get('name', 'unknown'),
            'domain': info.get('domain', '?'),
            'k': k, 'K': K,
            'p_value': p_value,
            'fold_enrichment': k / expected if expected > 0 else float('inf'),
            'genes': sorted(annotated_in_list),
        })

    results.sort(key=lambda r: r['p_value'])

    # Benjamini-Hochberg FDR
    m = len(results)
    for i, r in enumerate(results):
        r['fdr'] = min(r['p_value'] * m / (i + 1), 1.0)

    return results
```

```python
de_genes = ['TP53', 'BAX', 'CASP3', 'CASP9', 'BCL2', 'CDK2', 'CCNB1', 'CCND1', 'RB1', 'E2F1']
enrichment_results = go_enrichment(de_genes, GENE_ANNOTATIONS, GO_TERMS, background_size=20000)

print(f"{'GO ID':<14} {'Domain':<4} {'Name':<42} {'k':>3} {'p-value':>10} {'FDR':>10}")
for r in enrichment_results[:12]:
    print(f"{r['go_id']:<14} {r['domain']:<4} {r['name'][:40]:<42} {r['k']:>3} {r['p_value']:>10.2e} {r['fdr']:>10.2e}")
    print(f"{'':>14} Genes: {', '.join(r['genes'])}")
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing causes off-by-one errors.
- **Batch effects**: Always check for batch confounding before interpreting biological signal.
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously.
