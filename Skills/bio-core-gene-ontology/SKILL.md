---
name: bio-core-gene-ontology
description: package and adapt the example to match the actual API rather than retrying.
tool_type: python
primary_tool: NumPy
---

# Gene Ontology

## Complicated moments explained

- **The True Path Rule**: In GO, if a gene is annotated to a specific term, it is implicitly annotated to *all* its ancestor terms. For example, a protein annotated to 'protein kinase activity' (GO:0004672) is automatically annotated to 'catalytic activity', 'molecular function', etc. Enrichment analysis must propagate annotations up the DAG, or it will undercount genes in broad terms.
- **Background gene set matters enormously**: The hypergeometric test compares your gene list to a background set (all genes in the genome, all expressed genes in your experiment, or all genes in your database). Using the wrong background inflates or deflates enrichment scores. Use all expressed genes in your experiment as the background, not the whole genome.
- **Multiple testing correction is mandatory**: Testing 20,000 GO terms at α=0.05 would generate ~1,000 false positives by chance. Always apply Benjamini-Hochberg (FDR) correction. Report adjusted p-values (q-values or FDR), not raw p-values.
- **Gene ID mapping errors are common**: GO databases use Entrez Gene IDs or UniProt accessions; your analysis may use gene symbols or Ensembl IDs. Mapping errors (gene not found, multiple matches, deprecated IDs) silently reduce your gene list. Always report the number of genes successfully mapped.
- **ORA vs. GSEA**: Over-Representation Analysis (ORA) requires a binary gene list (in/out). GSEA uses a ranked list and is more powerful because it uses quantitative information (fold change, p-value). Use GSEA when you have a ranked list; use ORA when you have a hard threshold-based list.

## Environment check (run this first)

```python
# Environment check
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

print("Core imports ready (numpy, scipy, matplotlib).")
print()
print("GO ontology structure at a glance:")
print("  GO:0008150 - biological_process (BP)")
print("    GO:0009987 - cellular process")
print("    GO:0006950 - response to stress")
print("    GO:0006915 - apoptotic process  <-- example term")
print()
print("  GO:0003674 - molecular_function (MF)")
print("    GO:0003824 - catalytic activity")
print("    GO:0005488 - binding")
print()
print("  GO:0005575 - cellular_component (CC)")
print("    GO:0005623 - cell")
print("    GO:0005737 - cytoplasm")
print()
print("Evidence codes (most to least reliable):")
print("  EXP, IDA, IPI, IMP, IGI, IEP -> Experimental (gold standard)")
print("  ISS, ISO, ISA, ISM, IGC, IBA -> Computational (inferred)")
print("  IEA -> Inferred from Electronic Annotation (automated, lowest confidence)")
```python


### GO Annotations and Evidence Codes

A **GO annotation** links a gene product to a GO term, supported by an **evidence code** indicating how the annotation was determined.

```python
+----------+----------------------------------------------+----------+
| Code     | Meaning                                      | Quality  |
+----------+----------------------------------------------+----------+
| Experimental Evidence                                               |
+----------+----------------------------------------------+----------+
| EXP      | Inferred from Experiment                     | High     |
| IDA      | Inferred from Direct Assay                   | High     |
| IPI      | Inferred from Physical Interaction           | High     |
| IMP      | Inferred from Mutant Phenotype               | High     |
| IGI      | Inferred from Genetic Interaction            | High     |
| IEP      | Inferred from Expression Pattern             | High     |
+----------+----------------------------------------------+----------+
| Computational Evidence                                              |
+----------+----------------------------------------------+----------+
| ISS      | Inferred from Sequence Similarity            | Medium   |
| ISO      | Inferred from Sequence Orthology             | Medium   |
| ISA      | Inferred from Sequence Alignment             | Medium   |
| IBA      | Inferred from Biological Aspect of Ancestor  | Medium   |
+----------+----------------------------------------------+----------+
| Author / Curator Evidence                                           |
+----------+----------------------------------------------+----------+
| TAS      | Traceable Author Statement                   | Medium   |
| NAS      | Non-traceable Author Statement               | Low      |
+----------+----------------------------------------------+----------+
| Automatic                                                           |
+----------+----------------------------------------------+----------+
| IEA      | Inferred from Electronic Annotation          | Low      |
+----------+----------------------------------------------+----------+
```python

When filtering annotations, experimental evidence (EXP, IDA, IMP, ...) is the gold standard. IEA annotations are the most abundant but least reliable.

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
```python


## GO Enrichment Analysis

Given a list of interesting genes (e.g., differentially expressed genes from an RNA-seq experiment), we ask: **are any GO terms over-represented in this list compared to the background?**

### The Hypergeometric Test

The statistical model is the **hypergeometric distribution** (equivalent to Fisher's exact test for a 2x2 contingency table).

```python
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
```python

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
```python

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
```python

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
