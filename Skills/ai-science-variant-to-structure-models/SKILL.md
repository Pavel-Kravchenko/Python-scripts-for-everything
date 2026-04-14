---
name: ai-science-variant-to-structure-models
description: "From DNA Variants to Protein Structure: AlphaFold2, AlphaFold3, RoseTTAFold with NumPy"
tool_type: python
primary_tool: NumPy
---

# From DNA Variants to Protein Structure: AlphaFold2, AlphaFold3, RoseTTAFold

## Decision Framework

1. Prioritize variants from genomic models (splice/expression/regulatory)
2. Keep coding or coding-adjacent candidates
3. Map to protein regions (domain, interface, active-site proximity)
4. Run structural follow-up; assess confidence

Non-coding variants (splice, regulatory) rarely benefit from structure modeling. Save structural prediction for missense variants where mechanism depends on protein conformation.

**Confidence thresholds for decision-making:** require pLDDT ≥ 70 AND interface PAE ≤ 15 Å before using a structural prediction to guide hypothesis generation.

**AF3 access:** available via alphafoldserver.com (rate-limited) and open-source for non-commercial use — no pip install, run via provided scripts.

**RoseTTAFold2 vs AF2:** strong alternative for protein-protein interfaces; tighter Rosetta design integration; comparable performance — model choice depends on downstream workflow.

```python
import numpy as np
import pandas as pd

np.random.seed(23)
```

## Variant Triage Table

Combines splice (`max_ds`), regulatory expression (`expr_delta`), coding impact (`missense_prob`), and population rarity (`rarity_score`) into a composite structural-follow-up priority.

```python
variants = pd.DataFrame([
    {"var": "v1", "gene": "GENE1", "coding": 1, "max_ds": 0.12, "expr_delta": 0.45, "missense_prob": 0.80, "rarity_score": 0.90},
    {"var": "v2", "gene": "GENE2", "coding": 0, "max_ds": 0.91, "expr_delta": 0.15, "missense_prob": 0.05, "rarity_score": 0.70},
    {"var": "v3", "gene": "GENE3", "coding": 1, "max_ds": 0.20, "expr_delta": 0.30, "missense_prob": 0.65, "rarity_score": 0.85},
    {"var": "v4", "gene": "GENE4", "coding": 1, "max_ds": 0.05, "expr_delta": 0.10, "missense_prob": 0.40, "rarity_score": 0.60},
])

def structure_priority(row):
    if row["coding"] == 0:
        return 0.0
    return (
        0.45 * row["missense_prob"] +
        0.25 * abs(row["expr_delta"]) +
        0.20 * row["rarity_score"] +
        0.10 * row["max_ds"]
    )

variants["structure_priority"] = variants.apply(structure_priority, axis=1)
variants.sort_values("structure_priority", ascending=False)
```

## Model Selection

| Scenario | Prefer |
|---|---|
| Single-chain monomer, fast baseline | **AlphaFold2** |
| Complex with nucleic acids / ligands / modified residues | **AlphaFold3** |
| Rosetta-centric workflows and interface exploration | **RoseTTAFold** |

```python
def choose_structure_model(has_complex: bool, has_ligand_or_nucleic: bool, need_open_rosetta_workflow: bool) -> str:
    if has_ligand_or_nucleic:
        return "AlphaFold3"
    if need_open_rosetta_workflow:
        return "RoseTTAFold"
    if has_complex:
        return "AlphaFold3"
    return "AlphaFold2"

print(choose_structure_model(False, False, False))
print(choose_structure_model(True, False, False))
print(choose_structure_model(False, True, False))
print(choose_structure_model(False, False, True))
```

## Integrating Confidence with Variant Priority

Downweight high-priority variants landing in low-confidence structural regions.

```python
structure_results = pd.DataFrame([
    {"var": "v1", "mean_plddt": 84, "interface_pae": 6.0},
    {"var": "v3", "mean_plddt": 72, "interface_pae": 14.0},
    {"var": "v4", "mean_plddt": 58, "interface_pae": 22.0},
])

merged = variants.merge(structure_results, on="var", how="left")
merged["confidence_factor"] = (
    0.6 * (merged["mean_plddt"].fillna(50) / 100.0) +
    0.4 * (1.0 - np.clip(merged["interface_pae"].fillna(30) / 30.0, 0, 1))
)
merged["final_priority"] = merged["structure_priority"] * merged["confidence_factor"]
merged.sort_values("final_priority", ascending=False)[["var", "gene", "structure_priority", "confidence_factor", "final_priority"]]
```

## Key Points

- Structure prediction generates hypotheses, not final proof
- Low-confidence regions must not drive high-stakes decisions
- Pair structural predictions with orthogonal evidence (functional assays, literature, patient context)
- Use AF2/AF3/RoseTTAFold based on biological question and complex context
- Combine variant score and structural confidence into one interpretable ranking

## References

- [AlphaFold2 paper](https://www.nature.com/articles/s41586-021-03819-2)
- [AlphaFold3 paper](https://www.nature.com/articles/s41586-024-07487-w)
- [RoseTTAFold repository](https://github.com/RosettaCommons/RoseTTAFold)

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
