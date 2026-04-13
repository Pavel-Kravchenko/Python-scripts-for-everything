---
name: ai-science-variant-to-structure-models
description: "**Tier 5 — Modern AI for Science | Module 05 · Notebook 5**"
tool_type: python
source_notebook: "Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/05_variant_to_structure_models.ipynb"
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# From DNA Variants to Protein Structure: AlphaFold2, AlphaFold3, RoseTTAFold

*Source: Course notebook `Tier_5_Modern_AI_for_Science/05_Genomic_Foundation_Models/05_variant_to_structure_models.ipynb`*

# From DNA Variants to Protein Structure: AlphaFold2, AlphaFold3, RoseTTAFold

**Tier 5 — Modern AI for Science | Module 05 · Notebook 5**

*Prerequisites: Module 04 (AlphaFold & Protein Design), Notebook 3 (Splicing Models)*

---

**By the end of this notebook you will be able to:**
1. Build a practical genomics-to-structure triage workflow
2. Prioritize coding variants for structural follow-up
3. Select AF2 vs AF3 vs RoseTTAFold by use-case
4. Integrate sequence-level and structure-level evidence for interpretation

## Why this notebook matters

Computational variant interpretation is a multi-step triage pipeline: first you identify which variants exist (variant calling), then which are biologically plausible (genomic model scores), then which ones are worth spending experimental effort on (structural follow-up). This notebook covers the final stage: how to decide which variants deserve structural modeling, which model to use, and how to integrate structural confidence with genomic scores into an interpretable ranking.

## How to work through this notebook

1. The decision framework (Section 1) is the central concept — read it as a flowchart, not just text.
2. The triage table (Section 2) integrates scores from the previous four notebooks: splice (NB3), expression (NB2), missense probability, and rarity.
3. The model selection logic (Section 3) should be understood as heuristics, not rigid rules.
4. Section 4 shows how to downweight high-priority variants that land in low-confidence structural regions.

## Common sticking points

- **Coding vs non-coding**: non-coding variants (splice, regulatory) rarely benefit from structure modeling. Save structural prediction for missense variants where mechanism depends on protein conformation.
- **RoseTTAFold2 vs AF2**: RoseTTAFold2 (and RoseTTAFold-AA) are strong alternatives to AF2 for protein-protein interfaces and offer tighter integration with Rosetta design tools. Performance is generally comparable; model choice depends on downstream workflow.
- **AF3 access**: as of 2024–2025, AlphaFold3 is available through the AlphaFold Server (alphafoldserver.com) with rate limits, and as open-source code for non-commercial use. There is no pip install — you run it via the provided scripts.
- **Confidence thresholds for decision-making**: for clinical/experimental decision-making, require pLDDT ≥ 70 AND interface PAE ≤ 15 Å before using a structural prediction to guide hypothesis generation.

## 1. Decision Framework

A common path in variant interpretation:
1. Prioritize variants from genomic models (splice/expression/regulatory)
2. Keep coding or coding-adjacent candidates
3. Map to protein regions (domain, interface, active-site proximity)
4. Run structural follow-up and assess confidence

Not every variant needs structure modeling. Use it when mechanism plausibly involves protein conformation or interactions.

```python
import numpy as np
import pandas as pd

np.random.seed(23)
```

## 2. Variant Triage Table

We combine mock scores from multiple model families:
- splice (`max_ds`)
- regulatory expression (`expr_delta`)
- coding impact prior (`missense_prob`)
- population rarity (`rarity_score`)

Then we compute a composite structural-follow-up priority.

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

## 3. Which Structure Model to Use?

| Scenario | Prefer |
|---|---|
| Single-chain monomer, fast baseline | **AlphaFold2** |
| Complex with nucleic acids / ligands / modified residues | **AlphaFold3** |
| Open Rosetta-centric workflows and interface exploration | **RoseTTAFold** |

Model choice depends on biological context, not just benchmark averages.

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

## 4. Integrating Confidence with Variant Priority

After structure prediction, integrate confidence metrics (e.g., pLDDT/PAE-like proxies) to avoid overinterpreting weak models.

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

## 5. Practical Notes

- Structure prediction is strongest for generating hypotheses, not final proof.
- Keep uncertainty explicit; low-confidence regions should not drive high-stakes decisions.
- Pair structural predictions with orthogonal evidence (functional assays, literature, patient context).

## Optional CLI References (commented)

```python
# !docker run ... alphafold2
# !docker run ... alphafold3
# !python RoseTTAFold/run_e2e_ver.sh input.fa out_dir
```

These steps are intentionally commented to keep this notebook portable.

## Summary

- Use structure modeling selectively after genomic prioritization.
- Choose AF2/AF3/RoseTTAFold based on biological question and complex context.
- Combine variant score and structural confidence into one interpretable ranking.

## Source-backed Context

- AF2/AF3 and RoseTTAFold are best treated as complementary structural follow-up tools after genomic prioritization, not direct substitutes for genomic models.
- AF3 source documentation emphasizes complex-level prediction use cases.

## Validated Sources

Checked online during content expansion.

- [AlphaFold2 paper](https://www.nature.com/articles/s41586-021-03819-2)
- [AlphaFold3 paper](https://www.nature.com/articles/s41586-024-07487-w)
- [RoseTTAFold repository](https://github.com/RosettaCommons/RoseTTAFold)
