---
name: genomics-to-structure-triage
description: Route prioritized genomic variants into AlphaFold2/3 or RoseTTAFold structural follow-up.
tool_type: python
primary_tool: NumPy
---

# genomics-to-structure-triage

## Model Routing

| Scenario | Preferred Model |
|---|---|
| Monomer baseline | AlphaFold2 |
| Complex / nucleic acids / ligands | AlphaFold3 |
| Rosetta-centric open workflows | RoseTTAFold |

```python
def choose_structure_model(has_complex, has_ligand_or_nucleic, need_open_rosetta_workflow):
    if has_ligand_or_nucleic:
        return 'AlphaFold3'
    if need_open_rosetta_workflow:
        return 'RoseTTAFold'
    if has_complex:
        return 'AlphaFold3'
    return 'AlphaFold2'
```

## Priority Scoring

```python
import numpy as np

def structure_priority(coding, max_ds, expr_delta, missense_prob, rarity_score):
    """Only coding variants get structure follow-up."""
    if not coding:
        return 0.0
    return 0.45 * missense_prob + 0.25 * abs(expr_delta) + 0.20 * rarity_score + 0.10 * max_ds

def final_priority(structure_priority, mean_plddt, interface_pae):
    """Weight priority by structure confidence."""
    conf = 0.6 * (mean_plddt / 100.0) + 0.4 * (1.0 - np.clip(interface_pae / 30.0, 0, 1))
    return structure_priority * conf
```

## Pitfalls

- Non-coding variants should not go into structure pipelines
- Low-confidence regions (pLDDT < 50) should not be over-interpreted
- Always consider assay/clinical context alongside structural predictions

## Related Skills

- `alphafold-structure-prediction` -- AF2/AF3 usage details
- `genomic-foundation-models` -- upstream variant effect scoring
