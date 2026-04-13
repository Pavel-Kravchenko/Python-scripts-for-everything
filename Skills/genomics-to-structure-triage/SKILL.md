---
name: genomics-to-structure-triage
description: Workflow for routing prioritized genomic variants into AlphaFold2/3 or RoseTTAFold structural follow-up.
primary_tool: NumPy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# genomics-to-structure-triage

## When to Use

Use this skill when:
- Translating variant prioritization into structural hypothesis testing
- Choosing AF2 vs AF3 vs RoseTTAFold based on biological context
- Combining variant-level and structure-confidence signals into one rank

## Quick Reference

| Scenario | Preferred model |
|---|---|
| Monomer baseline | AlphaFold2 |
| Complex / nucleic acids / ligands | AlphaFold3 |
| Rosetta-centric open workflows | RoseTTAFold |

## Key Patterns

**Pattern 1: structure model routing**
```python
def choose_structure_model(has_complex, has_ligand_or_nucleic, need_open_rosetta_workflow):
    if has_ligand_or_nucleic:
        return 'AlphaFold3'
    if need_open_rosetta_workflow:
        return 'RoseTTAFold'
    if has_complex:
        return 'AlphaFold3'
    return 'AlphaFold2'
```python

**Pattern 2: final priority score**
```python
import numpy as np

def final_priority(structure_priority, mean_plddt, interface_pae):
    conf = 0.6 * (mean_plddt / 100.0) + 0.4 * (1.0 - np.clip(interface_pae / 30.0, 0, 1))
    return structure_priority * conf
```python

## Code Templates

### Coding-only structure priority
```python
def structure_priority(coding, max_ds, expr_delta, missense_prob, rarity_score):
    if not coding:
        return 0.0
    return (
        0.45 * missense_prob +
        0.25 * abs(expr_delta) +
        0.20 * rarity_score +
        0.10 * max_ds
    )
```python

## Common Pitfalls

- Sending non-coding variants directly into structure pipelines
- Overinterpreting low-confidence structural regions
- Ignoring assay/clinical context when ranking structural candidates

## Related Skills

- `alphafold-structure-prediction`
- `splicing-variant-models`
- `genomic-foundation-models`

