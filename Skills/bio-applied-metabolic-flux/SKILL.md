---
name: bio-applied-metabolic-flux
description: "Flux balance analysis and metabolic modeling with COBRApy. Use when predicting metabolic fluxes, simulating gene knockouts, or analyzing stoichiometric models."
tool_type: python
source_notebook: "Tier_3_Applied_Bioinformatics/36_Metabolomics/03_metabolic_flux.ipynb"
primary_tool: COBRApy
---

## Version Compatibility

Reference examples tested with: cobrapy 0.29+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# Metabolic Flux Analysis and Constraint-Based Modeling

*Source: Course notebook `Tier_3_Applied_Bioinformatics/36_Metabolomics/03_metabolic_flux.ipynb`*


**Tier 3 — Applied Bioinformatics | Module 36 · Notebook 3**

*Prerequisites: Notebook 2 (Metabolite Identification)*

---

**By the end of this notebook you will be able to:**
1. Explain stoichiometric metabolic network models (GEMs) and flux balance analysis (FBA)
2. Load a genome-scale metabolic model using COBRApy
3. Run FBA to predict growth rate and metabolic fluxes
4. Perform flux variability analysis (FVA) and gene knockout simulations
5. Integrate transcriptomic data into GEMs using GIMME / iMAT algorithms



**Key resources:**
- [COBRApy documentation](https://cobrapy.readthedocs.io/)
- [Recon3D model](https://www.vmh.life/)
- [BiGG Models database](http://bigg.ucsd.edu/)
- [AGORA2 microbial models](https://www.vmh.life/#microbes/overview)

## 1. Genome-Scale Metabolic Models

> Stoichiometric matrix S. Flux vector v. Steady-state assumption: S·v = 0. Reaction bounds (reversible/irreversible). Objective function (biomass maximization). Exchange reactions for nutrient uptake.

## 2. Flux Balance Analysis with COBRApy

> Load Recon3D or E. coli core model. Set medium constraints. Maximize biomass. Inspect flux distribution for central carbon metabolism pathways.

```python
# Example: FBA with COBRApy
# import cobra
# model = cobra.io.load_model('e_coli_core')
# solution = model.optimize()
# print(f'Objective value: {solution.objective_value:.4f}')
# print(solution.fluxes[solution.fluxes.abs() > 0.01])
```python

## 3. Gene Knockout and Drug Target Prediction

> Single gene deletion scanning. Synthetic lethality screen (two-gene knockouts). Identify essential genes. Compare to DepMap CRISPR screen data for validation.

## 4. Transcriptomics-Constrained FBA

> GIMME algorithm: penalize flux through low-expression reactions. iMAT: maximize agreement between predicted and measured gene expression. E-Flux method.

## 5. 13C Metabolic Flux Analysis (13C-MFA)

> Isotope labeling experiments (U-13C glucose). Measure mass isotopologue distributions. Fit fluxes to isotope data with INCA or OpenMebius. Resolve parallel pathways (glycolysis vs PPP).

## Common Pitfalls

- **Unbounded fluxes**: Always set upper bounds on exchange reactions; unbounded models give infinite growth
- **Infeasible models**: Missing biomass components or incorrect GPR rules cause infeasible solutions — check model.optimize().status
- **Gene knockout artifacts**: Single-gene knockouts may show no effect due to isozymes; check GPR rules before interpreting
