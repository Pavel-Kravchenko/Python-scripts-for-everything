---
name: bio-applied-metabolic-flux
description: Flux balance analysis and metabolic modeling with COBRApy. Use when predicting metabolic fluxes, simulating gene knockouts, or analyzing stoichiometric models.
tool_type: python
primary_tool: COBRApy
---

# Metabolic Flux Analysis and Constraint-Based Modeling

- [COBRApy documentation](https://cobrapy.readthedocs.io/)
- [Recon3D model](https://www.vmh.life/)
- [BiGG Models database](http://bigg.ucsd.edu/)

## Core Concepts

Stoichiometric matrix S, flux vector v, steady-state: S*v = 0. Reaction bounds define reversibility. Objective function typically maximizes biomass. Exchange reactions model nutrient uptake.

## FBA with COBRApy

Load model (Recon3D, E. coli core). Set medium constraints via exchange reaction bounds. `model.optimize()` maximizes biomass, returns flux distribution. Inspect fluxes for central carbon metabolism pathways.

## Gene Knockout and Drug Target Prediction

- Single gene deletion: `cobra.flux_analysis.single_gene_deletion(model)`
- Synthetic lethality: `cobra.flux_analysis.double_gene_deletion(model)`
- Compare essential genes to DepMap CRISPR screen data for validation

## Transcriptomics-Constrained FBA

- **GIMME**: penalize flux through low-expression reactions
- **iMAT**: maximize agreement between predicted and measured expression
- **E-Flux**: constrain reaction bounds proportional to expression

## 13C Metabolic Flux Analysis

Isotope labeling (U-13C glucose) + mass isotopologue distribution measurement. Fit fluxes with INCA or OpenMebius. Resolves parallel pathways (glycolysis vs PPP).

## Pitfalls

- **Unbounded fluxes**: Always set upper bounds on exchange reactions; unbounded models give infinite growth
- **Infeasible models**: Missing biomass components or incorrect GPR rules cause infeasible solutions -- check `model.optimize().status`
- **Gene knockout artifacts**: Single-gene knockouts may show no effect due to isozymes; check GPR rules before interpreting
