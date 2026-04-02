# 11. Modern Bioinformatics Workflows

This module covers contemporary tools and practices that are essential for modern bioinformatics work.

## Learning Objectives

After completing this module, you will be able to:
- Analyze single-cell RNA-seq data with Scanpy
- Build reproducible pipelines with Snakemake and Nextflow
- Write tests for bioinformatics code with pytest
- Set up continuous integration with GitHub Actions

## Notebooks

| # | Notebook | Topics |
|---|----------|--------|
| 1 | [Single-Cell Analysis with Scanpy](01_single_cell_scanpy.ipynb) | AnnData structure, QC, normalization, PCA/UMAP, clustering, marker genes |
| 2 | [Workflow Engines](02_workflow_engines.ipynb) | Snakemake rules, Nextflow processes/channels, parallelization, containers |
| 3 | [Testing & CI/CD](03_testing_cicd.ipynb) | pytest, fixtures, GitHub Actions, linting, coverage |

## Prerequisites

- Tier 2: Working with Biological Data
- Basic Python programming (Tier 1)
- Familiarity with command line

## Additional Resources

- [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Snakemake documentation](https://snakemake.readthedocs.io/)
- [Nextflow training](https://training.nextflow.io/)
- [nf-core pipelines](https://nf-co.re/)
- [pytest documentation](https://docs.pytest.org/)
