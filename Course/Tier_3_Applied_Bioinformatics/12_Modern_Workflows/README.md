# 12. Modern Bioinformatics Workflows

**Tier 3 – Applied Bioinformatics**

This module covers contemporary tools and practices that are essential for professional bioinformatics work. You will learn to analyze single-cell data, build reproducible pipelines with Snakemake and nf-core/Nextflow, and implement testing best practices.

---

## Learning Objectives

After completing this module, you will be able to:

- Analyze **single-cell RNA-seq** data using the standard Scanpy workflow
- Write **production-quality Snakefiles** with config files, wildcards, params, logs, benchmarks, conda environments, and containers
- Scale Snakemake to a **SLURM / AWS / GCP cluster** using executor plugins and profiles
- Use the **Nextflow DSL2** dataflow model (channels, processes, workflows)
- Run, configure, and customise **nf-core pipelines** (rnaseq, sarek, and more)
- Add nf-core **modules** to your own Nextflow pipeline
- Write **unit tests** for bioinformatics code with pytest
- Set up **continuous integration** with GitHub Actions

---

## Notebooks

| # | Notebook | Topics | Cells |
|---|----------|--------|-------|
| 1 | [Single-Cell Analysis with Scanpy](01_single_cell_scanpy.ipynb) | AnnData structure, QC metrics, normalization, PCA, UMAP, Leiden clustering, marker genes | ~45 |
| 2 | [Workflow Engines](02_workflow_engines.ipynb) | DAG model; Snakemake rules/wildcards/config/conda/containers/cluster; Nextflow DSL2 processes/channels; nf-core pipelines and modules | ~55 |
| 3 | [Testing and CI/CD](03_testing_cicd.ipynb) | pytest basics, fixtures, parametrized tests, GitHub Actions, code coverage, linting | ~40 |

---

## Prerequisites

- **Tier 2**: Core Bioinformatics (especially BioPython and file formats)
- **Tier 1**: Python fundamentals (functions, classes, file I/O)
- Basic command line experience

---

## Software Requirements

The following packages are required (included in `requirements.txt`):

```
scanpy>=1.9.0
anndata>=0.8.0
snakemake>=7.0.0
pytest>=7.0.0
pytest-cov>=4.0.0
```

---

## Additional Resources

### Single-Cell Analysis
- [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Best practices for single-cell analysis (Luecken & Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746)
- [scverse ecosystem](https://scverse.org/)

### Workflow Engines
- [Snakemake documentation](https://snakemake.readthedocs.io/)
- [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [Snakemake executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/)
- [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
- [Nextflow training](https://training.nextflow.io/)
- [Nextflow DSL2 documentation](https://www.nextflow.io/docs/latest/dsl2.html)
- [nf-core pipelines](https://nf-co.re/)
- [nf-core/configs (institutional profiles)](https://github.com/nf-core/configs)
- [nf-core/modules](https://github.com/nf-core/modules)

### Testing
- [pytest documentation](https://docs.pytest.org/)
- [GitHub Actions documentation](https://docs.github.com/en/actions)
- [Codecov for coverage reporting](https://codecov.io/)

---

[← Previous Module: Clinical Genomics](../11_Clinical_Genomics/) | [Back to Course Overview](../../README.md)
