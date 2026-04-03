# 2. Variant Calling and SNP Analysis

**Tier 3: Applied Bioinformatics**

The end-to-end pipeline for identifying and annotating genetic variants from sequencing data. Covers SNP and indel calling with GATK and bcftools, VCF format, variant annotation, and population genetics concepts.

## Topics Covered

- Variant types: SNPs, indels, structural variants, CNVs
- Variant calling pipeline: duplicate marking, BQSR, GATK HaplotypeCaller, bcftools, FreeBayes
- VCF format: header, INFO fields, FORMAT/genotype columns, FILTER
- Hard filtering vs. VQSR
- Variant annotation with SnpEff and VEP
- Population genetics: allele frequencies, Hardy-Weinberg equilibrium, LD
- GWAS concepts and ClinVar

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_variant_calling_and_snp_analysis.ipynb](01_variant_calling_and_snp_analysis.ipynb) | Variant calling pipeline, VCF format, annotation, and population genetics |

## Prerequisites

- Tier 3 Module 1: NGS Fundamentals

---

[<- Previous Module](../01_NGS_Fundamentals/) | [Back to Course Overview](../../README.md) | [Next Module ->](../03_RNA_seq_Analysis/)
