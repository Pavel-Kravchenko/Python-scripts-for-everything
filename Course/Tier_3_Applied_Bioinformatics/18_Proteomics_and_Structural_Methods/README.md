# 3.18 Proteomics and Structural Methods
**Tier 3: Applied Bioinformatics**

This module bridges experimental and computational biology: it covers mass spectrometry-based proteomics (from raw spectra to quantitative protein lists), computational protein engineering (conservation analysis, directed evolution, stability prediction), and the structural determination methods — X-ray crystallography, cryo-EM, and NMR — that underpin rational design. The programming exercises implement core algorithms from scratch: peptide mass calculation, in silico trypsin digestion, peptide mass fingerprinting, FDR estimation, and conservation scoring from a multiple sequence alignment.

## Topics Covered

- **Ionization methods**: MALDI vs. ESI — principles, charge state distributions, m/z formula
- **Mass analyzers**: TOF, Orbitrap, Q-TOF, ion trap, triple quadrupole — resolution and accuracy tradeoffs
- **MS/MS fragmentation**: b and y ion series, CID/HCD, de novo sequencing
- **Bottom-up proteomics workflow**: reduction, digestion, nano-LC, DDA data acquisition
- **In silico trypsin digestion**: cleavage rules, missed cleavages
- **Peptide mass fingerprinting**: mass tolerance matching, scoring
- **Database search engines**: Mascot, X!Tandem, MaxQuant, MSFragger — concepts and comparison
- **Target-decoy approach**: FDR estimation at PSM, peptide, and protein levels
- **Label-free quantification**: spectral counting, intensity-based LFQ, iBAQ
- **Isotope labeling**: SILAC (metabolic), TMT/iTRAQ (isobaric chemical tags)
- **DDA vs. DIA acquisition**: stochastic selection vs. comprehensive fragmentation
- **Quantitative data analysis**: median normalization, MA plots, volcano plots
- **Protein engineering strategies**: directed evolution, rational design, de novo design — historical development and tradeoffs (Suplatov/ФББ curriculum)
- **Conservation scoring**: Shannon entropy, mutability score from MSA
- **Stability prediction**: ΔΔG framework, FoldX concepts, scanning mutagenesis
- **X-ray crystallography**: Bragg's law, resolution, R-factor, phase problem, solution methods (MR, SIR, SAD)
- **Electron density maps**: 2Fo-Fc and Fo-Fc difference maps, Uppsala Electron Density Server (EDS) download
- **Resolution and density quality**: five resolution bands (1.0–3.5 Å) and what each reveals
- **Crystal lattices and unit cell**: 7 crystal systems, lattice parameters, CRYST1 record parsing (`parse_cryst1`)
- **Space groups and asymmetric units**: 230 space groups, 65 chiral groups compatible with proteins, symmetry operations (`apply_symmetry_operation`)
- **Cryo-EM**: single-particle analysis workflow, resolution revolution, comparison with X-ray
- **NMR spectroscopy**: chemical shifts, NOE restraints, structure calculation, dynamics, size limitations
- **Method selection**: decision framework for choosing X-ray / cryo-EM / NMR based on sample properties
- **PDB header interpretation**: EXPDTA, RESOLUTION, R-VALUE, CRYST1, SMTRY records

## Notebooks

| Notebook | Description |
|----------|-------------|
| [01_proteomics_and_structural_methods.ipynb](01_proteomics_and_structural_methods.ipynb) | Mass spectrometry proteomics, protein engineering computational design, and structural determination methods — theory, Python implementations, and 5 exercises |

## Prerequisites

- **2.07 Protein Structure** — PDB format, Bio.PDB, DSSP, Ramachandran plots (not repeated here)
- **3.09 Molecular Modeling and Docking** — force fields, MD simulation concepts (recommended; not required)

---
[← Previous: 3.17 Genome Assembly](../17_Genome_Assembly_and_Advanced_NGS/) | [Course Overview](../../README.md) | [Next: 3.19 GWAS →](../19_GWAS/)
