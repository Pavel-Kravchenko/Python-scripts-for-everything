---
name: metabolomics
description: "Metabolomics analysis workflow: LC-MS data processing, metabolite identification, statistical analysis, and pathway mapping. Use when conducting metabolomics experiments."
tool_type: python
primary_tool: COBRApy
---

# metabolomics

LC-MS metabolomics preprocessing, metabolite identification, and pathway enrichment analysis.

## Quick Reference

| Step | Tool | Notes |
|------|------|-------|
| Peak picking | XCMS / MZmine | centWave algorithm for LC-MS |
| Normalization | PQN, IS-based | Probabilistic Quotient Normalization |
| Exact mass match | HMDB / KEGG | ± 5 ppm tolerance |
| MS/MS matching | GNPS / MassBank | Cosine similarity threshold > 0.7 |
| Structure prediction | SIRIUS + CSI:FingerID | De novo from MS/MS |
| Statistics | limma / t-test | Log-transform before testing |
| Pathway enrichment | MetaboAnalyst MSEA | ORA or QEA |
| FBA modeling | COBRApy | `model.optimize()` |

## pyOpenMS: Load mzML Data

```python
from pyopenms import MSExperiment, MzMLFile
import numpy as np

exp = MSExperiment()
MzMLFile().load('sample.mzML', exp)

print(f'Number of spectra: {exp.size()}')

# Get MS1 spectra
ms1_spectra = [s for s in exp if s.getMSLevel() == 1]

# Extract first spectrum
spec = ms1_spectra[0]
mz, intensity = spec.get_peaks()
rt = spec.getRT()
print(f'RT: {rt:.2f} s, Peaks: {len(mz)}')
```

## Feature Table Processing (pandas)

```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load feature table (rowssamples, columnsfeatures)
features = pd.read_csv('feature_table.csv', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)

# Log transform (handle zeros)
features_log = np.log1p(features)

# Probabilistic Quotient Normalization (PQN)
def pqn_normalize(df):
    reference = df.median(axis=0)  # median sample as reference
    quotients = df.div(reference)
    quotient_medians = quotients.median(axis=1)
    return df.div(quotient_medians, axis=0)

features_norm = pqn_normalize(features_log)

# PCA
scaler = StandardScaler()
X_scaled = scaler.fit_transform(features_norm)
pca = PCA(n_components=2)
scores = pca.fit_transform(X_scaled)

# Plot
groups = metadata['condition']
for g in groups.unique():
    idx = groups == g
    plt.scatter(scores[idx, 0], scores[idx, 1], label=g, alpha=0.7)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(); plt.title('PCA of metabolomics data')
```

## Exact Mass Matching

```python
from rdkit.Chem import Descriptors, MolFromSmiles
import pandas as pd

# HMDB metabolite database (simplified)
def match_exact_mass(observed_mz, adduct='[M+H]+', ppm_tol=5):
    """Match observed m/z to HMDB metabolites within ppm tolerance."""
    adduct_mass = {'[M+H]+': 1.00728, '[M+Na]+': 22.98922,
                   '[M-H]-': -1.00728, '[M+NH4]+': 18.03437}
    neutral_mass = observed_mz - adduct_mass.get(adduct, 0)
    # Query HMDB CSV (pre-downloaded)
    hmdb = pd.read_csv('hmdb_metabolites.csv')  # monoisotopic_mass, name, formula
    ppm_error = abs((hmdb['monoisotopic_mass'] - neutral_mass) / neutral_mass * 1e6)
    matches = hmdb[ppm_error <= ppm_tol].copy()
    matches['ppm_error'] = ppm_error[ppm_error <= ppm_tol]
    return matches.sort_values('ppm_error')
```

## COBRApy Flux Balance Analysis

```python
import cobra
import cobra.test

# Load E. coli core model
model = cobra.test.create_test_model('textbook')

# Optimize (maximize biomass by default)
with model:
    solution = model.optimize()
    print(f'Growth rate: {solution.objective_value:.4f} h⁻¹')
    print(solution.fluxes[solution.fluxes.abs() > 0.1].sort_values())

# Gene knockout screen
essential_genes = []
for gene in model.genes:
    with model:
        gene.knock_out()
        sol = model.optimize()
        if sol.objective_value < 0.01:
            essential_genes.append(gene.id)
print(f'Essential genes: {len(essential_genes)}')

# Flux variability analysis
from cobra.flux_analysis import flux_variability_analysis
fva = flux_variability_analysis(model, fraction_of_optimum=0.9)
print(fva.sort_values('maximum', ascending=False).head(10))
```

## MSI Identification Levels

| Level | Evidence | Example |
|-------|----------|---------|
| 1 | Exact mass + MS/MS + RT vs authentic standard | Glucose confirmed by standard |
| 2 | Spectral match to library (no standard) | GNPS library hit (cosine > 0.7) |
| 3 | Putative: chemical class from MS/MS fragments | "Sphingolipid" from fragment pattern |
| 4 | Uncharacterized | "Feature 1234" |

## Pitfalls
- **In-source fragmentation**: fragments appear as MS1 features; check MS2 spectra
- **Adduct confusion**: `[M+Na]+` vs `[M+H]+` differ by 21.98 Da; annotate all adducts
- **Batch effects**: inject QC pools every 10 samples; LOESS-correct to QC
- **Missing values**: > 50% missing → likely noise; min/2 imputation for MCAR
- **Log transform before t-test**: raw intensities are log-normal distributed

## Module
Tier 3 · Module 36 (Metabolomics)
