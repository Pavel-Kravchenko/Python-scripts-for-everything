---
name: bio-applied-docking
description: "Molecular docking: ligand preparation, receptor setup, AutoDock Vina workflow, scoring functions, and binding pose analysis. Use when predicting protein-ligand interactions."
tool_type: python
primary_tool: NumPy
---

# Molecular Docking and Chemoinformatics

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
```

## Levels of Theory

| Level | Method | System Size | Use Case |
|-------|--------|-------------|----------|
| QM (ab initio/DFT) | Schrödinger equation | ~100–500 atoms | Bond breaking, charge distributions, FF parameterization |
| Semi-empirical (AM1/PM3) | QM + empirical params | ~1000 atoms | Quick estimates, QM/MM boundary |
| Molecular mechanics | Force fields | Whole proteins + solvent | MD, docking, free energy |
| Coarse-grained (MARTINI) | 4:1 mapping | Membranes, large assemblies | Long timescales |

**Tradeoff:** accuracy vs computational cost — always choose the simplest level that answers your question.

## Force Fields

Total energy: $E = E_{\text{bonds}} + E_{\text{angles}} + E_{\text{dihedrals}} + E_{\text{elec}} + E_{\text{vdW}}$

**Bonded terms:**
- Bond stretching: $E = \frac{1}{2} k_b (r - r_0)^2$
- Angle bending: $E = \frac{1}{2} k_\theta (\theta - \theta_0)^2$
- Dihedral: $E = \sum_n V_n [1 + \cos(n\phi - \gamma)]$

**Non-bonded:** Lennard-Jones (vdW) + Coulomb (electrostatic)

```python
def lennard_jones(r, epsilon=1.0, sigma=1.0):
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
```

### Common Force Fields

| Force Field | Typical Use | Key Features |
|------------|-------------|---------------|
| **AMBER** (ff14SB/ff19SB) | Proteins, nucleic acids | Widely used biomolecular MD |
| **CHARMM** (CGenFF) | Proteins, lipids, small molecules | Extensive lipid parameters |
| **OPLS-AA** | Proteins, organic molecules | Optimized for liquid-state thermodynamics |
| **GROMOS** | Proteins, lipids | United-atom (no aliphatic H); efficient |
| **MARTINI** | Coarse-grained | 4:1 heavy atom mapping; membranes |

**Never mix parameters from different force fields.**

## Energy Minimization

| Algorithm | Behavior | When to use |
|-----------|----------|-------------|
| Steepest Descent | Follows negative gradient; robust, slow near minimum | Initial clash removal |
| Conjugate Gradient | Uses prior steps to avoid oscillation | Refinement after SD |
| L-BFGS | Approximates Hessian; fastest on smooth surfaces | Final polish |

GROMACS typical: steepest descent until max force < 1000 kJ/mol/nm, then conjugate gradient.

```python
def rosenbrock(x, y, a=1, b=100):
    return (a - x)**2 + b * (y - x**2)**2

def rosenbrock_grad(x, y, a=1, b=100):
    dx = -2*(a - x) + b * 2 * (y - x**2) * (-2*x)
    dy = b * 2 * (y - x**2)
    return np.array([dx, dy])

def steepest_descent(x0, y0, lr=0.001, steps=500):
    path = [(x0, y0)]
    x, y = x0, y0
    for _ in range(steps):
        grad = rosenbrock_grad(x, y)
        x -= lr * grad[0]
        y -= lr * grad[1]
        path.append((x, y))
    return np.array(path)
```

## Homology Modeling

Proteins with similar sequences adopt similar structures. Homology modeling predicts target structure from a known template.

```
Target seq → Template search (BLAST/HHpred) → Alignment → Model building
→ Refinement (minimization/MD) → Quality check (Ramachandran/QMEAN/ProSA)
```

## Pitfalls

- **SMILES canonicalization** — different SMILES can encode the same molecule; always canonicalize before comparison
- **Stereochemistry** — ignoring chirality in fingerprints merges enantiomers with different bioactivity
- **Descriptor scaling** — molecular descriptors span different ranges; always standardize before ML
