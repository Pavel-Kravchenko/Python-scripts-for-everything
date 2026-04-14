---
name: bio-applied-molecular-modeling
description: Molecular Modeling with NumPy
tool_type: python
primary_tool: NumPy
---

# Molecular Modeling

**Libraries:** `rdkit`, `numpy`, `matplotlib`, `pandas`, `mdanalysis`
**External tools:** GROMACS, AutoDock Vina, SWISS-MODEL, Modeller, PyMOL

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
np.random.seed(42)
```

## Applications

| Application | What modeling provides | Example |
|------------|----------------------|--------|
| **Drug discovery** | Binding poses and affinities | Virtual screening |
| **Protein engineering** | Effects of mutations | Thermostable enzymes |
| **Mechanism understanding** | Conformational changes | Ion channel gating |
| **Structure prediction** | Models without crystal structure | Homology models |
| **Materials design** | Biomaterial properties | Self-assembling peptides |

**Levels of theory:** QM (electrons explicit, small systems) → MM/force fields (parametrized potentials, full proteins) → coarse-grained (sacrifices atomic detail for timescale). Always a tradeoff: accuracy vs. computational cost.

**QM methods:** Ab initio (Hartree-Fock, no empirical params; expensive), DFT (electron density, ~100–500 atoms), semi-empirical (AM1/PM3/PM7, quick estimates). QM needed for bond breaking/forming, charge distributions, FF parametrization, QM/MM hybrids.

```python
def lennard_jones(r, epsilon=1.0, sigma=1.0):
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

r = np.linspace(0.9, 3.0, 500)
r_min = 2**(1/6)  # equilibrium distance
```

## Force Fields

$$E_{\text{total}} = E_{\text{bonds}} + E_{\text{angles}} + E_{\text{dihedrals}} + E_{\text{electrostatic}} + E_{\text{vdW}}$$

**Bond stretching:** $E = \frac{1}{2} k_b (r - r_0)^2$

**Angle bending:** $E = \frac{1}{2} k_\theta (\theta - \theta_0)^2$

**Dihedral:** $E = \sum_n V_n [1 + \cos(n\phi - \gamma)]$

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

r = np.linspace(0.8, 2.0, 300)
r0, kb = 1.5, 300
axes[0, 0].plot(r, 0.5 * kb * (r - r0)**2, 'b-', lw=2)
axes[0, 0].axvline(x=r0, color='r', linestyle='--'); axes[0, 0].set_title('Bond Stretching')

theta = np.linspace(90, 150, 300)
axes[0, 1].plot(theta, 0.5 * 50 * (np.radians(theta) - np.radians(120))**2, 'g-', lw=2)
axes[0, 1].set_title('Angle Bending')

phi = np.linspace(-180, 180, 300)
E_dih = 0.5*(1+np.cos(np.radians(phi))) + 0.2*(1-np.cos(2*np.radians(phi))) + 0.8*(1+np.cos(3*np.radians(phi)))
axes[1, 0].plot(phi, E_dih, 'm-', lw=2); axes[1, 0].set_title('Dihedral Potential')

r_nb = np.linspace(2.5, 10.0, 300)
E_lj = 4*0.24*((3.4/r_nb)**12 - (3.4/r_nb)**6)
E_coul = 332.0637 * 0.5 * (-0.5) / r_nb
axes[1, 1].plot(r_nb, E_lj, 'b-', label='vdW'); axes[1, 1].plot(r_nb, E_coul, 'r-', label='Coulomb')
axes[1, 1].plot(r_nb, E_lj + E_coul, 'k--', label='Total'); axes[1, 1].set_ylim(-2, 2)
axes[1, 1].legend(fontsize=9); axes[1, 1].set_title('Non-bonded')
plt.tight_layout(); plt.show()
```

**Common force fields:**

| Force Field | Use | Notes |
|------------|-----|-------|
| AMBER | Proteins, nucleic acids | ff14SB, ff19SB |
| CHARMM | Proteins, lipids, small molecules | CGenFF for drug-like molecules |
| OPLS-AA | Proteins, organic molecules | Optimized for liquid-state thermodynamics |
| GROMOS | Proteins, lipids | United-atom (no aliphatic H); efficient |
| MARTINI | Coarse-grained | 4:1 heavy atom mapping; membranes |

Never mix parameters from different force fields.

## Energy Minimization

Remove steric clashes before simulation. **Steepest descent:** robust, slow near minimum — good initial relaxation. **Conjugate gradient:** uses prior steps, faster near minimum. **L-BFGS:** quasi-Newton, approximates Hessian, fastest for smooth surfaces. GROMACS: steepest descent until max force < 1000 kJ/mol/nm, then optionally conjugate gradient.

```python
def rosenbrock(x, y, a=1, b=100):
    return (a - x)**2 + b * (y - x**2)**2

def steepest_descent(x0, y0, lr=0.001, steps=3000):
    path = [(x0, y0)]; x, y = x0, y0
    for _ in range(steps):
        dx = -2*(1-x) + b*2*(y-x**2)*(-2*x)  # gradient x
        dy = 100*2*(y-x**2)                    # gradient y (b=100)
        x -= lr * (-2*(1-x) + 200*(y-x**2)*(-2*x))
        y -= lr * (200*(y-x**2))
        path.append((x, y))
    return np.array(path)
```

## Homology Modeling

Proteins with similar sequences adopt similar structures. Workflow:

```
Target sequence → Template search (BLAST/HHpred/SWISS-MODEL)
→ Target-template alignment (critical step)
→ Model building (copy coords, build loops, add side chains)
→ Refinement (energy minimization, MD relaxation)
→ Quality check (Ramachandran, QMEAN, ProSA)
```

## Pitfalls

- **Coordinate systems**: BED 0-based half-open; VCF/GFF 1-based inclusive — mixing causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) for thousands of features
