---
name: bio-applied-molecular-modeling
description: Molecular Modeling with NumPy
tool_type: python
primary_tool: NumPy
---

# Molecular Modeling

# Molecular Modeling and Docking

Molecular modeling bridges the gap between static structural data and dynamic biological reality. From predicting how a drug binds its target to understanding protein conformational changes, computational molecular modeling is indispensable in modern bioinformatics and drug discovery.


**Libraries:** `rdkit`, `numpy`, `matplotlib`, `pandas`, `mdanalysis` (concepts)
**External tools discussed:** GROMACS, AutoDock Vina, SWISS-MODEL, Modeller, PyMOL

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

%matplotlib inline
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 12
np.random.seed(42)
```

## Introduction to Molecular Modeling

### Why Model Biomolecules?

Experimental structure determination (X-ray crystallography, cryo-EM, NMR) provides snapshots, but biomolecules are dynamic. Molecular modeling lets us:

| Application | What modeling provides | Example |
|------------|----------------------|--------|
| **Drug discovery** | Predict binding poses and affinities | Virtual screening of millions of compounds |
| **Protein engineering** | Assess effects of mutations | Designing thermostable enzymes |
| **Mechanism understanding** | Visualize conformational changes | Ion channel gating, enzyme catalysis |
| **Structure prediction** | Build models when no crystal structure exists | Homology models for uncharacterized proteins |
| **Materials design** | Predict properties of biomaterials | Self-assembling peptide nanostructures |

### Levels of Theory

Molecular modeling spans multiple levels of approximation. The choice depends on the system size and the question being asked.


**Key insight:** There is always a tradeoff between accuracy and computational cost. Quantum mechanics treats electrons explicitly but is limited to small systems. Molecular mechanics uses parametrized potentials (force fields) and can handle entire proteins in explicit solvent. Coarse-grained models sacrifice atomic detail for access to longer timescales.

### Quantum Mechanics in Brief

Quantum mechanical methods solve (approximately) the Schrödinger equation:

$$\hat{H}\Psi = E\Psi$$

- **Ab initio methods** (Hartree-Fock, post-HF): No empirical parameters; computationally expensive.
- **Density Functional Theory (DFT):** Reformulates the problem in terms of electron density rather than the wavefunction. More efficient, widely used for systems of ~100–500 atoms.
- **Semi-empirical methods** (AM1, PM3, PM7): Use empirical parameters to speed up calculations. Useful for quick estimates.

QM is essential for:
- Chemical reactions (bond breaking/forming)
- Accurate charge distributions
- Parameterizing force fields
- QM/MM hybrid methods (quantum region embedded in a classical environment)

```python
# Visualization: the Lennard-Jones potential (a core non-bonded interaction)
# This is one of the fundamental functions in molecular mechanics

def lennard_jones(r, epsilon=1.0, sigma=1.0):
    """Compute Lennard-Jones potential: V(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)"""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

r = np.linspace(0.9, 3.0, 500)
V = lennard_jones(r)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(r, V, 'b-', linewidth=2)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Distance r (in units of σ)')
ax.set_ylabel('Potential Energy V(r) (in units of ε)')
ax.set_title('Lennard-Jones Potential')
ax.set_ylim(-1.5, 3)
ax.set_xlim(0.85, 3.0)

# Annotate equilibrium distance and well depth
r_min = 2**(1/6)  # equilibrium distance
ax.annotate(f'Equilibrium: r = 2^(1/6)σ ≈ {r_min:.3f}σ',
            xy=(r_min, -1), xytext=(1.8, -0.5),
            arrowprops=dict(arrowstyle='->', color='red'),
            fontsize=11, color='red')
ax.annotate('Repulsive\n(r⁻¹² dominates)', xy=(1.0, 2.0), fontsize=10, color='darkred')
ax.annotate('Attractive\n(r⁻⁶ dominates)', xy=(2.0, -0.3), fontsize=10, color='darkblue')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

## Force Fields and Molecular Mechanics

### The Force Field Concept

A **force field** is a mathematical function and associated parameter set that describes the potential energy of a molecular system. In molecular mechanics, the total energy is decomposed into bonded and non-bonded terms:

$$E_{\text{total}} = E_{\text{bonds}} + E_{\text{angles}} + E_{\text{dihedrals}} + E_{\text{electrostatic}} + E_{\text{van der Waals}}$$

### Bonded Interactions

**Bond stretching** (harmonic approximation):
$$E_{\text{bond}} = \frac{1}{2} k_b (r - r_0)^2$$

**Angle bending:**
$$E_{\text{angle}} = \frac{1}{2} k_\theta (\theta - \theta_0)^2$$

**Dihedral (torsional) angles:**
$$E_{\text{dihedral}} = \sum_n V_n [1 + \cos(n\phi - \gamma)]$$

Where:
- $k_b$, $k_\theta$ are force constants
- $r_0$, $\theta_0$ are equilibrium values
- $V_n$ is the torsional barrier height
- $n$ is the periodicity, $\gamma$ is the phase

```python
# Visualize the four main force field terms

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Bond stretching
r = np.linspace(0.8, 2.0, 300)
r0, kb = 1.5, 300  # equilibrium distance (Å), force constant (kcal/mol/Å²)
E_bond = 0.5 * kb * (r - r0)**2
axes[0, 0].plot(r, E_bond, 'b-', linewidth=2)
axes[0, 0].set_xlabel('Bond length r (Å)')
axes[0, 0].set_ylabel('Energy (kcal/mol)')
axes[0, 0].set_title('Bond Stretching (Harmonic)')
axes[0, 0].axvline(x=r0, color='r', linestyle='--', label=f'r₀ = {r0} Å')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Angle bending
theta = np.linspace(90, 150, 300)
theta0, ka = 120, 50
E_angle = 0.5 * ka * (np.radians(theta) - np.radians(theta0))**2
axes[0, 1].plot(theta, E_angle, 'g-', linewidth=2)
axes[0, 1].set_xlabel('Bond angle θ (degrees)')
axes[0, 1].set_ylabel('Energy (kcal/mol)')
axes[0, 1].set_title('Angle Bending (Harmonic)')
axes[0, 1].axvline(x=theta0, color='r', linestyle='--', label=f'θ₀ = {theta0}°')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Dihedral potential
phi = np.linspace(-180, 180, 300)
V1, V2, V3 = 0.5, 0.2, 0.8
E_dih = V1*(1 + np.cos(np.radians(phi))) + V2*(1 - np.cos(2*np.radians(phi))) + V3*(1 + np.cos(3*np.radians(phi)))
axes[1, 0].plot(phi, E_dih, 'm-', linewidth=2)
axes[1, 0].set_xlabel('Dihedral angle φ (degrees)')
axes[1, 0].set_ylabel('Energy (kcal/mol)')
axes[1, 0].set_title('Dihedral (Torsional) Potential')
axes[1, 0].grid(True, alpha=0.3)

# Non-bonded: Lennard-Jones + Coulomb
r_nb = np.linspace(2.5, 10.0, 300)
sigma, epsilon = 3.4, 0.24  # Å, kcal/mol (typical for carbon)
E_lj = 4 * epsilon * ((sigma/r_nb)**12 - (sigma/r_nb)**6)
# Coulomb with opposite charges
q1, q2 = 0.5, -0.5  # partial charges (e)
k_coul = 332.0637  # kcal·Å/(mol·e²)
E_coul = k_coul * q1 * q2 / r_nb
axes[1, 1].plot(r_nb, E_lj, 'b-', linewidth=2, label='van der Waals (LJ)')
axes[1, 1].plot(r_nb, E_coul, 'r-', linewidth=2, label='Electrostatic (Coulomb)')
axes[1, 1].plot(r_nb, E_lj + E_coul, 'k--', linewidth=2, label='Total non-bonded')
axes[1, 1].set_xlabel('Distance r (Å)')
axes[1, 1].set_ylabel('Energy (kcal/mol)')
axes[1, 1].set_title('Non-bonded Interactions')
axes[1, 1].set_ylim(-2, 2)
axes[1, 1].legend(fontsize=9)
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

### Common Force Fields

| Force Field | Developers | Typical Use | Key Features |
|------------|-----------|-------------|---------------|
| **AMBER** | Kollman group (UCSF) | Proteins, nucleic acids | Widely used for biomolecular MD; includes ff14SB, ff19SB |
| **CHARMM** | Karplus group (Harvard) | Proteins, lipids, small molecules | CGenFF for drug-like molecules; extensive lipid parameters |
| **OPLS-AA** | Jorgensen group (Yale) | Proteins, organic molecules | Optimized for liquid-state thermodynamics |
| **GROMOS** | van Gunsteren group (ETH) | Proteins, lipids | United-atom model (no aliphatic H); efficient |
| **MARTINI** | Marrink group (Groningen) | Coarse-grained simulations | 4:1 heavy atom mapping; membranes, large systems |

**Important:** A force field is only as good as its parametrization. Never mix parameters from different force fields. Always validate your system setup.

### Energy Minimization

Before any simulation, we need to remove bad contacts (steric clashes) by minimizing the potential energy. Common algorithms:

- **Steepest Descent:** Moves along the negative gradient. Robust but slow near the minimum. Good for initial relaxation.
- **Conjugate Gradient:** Uses information from previous steps to avoid oscillation. Faster convergence near the minimum.
- **L-BFGS:** Quasi-Newton method that approximates the Hessian. Fastest for smooth energy surfaces.

In GROMACS, energy minimization is typically performed with steepest descent until the maximum force drops below a threshold (e.g., 1000 kJ/mol/nm), then optionally refined with conjugate gradient.

```python
# Demonstrate steepest descent vs conjugate gradient on a 2D surface

def rosenbrock(x, y, a=1, b=100):
    """Rosenbrock function -- a classic optimization test surface."""
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

# Run steepest descent
path_sd = steepest_descent(-1.5, 2.0, lr=0.001, steps=3000)

# Plot
fig, ax = plt.subplots(figsize=(10, 7))
x_grid = np.linspace(-2, 2, 200)
y_grid = np.linspace(-1, 3, 200)
X, Y = np.meshgrid(x_grid, y_grid)
Z = rosenbrock(X, Y)

cs = ax.contour(X, Y, Z, levels=np.logspace(-0.5, 3.5, 20), cmap='viridis')
ax.plot(path_sd[:, 0], path_sd[:, 1], 'r.-', markersize=1, linewidth=0.5, alpha=0.7, label='Steepest descent path')
ax.plot(-1.5, 2.0, 'ro', markersize=10, label='Start')
ax.plot(1, 1, 'g*', markersize=15, label='Minimum (1, 1)')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Energy Minimization: Steepest Descent on Rosenbrock Surface')
ax.legend()
plt.colorbar(cs, label='Energy')
plt.tight_layout()
plt.show()

print(f"Final position: ({path_sd[-1, 0]:.4f}, {path_sd[-1, 1]:.4f})")
print(f"Final energy:    {rosenbrock(path_sd[-1, 0], path_sd[-1, 1]):.6f}")
print(f"True minimum:    (1.0, 1.0) with energy 0.0")
```

## Homology Modeling

### The Principle

Proteins with similar sequences tend to adopt similar 3D structures. **Homology modeling** (or comparative modeling) exploits this to predict the structure of a protein (target) using a known structure (template) as a scaffold.

The general workflow:

```python
Target sequence
       │
       ▼
┌─────────────────┐
│ Template search  │  ← BLAST, HHpred, SWISS-MODEL
│ (find homologs)  │
└────────┬────────┘
         ▼
┌─────────────────┐
│ Target-template  │  ← Sequence alignment (critical step!)
│ alignment        │
└────────┬────────┘
         ▼
┌─────────────────┐
│ Model building   │  ← Copy coordinates, build loops, add side chains
│                  │
└────────┬────────┘
         ▼
┌─────────────────┐
│ Refinement       │  ← Energy minimization, MD relaxation
└────────┬────────┘
         ▼
┌─────────────────┐
│ Quality check    │  ← Ramachandran, QMEAN, ProSA
└─────────────────┘
```

## Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
