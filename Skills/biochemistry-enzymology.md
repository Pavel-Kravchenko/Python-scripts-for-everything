---
name: biochemistry-enzymology
description: Enzyme kinetics fitting (Michaelis-Menten, Hill, inhibition models), EC number classification, stoichiometric matrices, and spectrophotometric assay analysis
---

# Biochemistry & Enzymology

## When to Use
- Fitting enzyme kinetics data (Michaelis-Menten, Hill equation) with confidence intervals
- Determining enzyme inhibition type from experimental velocity data
- Looking up or classifying enzymes by EC number
- Building stoichiometric matrices for basic metabolic pathway computation
- Protein thermostability analysis from activity-vs-temperature curves
- Processing enzyme assay data (progress curves, initial velocity extraction)

## Quick Reference

### Michaelis-Menten Equation and Linearizations
| Form | Equation | Plot axes | Slope / Intercept |
|------|----------|-----------|-------------------|
| Michaelis-Menten | `v = Vmax*[S] / (Km + [S])` | [S] vs v (hyperbola) | — |
| Lineweaver-Burk | `1/v = (Km/Vmax)*(1/[S]) + 1/Vmax` | 1/[S] vs 1/v | slope=Km/Vmax, y-int=1/Vmax |
| Eadie-Hofstee | `v = Vmax - Km*(v/[S])` | v/[S] vs v | slope=−Km, y-int=Vmax |
| Hanes-Woolf | `[S]/v = [S]/Vmax + Km/Vmax` | [S] vs [S]/v | slope=1/Vmax, x-int=−Km |

### Enzyme Inhibition Types
| Type | App. Km | App. Vmax | Lineweaver-Burk | Ki relation |
|------|---------|-----------|-----------------|-------------|
| Competitive | Km*(1+[I]/Ki) | unchanged | Intersect on y-axis | Ki=[I]/(alpha−1) |
| Uncompetitive | Km/(1+[I]/Ki') | Vmax/(1+[I]/Ki') | Parallel lines | Ki'=[I]/(alpha'−1) |
| Non-competitive (pure) | unchanged | Vmax/(1+[I]/Ki) | Intersect on x-axis | Ki=[I]/(Vmax/Vmax_app−1) |
| Mixed | Km*alpha/alpha' | Vmax/alpha' | Intersect left of y-axis | fit Ki and Ki' jointly |

### EC Number Classification
| EC Class | Name | Reaction | Example |
|----------|------|----------|---------|
| EC 1 | Oxidoreductases | oxidation/reduction | LDH (EC 1.1.1.27) |
| EC 2 | Transferases | group transfer | Hexokinase (EC 2.7.1.1) |
| EC 3 | Hydrolases | hydrolysis | Trypsin (EC 3.4.21.4) |
| EC 4 | Lyases | add/remove at double bonds | Aldolase (EC 4.1.2.13) |
| EC 5 | Isomerases | isomerization | PGI (EC 5.3.1.9) |
| EC 6 | Ligases | bond formation + ATP | Pyruvate carboxylase (EC 6.4.1.1) |

### Hill Coefficient Interpretation
| n | Cooperativity | Example |
|---|---------------|---------|
| n < 1 | Negative | Myoglobin subunits (theoretical) |
| n = 1 | None (hyperbolic) | Michaelis-Menten equivalent |
| 1 < n < 4 | Positive | Hemoglobin O2 (n ≈ 2.8) |

Hill equation: `v = Vmax * [S]^n / (K_half^n + [S]^n)`.  K_half = [S] at half-maximal response. n is empirical; it does not equal the number of binding sites.

### Common Enzyme Assay Types
| Type | Principle | Readout |
|------|-----------|---------|
| Continuous spectrophotometric | Monitor A at fixed λ over time | ΔA/min → velocity (e.g. NADH at 340 nm) |
| Coupled | Auxiliary enzyme links product to chromogenic reaction | ΔA/min |
| Fluorogenic | Fluorescent product released | ΔF/min (e.g. protease + AMC substrate) |
| Radiometric | Radiolabeled substrate; isolate product | DPM after TLC or filter binding |

### Key Conversions
| Quantity | Conversion |
|----------|------------|
| 1 U | 1 µmol substrate/min at defined conditions |
| 1 katal (kat) | 1 mol s⁻¹ = 6.0×10⁷ U |
| Specific activity | U / mg protein (increases with purity) |
| kcat | Vmax / [E]total (s⁻¹); requires active-site titration |

## Key Patterns

### Nonlinear Fitting (preferred over linearizations)
```python
from scipy.optimize import curve_fit
import numpy as np

def mm(S, Vmax, Km):
    return Vmax * S / (Km + S)

popt, pcov = curve_fit(mm, S, v, p0=[np.max(v), np.median(S)])
Vmax, Km = popt
Vmax_err, Km_err = np.sqrt(np.diag(pcov)) * 1.96  # 95% CI
```

### Inhibition Model Selection via AIC
```python
def aic(n_params, residuals):
    n = len(residuals)
    sse = np.sum(residuals**2)
    return n * np.log(sse / n) + 2 * n_params
# Fit each inhibition model, compute residuals, pick lowest AIC
```

### Beer-Lambert Law
`A = ε * c * l`  →  `c = A / (ε * l)`
NADH/NADPH: ε₃₄₀ = 6220 M⁻¹ cm⁻¹.  ΔA/min × (1/ε) × (1/l) × V_assay_L × 10⁶ = µmol/min.

## Code Templates

### `michaelis_menten_fit(substrate_conc, velocities)`
```python
import numpy as np
from scipy.optimize import curve_fit

def michaelis_menten_fit(substrate_conc, velocities):
    S, v = np.array(substrate_conc, float), np.array(velocities, float)
    def mm(S, Vmax, Km): return Vmax * S / (Km + S)
    popt, pcov = curve_fit(mm, S, v, p0=[np.max(v), np.median(S)], maxfev=5000)
    perr = np.sqrt(np.diag(pcov)) * 1.96
    Vmax, Km = popt
    return {"Vmax": Vmax, "Vmax_95ci": perr[0], "Km": Km, "Km_95ci": perr[1]}

# Hexokinase with glucose substrate
S = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]   # mM
v = [12.1, 21.3, 35.6, 58.9, 74.2, 85.0, 91.3]  # nmol/min/mg
res = michaelis_menten_fit(S, v)
print(f"Km = {res['Km']:.3f} ± {res['Km_95ci']:.3f} mM")
```

### `lineweaver_burk_plot(S, V)`
```python
import numpy as np, matplotlib.pyplot as plt
from scipy.stats import linregress

def lineweaver_burk_plot(S, V):
    inv_S, inv_V = 1.0/np.array(S, float), 1.0/np.array(V, float)
    slope, intercept, r, *_ = linregress(inv_S, inv_V)
    Vmax, Km = 1.0/intercept, slope/intercept
    plt.scatter(inv_S, inv_V)
    x = np.linspace(min(inv_S)*1.1, max(inv_S), 100)
    plt.plot(x, slope*x + intercept)
    plt.xlabel("1/[S] (mM⁻¹)"); plt.ylabel("1/v"); plt.show()
    return {"Km": Km, "Vmax": Vmax, "R2": r**2}
```

### `inhibition_analysis(S, V_ctrl, V_inh, I)`
```python
import numpy as np
from scipy.optimize import curve_fit

def inhibition_analysis(S, V_ctrl, V_inh, I):
    """Fit competitive/uncompetitive/noncompetitive; return best model by AIC."""
    S, V_ctrl, V_inh = [np.array(x, float) for x in (S, V_ctrl, V_inh)]

    def mm(S, Vmax, Km): return Vmax * S / (Km + S)
    (Vmax0, Km0), _ = curve_fit(mm, S, V_ctrl, p0=[max(V_ctrl), np.median(S)])

    models = {
        "competitive":    lambda S, Vmax, Km, Ki: Vmax*S/(Km*(1+I/Ki)+S),
        "uncompetitive":  lambda S, Vmax, Km, Ki: (Vmax/(1+I/Ki))*S/(Km/(1+I/Ki)+S),
        "noncompetitive": lambda S, Vmax, Km, Ki: (Vmax/(1+I/Ki))*S/(Km+S),
    }
    results = {}
    for name, model in models.items():
        try:
            popt, _ = curve_fit(model, S, V_inh, p0=[Vmax0, Km0, 1.0], maxfev=5000)
            res = V_inh - model(S, *popt)
            n = len(res)
            results[name] = {"params": popt, "AIC": n*np.log(np.sum(res**2)/n)+6}
        except RuntimeError:
            results[name] = {"params": None, "AIC": np.inf}
    results["best"] = min(results, key=lambda k: results[k]["AIC"])
    return results

# LDH inhibition by oxamate (competitive inhibitor)
S_mM  = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
v_ctrl = np.array([18.2, 31.0, 55.4, 72.1, 83.6, 91.0])
v_inh  = np.array([10.1, 18.5, 38.2, 55.6, 70.3, 82.1])
print(inhibition_analysis(S_mM, v_ctrl, v_inh, I=0.5)["best"])
```

### `hill_equation_fit(ligand_conc, response)`
```python
def hill_equation_fit(ligand_conc, response):
    L, Y = np.array(ligand_conc, float), np.array(response, float)
    def hill(L, Ymax, K_half, n): return Ymax * L**n / (K_half**n + L**n)
    popt, pcov = curve_fit(hill, L, Y, p0=[max(Y), np.median(L), 1.0], maxfev=5000)
    perr = np.sqrt(np.diag(pcov))
    Ymax, K_half, n = popt
    return {"Ymax": Ymax, "K_half": K_half, "n_hill": n, "n_err": perr[2]}

# Hemoglobin O2 binding
pO2 = [1, 2, 5, 10, 20, 40, 80, 100]    # mmHg
sat = [5, 10, 28, 52, 75, 90, 97, 98]    # % saturation
fit = hill_equation_fit(pO2, sat)
print(f"Hill n = {fit['n_hill']:.2f}, P50 = {fit['K_half']:.1f} mmHg")
```

### `ec_number_lookup(ec_number)` — KEGG REST
```python
import urllib.request

def ec_number_lookup(ec_number: str) -> dict[str, str]:
    """Query KEGG for enzyme entry. ec_number e.g. '1.1.1.27'."""
    url = f"https://rest.kegg.jp/get/ec:{ec_number}"
    with urllib.request.urlopen(url) as r:
        raw = r.read().decode()
    entry, current = {}, None
    for line in raw.splitlines():
        if line and line[0] != " ":
            current = line.split()[0]
            entry[current] = line[len(current):].strip()
        elif current:
            entry[current] += " " + line.strip()
    return entry

info = ec_number_lookup("1.1.1.27")  # LDH
print(info.get("NAME"), info.get("REACTION"))
```

### `stoichiometric_matrix(reactions)`
```python
import re

def stoichiometric_matrix(reactions: list[str]):
    """
    Build S matrix (compounds x reactions) from strings like '2 A + B -> C + D'.
    Substrates: negative coefficients. Products: positive.
    """
    def parse_side(text):
        terms = {}
        for token in re.split(r'\s*\+\s*', text.strip()):
            m = re.match(r'^(\d+)\s+(\S+)$', token.strip())
            coeff, name = (int(m.group(1)), m.group(2)) if m else (1, token.strip())
            terms[name] = coeff
        return terms

    compounds, rows = [], []
    for rxn in reactions:
        lhs, rhs = rxn.split("->")
        row = {n: -c for n, c in parse_side(lhs).items()}
        for n, c in parse_side(rhs).items():
            row[n] = row.get(n, 0) + c
        for n in row:
            if n not in compounds: compounds.append(n)
        rows.append(row)

    S = [[row.get(c, 0) for row in rows] for c in compounds]
    return S, compounds

rxns = ["Glucose + ATP -> G6P + ADP", "G6P -> F6P", "F6P + ATP -> FBP + ADP"]
S_mat, cpds = stoichiometric_matrix(rxns)
for c, row in zip(cpds, S_mat): print(f"{c:10s}: {row}")
```

### `beer_lambert(absorbance, epsilon, path_length)`
```python
def beer_lambert(absorbance: float, epsilon: float, path_length: float = 1.0) -> float:
    """Return molar concentration [M] from A = ε*c*l."""
    return absorbance / (epsilon * path_length)

# NADH at 340 nm (ε = 6220 M⁻¹ cm⁻¹)
c_M = beer_lambert(0.412, 6220)          # → 6.63e-5 M = 66.3 µM
delta_A_per_min = 0.045
v_umol_per_min = beer_lambert(delta_A_per_min, 6220) * 1e-3 * 1e6  # 1 mL assay
```

### `progress_curve_analysis(time, product)`
```python
from scipy.stats import linregress

def progress_curve_analysis(time, product, linear_fraction: float = 0.1):
    """Extract initial velocity from the linear phase of a progress curve."""
    t, p = np.array(time, float), np.array(product, float)
    cut = max(2, int(len(t) * linear_fraction))
    slope, _, r, _, se = linregress(t[:cut], p[:cut])
    return {"v0": slope, "R2": r**2, "cutoff_time": t[cut-1], "se": se}

# LDH: NADH produced over time (µM)
t_s    = [0, 10, 20, 30, 60, 90, 120, 180, 240]
nadh   = [0, 4.1, 8.3, 12.4, 22.0, 28.5, 32.1, 35.0, 36.2]
result = progress_curve_analysis(t_s, nadh, linear_fraction=0.15)
print(f"v0 = {result['v0']:.3f} µM/s  R² = {result['R2']:.4f}")
```

## Common Pitfalls
- **Substrate inhibition**: at very high [S] velocity drops; Michaelis-Menten underestimates Vmax. Use `v = Vmax*S/(Km + S*(1 + S/Ki_si))`.
- **Lineweaver-Burk distortion**: low-[S] points have high leverage and compounded error. Use only for visualization; always determine Km/Vmax by nonlinear fitting.
- **Ki ambiguity**: Ki (competitive) and Ki' (uncompetitive) are distinct constants; do not interchange them across inhibition model equations.
- **Hill n ≠ binding sites**: n is empirical. Hemoglobin has 4 O2 sites but n ≈ 2.8, not 4.
- **Non-linear progress curve region**: fitting a line to the full curve underestimates v0. Use only the first 5–15% of substrate consumed.
- **Product inhibition**: accumulating product slows the reaction independently of substrate depletion — visible as a faster-than-expected curve.
- **Wrong extinction coefficient**: ε is wavelength- and buffer-dependent. Always verify ε at your exact assay wavelength and conditions.
- **Unreported conditions**: Km and Vmax are condition-specific. Always report temperature, pH, ionic strength, and co-factor concentrations.
- **Poor p0 initialization**: `curve_fit` fails without good starting guesses. Use `p0=[max(v), np.median(S)]` for Michaelis-Menten.

## Related Skills
- `data-visualization-bio` — plotting kinetics curves (hyperbola, Lineweaver-Burk, residuals)
- `numpy-pandas-wrangling` — data manipulation for kinetics datasets and multi-condition comparisons
- `structural-bioinformatics` — active site analysis and PDB structure inspection
- `clinical-modeling-workflows` — molecular docking and drug binding builds on enzyme inhibition understanding
