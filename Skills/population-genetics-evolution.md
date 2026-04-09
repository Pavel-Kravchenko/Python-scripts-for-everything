---
name: population-genetics-evolution
description: Population-level processes and molecular evolution — genetic drift, selection, dN/dS, Tajima's D, Fst, LD, coalescent, and molecular clock
---

# Population Genetics & Molecular Evolution

## When to Use
- Simulating genetic drift and selection in finite populations
- Computing dN/dS (Ka/Ks) ratios to detect selection on coding sequences
- Estimating divergence times from sequence data with a molecular clock
- Testing for departures from neutrality (Tajima's D, McDonald-Kreitman)
- Computing population differentiation (Fst) from allele frequency data
- Measuring and modeling linkage disequilibrium (LD) decay
- Coalescent simulations of gene genealogies

## Quick Reference

### Wright-Fisher Model
| Parameter | Symbol | Interpretation |
|-----------|--------|----------------|
| Effective population size | Ne | Idealised size with same drift as real pop |
| Drift variance per gen | p(1−p) / 2Ne | Variance of allele freq change |
| Fixation prob (neutral) | p0 | Same as initial frequency |
| Fixation prob (selected) | (1 − e^(−2s)) / (1 − e^(−4Nes)) | Kimura formula |
| Expected fixation time | ≈ 4Ne generations | Neutral allele starting at 1/2N |

**Ne adjustments**: Ne < N due to bottlenecks, unequal sex ratio (Ne = 4·Nm·Nf / (Nm+Nf)), or variance in reproductive success. **Selection vs drift**: effective when |2Ne·s| >> 1; drift dominates when |2Ne·s| << 1.

### Selection Models
| Model | Fitness (AA, Aa, aa) | Equilibrium |
|-------|---------------------|-------------|
| Directional (additive) | 1+s, 1+s/2, 1 | Fixation of A; h = 0.5 |
| Directional (recessive) | 1+s, 1, 1 | Slow fixation; h = 0 |
| Balancing (overdominance) | 1, 1+s, 1 | Stable at p* = s2/(s1+s2) |
| Purifying | 1, 1−hs, 1−s | Deleterious allele removed |

### dN/dS (ω) Interpretation
| ω | Interpretation | Biological meaning |
|---|---------------|--------------------|
| < 1 | Purifying selection | Amino acid changes are deleterious |
| ≈ 1 | Neutral evolution | Amino acid changes tolerated |
| > 1 | Positive selection | Amino acid changes confer advantage |

**Caveats**: gene-wide ω < 1 can mask positive selection at a subset of sites; use site models or sliding windows. Avoid dN/dS when divergence > 1 substitution/site (saturation).

### Neutrality Tests
| Test | Null hypothesis | Positive value | Negative value |
|------|----------------|----------------|----------------|
| Tajima's D | Neutrality + constant size | Balancing sel / bottleneck | Purifying sel / expansion |
| Fu & Li's D | Same | Same | Same |
| McDonald-Kreitman | pN/pS = dN/dS | Excess divergence → adaptive | Excess polymorphism → deleterious |

### Fst Interpretation (Wright's scale)
| Fst | Differentiation |
|-----|----------------|
| 0 – 0.05 | Little / negligible |
| 0.05 – 0.15 | Moderate |
| 0.15 – 0.25 | Great |
| > 0.25 | Very great |

Fst = (HT − HS) / HT; HT = total heterozygosity, HS = mean within-subpopulation heterozygosity.

### Molecular Clock
```
divergence_time = jc_distance / (2 × substitution_rate)
```
Typical vertebrate mtDNA: ~2% per Myr pairwise. Nuclear rates ~10× slower. Requires fossil/biogeographic calibration.

### Linkage Disequilibrium Measures
| Measure | Formula | Range | Use |
|---------|---------|-------|-----|
| D | pAB − pA·pB | [−0.25, 0.25] | Raw covariance; scale-dependent |
| D' | D / Dmax | [−1, 1] | Normalised; |D'|=1 means no recombination observed |
| r² | D² / (pA·qA·pB·qB) | [0, 1] | Statistical correlation; best for GWAS power |

r² decays with distance: E[r²] ≈ 1 / (1 + 4Ne·c·d), c = recombination rate/bp, d = distance in bp.

## Key Patterns

### Wright-Fisher Drift
Each generation: allele count ~ Binomial(2N, p_current). Fixation (p=1) and loss (p=0) are absorbing states. Small Ne → fast fixation/loss; large Ne → slow drift. Many replicates reveal the distribution of outcomes.

### Nei-Gojobori dN/dS Method
1. Count synonymous (S) and non-synonymous (N) sites per codon by averaging over all single-step mutation paths
2. Count observed Sd (synonymous) and Nd (non-synonymous) differences
3. Correct for multiple hits: pS = Sd/S, pN = Nd/N, apply JC: dS = −0.75·ln(1 − 4pS/3)
4. ω = dN/dS

### Coalescent
Going backwards in time, two lineages coalesce with probability 1/(2Ne) per generation. E[TMRCA for k lineages] = 4Ne / k(k−1) generations. Gene trees can differ from species trees (incomplete lineage sorting, ILS) when Ne is large relative to the speciation interval.

### LD Decay
r² decays approximately as 1/(1 + C·d). Population-level recombination rate C = 4Ne·r. LD extends further in species with small Ne (domesticates, island populations) or after recent bottlenecks.

## Code Templates

### Wright-Fisher Simulation
```python
import numpy as np

def wright_fisher_sim(N: int, p0: float, generations: int,
                      n_replicates: int = 100) -> np.ndarray:
    """Returns (n_replicates, generations+1) allele frequency array."""
    freqs = np.zeros((n_replicates, generations + 1))
    freqs[:, 0] = p0
    for gen in range(1, generations + 1):
        freqs[:, gen] = np.random.binomial(2*N, freqs[:, gen-1]) / (2*N)
    return freqs

trajs = wright_fisher_sim(N=100, p0=0.2, generations=200, n_replicates=50)
print(f"Fixed: {(trajs[:,-1]==1).sum()}, Lost: {(trajs[:,-1]==0).sum()}")
```

### Selection Trajectory
```python
def selection_trajectory(N: int, p0: float, s: float,
                          h: float = 0.5, generations: int = 500) -> np.ndarray:
    """Allele frequency under selection + drift. h=0.5: additive, h=1: dominant."""
    freqs, p = np.zeros(generations + 1), p0
    freqs[0] = p0
    for gen in range(1, generations + 1):
        w_bar = p**2*(1+s) + 2*p*(1-p)*(1+h*s) + (1-p)**2
        p_sel = (p**2*(1+s) + p*(1-p)*(1+h*s)) / w_bar
        p = np.random.binomial(2*N, p_sel) / (2*N)
        freqs[gen] = p
        if p in (0, 1):
            freqs[gen:] = p; break
    return freqs
```

### Compute dN/dS (Nei-Gojobori)
```python
GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def _codon_sites(codon: str) -> tuple[float, float]:
    aa, syn, nonsyn = GENETIC_CODE.get(codon, '*'), 0.0, 0.0
    for i, base in enumerate(codon):
        for alt in 'ACGT':
            if alt == base: continue
            mut = codon[:i] + alt + codon[i+1:]
            mut_aa = GENETIC_CODE.get(mut, '*')
            if mut_aa == '*': continue
            if mut_aa == aa: syn += 1/3
            else: nonsyn += 1/3
    return syn, nonsyn

def compute_dnds(seq1: str, seq2: str) -> tuple[float, float, float]:
    """Returns (dN, dS, omega). omega=nan if dS=0."""
    S, N, Sd, Nd = 0.0, 0.0, 0.0, 0.0
    for i in range(0, len(seq1), 3):
        c1, c2 = seq1[i:i+3], seq2[i:i+3]
        if '-' in c1+c2: continue
        s1, n1 = _codon_sites(c1); s2, n2 = _codon_sites(c2)
        S += (s1+s2)/2; N += (n1+n2)/2
        if c1 != c2:
            if GENETIC_CODE.get(c1,'*') == GENETIC_CODE.get(c2,'*'): Sd += 1
            else: Nd += 1
    pS, pN = (Sd/S if S>0 else 0), (Nd/N if N>0 else 0)
    dS = -0.75*np.log(1-4*pS/3) if pS < 0.75 else float('inf')
    dN = -0.75*np.log(1-4*pN/3) if pN < 0.75 else float('inf')
    return dN, dS, dN/dS if dS > 0 else float('nan')

def sliding_window_dnds(seq1: str, seq2: str,
                         window: int = 30, step: int = 3) -> list[dict]:
    """dN/dS in overlapping codon windows (window/step in codons)."""
    n = len(seq1) // 3
    return [{'start': s, 'end': s+window,
              **dict(zip(('dN','dS','omega'), compute_dnds(seq1[s*3:(s+window)*3],
                                                            seq2[s*3:(s+window)*3])))}
            for s in range(0, n-window+1, step)]
```

### Tajima's D
```python
def tajimas_d(allele_freqs: list[float], n: int) -> float:
    """allele_freqs: derived freq per SNP; n: number of haploid sequences."""
    S = len(allele_freqs)
    if S == 0: return float('nan')
    pi = sum(2*f*(1-f)*n/(n-1) for f in allele_freqs)
    a1 = sum(1/i for i in range(1, n))
    a2 = sum(1/i**2 for i in range(1, n))
    theta_w = S / a1
    b1 = (n+1)/(3*(n-1)); b2 = 2*(n**2+n+3)/(9*n*(n-1))
    c1 = b1 - 1/a1; c2 = b2 - (n+2)/(a1*n) + a2/a1**2
    var_d = (c1/a1)*S + (c2/(a1**2+a2))*S*(S-1)
    return (pi - theta_w) / var_d**0.5

# 10 SNPs, 30 sequences
D = tajimas_d([0.3,0.7,0.5,0.4,0.6,0.2,0.8,0.45,0.55,0.35], n=30)
```

### Compute Fst
```python
def compute_fst(pop_freqs: list[list[float]]) -> float:
    """Mean Fst across biallelic loci. pop_freqs[pop][locus] = allele frequency."""
    vals = []
    for locus in range(len(pop_freqs[0])):
        ps = [pop[locus] for pop in pop_freqs]
        p_bar = sum(ps) / len(ps)
        if p_bar in (0, 1): continue
        H_T = 2*p_bar*(1-p_bar)
        H_S = sum(2*p*(1-p) for p in ps) / len(ps)
        vals.append((H_T - H_S) / H_T)
    return sum(vals)/len(vals) if vals else float('nan')

# 3 populations diverged across 5 SNPs
fst = compute_fst([[0.9,0.8,0.7,0.6,0.9],[0.5,0.4,0.5,0.5,0.5],[0.1,0.2,0.3,0.4,0.1]])
```

### Molecular Clock Divergence
```python
def molecular_clock_divergence(seq1: str, seq2: str,
                                rate: float) -> dict:
    """rate: substitutions/site/year. Returns p-distance, JC distance, divergence time."""
    p = sum(a!=b for a,b in zip(seq1,seq2) if '-' not in (a,b)) / len(seq1)
    if p >= 0.75: raise ValueError(f"Saturated sequences (p={p:.3f})")
    d = -0.75 * np.log(1 - 4*p/3)
    return {'p_distance': p, 'jc_distance': d, 'divergence_years': d/(2*rate)}
```

### LD Decay
```python
def ld_decay(genotypes: np.ndarray, positions: np.ndarray) -> tuple[list, list]:
    """Returns (distances_bp, r2_values) for all SNP pairs."""
    n = genotypes.shape[1]
    dists, r2s = [], []
    for i in range(n):
        for j in range(i+1, n):
            x, y = genotypes[:,i] - genotypes[:,i].mean(), genotypes[:,j] - genotypes[:,j].mean()
            denom = (x**2).sum() * (y**2).sum()
            if denom == 0: continue
            r2s.append((x*y).sum()**2 / denom)
            dists.append(abs(int(positions[j]) - int(positions[i])))
    return dists, r2s
```

### McDonald-Kreitman Test
```python
from scipy.stats import fisher_exact

def mcdonald_kreitman_test(Pn: int, Ps: int, Dn: int, Ds: int) -> dict:
    """
    Pn/Ps: within-species non-syn/syn polymorphisms.
    Dn/Ds: between-species non-syn/syn fixed differences.
    NI < 1 and significant p → positive selection.
    """
    _, p_value = fisher_exact([[Pn, Ps], [Dn, Ds]])
    ni = (Pn/Ps)/(Dn/Ds) if Ps>0 and Dn>0 and Ds>0 else float('nan')
    return {'NI': ni, 'alpha': 1-ni, 'p_value': p_value}

# Drosophila Adh (classic example): NI < 1 → adaptive divergence
result = mcdonald_kreitman_test(Pn=4, Ps=19, Dn=7, Ds=17)
```

## Common Pitfalls

- **Drift mimics selection**: large allele freq swings arise by chance alone in small Ne; verify 2Ne·|s| >> 1 before claiming a selective sweep
- **dN/dS averaging hides positive selection**: gene-wide ω < 1 can mask adaptive sites; use PAML site models (M8) or sliding windows
- **Tajima's D confounded by demography**: expansion → negative D (mimics purifying selection); bottleneck → positive D (mimics balancing selection); always compare to a demographic null
- **Fst inflated by rare variants and marker type**: microsatellites give higher Fst than SNPs; use Jost's D or G_ST for highly polymorphic markers
- **Molecular clock violations**: generation-time and metabolic-rate effects cause among-lineage rate heterogeneity; use relaxed-clock Bayesian methods (BEAST) when life histories differ substantially
- **dN/dS on saturated sequences**: synonymous sites saturate first; dS > 2 is unreliable; apply only to sequences with p-distance < 0.5
- **Effective vs census size**: Ne << N under bottlenecks, skewed sex ratios, or high reproductive variance; demographic events leave signatures for ~4Ne generations
- **Spurious LD from population structure**: admixture between populations with different allele frequencies creates LD between unlinked loci; check for stratification before interpreting LD patterns

## Related Skills
- `phylogenetics-evolution` — tree building, bootstrap, comparative genomics (this skill adds within/between-population processes)
- `ngs-variant-calling` — HWE test and VCF variant data feeding into population analyses
- `genetics-computational` — genetic code tables and mutation types underlying dN/dS
- `biostatistics-r` / `probability-statistics-python` — Fisher exact, chi-square, and regression for neutrality and LD tests
