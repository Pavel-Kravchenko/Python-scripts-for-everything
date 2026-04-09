---
name: genetics-computational
description: Computational genetics — genetic code analysis, codon usage/CAI, restriction digests, genetic mapping from cross data, Ts/Tv ratios, ORF finding, and mutation spectra
---

# Computational Genetics

## When to Use
- Working with the genetic code (codon lookups, translation verification)
- Codon usage analysis and optimization (CAI, RSCU)
- Restriction enzyme analysis beyond basic cut-site mapping
- Genetic mapping from two-point or three-point cross data
- Mutation analysis (Ts/Tv ratios, mutation spectra, C>T bias)
- ORF finding and annotation across all six reading frames

## Quick Reference

### Standard Genetic Code — Amino Acid to Codons
| AA | Codons | AA | Codons |
|----|--------|----|--------|
| F (Phe) | TTT TTC | L (Leu) | TTA TTG CTT CTC CTA CTG |
| I (Ile) | ATT ATC ATA | M (Met) | ATG |
| V (Val) | GTT GTC GTA GTG | S (Ser) | TCT TCC TCA TCG AGT AGC |
| P (Pro) | CCT CCC CCA CCG | T (Thr) | ACT ACC ACA ACG |
| A (Ala) | GCT GCC GCA GCG | Y (Tyr) | TAT TAC |
| H (His) | CAT CAC | Q (Gln) | CAA CAG |
| N (Asn) | AAT AAC | K (Lys) | AAA AAG |
| D (Asp) | GAT GAC | E (Glu) | GAA GAG |
| C (Cys) | TGT TGC | W (Trp) | TGG |
| R (Arg) | CGT CGC CGA CGG AGA AGG | G (Gly) | GGT GGC GGA GGG |
| * (Stop) | TAA TAG TGA | | |

### Codon Degeneracy Classes
| Class | Meaning | Examples |
|-------|---------|---------|
| 0-fold (non-degenerate) | All changes alter AA | ATG (Met), TGG (Trp) |
| 2-fold | Transitions silent; transversions change AA | AAA/AAG (Lys), CAA/CAG (Gln) |
| 4-fold | All third-position changes silent | GCN (Ala), GGN (Gly), CCN (Pro), GTN (Val) |
| 6-fold | Mix of 2-fold + 4-fold synonymy | Leu (6 codons), Arg (6), Ser (6) |

### Common Restriction Enzymes — Cut Sites and Overhang Types
| Enzyme | Recognition Site | Overhang |
|--------|-----------------|---------|
| EcoRI | G^AATTC | 5′ 4 nt |
| BamHI | G^GATCC | 5′ 4 nt |
| HindIII | A^AGCTT | 5′ 4 nt |
| NotI | GC^GGCCGC | 5′ 4 nt |
| XhoI | C^TCGAG | 5′ 4 nt |
| NdeI | CA^TATG | 5′ 2 nt |
| SmaI | CCC^GGG | Blunt |
| EcoRV | GAT^ATC | Blunt |
| KpnI | GGTAC^C | 3′ 4 nt |
| SacI | GAGCT^C | 3′ 4 nt |

For basic cut-site finding see `restriction_map` in `python-collections-regex`. `virtual_digest` below extends it with fragment sizes.

### Recombination Frequency → Map Distance
| RF (%) | Kosambi cM | Notes |
|--------|-----------|-------|
| ≤10 | ≈ RF × 100 | Near-linear; interference negligible |
| 20 | 22.2 | Correction begins to matter |
| 30 | 35.3 | Significant correction |
| 50 | ∞ | Independent assortment (unlinked) |

`d (Morgans) = 0.25 × ln[(1 + 2r) / (1 − 2r)]`, where r = RF as a fraction.

### Mutation Type Classification — Ts/Tv Matrix
| From \ To | A | G | C | T |
|-----------|---|---|---|---|
| A | — | Ts | Tv | Tv |
| G | Ts | — | Tv | Tv |
| C | Tv | Tv | — | Ts |
| T | Tv | Tv | Ts | — |

Transitions (Ts): A↔G, C↔T. Transversions (Tv): purine↔pyrimidine.
Expected Ts/Tv ~0.5 under neutral model; observed ~2–3 in coding regions (CpG hypermutation).

## Key Patterns

### Codon Extraction and Reading Frames
```python
# Full CDS — partial codons at end discarded
def extract_codons(cds: str) -> list[str]:
    return [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]

# All 3 forward frames
def all_frames(seq: str) -> list[list[str]]:
    return [[seq[i:i+3] for i in range(f, len(seq) - 2, 3)] for f in range(3)]
```

### RSCU and CAI Concepts
- RSCU(codon) = observed_count / (expected if uniform among synonymous codons); 1.0 = no bias
- CAI = geometric mean of w(codon) = RSCU_i / max_RSCU_in_synonymous_family, over all codons
- Reference set must be highly expressed genes (ribosomal proteins, housekeeping); never mix all genes

### Three-Point Cross: Gene Order Determination
1. Sort all 8 classes by count — two most frequent = parentals, two least = DCO.
2. Which gene flipped in the DCO class relative to the parental arrangement is the middle gene.
3. RF(interval) = (single-CO in that interval + DCO) / total.
4. Interference = 1 − CoC, where CoC = observed DCO / (RF1 × RF2 × total).

### ORF Finding Strategy
- Search all 6 frames: 3 forward on the sequence, 3 reverse on the reverse complement.
- ORF = ATG to first in-frame stop codon (TAA, TAG, TGA); filter by minimum length (default 100 nt).
- Nested ORFs share a stop codon — report the longest.

## Code Templates

### `codon_usage_table(seq)` — codon frequencies and RSCU
```python
from collections import Counter, defaultdict

GENETIC_CODE = {  # see full table in python-collections-regex
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def codon_usage_table(cds: str) -> dict[str, dict]:
    codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    counts = Counter(c for c in codons if GENETIC_CODE.get(c, '*') != '*')

    aa_to_codons = defaultdict(list)
    for codon, aa in GENETIC_CODE.items():
        if aa != '*':
            aa_to_codons[aa].append(codon)

    total = sum(counts.values())
    result = {}
    for codon, count in counts.items():
        aa = GENETIC_CODE[codon]
        aa_total = sum(counts[s] for s in aa_to_codons[aa])
        expected = aa_total / len(aa_to_codons[aa])
        result[codon] = {
            'aa': aa, 'count': count,
            'freq_per_1000': count / total * 1000,
            'rscu': count / expected if expected > 0 else 0.0,
        }
    return result
```

### `codon_adaptation_index(seq, ref_table)` — CAI
```python
import math

def codon_adaptation_index(cds: str, ref_table: dict[str, dict]) -> float:
    """CAI = geometric mean of w(codon); ref_table from codon_usage_table on high-expression genes."""
    aa_max = defaultdict(float)
    for codon, info in ref_table.items():
        aa_max[info['aa']] = max(aa_max[info['aa']], info['rscu'])

    codons = [cds[i:i+3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    log_sum, n = 0.0, 0
    for codon in codons:
        aa = GENETIC_CODE.get(codon)
        if not aa or aa == '*':
            continue
        rscu = ref_table.get(codon, {}).get('rscu', 0.0)
        w = rscu / aa_max[aa] if aa_max[aa] > 0 else 0.0
        if w > 0:
            log_sum += math.log(w); n += 1
    return math.exp(log_sum / n) if n > 0 else 0.0
```

CAI range 0–1; values ≥ 0.8 indicate high codon optimization. Codons with w = 0 are excluded from the geometric mean.

### `virtual_digest(seq, enzymes)` — fragments and overhang types
```python
import re

# Extends restriction_map from python-collections-regex with fragment sizes.
ENZYME_CUT_INFO = {  # (cut_offset, overhang_len, overhang_type)
    'EcoRI':(1,4,"5'"), 'BamHI':(1,4,"5'"), 'HindIII':(1,4,"5'"),
    'NotI':(2,4,"5'"),  'XhoI':(1,4,"5'"),  'NdeI':(2,2,"5'"),
}

def virtual_digest(seq: str, enzymes: dict[str, str]) -> dict:
    results = {}
    for name, site in enzymes.items():
        cuts = sorted(m.start() for m in re.finditer(site, seq))
        if not cuts:
            continue
        boundaries = [0] + cuts + [len(seq)]
        fragments = sorted(
            [boundaries[i+1] - boundaries[i] for i in range(len(boundaries)-1)],
            reverse=True
        )
        results[name] = {
            'cut_positions': cuts, 'n_cuts': len(cuts),
            'fragment_sizes': fragments,
            'overhang_type': ENZYME_CUT_INFO.get(name, (None,None,'unknown'))[2],
        }
    return results
```

### `three_point_cross(data)` — gene order and map distances
```python
def three_point_cross(data: dict[str, int]) -> dict:
    """data: genotype string -> count. Parental = most frequent, DCO = least frequent."""
    import math
    classes = sorted(data.items(), key=lambda x: x[1], reverse=True)
    total = sum(data.values())
    dco_count  = classes[-1][1] + classes[-2][1]
    sco1_count = classes[2][1]  + classes[3][1]
    sco2_count = classes[4][1]  + classes[5][1]
    rf1 = (sco1_count + dco_count) / total
    rf2 = (sco2_count + dco_count) / total
    kosambi = lambda r: 25 * math.log((1 + 2*r) / (1 - 2*r))  # cM
    expected_dco = rf1 * rf2 * total
    coc = dco_count / expected_dco if expected_dco > 0 else float('nan')
    return {
        'parental': [classes[0][0], classes[1][0]],
        'dco': [classes[-1][0], classes[-2][0]],
        'rf1_pct': round(rf1*100, 2), 'rf2_pct': round(rf2*100, 2),
        'map_d1_cM': round(kosambi(min(rf1, 0.499)), 2),
        'map_d2_cM': round(kosambi(min(rf2, 0.499)), 2),
        'coc': round(coc, 3), 'interference': round(1-coc, 3),
    }
```

### `ts_tv_ratio(seq1, seq2)` — Ts/Tv from aligned sequences
```python
def ts_tv_ratio(seq1: str, seq2: str) -> tuple[int, int, float]:
    """seq1 and seq2 must be equal-length aligned strings (gaps as '-')."""
    TRANSITIONS = {('A','G'),('G','A'),('C','T'),('T','C')}
    ts = tv = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a == b or '-' in (a, b) or 'N' in (a, b):
            continue
        if (a, b) in TRANSITIONS:
            ts += 1
        else:
            tv += 1
    return ts, tv, round(ts / tv, 4) if tv > 0 else float('inf')
```

### `find_orfs(seq, min_length=100)` — all 6 reading frames
```python
import re

RC_TABLE = str.maketrans("ATGC", "TACG")

def find_orfs(seq: str, min_length: int = 100) -> list[dict]:
    seq = seq.upper()
    rc  = seq.translate(RC_TABLE)[::-1]
    orfs = []
    for strand, template in (('+', seq), ('-', rc)):
        for frame in range(3):
            for m in re.finditer(r'ATG', template[frame:]):
                start = m.start() + frame
                region = template[start:]
                stop = next(
                    (i+3 for i in range(0, len(region)-2, 3)
                     if region[i:i+3] in ('TAA','TAG','TGA')), None
                )
                if stop and stop >= min_length:
                    orfs.append({'start': start, 'end': start+stop,
                                 'strand': strand, 'frame': frame+1,
                                 'length': stop, 'seq': region[:stop]})
    return sorted(orfs, key=lambda x: x['length'], reverse=True)
```

### `mutation_spectrum(variants)` — classify C>T, C>A, etc.
```python
from collections import Counter

def mutation_spectrum(variants: list[tuple[str, str]]) -> dict[str, int]:
    """Normalise to pyrimidine reference; variants = list of (ref, alt) single bases."""
    COMP = {'A':'T','T':'A','G':'C','C':'G'}
    counts: Counter = Counter()
    for ref, alt in variants:
        ref, alt = ref.upper(), alt.upper()
        if len(ref) != 1 or len(alt) != 1 or ref == alt or ref not in 'ACGT':
            continue
        if ref in ('A', 'G'):      # normalise to pyrimidine ref
            ref, alt = COMP[ref], COMP[alt]
        counts[f'{ref}>{alt}'] += 1
    return dict(counts)
```

Six canonical types: C>A, C>G, C>T, T>A, T>C, T>G. C>T at CpG is the most frequent somatic mutation in human cancer.

### `gc_at_codon_positions(seq)` — GC at 1st, 2nd, 3rd positions
```python
def gc_at_codon_positions(cds: str) -> dict[str, float]:
    pos = {1:[], 2:[], 3:[]}
    for i in range(0, len(cds)-2, 3):
        for j, b in enumerate(cds[i:i+3].upper(), 1):
            pos[j].append(b)
    gc = lambda bases: sum(1 for b in bases if b in 'GC') / len(bases) if bases else 0.0
    return {f'GC{k}': round(gc(v), 4) for k, v in pos.items()}
```

GC3 is highly variable across organisms; used as a proxy for mutational bias. GC1/GC2 are constrained by amino-acid composition.

## Common Pitfalls

| Pitfall | Fix |
|---------|-----|
| Translating partial codons | Truncate CDS to `len - len % 3` before splitting |
| Including stop codons in RSCU | Filter out TAA, TAG, TGA before counting |
| Using genomic sequence for codon usage | Use CDS only (exons in-frame from ATG) |
| Forgetting reverse complement in ORF search | Always search both strands; `find_orfs` handles this |
| Overlapping restriction sites | `re.finditer` finds all non-overlapping; add `(?=...)` lookahead for overlapping |
| Observed DCO > expected (CoC > 1) | Negative interference — rare; verify class assignments |
| RF ≈ 0.5 produces infinite Kosambi distance | Cap r at 0.499 before applying mapping function |
| HWE testing with small n | `hwe_test` in `ngs-variant-calling`; n < 50 yields unreliable chi-squared p-values |
| Codons absent from ref set give w = 0 | These are silently excluded from CAI geometric mean — document it |
| Mis-assigning parental vs DCO classes | Always sort all 8 class frequencies before role assignment |

## Related Skills
- `python-core-bio` — DNA/RNA string manipulation, reverse complement
- `python-collections-regex` — `restriction_map`, `RESTRICTION_ENZYMES`, `GENETIC_CODE` dict
- `ngs-variant-calling` — `hwe_test` base implementation; VCF parsing for mutation spectrum input
- `biopython-databases` — fetching CDS sequences from NCBI/Ensembl for codon usage
- `sequence-alignment` — aligned sequences required by `ts_tv_ratio`
- `phylogenetics-evolution` — dN/dS builds on Ts/Tv ratios and codon degeneracy classes
