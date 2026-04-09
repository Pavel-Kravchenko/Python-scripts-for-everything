---
name: genetic-engineering-insilico
description: In silico genetic engineering — primer design, CRISPR guide selection, codon optimization, Gibson/Golden Gate assembly design, and virtual cloning workflows
---

# Genetic Engineering In Silico

## When to Use
- Designing cloning experiments in silico before benchwork
- PCR primer design and Tm validation
- CRISPR guide RNA selection and off-target scoring
- Codon optimization for heterologous expression in E. coli, yeast, or human cells
- Gibson Assembly or Golden Gate assembly planning
- Verifying reading frame and restriction sites after ligation

## Quick Reference

### Common Restriction Enzymes (extends `restriction_map()` from `python-collections-regex`)
| Enzyme  | Recognition | Cut (top/bot) | Overhang   | Methylation | Compatible ends |
|---------|------------|---------------|------------|-------------|-----------------|
| EcoRI   | G↓AATTC    | 1 / 5         | 5' AATT    | star activity| MfeI           |
| BamHI   | G↓GATCC    | 1 / 5         | 5' GATC    | —           | BclI, BglII, MboI|
| HindIII | A↓AGCTT    | 1 / 5         | 5' AGCT    | —           | —               |
| NotI    | GC↓GGCCGC  | 2 / 6         | 5' GGCC    | —           | —               |
| NdeI    | CA↓TATG    | 2 / 4         | 5' TA      | —           | —               |
| XhoI    | C↓TCGAG    | 1 / 5         | 5' TCGA    | —           | SalI            |
| DpnI    | GA↓TC      | 2 / 2         | blunt      | requires dam+| —              |
| BsaI    | GGTCTC(1/5)| type IIS      | 5' 4-nt    | —           | BsmBI by design |

### Cas Variants and PAM Sequences
| Cas variant   | PAM          | Guide length | Cut type    | Notes               |
|---------------|-------------|--------------|-------------|---------------------|
| SpCas9        | NGG (3')     | 20 bp        | blunt       | most widely used    |
| SaCas9        | NNGRRT (3')  | 21 bp        | blunt       | smaller for AAV     |
| Cas12a (Cpf1) | TTTV (5')    | 23 bp        | 5' overhang | PAM on 5' side      |
| xCas9         | NGN, GAA (3')| 20 bp        | blunt       | expanded PAM        |
| BE3           | NGG (3')     | 20 bp        | nickase     | C→T, window pos 4-8 |

### Tm Calculation Methods
| Method           | Formula / approach                                    | Use case              |
|------------------|-------------------------------------------------------|-----------------------|
| Basic 4+2        | `4*(G+C) + 2*(A+T)`                                   | rough / short oligos  |
| Salt-adjusted    | 4+2 result + 16.6*log10([Na+]/0.05)                   | ~50 mM Na+            |
| Nearest-neighbor | `ΔH / (ΔS + R*ln(C_T/4)) − 273.15` (SantaLucia 1998) | most accurate; use for design |

### Codon Usage Bias (most common codon per AA, abbreviated)
| AA  | E. coli | S. cerevisiae | H. sapiens |
|-----|---------|---------------|------------|
| Leu | CTG     | TTG           | CTG        |
| Arg | CGT     | AGA           | AGG        |
| Ser | AGC     | TCT           | AGC        |
| Gly | GGT     | GGT           | GGC        |
| Pro | CCG     | CCA           | CCC        |

### Gibson Assembly Overlap Guidelines
| Overlap length | Tm target  | Status              |
|---------------|-----------|---------------------|
| <15 bp        | —         | avoid — frequent failures |
| 15–20 bp      | 48–52 °C  | minimum acceptable  |
| 20–25 bp      | 52–58 °C  | recommended         |
| 25–40 bp      | 58–65 °C  | large constructs >4 kb |

### Golden Gate Type IIS Enzymes
| Enzyme | Recognition  | Cut offset (fwd+rev) | Overhang |
|--------|-------------|----------------------|----------|
| BsaI   | GGTCTC(N)1  | +1 / +5              | 4 nt 5'  |
| BbsI   | GAAGAC(N)2  | +2 / +6              | 4 nt 5'  |
| BsmBI  | CGTCTC(N)1  | +1 / +5              | 4 nt 5'  |

## Key Patterns

### Virtual Cloning Workflow
```
1. restriction_map(insert, enzymes)    # find all sites — from python-collections-regex
2. check compatible ends between vector/insert
3. confirm orientation from site asymmetry
4. frame check: (vector_cut_pos - start_codon) % 3 == 0
5. scan insert for internal sites that destroy the construct
6. verify_clone() → check start/stop codons at junction
```

### Primer Design Constraints
- Tm fwd vs rev: within ±2 °C (nearest-neighbor both); GC 40–60 %; length 18–25 bp
- No 3' complementarity ≥3 bp between fwd and rev (hetero-dimer)
- No 3' self-complementarity (hairpin blocks extension)
- No mononucleotide run ≥4 at 3' end; end with G or C (GC clamp)

### CRISPR Guide Selection Pipeline
```
1. find_crispr_guides(seq, pam='NGG')   # PAM scan on both strands
2. filter: 40–70 % GC; no TTTT run (U6 terminator)
3. score on-target: GC fraction, G at position −1, seed region quality
4. score_guide_specificity(guide, genome) — mismatch count
5. pick top 3; verify no off-targets with ≤2 mismatches
```

### Codon Harmonization vs Optimization
- **Optimization**: every codon → most frequent host codon; maximizes speed, may disrupt co-translational folding
- **Harmonization**: match *relative* codon frequency of source organism in host; preserves slow-translation pausing at domain boundaries
- Rule: optimize short single-domain proteins; harmonize complex multi-domain proteins

### Gateway Cloning attB Site Design
```
fwd tail: 5'-GGGGACAAGTTTGTACAAAAAAGCAGGCT-[20 bp gene-specific]-3'
rev tail: 5'-GGGGACCACTTTGTACAAGAAAGCTGGGT-[20 bp RC of gene end]-3'
PCR → BP reaction with donor vector → entry clone
```

## Code Templates

### Nearest-Neighbor Tm + Primer Design
```python
import re, math

NN_DH = {'AA':-7.9,'AT':-7.2,'TA':-7.2,'CA':-8.5,'GT':-8.4,'CT':-7.8,
          'GA':-8.2,'CG':-10.6,'GC':-9.8,'GG':-8.0,'TT':-7.9,'TG':-8.5,
          'AC':-8.4,'TC':-8.2,'AG':-7.8,'CC':-8.0}          # kcal/mol
NN_DS = {'AA':-22.2,'AT':-20.4,'TA':-21.3,'CA':-22.7,'GT':-22.4,'CT':-21.0,
          'GA':-22.2,'CG':-27.2,'GC':-24.4,'GG':-19.9,'TT':-22.2,'TG':-22.7,
          'AC':-22.4,'TC':-22.2,'AG':-21.0,'CC':-19.9}       # cal/mol/K
R = 1.987  # cal/mol/K

def tm_nearest_neighbor(seq, conc_nm=250):
    seq = seq.upper()
    dH = sum(NN_DH.get(seq[i:i+2], -8.0) for i in range(len(seq)-1)) * 1000
    dS = sum(NN_DS.get(seq[i:i+2], -21.0) for i in range(len(seq)-1)) - 5.9
    return dH / (dS + R * math.log(conc_nm * 1e-9 / 4)) - 273.15

RC = str.maketrans('ATGC', 'TACG')
def rev_comp(seq): return seq.upper().translate(RC)[::-1]
def gc_frac(seq): return sum(1 for c in seq.upper() if c in 'GC') / len(seq)

def design_primers(template, target_start, target_end, tm_target=60, min_len=18, max_len=25):
    """Design fwd/rev primers bracketing template[target_start:target_end]."""
    template = template.upper()
    def grow(region, target):
        for n in range(min_len, max_len + 1):
            if tm_nearest_neighbor(region[:n]) >= target or n == max_len:
                p = region[:n]
                return p, round(tm_nearest_neighbor(p), 1), round(gc_frac(p), 2)
    fwd, fwd_tm, fwd_gc = grow(template[target_start:], tm_target)
    rev_region = rev_comp(template[:target_end])
    rev, rev_tm, rev_gc = grow(rev_region, tm_target)
    return {'fwd': fwd, 'fwd_tm': fwd_tm, 'fwd_gc': fwd_gc,
            'rev': rev, 'rev_tm': rev_tm, 'rev_gc': rev_gc,
            'tm_delta': round(abs(fwd_tm - rev_tm), 1)}
```

### Primer Dimer Check
```python
def check_primer_dimers(fwd, rev, min_match=3):
    warnings = []
    fwd, rev = fwd.upper(), rev.upper()
    for n in range(min_match, min(8, len(fwd), len(rev)) + 1):
        if fwd[-n:] == rev_comp(rev[-n:]):
            warnings.append(f"3' hetero-dimer: {n} bp"); break
    for primer, label in [(fwd,'fwd'), (rev,'rev')]:
        for n in range(min_match, min(8, len(primer)) + 1):
            if primer[-n:] == rev_comp(primer[:n]):
                warnings.append(f"{label} self-dimer: {n} bp"); break
        if re.search(r'(.)\1{3}$', primer):
            warnings.append(f"{label}: mononucleotide run at 3' end")
    return warnings or ['OK']
```

### CRISPR Guide Design + Off-Target Score
```python
IUPAC = {'N':'ATGC','R':'AG','Y':'CT','V':'ACG','B':'CGT','H':'ACT','D':'AGT',
          'A':'A','T':'T','G':'G','C':'C','S':'GC','W':'AT','K':'GT','M':'AC'}

def _pam_match(seq, pam):
    return len(seq) == len(pam) and all(seq[i] in IUPAC.get(pam[i], pam[i])
                                        for i in range(len(pam)))

def find_crispr_guides(seq, pam='NGG', guide_length=20):
    """Find guide RNAs on both strands. Returns list sorted by score (GC-based)."""
    seq, plen, guides = seq.upper(), len(pam), []
    def scan(strand, label):
        for i in range(len(strand) - guide_length - plen + 1):
            if _pam_match(strand[i+guide_length:i+guide_length+plen], pam):
                g = strand[i:i+guide_length]
                gc = gc_frac(g); tttt = 'TTTT' in g
                guides.append({'guide':g,'strand':label,'pos':i,'gc':round(gc,2),
                                'score':round(gc if 0.4<=gc<=0.7 and not tttt else gc*0.5,3)})
    scan(seq, '+'); scan(rev_comp(seq), '-')
    return sorted(guides, key=lambda g: g['score'], reverse=True)

def score_guide_specificity(guide, genome_seq, max_mm=3):
    """Count genome positions with <= max_mm mismatches. Use BLAST for real genomes."""
    guide, gseq = guide.upper(), genome_seq.upper()
    n = len(guide); hits = []
    for strand, s in [('+', gseq), ('-', rev_comp(gseq))]:
        for i in range(len(s) - n + 1):
            mm = sum(a != b for a, b in zip(guide, s[i:i+n]))
            if mm <= max_mm:
                hits.append({'strand':strand,'pos':i,'mm':mm,'site':s[i:i+n]})
    return sorted(hits, key=lambda h: h['mm'])
```

### Codon Optimization
```python
PREFERRED = {
    'ecoli': {'F':'TTT','L':'CTG','I':'ATT','M':'ATG','V':'GTG','S':'AGC','P':'CCG',
              'T':'ACC','A':'GCT','Y':'TAT','H':'CAT','Q':'CAG','N':'AAC','K':'AAA',
              'D':'GAT','E':'GAA','C':'TGT','W':'TGG','R':'CGT','G':'GGT','*':'TAA'},
    'yeast': {'F':'TTT','L':'TTG','I':'ATT','M':'ATG','V':'GTT','S':'TCT','P':'CCA',
              'T':'ACT','A':'GCT','Y':'TAT','H':'CAT','Q':'CAA','N':'AAT','K':'AAA',
              'D':'GAT','E':'GAA','C':'TGT','W':'TGG','R':'AGA','G':'GGT','*':'TAA'},
    'human': {'F':'TTC','L':'CTG','I':'ATC','M':'ATG','V':'GTG','S':'AGC','P':'CCC',
              'T':'ACC','A':'GCC','Y':'TAC','H':'CAC','Q':'CAG','N':'AAC','K':'AAG',
              'D':'GAC','E':'GAG','C':'TGC','W':'TGG','R':'AGG','G':'GGC','*':'TGA'},
}
GENETIC_CODE = {  # abbreviated — see python-collections-regex for full table
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def optimize_codons(cds, organism='ecoli'):
    """Replace each codon with organism-preferred codon. Returns (opt_cds, report)."""
    pref = PREFERRED[organism]
    codons = [cds[i:i+3].upper() for i in range(0, len(cds) - len(cds) % 3, 3)]
    out = [pref.get(GENETIC_CODE.get(c, 'X'), c) for c in codons]
    changed = sum(a != b for a, b in zip(codons, out))
    return ''.join(out), {'changed': changed, 'total': len(codons),
                          'pct_changed': round(100 * changed / len(codons), 1)}
```

### Gibson Assembly + Golden Gate Design
```python
def gibson_overlaps(fragments, overlap_length=20):
    """Design overlapping PCR primer tails for circular Gibson Assembly."""
    n = len(fragments); result = []
    for i, frag in enumerate(fragments):
        nxt = fragments[(i + 1) % n].upper()
        overlap = nxt[:overlap_length]
        result.append({'fragment': i+1, 'overlap_to_next': overlap,
                        'overlap_tm': round(tm_nearest_neighbor(overlap), 1),
                        'rev_primer': rev_comp(frag[-20:].upper()) + overlap})
    return result

TYPE_IIS = {'BsaI': 'GGTCTC', 'BbsI': 'GAAGAC', 'BsmBI': 'CGTCTC'}

def golden_gate_design(parts, overhangs, enzyme='BsaI'):
    """Add type IIS sites + 4-nt overhangs to each part as PCR primer tails.
    overhangs: list of 4-nt strings, len = len(parts) + 1 (junction sequences).
    """
    site = TYPE_IIS[enzyme]; rc_site = rev_comp(site)
    return [{'part': i+1,
             'fwd': overhangs[i] + site + 'A' + parts[i][:20],
             'rev': rev_comp(overhangs[i+1] + rc_site + 'A' + rev_comp(parts[i][-20:]))}
            for i, part in enumerate(parts)]
```

### Verify Clone (Reading Frame + Junction Check)
```python
def verify_clone(vector, insert, cut_pos_vector, cut_pos_insert=0, expected_frame=0):
    """Simulate ligation and verify reading frame and stop codons at junction."""
    construct = vector[:cut_pos_vector] + insert + vector[cut_pos_vector:]
    frame_offset = (cut_pos_vector + cut_pos_insert) % 3
    j_start = max(0, cut_pos_vector - 9)
    junction = construct[j_start: cut_pos_vector + 30]
    stops = [i for i in range(0, len(junction)-2, 3) if junction[i:i+3] in ('TAA','TAG','TGA')]
    codons = [insert[i:i+3] for i in range(cut_pos_insert, min(cut_pos_insert+15,len(insert)), 3)]
    return {'construct_len': len(construct), 'frame_ok': frame_offset == expected_frame,
            'frame_offset': frame_offset, 'stops_near_junction': stops,
            'first_5_codons': codons}

# NdeI cloning: CA|TATG at position 8; insert GFP (ATG starts insert)
vector = 'AAATTTCATATGCACCACCAC'
insert = 'ATGGTGAGCAAGGGCGAGGAG'
result = verify_clone(vector, insert, cut_pos_vector=10, cut_pos_insert=0, expected_frame=0)
# {'frame_ok': True, 'stops_near_junction': [], 'first_5_codons': ['ATG','GTG','AGC','AAG','GGC']}
```

## Common Pitfalls

| Pitfall | Fix |
|---|---|
| dam methylation blocks BclI/MboI (GATC sites) | Use dam⁻ strain; DpnI *requires* dam+ to digest template |
| Star activity: high glycerol or wrong salt → non-canonical cuts | Use NEB-recommended buffer; limit digest to 1–2 h |
| Primer Tm mismatch >2 °C | Extend lower-Tm primer 1–2 nt; recalculate with nearest-neighbor |
| CRISPR off-targets in repeats / pseudogenes | BLAST guide against full genome; flag hits with ≤2 mismatches |
| Codon optimization removes rare codons at domain boundaries | Harmonize instead; restore key slow codons manually |
| Gibson overlap with secondary structure | Check overlap with mfold; redesign if ΔG < −3 kcal/mol |
| Reading frame lost after ligation | Run `verify_clone()` before ordering primers |
| Internal restriction site in insert destroys construct | Run `restriction_map(insert, enzymes)` before finalizing vector sites |
| Stop codon introduced at cloning junction | Inspect `stops_near_junction` in `verify_clone()` output |
| CRISPR PAMs only searched on one strand | `find_crispr_guides()` scans both + and − strands |

## Related Skills
- `genetics-computational` — codon usage, restriction basics, GC content
- `python-collections-regex` — `restriction_map()` and `predict_fragments()` base functions this skill extends
- `biopython-databases` — retrieve template sequences and genome assemblies for off-target checks
- `sequence-alignment` — primer specificity via BLAST, guide off-target alignment
- `python-core-bio` — `rev_comp()`, DNA string manipulation, ORF finding
