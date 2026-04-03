---
name: sequence-alignment
description: Pairwise alignment (Needleman-Wunsch, Smith-Waterman), BLOSUM/PAM matrices, BLAST searching, and multiple sequence alignment
---

# Sequence Alignment & BLAST

## When to Use
- Pairwise alignment: comparing two sequences for homology, finding shared domains
- Global (NW): both sequences expected to align end-to-end (full orthologs, similar length)
- Local (SW): finding best sub-region match (domain search, short read vs reference)
- BLAST: database search (too slow to run SW against all sequences)
- MSA: 3+ sequences to find conserved positions, build phylogenies, design primers

---

## Quick Reference

### Algorithm Comparison

| Feature | Needleman-Wunsch | Smith-Waterman |
|---|---|---|
| Scope | Global (end-to-end) | Local (sub-region) |
| Init | `F[i,0] = gap*i`, `F[0,j] = gap*j` | All zeros |
| Recurrence | max(diag, up, left) | max(**0**, diag, up, left) |
| Traceback start | Bottom-right | Cell with max score |
| Traceback end | Top-left (0,0) | First cell with score 0 |
| Complexity | O(mn) time and space | O(mn) time and space |

### DP Recurrences

**Needleman-Wunsch (global):**
```
F(i,j) = max{
  F(i-1,j-1) + s(xi, yj)   # match/mismatch
  F(i-1,j) - d              # gap in seq1
  F(i,j-1) - d              # gap in seq2
}
```

**Smith-Waterman (local)** — same but with zero floor:
```
F(i,j) = max{
  0,
  F(i-1,j-1) + s(xi, yj),
  F(i-1,j) - d,
  F(i,j-1) - d
}
```

### Gap Penalty Models

**Linear:** `cost(k) = d * k` — every gap position costs the same.

**Affine:** `cost(k) = d + (k-1) * e` where `d` = open penalty, `e` = extend penalty (`e < d`).
- Affine is biologically realistic: one mutational event creates a long gap, so extension should be cheap.
- Typical: `d = 10, e = 0.5` (BLAST default for proteins: open=-11, extend=-1).

### Substitution Matrix Selection

| Matrix | Best for | Notes |
|---|---|---|
| **BLOSUM62** | Standard protein alignment (~62% identity blocks) | BLAST default |
| **BLOSUM80** | Closely related sequences (>80% identity) | High specificity |
| **BLOSUM45** | Distantly related sequences (<45% identity) | More permissive |
| **PAM40** | Closely related proteins (~85% identity) | Legacy |
| **PAM250** | Distantly related proteins | Legacy, less common now |

**PAM vs BLOSUM:** PAM numbers increase with evolutionary distance; BLOSUM numbers increase with sequence similarity (inverse relationship). BLOSUM is derived directly from conserved blocks; prefer BLOSUM for most work.

Formula: `S(a,b) = log2(q_ab / (p_a * p_b))` — positive = conservative substitution, negative = disruptive.

### Identity Thresholds (Proteins)

| %Identity | Interpretation |
|---|---|
| >30% | Almost certainly homologous |
| 20–30% | Twilight zone — may or may not be homologous |
| <20% | Midnight zone — cannot confirm by identity alone |

---

## BLAST

### Five Variants

| Program | Query | Database | Translation | Use |
|---|---|---|---|---|
| **blastn** | Nucleotide | Nucleotide | None | Gene finding, same/close species |
| **blastp** | Protein | Protein | None | Protein homology |
| **blastx** | Nucleotide | Protein | Query (6 frames) | Annotate ESTs/mRNA |
| **tblastn** | Protein | Nucleotide | DB (6 frames) | Find unannotated genes |
| **tblastx** | Nucleotide | Nucleotide | Both (6 frames) | Compare unannotated genomes (slow) |

**Nucleotide sub-variants (blastn):**
- `megablast` (W=28): >95% identity, same species
- `blastn` (W=11): cross-species
- `discontiguous megablast`: ~80% identity, sensitive

**PSI-BLAST:** Iterative PSSM-based protein search for remote homologs.
1. Standard blastp → collect hits → build PSSM → re-search with PSSM → repeat until convergence.
- Use when blastp finds few/no hits; protein belongs to large diverse family.
- Risk: profile corruption from false-positive hits — use inclusion E-value ≤ 0.001.

### Seed-and-Extend Mechanism

1. Break query into overlapping words of length W
2. Generate neighborhood words scoring ≥ T against substitution matrix
3. Look up words in pre-indexed database → High-scoring Segment Pairs (HSPs)
4. Ungapped extension; if score above threshold → gapped alignment (SW-like)
5. Compute E-value; report significant hits

### E-value and Bit Score

**E-value:** `E = K * m * n * e^(-λS)` = expected random hits at this score in database of size n.
- `m` = query length, `n` = database length, `S` = raw score, `K/λ` = Karlin-Altschul params
- E-value is **not** a probability; it is database-size dependent.

**Bit score:** `S' = (λS - ln K) / ln 2` — normalized, database-independent.
- Bit score 50 means 2^50× more likely than random.

| E-value | Interpretation |
|---|---|
| < 1e-50 | Virtually identical |
| 1e-50 to 1e-10 | Clear homolog, same family |
| 1e-10 to 1e-5 | Probable distant homolog |
| 1e-5 to 0.01 | Possible homology — verify further |
| > 1 | Likely not significant |

### Parameter Tuning

| Parameter | Default | When to change |
|---|---|---|
| E-value | 10 | Lower to 0.001 for stringency; higher for remote homologs |
| Word size W | 3 (prot) / 11 (DNA) | Decrease for sensitivity; increase for speed |
| Matrix | BLOSUM62 | BLOSUM45 distant; BLOSUM80 close |
| Gap open/extend | -11/-1 | Experiment if unreasonable gaps appear |
| Low-complexity filter (SEG/DUST) | ON | Only disable if query has genuine low-complexity regions |
| Composition-based stats | ON | Reduces bias from amino acid composition |

### Database Selection

| Database | Quality | Use |
|---|---|---|
| **swissprot** | High (curated) | Start here for proteins — fast, reliable annotations |
| **nr** | Mixed | Comprehensive when swissprot misses hits |
| **refseq_protein** | High | Reference annotations |
| **pdb** | High | Structural homologs |
| **nt** | Mixed | Nucleotide searches |

---

## Multiple Sequence Alignment (MSA)

### Why Not Exact DP?
Exact MSA for s sequences of length n requires O(n^s) time — NP-hard. For 3 sequences, each cell has 7 predecessors (2^s - 1 in general). Practical tools use heuristics.

### Progressive Alignment Pipeline
1. All-vs-all pairwise distances (k-mer or alignment-based)
2. Build guide tree (UPGMA or NJ)
3. Align most similar pair first, then progressively merge along the tree
4. Align-to-alignment uses **profile columns** (frequency vectors): score = average of all pairwise substitution scores across column residues
5. Optional iterative refinement

**Limitation:** Errors introduced early propagate — guide tree quality matters.

### MSA Tools

| Tool | Algorithm | Speed | Accuracy | Best For |
|---|---|---|---|---|
| **ClustalW** | Progressive (NJ) | Slow | Moderate | Small datasets, legacy |
| **Clustal Omega** | HMM profile-profile + mBed | Fast | Good | Large datasets (>2000 seqs) |
| **MUSCLE** | Progressive + iterative | Fast | Good | Medium datasets |
| **MAFFT** | FFT-based + iterative | Very fast | Very good | Large/diverse datasets |
| **T-Coffee** | Consistency-based | Slow | Excellent | Small, high-accuracy needs |
| **PROBCONS** | Probabilistic consistency | Slow | Excellent | Benchmark quality |

**Rule:** Default to `mafft --auto`; use T-Coffee/PROBCONS for small critical alignments.

### Conservation Scoring
- **Conservation fraction:** `max_freq_at_position / n_seqs`
- **Information content (bits):** `IC(pos) = log2(alphabet_size) - entropy(pos)` where entropy = -Σ f_i * log2(f_i)
  - Max IC = log2(20) ≈ 4.32 bits (protein), log2(4) = 2 bits (DNA)
  - IC ≈ 0: fully variable; IC = max: fully conserved
- **Sum-of-pairs score:** sum of all pairwise substitution scores per column — overall MSA quality metric

---

## Code Templates

### Pairwise Alignment with BioPython

```python
from Bio.Align import PairwiseAligner, substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

aligner = PairwiseAligner()
aligner.mode = "global"          # or "local"
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

alignments = aligner.align(seq1, seq2)
best = alignments[0]
print(f"Score: {best.score}")
print(best)
```

### Load Substitution Matrix

```python
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")
score = blosum62["L", "I"]  # lookup specific pair
```

### BLAST via BioPython (remote)

```python
from Bio.Blast import NCBIWWW, NCBIXML

result_handle = NCBIWWW.qblast(
    program="blastp",
    database="swissprot",
    sequence=query_seq,
    expect=0.001,
    hitlist_size=20,
)
blast_record = NCBIXML.read(result_handle)

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        pct_id = 100 * hsp.identities / hsp.align_length
        print(f"{alignment.title[:60]}: E={hsp.expect:.1e}, id={pct_id:.1f}%")
```

### Local BLAST (command-line)

```bash
makeblastdb -in proteins.fasta -dbtype prot -out mydb -parse_seqids
blastp -query query.fasta -db mydb -out results.txt \
       -outfmt "6 qseqid sseqid pident evalue bitscore stitle" \
       -evalue 0.001 -num_threads 4
```

```python
from Bio.Blast.Applications import NcbiblastpCommandline
cmd = NcbiblastpCommandline(query="q.fasta", db="mydb",
                             evalue=0.001, outfmt=5, out="out.xml", num_threads=4)
cmd()
```

### Parse Tabular BLAST Output (outfmt 6)

```python
import pandas as pd
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv("results.txt", sep="\t", names=cols)
df_filtered = df[(df.evalue <= 1e-5) & (df.pident >= 25) & (df.bitscore >= 50)]
```

### BioPython AlignIO

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read alignment
aln = AlignIO.read("alignment.fasta", "fasta")
print(f"{len(aln)} seqs, {aln.get_alignment_length()} columns")

# Write in different format
AlignIO.write(aln, "alignment.clustal", "clustal")

# Access columns
col_5 = aln[:, 5]         # column 5 as string
subregion = aln[:, 10:20] # columns 10-19

# Build programmatically
records = [SeqRecord(Seq("MVLS-PADKTNVK"), id="seq1"),
           SeqRecord(Seq("MVLSGEDKS-IVK"), id="seq2")]
aln = MultipleSeqAlignment(records)
```

### Conservation Analysis

```python
from collections import Counter
import math

def conservation_at_pos(column):
    """column: list of residue chars (gaps allowed)"""
    counts = Counter(c for c in column if c != '-')
    if not counts:
        return 0.0
    return max(counts.values()) / len(column)

def information_content(column, alphabet_size=20):
    counts = Counter(c for c in column if c != '-')
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = -sum((c/total) * math.log2(c/total) for c in counts.values())
    return math.log2(alphabet_size) - entropy

# Apply across alignment
aln_strings = [str(r.seq) for r in aln]
n_cols = len(aln_strings[0])
cons_scores = [conservation_at_pos([s[i] for s in aln_strings]) for i in range(n_cols)]
ic_scores   = [information_content([s[i] for s in aln_strings]) for i in range(n_cols)]
```

### Run MSA Tool (subprocess)

```python
import subprocess, shutil

def run_mafft(input_fasta, output_fasta):
    if not shutil.which("mafft"):
        raise RuntimeError("mafft not in PATH — install: conda install -c bioconda mafft")
    result = subprocess.run(["mafft", "--auto", input_fasta],
                            capture_output=True, text=True, timeout=120)
    with open(output_fasta, "w") as f:
        f.write(result.stdout)
```

---

## Common Pitfalls

- **Wrong BLAST program:** Using blastn when the DB is protein; use blastx to translate a nucleotide query.
- **Default E-value = 10 is permissive by design:** Many returned hits will be non-significant — always filter.
- **E-value without database context:** Same alignment has a worse E-value in a larger database; always report both.
- **Percent identity alone is insufficient:** 90% identity over 10 residues is meaningless; coverage matters too.
- **Low-complexity filtering off:** Poly-X tails, proline-rich regions → mass false positives.
- **PSI-BLAST profile corruption:** One false positive in early iteration poisons the PSSM — use strict inclusion threshold (e.g., 0.001).
- **Progressive alignment error propagation:** Errors from the first pair persist — guide tree quality matters; consider iterative refinement.
- **Not trimming poorly aligned MSA regions:** Downstream phylogenetics and conservation analyses are distorted by noisy ends/insertions.
- **NCBI rate limits:** Max one remote BLAST request per 10 seconds; use local BLAST for high throughput.

---

## Related Skills
- `phylogenetics-evolution` — building trees from MSA; UPGMA, NJ, maximum likelihood
- `biopython-databases` — SeqRecord, SeqIO, NCBI Entrez, UniProt, PDB access
- `string-algorithms` — exact/approximate string matching underlying seed-and-extend
