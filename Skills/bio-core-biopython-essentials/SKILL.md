---
name: bio-core-biopython-essentials
description: BioPython essentials — Seq/SeqRecord objects, SeqIO file I/O, Entrez API, pairwise alignment with PairwiseAligner.
tool_type: python
primary_tool: Python
---

# BioPython Essentials

## Pitfalls

- **Seq vs str**: `Seq` behaves like a string (slicing, `len()`, iteration) but adds `.complement()`, `.reverse_complement()`, `.transcribe()`, `.translate()`. In BioPython >=1.79, `Seq` and `str` interoperate freely. Pre-1.79 required explicit `str()` conversion.
- **`translate()` stop codons**: `dna.translate()` represents stops as `*`. Use `dna.translate(to_stop=True)` to stop at the first stop and omit `*` — always use this for CDS protein extraction.
- **`SeqIO.parse()` vs `SeqIO.read()`**: `parse()` returns an *iterator* — always safe for multi-record files. `read()` expects exactly one record and raises an error otherwise.
- **Entrez rate limits**: Always call `handle.close()`. NCBI allows 3 req/s without an API key. Set `Entrez.api_key = "yourkey"` for batch work (10 req/s).
- **Genetic code table**: Mitochondrial code (table=2) uses TGA as Trp, not stop. Bacterial code (table=11) differs from standard (table=1). Always specify `table=` for non-standard organisms.

## Seq Object API

```python
from Bio.Seq import Seq, MutableSeq
from Bio.SeqUtils import gc_fraction, molecular_weight

dna = Seq("ATGCGATCGATCGTAA")
dna.complement()           # 3'→5' complement
dna.reverse_complement()   # 5'→3' reverse complement
dna.transcribe()           # DNA → RNA (T→U)
dna.translate()            # all codons, stops as *
dna.translate(to_stop=True)  # stop at first stop, omit *
dna.translate(table=11)    # bacterial genetic code

gc_fraction(dna)           # 0.0–1.0
molecular_weight(dna)      # Da
molecular_weight(protein, seq_type='protein')

# Point mutation
mut = MutableSeq("ATGCGATCG")
mut[3] = 'T'   # G→T at position 3
```

## SeqRecord

```python
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

record = SeqRecord(
    Seq("ATGCGATCGATCGTAA"),
    id="BRCA1_001", name="BRCA1",
    description="BRCA1 partial CDS"
)
record.annotations["organism"] = "Homo sapiens"

# Add feature
cds = SeqFeature(FeatureLocation(0, 15), type="CDS",
                 qualifiers={"gene": ["BRCA1"]})
record.features.append(cds)

# Extract feature sequence and translate
cds_seq = cds.location.extract(record.seq)
protein = cds_seq.translate(to_stop=True)

# Per-letter annotations (e.g. FASTQ quality)
record.letter_annotations["phred_quality"] = [30, 30, 28, 35, ...]
```

## SeqIO

```python
from Bio import SeqIO

# Read multi-record file (always safe)
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(record.id, len(record))

# Read exactly one record
record = SeqIO.read("single.gb", "genbank")

# Load all into dict keyed by ID
records_dict = SeqIO.to_dict(SeqIO.parse("seqs.fasta", "fasta"))

# Write records
SeqIO.write(records, "output.fasta", "fasta")

# Format conversion (one-liner)
SeqIO.convert("reads.fastq", "fastq", "reads.fasta", "fasta")

# Quality filtering from FASTQ
good = [r for r in SeqIO.parse("reads.fastq", "fastq")
        if sum(r.letter_annotations["phred_quality"]) / len(r) >= 25]
SeqIO.write(good, "filtered.fastq", "fastq")
```

**Supported formats**: `fasta`, `fastq`, `genbank` (or `gb`), `embl`, `stockholm`, `clustal`, `phylip`

## Bio.Entrez

```python
from Bio import Entrez
Entrez.email = "you@example.com"
# Entrez.api_key = "yourkey"  # 10 req/s vs 3 req/s

# Search
handle = Entrez.esearch(db="nucleotide", term="BRCA1[gene] AND Homo sapiens[organism]")
record = Entrez.read(handle); handle.close()
ids = record["IdList"]

# Fetch sequences
handle = Entrez.efetch(db="nucleotide", id=ids[:5], rettype="fasta", retmode="text")
for rec in SeqIO.parse(handle, "fasta"):
    print(rec.id, len(rec))
handle.close()

# Fetch GenBank with features
handle = Entrez.efetch(db="nucleotide", id="NM_007294", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank"); handle.close()

# Cross-database links (nucleotide → protein)
handle = Entrez.elink(dbfrom="nucleotide", db="protein", id=ids[0])
link_record = Entrez.read(handle); handle.close()
```

## Pairwise Alignment

```python
from Bio.Align import PairwiseAligner, substitution_matrices

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

alignments = aligner.align(seq1, seq2)
best = alignments[0]
print(best.score)
print(best)   # formatted alignment

# Local alignment (Smith-Waterman)
aligner.mode = 'local'
```
