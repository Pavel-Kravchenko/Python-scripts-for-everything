---
name: biopython-databases
description: BioPython Seq/SeqRecord/SeqIO, NCBI Entrez API, UniProt queries, and biological format conversion
primary_tool: Python
---

# BioPython  Biological Databases

## When to Use
- Fetching sequences, annotations, or cross-references from NCBI, UniProt, PDB, or Ensembl
- Reading/writing/converting FASTA, GenBank, EMBL, FASTQ files
- Manipulating sequences (complement, translate, GC%, MW, ORF finding)
- Batch downloading gene/protein records programmatically

## Quick Reference

### NCBI Accession Prefixes
| Prefix | Type | DB |
|--------|------|----|
| NM_ | curated mRNA | RefSeq |
| NR_ | curated ncRNA | RefSeq |
| NP_ | curated protein | RefSeq |
| NC_ | complete chromosome | RefSeq |
| XM_/XP_ | predicted mRNA/protein | RefSeq |
| No prefix (e.g. AY123456) | submitted sequence | GenBank |
| SRP/SRS/SRX/SRR | study/sample/experiment/run | SRA |
| GSE/GSM/GPL | series/sample/platform | GEO |

### UniProt Accession Format
- SwissProt (reviewed): 6-char alphanumeric, e.g. `P01308`, `P68871`
- TrEMBL (unreviewed): same format, much larger volume
- Search filter: `reviewed:true` for SwissProt only

### SeqIO Format Strings
| Format | Read | Write | Notes |
|--------|------|-------|-------|
| `"fasta"` | yes | yes | no quality or annotation |
| `"genbank"` / `"gb"` | yes | yes | full annotations + features |
| `"fastq"` | yes | yes | per-base quality scores |
| `"embl"` | yes | yes | EBI equivalent of GenBank |

## Key Patterns

### Setup
```python
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO, Entrez, AlignIO, Align
from Bio.SeqUtils import gc_fraction, molecular_weight

Entrez.email = "your.email@example.com"  # NCBI requires this — set before any E-utility call
```

### Seq Object
```python
dna = Seq("ATGCGATCGATCGTAA")
dna.complement()           # 3'->5' complement
dna.reverse_complement()   # 5'->3' complement of the other strand
dna.transcribe()           # DNA -> RNA (T -> U)
rna.back_transcribe()      # RNA -> DNA
dna.translate()            # returns protein, * = stop codon
dna.translate(to_stop=True)  # stops before stop codon (most common)
dna.translate(table=2)     # mitochondrial genetic code

gc_fraction(dna) * 100     # GC percentage
molecular_weight(dna)      # DNA MW in Da
molecular_weight(protein, seq_type="protein")

# Mutable (in-place editing)
m = MutableSeq("ATGCGATCG")
m[3] = "T"  # point mutation
```

### SeqRecord Object
```python
record = SeqRecord(
    Seq("ATGCGATCG"),
    id="GENE_001",
    name="GENE",
    description="Description string"
)
record.annotations["organism"] = "Homo sapiens"
record.annotations["molecule_type"] = "DNA"  # required for GenBank output

# Add feature
feat = SeqFeature(FeatureLocation(0, 27), type="CDS",
                  qualifiers={"gene": ["BRCA1"], "product": ["BRCA1 protein"]})
record.features.append(feat)

# Per-letter annotations (e.g. FASTQ quality)
record.letter_annotations["phred_quality"] = [30, 28, 35, ...]

# Extract subsequence from feature
cds_seq = feature.location.extract(record.seq)
```

### SeqIO: Reading & Writing
```python
# Read multiple records (returns iterator — wrap in list() to reuse)
records = list(SeqIO.parse("file.fasta", "fasta"))
records = list(SeqIO.parse("file.gb", "genbank"))

# Read exactly one record
record = SeqIO.read("single.gb", "genbank")

# Write
SeqIO.write(records, "out.fasta", "fasta")

# Load into dict keyed by record.id
d = SeqIO.to_dict(SeqIO.parse("file.fasta", "fasta"))
d["INS_HUMAN"].seq

# One-liner format conversion (returns count)
count = SeqIO.convert("reads.fastq", "fastq", "reads.fasta", "fasta")
```

## Code Templates

### NCBI esearch + efetch (single record)
```python
handle = Entrez.esearch(
    db="nucleotide",
    term="insulin[Gene] AND Homo sapiens[Organism] AND mRNA[Filter] AND RefSeq[Filter]",
    retmax=10
)
results = Entrez.read(handle)
handle.close()
ids = results["IdList"]   # list of NCBI UIDs

handle = Entrez.efetch(db="nucleotide", id="NM_000207.3",
                       rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
```

### Batch fetch (comma-separated IDs)
```python
accessions = ["NM_000518.5", "NM_000207.3", "NM_000546.6"]
handle = Entrez.efetch(db="nucleotide", id=",".join(accessions),
                       rettype="fasta", retmode="text")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
SeqIO.write(records, "batch.fasta", "fasta")
```

### elink: cross-database navigation
```python
# Nucleotide - Protein
handle = Entrez.elink(dbfrom="nucleotide", db="protein", id="NM_000207.3")
link_results = Entrez.read(handle)
handle.close()
protein_ids = [link["Id"] for linkset in link_results
               for linkdb in linkset["LinkSetDb"]
               for link in linkdb["Link"]]

# Gene - PubMed (gene ID 7157  TP53)
handle = Entrez.elink(dbfrom="gene", db="pubmed", id="7157")
```

### Extract CDS and translate from GenBank record
```python
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3",
                       rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.location.extract(record.seq)
        protein = cds_seq.translate(to_stop=True)
        print(feature.qualifiers.get("product", ["?"])[0], protein)
```

### UniProt REST API
```python
import urllib.request, json, urllib.parse

def fetch_uniprot(accession, fmt="json"):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.{fmt}"
    with urllib.request.urlopen(url) as r:
        return json.loads(r.read()) if fmt == "json" else r.read().decode()

entry = fetch_uniprot("P01308")          # human insulin, JSON
fasta  = fetch_uniprot("P01308", "fasta")

def search_uniprot(query, limit=5):
    q = urllib.parse.quote(query)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&size={limit}&format=json"
    with urllib.request.urlopen(url) as r:
        return json.loads(r.read())

results = search_uniprot("insulin AND reviewed:true AND organism_id:9606")
```

### PDB REST API
```python
# Metadata
def get_pdb_info(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    with urllib.request.urlopen(url) as r:
        return json.loads(r.read())

# Download structure file
def download_pdb(pdb_id, fmt="pdb"):
    base = "https://files.rcsb.org/download"
    url = f"{base}/{pdb_id}.{fmt}"          # fmt: "pdb" or "cif"
    urllib.request.urlretrieve(url, f"{pdb_id}.{fmt}")

# Full-text search
def search_pdb(query_text, max_results=5):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = json.dumps({
        "query": {"type": "terminal", "service": "full_text",
                  "parameters": {"value": query_text}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_results}}
    }).encode()
    req = urllib.request.Request(url, data=payload,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())
```

### Ensembl REST API
```python
def ensembl_get(endpoint):
    url = f"https://rest.ensembl.org{endpoint}"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())

# Gene info by Ensembl ID
info = ensembl_get("/lookup/id/ENSG00000254647?expand=1")

# Sequence (seq_type: genomic  cds  cdna  protein)
def ensembl_sequence(ensembl_id, seq_type="cds"):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?type={seq_type}"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())

# Orthologs
def get_orthologs(gene_id, target_species):
    ep = f"/homology/id/{gene_id}?target_species={target_species}&type=orthologues"
    return ensembl_get(ep)
```

### Pairwise alignment
```python
aligner = Align.PairwiseAligner()
aligner.mode = "global"   # Needleman-Wunsch; "local" = Smith-Waterman
alignments = aligner.align(seq1, seq2)
print(alignments[0].score)
print(alignments[0])
```

### ORF finder (all 6 frames)
```python
def find_orfs(seq, min_length=100):
    orfs = []
    for strand, s in [('+', seq), ('-', seq.reverse_complement())]:
        for frame in range(3):
            for i in range(frame, len(s) - 2, 3):
                if s[i:i+3] == 'ATG':
                    for j in range(i, len(s) - 2, 3):
                        if str(s[j:j+3]) in ('TAA','TAG','TGA'):
                            if j - i >= min_length:
                                orfs.append((strand, frame, i, j+3, s[i:j+3]))
                            break
    return orfs
```

## Pitfalls
- **Always set `Entrez.email`** before any E-utility call; NCBI will block unidentified clients.
- `SeqIO.parse()` returns an **iterator**, not a list — iterating it twice yields nothing. Wrap in `list()` if you need to reuse.
- `SeqIO.read()` raises if the file has 0 or >1 records — use `parse()` for multi-record files.
- `record.annotations["molecule_type"]` is **required** when writing GenBank format; omitting it raises `ValueError`.
- `translate()` without `to_stop=True` includes the stop codon as `*`; `translate(to_stop=True)` is almost always what you want for a clean protein sequence.
- `efetch` with `rettype="gb"` + `SeqIO.read(..., "genbank")` — use `"genbank"` not `"gb"` as the SeqIO format string.
- Rate-limit NCBI calls: max 3/second without API key, 10/second with one. Add `time.sleep(0.4)` in loops.
- UniProt REST base URL changed in 2022: use `rest.uniprot.org` (not `www.uniprot.org/uniprot/`).

## Related Skills
- `sequence-alignment` — pairwise/multiple alignment, scoring matrices
- `ngs-variant-calling` — BLAST via BioPython (`Bio.Blast.NCBIWWW`) and QC pipelines
- `rnaseq` — RNA-seq differential expression, DESeq2/edgeR
