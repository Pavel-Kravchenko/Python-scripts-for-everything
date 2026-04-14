---
name: biopython-databases
description: BioPython Seq/SeqRecord/SeqIO, NCBI Entrez API, UniProt queries, and biological format conversion
tool_type: python
primary_tool: Python
---

# BioPython — Biological Databases

**Use when:** fetching sequences/annotations from NCBI/UniProt/PDB/Ensembl; reading/writing FASTA, GenBank, EMBL, FASTQ; manipulating sequences (complement, translate, GC%, ORF finding); batch downloading records programmatically.

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

Entrez.email = "your.email@example.com"  # required before any E-utility call
```

### Seq Object
```python
dna = Seq("ATGCGATCGATCGTAA")
dna.complement()           # 3'->5' complement
dna.reverse_complement()   # 5'->3' complement of the other strand
dna.transcribe()           # DNA -> RNA (T -> U)
dna.translate()            # returns protein, * = stop codon
dna.translate(to_stop=True)  # stops before stop codon (most common)
dna.translate(table=2)     # mitochondrial genetic code

gc_fraction(dna) * 100     # GC percentage
molecular_weight(dna)      # DNA MW in Da
molecular_weight(protein, seq_type="protein")

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

feat = SeqFeature(FeatureLocation(0, 27), type="CDS",
                  qualifiers={"gene": ["BRCA1"], "product": ["BRCA1 protein"]})
record.features.append(feat)

record.letter_annotations["phred_quality"] = [30, 28, 35, ...]
cds_seq = feature.location.extract(record.seq)
```

### SeqIO: Reading & Writing
```python
records = list(SeqIO.parse("file.fasta", "fasta"))   # iterator -> wrap in list() to reuse
record  = SeqIO.read("single.gb", "genbank")          # exactly one record

SeqIO.write(records, "out.fasta", "fasta")

d = SeqIO.to_dict(SeqIO.parse("file.fasta", "fasta"))
d["INS_HUMAN"].seq

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
results = Entrez.read(handle); handle.close()
ids = results["IdList"]

handle = Entrez.efetch(db="nucleotide", id="NM_000207.3", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank"); handle.close()
```

### Batch fetch
```python
accessions = ["NM_000518.5", "NM_000207.3", "NM_000546.6"]
handle = Entrez.efetch(db="nucleotide", id=",".join(accessions), rettype="fasta", retmode="text")
records = list(SeqIO.parse(handle, "fasta")); handle.close()
SeqIO.write(records, "batch.fasta", "fasta")
```

### elink: cross-database navigation
```python
handle = Entrez.elink(dbfrom="nucleotide", db="protein", id="NM_000207.3")
link_results = Entrez.read(handle); handle.close()
protein_ids = [link["Id"] for linkset in link_results
               for linkdb in linkset["LinkSetDb"]
               for link in linkdb["Link"]]

handle = Entrez.elink(dbfrom="gene", db="pubmed", id="7157")  # TP53
```

### Extract CDS and translate from GenBank record
```python
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank"); handle.close()

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

entry = fetch_uniprot("P01308")           # human insulin, JSON
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
def get_pdb_info(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    with urllib.request.urlopen(url) as r:
        return json.loads(r.read())

def download_pdb(pdb_id, fmt="pdb"):
    url = f"https://files.rcsb.org/download/{pdb_id}.{fmt}"
    urllib.request.urlretrieve(url, f"{pdb_id}.{fmt}")

def search_pdb(query_text, max_results=5):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = json.dumps({
        "query": {"type": "terminal", "service": "full_text",
                  "parameters": {"value": query_text}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_results}}
    }).encode()
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"})
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

info = ensembl_get("/lookup/id/ENSG00000254647?expand=1")

def ensembl_sequence(ensembl_id, seq_type="cds"):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?type={seq_type}"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())

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
- **Always set `Entrez.email`** before any E-utility call; NCBI blocks unidentified clients.
- `SeqIO.parse()` returns an **iterator** — iterating twice yields nothing. Wrap in `list()` to reuse.
- `SeqIO.read()` raises if file has 0 or >1 records — use `parse()` for multi-record files.
- `record.annotations["molecule_type"]` is **required** for GenBank output; omitting it raises `ValueError`.
- Use `translate(to_stop=True)` for clean protein sequences (omits trailing `*`).
- `efetch` with `rettype="gb"` + SeqIO: use `"genbank"` not `"gb"` as the format string.
- Rate-limit NCBI: max 3/sec without API key, 10/sec with. Add `time.sleep(0.4)` in loops.
- UniProt REST base URL changed in 2022: use `rest.uniprot.org` (not `www.uniprot.org/uniprot/`).
