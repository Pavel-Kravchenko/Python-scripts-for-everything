---
name: bio-core-biological-databases
description: Programmatic access to NCBI (Entrez), UniProt, GEO, and SRA — accession types, API patterns, rate limits, and cross-database linking with BioPython
tool_type: python
primary_tool: Python
---

# Biological Databases

## Database Selector

| Need | Database | Preferred subset |
|------|----------|-----------------|
| Gene info + orthologs | NCBI Gene | — |
| mRNA/protein sequence | NCBI Nucleotide / Protein | RefSeq (`NM_`, `NP_`, `NC_`) |
| Protein function/structure | UniProt | SwissProt (`reviewed:true`) |
| Raw sequencing reads | SRA | SRR accession |
| Processed expression data | GEO | GSE (series), GSM (sample) |
| Literature | PubMed | — |

## NCBI Accession Prefixes

| Prefix | Type | Curated? |
|--------|------|----------|
| `NM_` | curated mRNA | Yes (RefSeq) |
| `NP_` | curated protein | Yes (RefSeq) |
| `NC_` | chromosome/complete genome | Yes (RefSeq) |
| `XM_`/`XP_` | predicted/model | No |
| `NR_` | non-coding RNA | Yes (RefSeq) |

GEO accessions: `GSE` (series/experiment), `GSM` (sample), `GPL` (platform), `GDS` (curated dataset).  
SRA hierarchy: Study (SRP) > Sample (SRS) > Experiment (SRX) > Run (SRR).

## Entrez API (BioPython)

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"   # required by NCBI

# Search
handle = Entrez.esearch(
    db="nucleotide",
    term="insulin[Gene] AND Homo sapiens[Organism] AND RefSeq[Filter]",
    retmax=10
)
results = Entrez.read(handle); handle.close()
# results["IdList"], results["Count"]

# Fetch sequence (FASTA)
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3",
                       rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta"); handle.close()

# Fetch full annotations (GenBank)
handle = Entrez.efetch(db="nucleotide", id="NM_000207.3",
                       rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank"); handle.close()

# Extract CDS and translate
for feat in record.features:
    if feat.type == "CDS":
        cds_seq = feat.location.extract(record.seq)
        protein = cds_seq.translate(to_stop=True)
        annotated = feat.qualifiers.get("translation", [""])[0]

# Cross-database link
handle = Entrez.elink(dbfrom="nucleotide", db="protein", id="NM_000207.3")
link_results = Entrez.read(handle); handle.close()

# Batch fetch (comma-separated IDs — faster than looping)
handle = Entrez.efetch(db="nucleotide",
                       id="NM_000207.3,NM_001301717.2",
                       rettype="fasta", retmode="text")
```

### Entrez Search Syntax

```
insulin[Gene]                    # field tag
Homo sapiens[Organism]
mRNA[Filter] AND RefSeq[Filter]
2020:2024[PDAT]                  # date range
CRISPR AND (cancer OR tumor)
```

## UniProt REST API

```python
import urllib.request, json

def fetch_uniprot(accession, fmt="json"):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.{fmt}"
    with urllib.request.urlopen(url) as r:
        return json.loads(r.read()) if fmt == "json" else r.read().decode()

# SwissProt entry (human insulin)
insulin = fetch_uniprot("P01308")
# insulin["proteinDescription"]["recommendedName"]["fullName"]["value"]
# insulin["sequence"]["value"]
# insulin["genes"][0]["geneName"]["value"]

# Extract PDB cross-references
pdb_refs = [r for r in insulin.get("uniProtKBCrossReferences", [])
            if r["database"] == "PDB"]

# Search SwissProt only
url = "https://rest.uniprot.org/uniprotkb/search?query=insulin+AND+reviewed:true&format=json"
```

## Pitfalls

- **GenBank vs RefSeq**: GenBank has all submitted sequences (redundant, unreviewed); RefSeq is curated. For analyses, always use RefSeq accessions
- **SwissProt vs TrEMBL**: SwissProt (~570K entries, manually curated) vs TrEMBL (~250M, auto-annotated); filter with `reviewed:true` for reliable annotations
- **Rate limits**: NCBI allows 3 req/s without API key, 10/s with one (`Entrez.api_key = "..."`); batch with comma-separated IDs, never loop individual fetches
- **GI numbers deprecated**: use accession numbers (`NM_000207.3`), not numeric GI IDs
- **Coordinate systems**: GenBank features use 1-based inclusive; BioPython converts to 0-based half-open internally — use `feat.location.extract(record.seq)` rather than manual slicing
- **`efetch` handle must be closed**: always `handle.close()` or use `with` — unclosed handles exhaust NCBI connections
