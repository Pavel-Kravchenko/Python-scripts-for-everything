---
name: bio-core-biological-databases
description: "Modern bioinformatics depends on large, well-curated public databases. This notebook covers the major databases you will use daily as a bioinformatician: where biological data lives, what each databas"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/01_Biological_Databases/01_biological_databases.ipynb"
---

# Biological Databases

*Source: Course notebook `Tier_2_Core_Bioinformatics/01_Biological_Databases/01_biological_databases.ipynb`*

# Biological Databases

Modern bioinformatics depends on large, well-curated public databases. This notebook covers the major databases you will use daily as a bioinformatician: where biological data lives, what each database contains, and how to access it programmatically.

## Learning Objectives

- Navigate the NCBI suite of databases (GenBank, RefSeq, Gene, SRA, GEO)
- Understand UniProt (SwissProt vs TrEMBL) and search for protein information
- Use the PDB to find and evaluate protein structures
- Access Ensembl genome browser and BioMart
- Retrieve data programmatically using Bio.Entrez and REST APIs

**Prerequisites:** Basic Python, familiarity with DNA/RNA/protein sequences

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **GenBank vs. RefSeq**: GenBank accepts all submitted sequences (may be redundant, unreviewed, or contain sequencing errors). RefSeq is NCBI's curated, non-redundant reference set. For most analyses, use RefSeq. Identify RefSeq records by their prefix: `NM_` (curated mRNA), `NP_` (curated protein), `NC_` (chromosome/complete genome), `XM_`/`XP_` (predicted/model).
- **UniProt/SwissProt vs. TrEMBL**: SwissProt (~570,000 entries) is manually curated with experimental annotations. TrEMBL (~250 million entries) is computationally annotated. For well-studied proteins, always prefer SwissProt. Use the filter `reviewed:true` in UniProt search to restrict to SwissProt.
- **E-utilities rate limits**: NCBI allows 3 requests/second without an API key, 10/second with one. Always set `Entrez.email` (NCBI's requirement) and use `Entrez.sleep_between_tries = True`. For large batch downloads, use `efetch` with comma-separated IDs rather than looping.
- **Accession vs. GI vs. UID**: Modern NCBI uses accession numbers (e.g., NM_000207.3) as the stable identifier. Numeric GI numbers are being deprecated. When retrieving records programmatically, use accession numbers.

## Environment check (run this first)

```python
# Environment check
from Bio import Entrez, SeqIO
import urllib.request
import json

Entrez.email = "your.email@example.com"  # Always set this

# Quick connectivity test (does not consume rate limits)
try:
    handle = Entrez.einfo()
    record = Entrez.read(handle)
    handle.close()
    db_count = len(record["DbList"])
    print(f"BioPython Entrez: connected to NCBI. {db_count} databases available.")
except Exception as e:
    print(f"NCBI connection failed: {e}")
    print("Check your internet connection. You can still read the notebook offline.")

print("\nImports ready: Bio.Entrez, Bio.SeqIO, urllib.request, json")
```

### 1.3 NCBI Gene

The Gene database provides a **gene-centric** view that links to all associated data:
- RefSeq transcripts and proteins
- Genomic context and neighboring genes
- Orthologs and homologs
- Expression data, pathways, and phenotypes
- Literature references

Each gene has a unique Gene ID (e.g., Gene ID 672 = BRCA1).

URL pattern: `https://www.ncbi.nlm.nih.gov/gene/<GeneID>`

### 1.4 SRA (Sequence Read Archive)

SRA stores **raw sequencing reads** from next-generation sequencing experiments.

- Largest publicly available repository of sequencing data
- Contains reads from Illumina, PacBio, Oxford Nanopore, etc.
- Organized by: Study (SRP) > Sample (SRS) > Experiment (SRX) > Run (SRR)
- Download with `fastq-dump` or `fasterq-dump` from the SRA Toolkit

**Example:** SRR1234567 is a specific sequencing run you can download.

### 1.5 GEO (Gene Expression Omnibus)

GEO stores **processed** gene expression data:
- Microarray data
- RNA-seq quantifications
- ChIP-seq peak calls
- Methylation arrays

GEO accession types:
- **GSE** (Series): a complete experiment with multiple samples
- **GSM** (Sample): data from one biological sample
- **GPL** (Platform): the array/technology description
- **GDS** (DataSet): curated, analysis-ready datasets

---
## 2. Programmatic Access with Bio.Entrez

The Entrez Programming Utilities (E-utilities) let you search and retrieve data from all NCBI databases via a single API. BioPython wraps these in `Bio.Entrez`.

The three main functions:

| Function | Purpose |
|----------|---------|
| `Entrez.einfo()` | Get database information |
| `Entrez.esearch()` | Search a database |
| `Entrez.efetch()` | Download records |

```python
# List all available NCBI databases
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

db_list = record["DbList"]
print(f"NCBI has {len(db_list)} databases:\n")

# Display in columns
for i in range(0, len(db_list), 5):
    row = db_list[i:i+5]
    print("  ".join(f"{db:<18}" for db in row))
```

```python
# Get details about a specific database
handle = Entrez.einfo(db="nucleotide")
record = Entrez.read(handle)
handle.close()

info = record["DbInfo"]
print(f"Database: {info['DbName']}")
print(f"Description: {info['Description']}")
print(f"Total records: {info['Count']}")
print(f"Last updated: {info['LastUpdate']}")
print(f"\nSearchable fields ({len(info['FieldList'])}):\n")
for field in info['FieldList'][:10]:
    print(f"  [{field['Name']}] - {field['FullName']}: {field['Description']}")
```

### 2.1 Searching NCBI with esearch()

The `esearch()` function searches a database and returns matching record IDs.

**Search syntax tips:**
- `insulin[Gene]` -- search in the Gene field
- `Homo sapiens[Organism]` -- search by organism
- `AND`, `OR`, `NOT` -- Boolean operators
- `RefSeq[Filter]` -- limit to RefSeq entries
- `2020:2024[PDAT]` -- publication date range

```python
# Search for human insulin mRNA sequences in RefSeq
handle = Entrez.esearch(
    db="nucleotide",
    term="insulin[Gene] AND Homo sapiens[Organism] AND mRNA[Filter] AND RefSeq[Filter]",
    retmax=10
)
results = Entrez.read(handle)
handle.close()

print(f"Total matches: {results['Count']}")
print(f"IDs returned: {results['IdList']}")
print(f"\nNote: retmax={10}, so up to 10 IDs are returned.")
print(f"Set retmax higher to get more IDs.")
```

```python
# Search PubMed for recent bioinformatics papers
handle = Entrez.esearch(
    db="pubmed",
    term="CRISPR AND bioinformatics AND 2023:2024[PDAT]",
    retmax=5
)
results = Entrez.read(handle)
handle.close()

print(f"Total CRISPR + bioinformatics papers (2023-2024): {results['Count']}")
print(f"First 5 PubMed IDs: {results['IdList']}")
```

### 2.2 Fetching Records with efetch()

Once you have IDs from `esearch()`, use `efetch()` to download actual records.

Key parameters:
- `db` -- which database
- `id` -- accession or UID (can be a comma-separated list)
- `rettype` -- return type: `"fasta"`, `"gb"` (GenBank), `"gp"` (GenPept)
- `retmode` -- return mode: `"text"` or `"xml"`

```python
# Fetch a sequence in FASTA format
handle = Entrez.efetch(
    db="nucleotide",
    id="NM_000207.3",   # Human insulin mRNA
    rettype="fasta",
    retmode="text"
)
record = SeqIO.read(handle, "fasta")
handle.close()

print(f"ID: {record.id}")
print(f"Description: {record.description}")
print(f"Length: {len(record)} bp")
print(f"First 60 nt: {record.seq[:60]}")
```

```python
# Fetch in GenBank format to get full annotations
handle = Entrez.efetch(
    db="nucleotide",
    id="NM_000207.3",
    rettype="gb",
    retmode="text"
)
record = SeqIO.read(handle, "genbank")
handle.close()

print(f"ID: {record.id}")
print(f"Name: {record.name}")
print(f"Description: {record.description}")
print(f"Sequence length: {len(record)} bp")
print(f"\nAnnotations:")
for key in ['organism', 'taxonomy', 'keywords']:
    if key in record.annotations:
        print(f"  {key}: {record.annotations[key]}")

print(f"\nFeatures ({len(record.features)}):")
for feat in record.features:
    print(f"  {feat.type:12s} {feat.location}")
```

```python
# Extract CDS from the GenBank record and translate it
for feature in record.features:
    if feature.type == "CDS":
        # Get the CDS nucleotide sequence
        cds_seq = feature.location.extract(record.seq)
        print(f"CDS location: {feature.location}")
        print(f"CDS length: {len(cds_seq)} bp")
        
        # Get the annotated translation (from the GenBank file)
        if "translation" in feature.qualifiers:
            annotated_protein = feature.qualifiers["translation"][0]
            print(f"\nAnnotated protein ({len(annotated_protein)} aa):")
            print(f"  {annotated_protein[:60]}...")
        
        # Translate it ourselves for verification
        our_protein = cds_seq.translate(to_stop=True)
        print(f"\nOur translation ({len(our_protein)} aa):")
        print(f"  {str(our_protein)[:60]}...")
        
        # Gene name
        if "gene" in feature.qualifiers:
            print(f"\nGene: {feature.qualifiers['gene'][0]}")
        break
```

```python
# Fetch a PubMed abstract
handle = Entrez.efetch(
    db="pubmed",
    id="27886196",
    rettype="abstract",
    retmode="text"
)
abstract_text = handle.read()
handle.close()

print("PubMed record:")
print(abstract_text[:1000])
```

### 2.3 Linking Between Databases with elink()

NCBI databases are cross-linked. Use `Entrez.elink()` to find related records across databases.

```python
# Find the protein record linked to our nucleotide record
handle = Entrez.elink(
    dbfrom="nucleotide",
    db="protein",
    id="NM_000207.3"
)
link_results = Entrez.read(handle)
handle.close()

# Extract linked protein IDs
for linkset in link_results:
    for link_db in linkset["LinkSetDb"]:
        print(f"Link name: {link_db['LinkName']}")
        for link in link_db["Link"][:5]:
            print(f"  Protein ID: {link['Id']}")
```

---
## 3. UniProt: The Protein Knowledgebase

UniProt is the most comprehensive protein database, organized into two sections:

| Section | Curation | Size | Use when you need |
|---------|----------|------|-------------------|
| **SwissProt** | Manually reviewed | ~570,000 entries | Reliable, well-annotated protein info |
| **TrEMBL** | Automatically annotated | ~250,000,000 entries | Broader coverage, less curated |

**What UniProt provides per protein entry:**
- Function description
- Subcellular location
- Post-translational modifications
- Protein domains and families
- Disease associations
- 3D structure cross-references
- Gene Ontology (GO) annotations
- Sequence variants and isoforms

```python
# Access UniProt via their REST API
def fetch_uniprot(accession, fmt="json"):
    """Fetch a UniProt entry by accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.{fmt}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        if fmt == "json":
            return json.loads(response.read().decode())
        return response.read().decode()

# Fetch human insulin (SwissProt entry)
insulin = fetch_uniprot("P01308")

print(f"Protein: {insulin['proteinDescription']['recommendedName']['fullName']['value']}")
print(f"Gene: {insulin['genes'][0]['geneName']['value']}")
print(f"Organism: {insulin['organism']['scientificName']}")
print(f"Length: {insulin['sequence']['length']} aa")
print(f"\nSequence:")
seq = insulin['sequence']['value']
for i in range(0, len(seq), 60):
    print(f"  {seq[i:i+60]}")
```

```python
# Explore protein function and annotations
print("Function:")
for comment in insulin.get('comments', []):
    if comment['commentType'] == 'FUNCTION':
        for text in comment['texts']:
            print(f"  {text['value'][:200]}")

print("\nSubcellular location:")
for comment in insulin.get('comments', []):
    if comment['commentType'] == 'SUBCELLULAR LOCATION':
        for subloc in comment.get('subcellularLocations', []):
            loc = subloc.get('location', {}).get('value', 'N/A')
            print(f"  {loc}")

print("\nCross-references to PDB structures:")
pdb_refs = [ref for ref in insulin.get('uniProtKBCrossReferences', [])
            if ref['database'] == 'PDB']
for ref in pdb_refs[:5]:
    props = {p['key']: p['value'] for p in ref.get('properties', [])}
    print(f"  {ref['id']} - Method: {props.get('Method', 'N/A')}, "
          f"Resolution: {props.get('Resolution', 'N/A')}")
```
