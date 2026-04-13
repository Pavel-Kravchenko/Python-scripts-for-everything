---
name: bio-core-biopython-essentials
description: "BioPython is the most widely used Python library for bioinformatics. This notebook covers its core modules through practical, hands-on examples -- from sequence manipulation to fetching data from NCBI"
tool_type: python
source_notebook: "Tier_2_Core_Bioinformatics/02_BioPython_Essentials/01_biopython_essentials.ipynb"
primary_tool: Python
---

## Version Compatibility

Reference examples tested with: biopython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.


# BioPython Essentials

*Source: Course notebook `Tier_2_Core_Bioinformatics/02_BioPython_Essentials/01_biopython_essentials.ipynb`*


BioPython is the most widely used Python library for bioinformatics. This notebook covers its core modules through practical, hands-on examples -- from sequence manipulation to fetching data from NCBI and building a complete gene-to-protein analysis workflow.

## Learning Objectives

- Create and manipulate `Seq` objects (complement, reverse complement, translate, transcribe)
- Work with `SeqRecord` objects (annotations, features, letter_annotations)
- Read and write biological file formats with `SeqIO` (FASTA, GenBank, FASTQ)
- Search and fetch data from NCBI using `Bio.Entrez`
- Parse sequence alignments with `AlignIO`
- Build a complete workflow: from fetching a gene to analyzing its protein product

**Prerequisites:** Basic Python, understanding of DNA/RNA/protein central dogma

## How to use this notebook

1. Run cells top-to-bottom in order — later cells depend on earlier ones
2. Run the environment check cell first to confirm all imports work
3. Read the explanatory text before each code cell — the context matters
4. The exercises at the end are designed to be done after reading each section
5. If a code cell requires internet access (Entrez, PDB download), it is marked — these can be skipped if offline

## Complicated moments explained

- **Seq vs. str**: A `Seq` object behaves like a Python string (supports slicing, len(), iteration) but adds biological methods: `.complement()`, `.reverse_complement()`, `.transcribe()`, `.translate()`. In modern BioPython (>=1.79), `Seq` and `str` can be freely combined in most contexts. Pre-1.79, you had to call `str()` explicitly.
- **translate() and stop codons**: `dna.translate()` translates all codons, representing stop codons as `*`. `dna.translate(to_stop=True)` stops at the first stop codon and omits the `*`. For extracting the protein product of a CDS, always use `to_stop=True`.
- **SeqIO.parse() vs SeqIO.read()**: `parse()` returns an *iterator* of SeqRecord objects — use it for multi-record files (always safe). `read()` returns a single SeqRecord — use it only when the file is guaranteed to contain exactly one record (raises an error otherwise).
- **Entrez rate limits**: Always call `handle.close()` after reading to release the connection. NCBI allows 3 requests/second without an API key. For batch work, get a free NCBI API key and set `Entrez.api_key = "yourkey"`.

## Environment check (run this first)

```python
# Environment check
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO, Entrez, AlignIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import gc_fraction, molecular_weight

print(f"BioPython version: {Bio.__version__}")

# Quick sanity checks
dna = Seq("ATGAAACCCGGGTAA")
protein = dna.translate(to_stop=True)
print(f"Seq object: {dna} -> translate(to_stop=True) -> {protein}")
print(f"GC content: {gc_fraction(dna)*100:.1f}%")

blosum62 = substitution_matrices.load("BLOSUM62")
print(f"BLOSUM62 loaded: score(L,I)={blosum62['L','I']}, score(K,D)={blosum62['K','D']}")
print("\nAll imports ready. Proceed to Section 1.")
```python

```python
dna = Seq("ATGCGATCGATCGTAA")

print(f"Original:           5'-{dna}-3'")
print(f"Complement:         3'-{dna.complement()}-5'")
print(f"Reverse complement: 5'-{dna.reverse_complement()}-3'")
print()

# Transcription: DNA -> RNA (replaces T with U)
rna = dna.transcribe()
print(f"DNA: {dna}")
print(f"RNA: {rna}")

# Back-transcription: RNA -> DNA
print(f"Back to DNA: {rna.back_transcribe()}")
print(f"Round-trip OK: {rna.back_transcribe() == dna}")
```python

### 1.2 Translation: DNA/RNA to Protein

Each triplet of nucleotides (codon) codes for one amino acid. The `*` character represents a stop codon.

```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# Full translation (includes stop codon as *)
print(f"Full:        {coding_dna.translate()}")

# Stop at first stop codon (most common usage)
print(f"to_stop:     {coding_dna.translate(to_stop=True)}")
print()

# Different genetic codes
dna2 = Seq("ATGTGAATGGAA")
print(f"Codons:         ATG TGA ATG GAA")
print(f"Standard (1):   {dna2.translate(table=1)}  (TGA = stop)")
print(f"Mito (2):       {dna2.translate(table=2)}  (TGA = Trp)")
print(f"Bacterial (11): {dna2.translate(table=11)}  (TGA = stop)")
```python

```python
# GC content and molecular weight
dna = Seq("ATGCGATCGATCGATCG")
print(f"GC content: {gc_fraction(dna) * 100:.1f}%")
print(f"DNA MW: {molecular_weight(dna):.1f} Da")

protein = Seq("MKPG")
print(f"Protein MW: {molecular_weight(protein, seq_type='protein'):.1f} Da")

# Mutable sequences (for in-place editing)
mutable = MutableSeq("ATGCGATCG")
mutable[3] = "T"  # G -> T point mutation
print(f"\nMutated: {mutable} (position 3: G->T)")
```python

---
## 2. SeqRecord Objects

`SeqRecord` wraps a `Seq` object and adds metadata: ID, name, description, annotations, features, and per-letter annotations. This is the main data structure you get when reading files with `SeqIO`.

```python
# Creating a SeqRecord with annotations
record = SeqRecord(
    Seq("ATGCGATCGATCGATCGATCGATCGTAA"),
    id="BRCA1_001",
    name="BRCA1",
    description="Breast cancer type 1 susceptibility protein, partial CDS"
)

# Record-level annotations
record.annotations["molecule_type"] = "DNA"
record.annotations["organism"] = "Homo sapiens"
record.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Mammalia"]

print(f"ID: {record.id}  |  Name: {record.name}")
print(f"Description: {record.description}")
print(f"Sequence: {record.seq}  ({len(record)} bp)")
print(f"Annotations: {dict(record.annotations)}")
```python

```python
# Adding features (annotated regions on the sequence)
gene_feature = SeqFeature(
    FeatureLocation(start=0, end=28), type="gene",
    qualifiers={"gene": ["BRCA1"]}
)
cds_feature = SeqFeature(
    FeatureLocation(start=0, end=27), type="CDS",
    qualifiers={"gene": ["BRCA1"], "product": ["BRCA1 protein"]}
)
record.features = [gene_feature, cds_feature]

for feat in record.features:
    print(f"{feat.type:6s} {feat.location}  qualifiers={dict(feat.qualifiers)}")

# Extract the nucleotide sequence of a feature
cds_seq = cds_feature.location.extract(record.seq)
print(f"\nCDS protein: {cds_seq.translate(to_stop=True)}")
```python

```python
# Per-letter annotations (e.g., quality scores in FASTQ)
short_record = SeqRecord(Seq("ATGCGATCG"), id="read001")
short_record.letter_annotations["phred_quality"] = [30, 30, 28, 35, 35, 33, 30, 25, 20]

for base, qual in zip(short_record.seq, short_record.letter_annotations["phred_quality"]):
    print(f"  {base} Q={qual}", end="")
print(f"\nMean quality: {sum(short_record.letter_annotations['phred_quality']) / len(short_record):.1f}")
```python

---
## 3. SeqIO: Reading and Writing Sequence Files

`SeqIO` provides a uniform interface for many biological file formats:
- `SeqIO.parse(file, fmt)` -- read multiple records (iterator)
- `SeqIO.read(file, fmt)` -- read exactly one record
- `SeqIO.write(records, file, fmt)` -- write records
- `SeqIO.to_dict(...)` -- load into a dictionary keyed by ID
- `SeqIO.convert(in_file, in_fmt, out_file, out_fmt)` -- format conversion

```python
# Create a sample FASTA file
fasta_content = """>INS_HUMAN Insulin [Homo sapiens]
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT
RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
>INS_MOUSE Insulin [Mus musculus]
MALWTRLLPLLALLALWAPAPAHAFVNQHLCGSHLVEALYLVCGERGFFYTPKS
RREVEDPQVAQLELGGGPGADDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN
>INS_RAT Insulin [Rattus norvegicus]
MALWMHLLTVLALLALWGPNSVQAFVKQHLCGSHLVEALYLVCGERGFFYTPKS
RREVEDPQVEQLELGGSPGDLQTLALEVARQKRGIVDQCCTSICSLYQLENYCN
"""
with open("insulin_seqs.fasta", "w") as f:
    f.write(fasta_content)

# Read and display
for record in SeqIO.parse("insulin_seqs.fasta", "fasta"):
    print(f"{record.id:15s} {len(record):4d} aa  {record.description}")
```python

```python
# Load into dictionary, write new file
records_dict = SeqIO.to_dict(SeqIO.parse("insulin_seqs.fasta", "fasta"))
print(f"Loaded IDs: {list(records_dict.keys())}")
print(f"Human insulin: {records_dict['INS_HUMAN'].seq[:40]}...")

# Write selected records
new_records = [
    SeqRecord(Seq("ATGCGATCGATCG"), id="seq1", description="Test sequence 1"),
    SeqRecord(Seq("GCTAGCTAGCTAG"), id="seq2", description="Test sequence 2"),
]
count = SeqIO.write(new_records, "output.fasta", "fasta")
print(f"\nWrote {count} sequences to output.fasta")
```python

### 3.2 GenBank Files

```python
# Create a sample GenBank file
genbank_content = """LOCUS       EXAMPLE_GENE          120 bp    DNA     linear   PRI 01-JAN-2024
DEFINITION  Homo sapiens example gene, complete CDS.
ACCESSION   EX000001
VERSION     EX000001.1
KEYWORDS    .
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata;
            Mammalia; Primates; Haplorrhini; Catarrhini; Hominidae;
            Homo.
FEATURES             Location/Qualifiers
     source          1..120
                     /organism="Homo sapiens"
                     /mol_type="mRNA"
     gene            1..120
                     /gene="EXAMPLE"
     CDS             11..100
                     /gene="EXAMPLE"
                     /codon_start=1
                     /product="example protein"
                     /protein_id="EXP00001.1"
                     /translation="MAISRVTLGERILKLCHELY"
     5'UTR           1..10
                     /gene="EXAMPLE"
     3'UTR           101..120
                     /gene="EXAMPLE"
ORIGIN
        1 cttcacagcg atggctatca gccgtgtgac actgggcgag cgtatcctga agctgtgcca
       61 tgagctctat tagcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga
//
"""
with open("example.gb", "w") as f:
    f.write(genbank_content)

# Read and explore
record = SeqIO.read("example.gb", "genbank")
print(f"ID: {record.id}  Name: {record.name}")
print(f"Description: {record.description}")
print(f"Length: {len(record)} bp")
print(f"Organism: {record.annotations.get('organism', 'N/A')}")

print(f"\nFeatures ({len(record.features)}):")
for feat in record.features:
    gene = feat.qualifiers.get('gene', [''])[0]
    print(f"  {feat.type:10s} {str(feat.location):20s} {gene}")
```python

```python
# Extract CDS and translate
for feat in record.features:
    if feat.type == "CDS":
        cds_seq = feat.location.extract(record.seq)
        annotated = feat.qualifiers.get("translation", [""])[0]
        our_translation = cds_seq.translate(to_stop=True)
        
        print(f"Gene: {feat.qualifiers.get('gene', ['?'])[0]}")
        print(f"CDS location: {feat.location}")
        print(f"Annotated protein: {annotated}")
        print(f"Our translation:   {our_translation}")
        break
```python

```python
# Create sample FASTQ, read, and filter by quality
fastq_content = """@read001 Example read 1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIHHHHHGGGGFFFFEEEEDDDDDCCCCBBBB
@read002 Example read 2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
+
IIIIIIIIIHHHHHHHGGGGGFFFFFEEEEEDDDDDCC
@read003 Low quality read
TTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAAC
+
BBBBAAA@@@@?????>>>>>>======<<<<<<;;;;;
"""
with open("reads.fastq", "w") as f:
    f.write(fastq_content)

# Read and show quality stats
for record in SeqIO.parse("reads.fastq", "fastq"):
    quals = record.letter_annotations["phred_quality"]
    mean_q = sum(quals) / len(quals)
    print(f"{record.id:10s} len={len(record):3d}  mean_Q={mean_q:.1f}  min_Q={min(quals)}")

# Quality filtering
good = [r for r in SeqIO.parse("reads.fastq", "fastq")
        if sum(r.letter_annotations["phred_quality"]) / len(r) >= 25]
SeqIO.write(good, "filtered.fastq", "fastq")
print(f"\nKept {len(good)}/3 reads with mean Q >= 25")
```python

```python
# Format conversion
# GenBank -> FASTA
records = list(SeqIO.parse("example.gb", "genbank"))
SeqIO.write(records, "from_genbank.fasta", "fasta")

# FASTQ -> FASTA (one-liner, note: quality info is lost)
count = SeqIO.convert("reads.fastq", "fastq", "reads.fasta", "fasta")
print(f"Converted {count} FASTQ records to FASTA")

with open("reads.fasta") as f:
    print(f.read())
```python

---
## 4. Bio.Entrez: Accessing NCBI Databases

The `Entrez` module provides programmatic access to all NCBI databases. The key functions:
- `Entrez.esearch(db, term)` -- search, returns IDs
- `Entrez.efetch(db, id, rettype, retmode)` -- download records
- `Entrez.elink(dbfrom, db, id)` -- find cross-database links

## Common Pitfalls

- **Coordinate systems**: BED uses 0-based half-open; VCF/GFF use 1-based inclusive — mixing them causes off-by-one errors
- **Batch effects**: Always check for batch confounding before interpreting biological signal
- **Multiple testing**: Apply FDR correction (Benjamini-Hochberg) when testing thousands of features simultaneously
