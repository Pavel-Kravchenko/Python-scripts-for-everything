# 🔬 Module 7: Structural Bioinformatics

---

## 📖 Overview

This module covers structural bioinformatics based on materials from the Kodomo program at Moscow State University.
Topics include X-ray crystallography, protein structure analysis, sequence databases, NGS analysis, phylogenetics, motifs/domains, and functional annotation.

---

## 📂 Structure

```
07_Kodomo_Structural_Bioinformatics/
│
├── 01_Electron_Density/           # X-ray crystallography & electron density
│   └── 01_Visualizing_Electron_Density.ipynb
│
├── 02_Crystallography/            # Crystal lattices & space groups
│   └── 01_Crystal_Lattices.ipynb
│
├── 03_PDB_Analysis/               # PDB file format & parsing
│   └── 01_PDB_File_Format.ipynb
│
├── 04_Secondary_Structure/        # DSSP, Ramachandran plots
│   └── 01_Secondary_Structure.ipynb
│
├── 05_Structural_Alignment/       # RMSD, Kabsch algorithm, TM-score
│   └── 01_Superposition_RMSD.ipynb
│
├── 06_Molecular_Surfaces/         # SASA, molecular surfaces
│   └── 01_Surface_Analysis.ipynb
│
├── 07_Sequence_Databases/         # GenBank, FASTA, BioPython Entrez
│   └── 01_GenBank_Databases.ipynb
│
├── 08_BLAST_Search/               # BLAST algorithm & usage
│   └── 01_BLAST_Fundamentals.ipynb
│
├── 09_NGS_Analysis/               # FASTQ, quality control, VCF
│   └── 01_NGS_Fundamentals.ipynb
│
├── 10_Nucleic_Acids/              # DNA/RNA structure, sequencing, SNPs
│   ├── 01_DNA_RNA_Structure.ipynb
│   ├── 02_Sanger_Sequencing.ipynb
│   └── 03_SNP_Analysis.ipynb
│
├── 11_Phylogenetics/              # Phylogenetic trees, Neighbor-Joining
│   └── 01_Phylogenetic_Trees.ipynb
│
├── 12_Motifs_Domains/             # PWM, sequence motifs, Pfam
│   ├── 01_PWM_Sequence_Motifs.ipynb
│   └── 02_Protein_Domains_Pfam.ipynb
│
├── 13_Functional_Annotation/      # KEGG, GO, transmembrane proteins
│   ├── 01_KEGG_Pathways.ipynb
│   ├── 02_Gene_Ontology.ipynb
│   └── 03_Transmembrane_Proteins.ipynb
│
└── images/                        # Supporting images
```

---

## 🎯 Topics Covered

| # | Topic | Notebook | Key Concepts |
|---|-------|----------|--------------|
| 1 | Electron Density | `01_Visualizing_Electron_Density.ipynb` | 2Fo-Fc maps, Fo-Fc maps, resolution |
| 2 | Crystallography | `01_Crystal_Lattices.ipynb` | Unit cells, 230 space groups, asymmetric units |
| 3 | PDB Analysis | `01_PDB_File_Format.ipynb` | ATOM records, PDBParser class, headers |
| 4 | Secondary Structure | `01_Secondary_Structure.ipynb` | DSSP codes, φ/ψ angles, Ramachandran |
| 5 | Structural Alignment | `01_Superposition_RMSD.ipynb` | RMSD, Kabsch algorithm, TM-score |
| 6 | Molecular Surfaces | `01_Surface_Analysis.ipynb` | VdW, SAS, SES, SASA calculation |
| 7 | Sequence Databases | `01_GenBank_Databases.ipynb` | INSDC, FASTA, GenBank, Entrez |
| 8 | BLAST Search | `01_BLAST_Fundamentals.ipynb` | Seeding, E-values, blastp/blastn |
| 9 | NGS Analysis | `01_NGS_Fundamentals.ipynb` | FASTQ, Phred scores, VCF format |
| 10 | Nucleic Acids | `01_DNA_RNA_Structure.ipynb` | DNA/RNA chemistry, A/B forms, secondary structure |
| 10b | Sanger Sequencing | `02_Sanger_Sequencing.ipynb` | Dideoxy method, chromatograms |
| 10c | SNP Analysis | `03_SNP_Analysis.ipynb` | VCF format, variant annotation, dbSNP |
| 11 | Phylogenetics | `01_Phylogenetic_Trees.ipynb` | Neighbor-Joining, bootstrap, Newick format |
| 12 | Motifs | `01_PWM_Sequence_Motifs.ipynb` | PFM, PWM, sequence logos, JASPAR |
| 12b | Domains | `02_Protein_Domains_Pfam.ipynb` | Pfam, HMM profiles, InterPro |
| 13 | KEGG | `01_KEGG_Pathways.ipynb` | Pathway enrichment, KEGG IDs |
| 13b | Gene Ontology | `02_Gene_Ontology.ipynb` | GO terms, DAG structure, enrichment |
| 13c | Transmembrane | `03_Transmembrane_Proteins.ipynb` | Hydropathy, TM helix prediction |

---

## 🛠️ Tools & Libraries

| Tool | Purpose |
|------|---------|
| **BioPython** | Sequence analysis, Entrez queries, phylogenetics |
| **NumPy** | Numerical calculations (RMSD, Kabsch, PWM) |
| **Jmol/PyMOL** | 3D molecular visualization |
| **BLAST+** | Sequence similarity search |
| **FastQC** | NGS quality control |
| **JASPAR** | Transcription factor binding motifs |
| **Pfam/InterPro** | Protein domain databases |
| **KEGG/GO** | Functional annotation databases |

---

## 💻 Code Highlights

### PDB Parser (from notebook 03)

```python
class PDBParser:
    def parse(self, filename):
        for line in open(filename):
            if line.startswith('ATOM'):
                atom = {
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resnum': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54])
                }
```

### RMSD Calculation (from notebook 05)

```python
def rmsd(coords1, coords2):
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
```

### Phred Score (from notebook 09)

```python
def phred_to_probability(phred):
    return 10 ** (-phred / 10)
```

---

## 📚 Resources

- [RCSB PDB](https://www.rcsb.org/) - Protein Data Bank
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) - Sequence search
- [UniProt](https://www.uniprot.org/) - Protein sequences
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - QC tool
- [JASPAR](https://jaspar.elixir.no/) - TF binding profiles
- [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/) - Protein families
- [KEGG](https://www.genome.jp/kegg/) - Pathways database
- [Gene Ontology](http://geneontology.org/) - GO terms
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) - SNP database

---

## 🎓 Learning Path

1. Start with **Electron Density** and **Crystallography** for structure determination basics
2. Learn **PDB Analysis** to work with structure files
3. Study **Secondary Structure** and **Structural Alignment** for protein comparison
4. Explore **Molecular Surfaces** for interaction analysis
5. Use **Sequence Databases** and **BLAST** for homology searches
6. Study **NGS Analysis** for sequencing data
7. Learn **Nucleic Acids** and **Phylogenetics** for evolutionary analysis
8. Master **Motifs/Domains** for sequence pattern recognition
9. Finish with **Functional Annotation** (KEGG, GO, Transmembrane) for biological interpretation
