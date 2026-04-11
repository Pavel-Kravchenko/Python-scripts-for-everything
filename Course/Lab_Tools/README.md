# Lab Tools

Utility scripts and templates for bioinformatics coursework, organized by category.

## Structure

```
Lab_Tools/
├── Alignment/
│   └── alignment_stats.py   # MSA identity/similarity calculator
└── Visualization/
    └── protein_presentation.spt  # Jmol protein presentation script
```

> **Note:** Rosalind problem-solving tools have been moved to `Assignments/Rosalind/`.

---

## Alignment

### `Alignment/alignment_stats.py`

Calculate identity and similarity percentages for multiple sequence alignments.

```bash
python Alignment/alignment_stats.py -i alignment.fasta -m ../Data/EBLOSUM62.txt
```

**Output:**
- Identical positions (all sequences have same residue)
- Similar positions (all pairwise BLOSUM scores positive)

---

## Visualization

### `Visualization/protein_presentation.spt`

Jmol script template for protein structure presentations.

**Usage in Jmol:**
1. Open Jmol application
2. Load your PDB file: `File → Open`
3. Run script: `script "protein_presentation.spt"` in console
4. Press any key to advance through scenes

**Scenes included:**
1. Chain overview
2. Secondary structure coloring
3. Ligand visualization
4. Active site detail with H-bonds
5. Alpha helix with backbone H-bonds
6. Beta sheet with inter-strand H-bonds
7. Van der Waals surface
8. Hydrophobic residue highlighting
9. Spinning final view

**HTML with Jmol Applet template:**
```html
<!DOCTYPE html>
<html>
<head>
    <title>Protein Viewer</title>
    <script src="JSmol.min.js"></script>
</head>
<body>
    <div id="jmol"></div>
    <button onclick="Jmol.script(jmolApplet, 'cartoon only; color structure')">
        Secondary Structure
    </button>
    <script>
        var jmolApplet = Jmol.getApplet("jmolApplet", {
            width: 600,
            height: 400,
            script: "load protein.pdb"
        });
    </script>
</body>
</html>
```

---

## Resources

- **Jmol Wiki**: https://wiki.jmol.org
- **BioPython**: https://biopython.org
- **BLOSUM Matrices**: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/
