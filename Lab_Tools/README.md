# Lab Tools

Utility scripts and templates for bioinformatics work.

## Scripts

### `rosalind.py`
Python wrapper for working with Rosalind bioinformatics problems (rosalind.info).

**Features:**
- List available problems from any track (Stronghold, Python Village, Armory, etc.)
- Fetch problem details including sample input/output
- Test solutions against sample data
- Generate Jupyter notebook templates for problems

**Quick Start:**
```python
from Lab_Tools.rosalind import Rosalind

ros = Rosalind()

# List first 10 problems
for p in ros.list_problems()[:10]:
    print(f"{p['id']} - {p['title']}")

# Get problem details
dna = ros.get_problem("DNA")
print(dna.sample_input)
print(dna.sample_output)

# Test your solution
def solve_dna(data):
    s = data.strip()
    return f"{s.count('A')} {s.count('C')} {s.count('G')} {s.count('T')}"

ros.test_solution("DNA", solve_dna)
# ✓ Solution passed for DNA!

# Generate a notebook template
ros.create_notebook("RNA", output_dir="./problems")
```

**Available Tracks:**
- `stronghold` - Bioinformatics Stronghold (default)
- `python-village` - Python basics
- `armory` - Bioinformatics tools
- `textbook` - Textbook Track
- `algorithmic` - Algorithmic Heights

### `alignment_stats.py`
Calculate identity and similarity percentages for multiple sequence alignments.

```bash
python alignment_stats.py -i alignment.fasta -m ../Data/EBLOSUM62.txt
```

**Output:**
- Identical positions (all sequences have same residue)
- Similar positions (all pairwise BLOSUM scores positive)

### `protein_presentation.spt`
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

## Templates

### HTML with Jmol Applet
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

## Resources

- **Rosalind**: https://rosalind.info
- **Jmol Wiki**: https://wiki.jmol.org
- **BioPython**: https://biopython.org
- **BLOSUM Matrices**: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/
