# Lab Tools

Utility scripts and templates for bioinformatics work.

## Scripts

### `rosalind.py`
Python wrapper for working with Rosalind bioinformatics problems (rosalind.info).

#### Basic Usage
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

#### 🎓 AI-Assisted Learning Mode (Copilot Compatible)

The tutor module helps you **learn** without giving away answers:

```python
from Lab_Tools.rosalind import RosalindTutor, copilot_help, learning_mode

# Get learning guidance (concepts, not solutions)
print(learning_mode("GC"))

# Get progressive hints
tutor = RosalindTutor()
print(tutor.get_hint("GC", level=1))  # Vague hint
print(tutor.get_hint("GC", level=2))  # More specific
print(tutor.get_hint("GC", level=3))  # Most direct (still no answer!)

# Generate AI prompt for help (Socratic, no spoilers)
prompt = copilot_help("GC", stuck_on="parsing FASTA format")
print(prompt)  # Paste this to Copilot for guided learning

# Create a learning notebook with built-in scaffolds
tutor.create_tutor_notebook("DNA", output_dir="./learn")
```

#### AI Prompt Types

**1. Learning Assistant Prompt** - For when you're stuck
```python
prompt = tutor.copilot_prompt("PROT", stuck_on="codon translation")
```

**2. Code Review Prompt** - Socratic feedback on your code
```python
prompt = tutor.review_prompt("DNA", code=my_solution_code)
```

**3. Debug Assistance Prompt** - Teaches debugging skills
```python
prompt = tutor.debug_prompt("GC", code=my_code, error="IndexError...")
```

**4. Concept Explanation Prompt** - Understand underlying concepts
```python
prompt = tutor.explain_concept("dynamic programming")
```

#### Learning Path
```python
# Get a recommended sequence of problems
path = tutor.learning_path("DNA")
for p in path:
    print(f"{p['id']} ({'⭐'*p['difficulty']}) - {p['title']}")
```

#### Available Tracks
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
