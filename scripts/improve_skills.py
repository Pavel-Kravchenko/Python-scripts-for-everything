#!/usr/bin/env python3
"""Auto-fix common quality issues in skills.

Usage:
    python3 scripts/improve_skills.py --dry-run          # preview changes
    python3 scripts/improve_skills.py --apply             # apply all fixes
    python3 scripts/improve_skills.py --apply --fix version_compat,primary_tool
"""

import argparse
import csv
import os
import re
from pathlib import Path

SKILLS_DIR = Path(__file__).resolve().parent.parent / "Skills"
REPORTS_DIR = Path(__file__).resolve().parent.parent / "reports"

# ── Package detection ────────────────────────────────────────────────

IMPORT_TO_VERSION = {
    "numpy": "numpy 1.26+", "np": "numpy 1.26+",
    "pandas": "pandas 2.1+", "pd": "pandas 2.1+",
    "matplotlib": "matplotlib 3.8+", "plt": "matplotlib 3.8+",
    "seaborn": "seaborn 0.13+", "sns": "seaborn 0.13+",
    "scanpy": "scanpy 1.10+", "sc": "scanpy 1.10+",
    "anndata": "anndata 0.10+",
    "sklearn": "scikit-learn 1.4+",
    "scipy": "scipy 1.12+",
    "Bio": "biopython 1.83+",
    "networkx": "networkx 3.2+", "nx": "networkx 3.2+",
    "torch": "pytorch 2.2+",
    "transformers": "transformers 4.38+",
    "statsmodels": "statsmodels 0.14+",
    "lifelines": "lifelines 0.28+",
    "rdkit": "rdkit 2024.03+",
    "pysam": "pysam 0.22+",
    "pydeseq2": "pydeseq2 0.4+",
    "muon": "muon 0.1+",
    "episcanpy": "episcanpy 0.4+",
    "scvelo": "scvelo 0.3+",
    "cobra": "cobrapy 0.29+",
}

TOPIC_TO_TOOL = {
    "scanpy": "scanpy", "scrna": "scanpy", "single-cell": "scanpy",
    "seurat": "Seurat", "deseq": "DESeq2", "edger": "edgeR",
    "bismark": "Bismark", "star": "STAR", "salmon": "salmon",
    "samtools": "samtools", "bcftools": "bcftools", "gatk": "GATK",
    "bwa": "BWA", "bowtie": "Bowtie2", "mageck": "MAGeCK",
    "blast": "BLAST+", "macs": "MACS2/3", "picard": "Picard",
    "rdkit": "RDKit", "networkx": "NetworkX",
    "numpy": "NumPy", "pandas": "Pandas",
    "matplotlib": "Matplotlib", "pytorch": "PyTorch",
    "transformers": "HuggingFace Transformers",
    "cobra": "COBRApy", "kraken": "Kraken2",
    "metaphlan": "MetaPhlAn", "qiime": "QIIME2",
}

VERSION_COMPAT_TEMPLATE = """## Version Compatibility

Reference examples tested with: {packages}

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.
"""

# ── Fix functions ────────────────────────────────────────────────────

def detect_packages(content: str) -> list[str]:
    """Detect imported packages from code blocks."""
    found = set()
    for m in re.finditer(r"(?:import|from)\s+(\w+)", content):
        pkg = m.group(1)
        if pkg in IMPORT_TO_VERSION:
            found.add(IMPORT_TO_VERSION[pkg])
    return sorted(found) if found else ["Python 3.10+"]


def detect_primary_tool(name: str, content: str) -> str:
    """Detect primary tool from skill name and content."""
    name_lower = name.lower()
    for keyword, tool in TOPIC_TO_TOOL.items():
        if keyword in name_lower:
            return tool
    # Fallback: most-imported package
    import_counts = {}
    for m in re.finditer(r"(?:import|from)\s+(\w+)", content):
        pkg = m.group(1)
        if pkg in TOPIC_TO_TOOL:
            import_counts[pkg] = import_counts.get(pkg, 0) + 1
    if import_counts:
        top = max(import_counts, key=import_counts.get)
        return TOPIC_TO_TOOL[top]
    return "Python"


def fix_version_compat(content: str) -> str:
    """Add Version Compatibility section after frontmatter if missing."""
    if re.search(r"##\s*Version Compatibility", content):
        return content
    packages = detect_packages(content)
    section = VERSION_COMPAT_TEMPLATE.format(packages=", ".join(packages))
    # Insert after frontmatter
    m = re.match(r"(^---\n.*?\n---\n)", content, re.DOTALL)
    if m:
        return content[:m.end()] + "\n" + section + "\n" + content[m.end():]
    return section + "\n" + content


def fix_primary_tool(content: str, name: str) -> str:
    """Add primary_tool to frontmatter if missing."""
    m = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if not m:
        return content
    fm_text = m.group(1)
    if "primary_tool:" in fm_text:
        return content
    tool = detect_primary_tool(name, content)
    new_fm = fm_text + f"\nprimary_tool: {tool}"
    return content[:m.start(1)] + new_fm + content[m.end(1):]


def fix_duplicate_h1(content: str) -> str:
    """Remove duplicate H1 headers (keep the first occurrence)."""
    lines = content.split("\n")
    seen_h1 = None
    result = []
    for line in lines:
        if line.startswith("# ") and not line.startswith("## "):
            if seen_h1 is None:
                seen_h1 = line
                result.append(line)
            elif line == seen_h1:
                continue  # Skip duplicate
            else:
                result.append(line)  # Different H1, keep it
        else:
            result.append(line)
    return "\n".join(result)


def fix_add_frontmatter(content: str, name: str) -> str:
    """Add YAML frontmatter to files that lack it."""
    if content.startswith("---"):
        return content
    # Extract title and first paragraph for description
    title_m = re.search(r"^#\s+(.+)", content, re.MULTILINE)
    title = title_m.group(1).strip() if title_m else name
    # Find first substantive paragraph
    desc = title
    for para in re.split(r"\n\n+", content)[1:]:
        para = para.strip()
        if len(para) > 30 and not para.startswith("#") and not para.startswith("|"):
            desc = para.replace("\n", " ")[:200]
            break
    tool_type = "python"
    if re.search(r"```r\b", content):
        tool_type = "mixed"
    elif re.search(r"```bash", content):
        tool_type = "cli"
    tool = detect_primary_tool(name, content)
    fm = f"""---
name: {name}
description: "{desc}"
tool_type: {tool_type}
primary_tool: {tool}
---

"""
    return fm + content


# ── Orchestrator ─────────────────────────────────────────────────────

ALL_FIXES = {
    "version_compat": ("Add Version Compatibility section", fix_version_compat),
    "primary_tool": ("Add primary_tool to frontmatter", None),  # needs name
    "duplicate_h1": ("Remove duplicate H1 headers", fix_duplicate_h1),
    "add_frontmatter": ("Add missing frontmatter", None),  # needs name
}

def apply_fixes(skills_dir: Path, fixes: list[str], dry_run: bool) -> dict:
    stats = {"checked": 0, "modified": 0, "fixes_applied": 0}
    details = []

    for item in sorted(os.listdir(skills_dir)):
        p = skills_dir / item
        # Find the skill file
        if p.is_dir() and (p / "SKILL.md").exists():
            skill_path = p / "SKILL.md"
            name = item
        elif p.is_file() and p.suffix == ".md" and item != "README.md":
            skill_path = p
            name = item.replace(".md", "")
        else:
            continue

        stats["checked"] += 1
        original = skill_path.read_text(errors="replace")
        content = original
        applied = []

        if "add_frontmatter" in fixes and not content.startswith("---"):
            content = fix_add_frontmatter(content, name)
            applied.append("add_frontmatter")

        if "version_compat" in fixes:
            new = fix_version_compat(content)
            if new != content:
                content = new
                applied.append("version_compat")

        if "primary_tool" in fixes:
            new = fix_primary_tool(content, name)
            if new != content:
                content = new
                applied.append("primary_tool")

        if "duplicate_h1" in fixes:
            new = fix_duplicate_h1(content)
            if new != content:
                content = new
                applied.append("duplicate_h1")

        if applied:
            stats["modified"] += 1
            stats["fixes_applied"] += len(applied)
            if dry_run:
                details.append(f"  [DRY-RUN] {name}: {', '.join(applied)}")
            else:
                skill_path.write_text(content)
                details.append(f"  [FIXED] {name}: {', '.join(applied)}")

    return stats, details


def main():
    ap = argparse.ArgumentParser(description="Auto-fix common skill quality issues")
    ap.add_argument("--skills-dir", default=str(SKILLS_DIR))
    ap.add_argument("--dry-run", action="store_true", help="Preview fixes without applying")
    ap.add_argument("--apply", action="store_true", help="Apply fixes")
    ap.add_argument("--fix", default="version_compat,primary_tool,duplicate_h1,add_frontmatter",
                    help="Comma-separated list of fixes to apply")
    args = ap.parse_args()

    if not args.dry_run and not args.apply:
        print("Specify --dry-run or --apply")
        return

    fixes = [f.strip() for f in args.fix.split(",")]
    stats, details = apply_fixes(Path(args.skills_dir), fixes, dry_run=args.dry_run)

    print(f"\n{'='*60}")
    mode = "DRY RUN" if args.dry_run else "APPLIED"
    print(f"  SKILL IMPROVEMENT — {mode}")
    print(f"{'='*60}")
    print(f"  Skills checked:  {stats['checked']}")
    print(f"  Skills modified: {stats['modified']}")
    print(f"  Fixes applied:   {stats['fixes_applied']}")
    print()
    for d in details[:30]:
        print(d)
    if len(details) > 30:
        print(f"  ... and {len(details) - 30} more")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
