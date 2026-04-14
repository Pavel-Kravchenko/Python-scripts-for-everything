#!/usr/bin/env python3
"""Compress skills from course-notebook dumps into practical reference cards.

Aggressively strips tutorial prose, demo code, boilerplate sections, and
course metadata. Keeps: pitfalls, real code patterns, API tables, recipes.

Usage:
    python3 scripts/compress_skills.py                    # apply
    python3 scripts/compress_skills.py --dry-run          # preview stats
"""

import argparse
import os
import re
from pathlib import Path

SKILLS_DIR = Path(__file__).resolve().parent.parent / "Skills"

# ── Section-level removal (H2/H3 heading text, case-insensitive) ────

REMOVE_SECTION_PATTERNS = [
    r"version\s+compatibility",
    r"learning\s+objectives?",
    r"how\s+to\s+use\s+this\s+(?:notebook|module)",
    r"why\s+(?:this\s+)?(?:notebook|module|single.cell|rna.seq|self.balancing)\s+matters",
    r"why\s+\w+\??\s*$",                         # "Why RNA-seq?" etc
    r"environment\s+setup",
    r"what\s+you(?:'ll|\s+will)\s+learn",
    r"prerequisites?\s*$",
    r"notebook\s+setup",
    r"setup\s+and\s+imports?\s*$",
    r"exercises?\s*$",                             # Exercise sections
    r"practice\s+exercises?\s*$",
    r"summary\s+and\s+next\s+steps?\s*$",
    r"next\s+steps?\s*$",
    r"what(?:'s)?\s+next\s*$",
    r"recap\s*$",
    r"review\s+questions?\s*$",
    r"further\s+reading\s*$",
    r"additional\s+resources?\s*$",
    r"wrap.?up\s*$",
    r"conclusion\s*$",
]

# ── Line-level removal patterns ─────────────────────────────────────

REMOVE_LINE_PATTERNS = [
    r"^\*Source:\s+Course\s+notebook",
    r"^\*\*Navigation:\*\*",
    r"^Split\s+from\s+`.*`\s+to\s+keep",
    r"^\*\*Prerequisites?:\*\*",
    r"^\*\*Estimated\s+time:\*\*",
    r"^\*\*Tier\s+\d+\s+[—–-]",
    r"^Tier\s+\d+\s+[—–-]",
    r"^---\s*$",
    r"^\*Prerequisites?:",
    r"^By\s+the\s+end\s+of\s+this\s+(?:notebook|module)",
    r"^\*\*By\s+the\s+end\s+of\s+this\s+(?:notebook|module)",
    r"^\*\*Key\s+resources?:\*\*\s*$",
]

# ── Tutorial paragraph detection ─────────────────────────────────────

TUTORIAL_PARA_PATTERNS = [
    # "A variable is a name that..." / "The GC content is..."
    r"^(?:A|An|The)\s+\*?\*?\w+\*?\*?\s+(?:is|are|was|were|refers?|represents?|means?)\s",
    # "In Python, ..." / "In this section..."
    r"^In\s+(?:Python|this|general|practice|bioinformatics|the)",
    # "Think of it as..."
    r"^Think\s+of\s+",
    # "Every/Each value in..."
    r"^(?:Every|Each)\s+\w+\s+(?:in|has|is)",
    # "Before we..." / "Now that you..." / "Once you..."
    r"^(?:Before|Now\s+that|Once)\s+(?:we|you)",
    # "Let's..." / "Let us..."
    r"^Let(?:'s|\s+us)\s+",
    # "This is the/a..." / "This section..."
    r"^This\s+(?:is|section|module|notebook|chapter|approach)",
    # "Remember that..." / "Note that..."
    r"^(?:Remember|Recall|Note)\s+that\s+",
    # "We will..." / "We can..."/ "We need..."
    r"^We\s+(?:will|can|need|have|use|start|begin|now|first)",
    # "You will..." / "You can..."
    r"^You\s+(?:will|can|should|need|may)\s+(?:often|also|use|see|notice|find|need|want|have|learn)",
    # "Here we..." / "Here's..."
    r"^Here(?:\s+we|'s)\s+",
    # "For example,..." used as filler
    r"^For\s+example,\s+(?:if|when|consider|suppose|imagine)",
    # Contains "this notebook" or "this module"
    r"this\s+(?:notebook|module|exercise|lab|tutorial)",
    # "Good/Bad variable names" demo labels
    r"^#\s+(?:Good|Bad)\s+\w+\s+names",
    # "Run cells..." / "Try running..."
    r"^(?:Run|Try)\s+(?:cells?|running|the\s+cell)",
]

# ── Frontmatter ──────────────────────────────────────────────────────

REMOVE_FM_KEYS = {"source_notebook"}


def parse_frontmatter(content: str) -> tuple[dict, str]:
    m = re.match(r"^---\n(.*?)\n---\n?", content, re.DOTALL)
    if not m:
        return {}, content
    fm = {}
    for line in m.group(1).splitlines():
        if ":" in line:
            key, _, val = line.partition(":")
            fm[key.strip()] = val.strip().strip('"').strip("'")
    return fm, content[m.end():]


def serialize_frontmatter(fm: dict) -> str:
    lines = ["---"]
    for k, v in fm.items():
        if any(c in v for c in ':"{}[]'):
            lines.append(f'{k}: "{v}"')
        else:
            lines.append(f"{k}: {v}")
    lines.append("---")
    return "\n".join(lines)


# ── Description improvement ──────────────────────────────────────────

def improve_description(fm: dict, body: str, name: str) -> str:
    desc = fm.get("description", "")

    is_course_desc = any(p in desc.lower() for p in [
        "tier ", "module ", "by the end of", "learning objective",
        "you will be able to", "split from", "1. ", "2. ", "3. ",
        "reference examples tested",
    ])
    # Also bad: description is just the skill name repeated
    if desc.lower().replace("-", " ").replace("_", " ") == name.replace("-", " "):
        is_course_desc = True

    if not is_course_desc and len(desc) > 30 and "notebook" not in desc.lower():
        return desc

    # Extract H1 title (before any stripping)
    h1_match = re.search(r"^#\s+(.+)", body, re.MULTILINE)
    h1 = h1_match.group(1).strip() if h1_match else ""
    # Clean emoji and module prefix
    h1 = re.sub(r"[^\w\s:,/().+—–\-']", "", h1).strip()
    h1 = re.sub(r"^Module\s+\d+\.?\d*[a-z]?:\s*", "", h1).strip()

    # Find first non-boilerplate paragraph
    first_para = ""
    in_code = False
    skip_phrases = [
        "Source:", "**Tier", "Navigation", "Split from", "Prerequisites",
        "Reference examples tested with", "Before using code patterns",
        "If code throws ImportError",
    ]
    for line in body.splitlines():
        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("|") or s.startswith("-"):
            continue
        if any(skip in s for skip in skip_phrases):
            continue
        if any(re.search(p, s, re.IGNORECASE) for p in TUTORIAL_PARA_PATTERNS[:5]):
            continue
        if len(s) > 40:
            first_para = re.sub(r"\*\*?|`", "", s)
            break

    tool = fm.get("primary_tool", "")
    tool_suffix = f" with {tool}" if tool and tool not in ("Python", "") else ""

    if h1 and len(h1) > 15:
        desc = h1 + tool_suffix
    elif first_para:
        desc = first_para[:200]
    else:
        desc = name.replace("-", " ").replace("bio ", "bioinformatics ").title() + tool_suffix

    if len(desc) > 200:
        desc = desc[:197] + "..."
    return desc


# ── Body transformations ─────────────────────────────────────────────

def remove_sections(body: str) -> str:
    """Remove entire H2/H3 sections matching patterns."""
    lines = body.splitlines()
    result = []
    skipping = False
    skip_level = 0

    for line in lines:
        heading_match = re.match(r"^(#{2,3})\s+(.+)", line)
        if heading_match:
            level = len(heading_match.group(1))
            text = heading_match.group(2).strip()
            should_remove = any(
                re.search(pat, text, re.IGNORECASE)
                for pat in REMOVE_SECTION_PATTERNS
            )
            if should_remove:
                skipping = True
                skip_level = level
                continue
            elif skipping and level <= skip_level:
                skipping = False

        if line.startswith("# ") and not line.startswith("## "):
            skipping = False

        if not skipping:
            result.append(line)

    return "\n".join(result)


def remove_lines(body: str) -> str:
    """Remove individual matching lines + objective lists."""
    lines = body.splitlines()
    result = []
    in_objectives_list = False

    for line in lines:
        s = line.strip()

        if any(re.search(pat, s, re.IGNORECASE) for pat in REMOVE_LINE_PATTERNS):
            if re.match(r"^(?:\*\*)?By\s+the\s+end", s, re.IGNORECASE):
                in_objectives_list = True
            continue

        if in_objectives_list:
            if re.match(r"^\d+\.\s", s) or s == "":
                continue
            in_objectives_list = False

        result.append(line)

    return "\n".join(result)


def remove_tutorial_paragraphs(body: str) -> str:
    """Remove paragraphs that are tutorial explanations (not reference material)."""
    lines = body.splitlines()
    result = []
    in_code = False

    i = 0
    while i < len(lines):
        line = lines[i]

        # Track code blocks — never touch content inside them
        if line.strip().startswith("```"):
            in_code = not in_code
            result.append(line)
            i += 1
            continue

        if in_code:
            result.append(line)
            i += 1
            continue

        s = line.strip()

        # Keep: headings, blank lines, tables, bullet lists
        if not s or s.startswith("#") or s.startswith("|") or s.startswith("-") or s.startswith("*"):
            # But check bullet lists for tutorial content too
            if s.startswith("- ") or s.startswith("* "):
                # Keep bullet lists in pitfall sections (they contain **)
                result.append(line)
            elif s.startswith("*") and not s.startswith("**"):
                # Italic text — often "Source:" or navigation, skip these
                if "source" in s.lower() or "navigation" in s.lower():
                    i += 1
                    continue
                result.append(line)
            else:
                result.append(line)
            i += 1
            continue

        # Check if this line starts a tutorial paragraph
        is_tutorial = any(
            re.search(pat, s, re.IGNORECASE)
            for pat in TUTORIAL_PARA_PATTERNS
        )

        if is_tutorial:
            # Skip this line and any continuation lines (non-blank, non-heading, non-code)
            i += 1
            while i < len(lines):
                next_s = lines[i].strip()
                if not next_s or next_s.startswith("#") or next_s.startswith("```"):
                    break
                if next_s.startswith("- ") or next_s.startswith("* ") or next_s.startswith("|"):
                    break
                i += 1
            continue

        result.append(line)
        i += 1

    return "\n".join(result)


def remove_demo_code_blocks(body: str) -> str:
    """Remove code blocks that only contain print()/assignment demos."""
    lines = body.splitlines()
    result = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Detect start of code block
        if line.strip().startswith("```"):
            # Collect the entire code block
            block = [line]
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                block.append(lines[i])
                i += 1
            if i < len(lines):
                block.append(lines[i])  # closing ```
                i += 1

            # Analyze block content (skip opening/closing ```)
            code_lines = [l.strip() for l in block[1:-1] if l.strip() and not l.strip().startswith("#")]

            if not code_lines:
                # Empty code block — skip
                continue

            # Check if this is a "diagram" block (ASCII art, not real code)
            is_diagram = all(
                not any(kw in l for kw in ["import", "def ", "class ", "return", "for ", "if ", "while ", "with ", "try:", "except", "raise", "yield", "async ", "await "])
                and not re.match(r"^\w+\s*[=(]", l)
                for l in code_lines
            ) and any(
                any(c in l for c in ["──", "→", "---", "...", "^^^", "├", "└", "│"])
                for l in code_lines
            )

            if is_diagram and len(code_lines) < 15:
                continue

            # Check if purely demo: only print(), simple assignments, comments
            is_demo = len(code_lines) <= 12 and all(
                l.startswith("print(") or l.startswith("print (") or
                re.match(r"^\w+\s*=\s*(?:\d|['\"]|True|False|None|\[|\(|{)", l) or
                l.startswith("# ") or l == "" or
                re.match(r"^(?:type|len|id)\(", l) or
                re.match(r"^f?['\"]", l) or
                l.startswith("for v in")
                for l in code_lines
            )

            # Don't remove if it contains f-string formatting patterns useful as reference
            has_fstring_ref = any(":.1f" in l or ":.2e" in l or ":," in l for l in code_lines)

            if is_demo and not has_fstring_ref:
                continue

            # Keep the block
            result.extend(block)
        else:
            result.append(line)
            i += 1

    return "\n".join(result)


def clean_h1(body: str) -> str:
    """Clean H1 headers: remove emoji, module prefix, dedup."""
    lines = body.splitlines()
    result = []
    seen_h1 = False

    for line in lines:
        if line.startswith("# ") and not line.startswith("## "):
            # Remove emoji (keep alphanumeric, whitespace, punctuation)
            cleaned = re.sub(r"[^\w\s:,/().+—–\-#*`']", "", line).strip()
            # Remove "Module X.Y:" prefix
            cleaned = re.sub(r"^#\s+Module\s+\d+\.?\d*[a-z]?:\s*", "# ", cleaned)
            cleaned = re.sub(r"#\s+[:\-—–]+\s*", "# ", cleaned)

            if not seen_h1:
                seen_h1 = True
                result.append(cleaned)
            else:
                # Second H1 — only keep if substantially different from first
                if cleaned != result[0] if result else True:
                    result.append(cleaned)
        else:
            result.append(line)

    return "\n".join(result)


def strip_section_numbers(body: str) -> str:
    """Remove numbering from section headers: '## 1. Variables' → '## Variables'."""
    return re.sub(
        r"^(#{2,4})\s+\d+\.?\d*\.?\s+",
        r"\1 ",
        body,
        flags=re.MULTILINE,
    )


def rename_sections(body: str) -> str:
    renames = {
        r"common\s+stumbling\s+points": "Pitfalls",
        r"common\s+(?:\w+\s+)?pitfalls": "Pitfalls",
        r"key\s+resources?": "References",
        r"common\s+mistakes?\s+and\s+gotchas?": "Pitfalls",
    }
    lines = body.splitlines()
    result = []
    for line in lines:
        m = re.match(r"^(#{2,3})\s+(.+)", line)
        if m:
            level, text = m.group(1), m.group(2).strip()
            for pat, replacement in renames.items():
                if re.search(pat, text, re.IGNORECASE):
                    line = f"{level} {replacement}"
                    break
        result.append(line)
    return "\n".join(result)


def remove_key_resources(body: str) -> str:
    """Remove 'Key resources:' inline sections with URL lists."""
    lines = body.splitlines()
    result = []
    in_url_list = False

    for line in lines:
        s = line.strip()
        if re.match(r"^\*\*Key\s+resources?:\*\*\s*$", s, re.IGNORECASE):
            in_url_list = True
            continue
        if in_url_list:
            if s.startswith("- [") or s.startswith("- *") or s.startswith("- http") or s == "":
                continue
            in_url_list = False
        result.append(line)

    return "\n".join(result)


def collapse_whitespace(body: str) -> str:
    body = re.sub(r"\n{4,}", "\n\n\n", body)
    # Also collapse: heading followed by 2+ blank lines
    body = re.sub(r"(^#{1,4}\s+.+\n)\n{2,}", r"\1\n", body, flags=re.MULTILINE)
    lines = [l.rstrip() for l in body.splitlines()]
    while lines and lines[0] == "":
        lines.pop(0)
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def remove_duplicate_h1(body: str) -> str:
    """If there are multiple H1s and they're similar, keep only the first."""
    lines = body.splitlines()
    h1_indices = [i for i, l in enumerate(lines) if l.startswith("# ") and not l.startswith("## ")]
    if len(h1_indices) <= 1:
        return body

    # Keep the first H1, remove subsequent ones that are just variations
    first_h1 = lines[h1_indices[0]].lower().replace(" ", "")
    to_remove = set()
    for idx in h1_indices[1:]:
        other = lines[idx].lower().replace(" ", "")
        # If >50% word overlap, it's a duplicate
        first_words = set(re.findall(r"\w+", first_h1))
        other_words = set(re.findall(r"\w+", lines[idx].lower()))
        if first_words and other_words:
            overlap = len(first_words & other_words) / max(len(first_words), len(other_words))
            if overlap > 0.4:
                to_remove.add(idx)

    return "\n".join(l for i, l in enumerate(lines) if i not in to_remove)


# ── Main pipeline ────────────────────────────────────────────────────

def compress_skill(content: str, name: str) -> str:
    fm, body = parse_frontmatter(content)

    for key in REMOVE_FM_KEYS:
        fm.pop(key, None)

    # Improve description BEFORE stripping body (so we have full content to extract from)
    fm["description"] = improve_description(fm, body, name)

    # Apply transformations in order (most structural first)
    body = remove_sections(body)
    body = remove_lines(body)
    body = remove_key_resources(body)
    body = remove_tutorial_paragraphs(body)
    body = remove_demo_code_blocks(body)
    body = clean_h1(body)
    body = remove_duplicate_h1(body)
    body = strip_section_numbers(body)
    body = rename_sections(body)
    body = collapse_whitespace(body)

    return serialize_frontmatter(fm) + "\n\n" + body


# ── Orchestrator ─────────────────────────────────────────────────────

def compress_all(skills_dir: Path, dry_run: bool) -> dict:
    stats = {"checked": 0, "compressed": 0, "lines_before": 0, "lines_after": 0, "details": []}

    for item in sorted(os.listdir(skills_dir)):
        p = skills_dir / item
        skill_md = p / "SKILL.md"
        if not (p.is_dir() and skill_md.exists()):
            continue

        stats["checked"] += 1
        original = skill_md.read_text(errors="replace")
        before = original.count("\n")
        stats["lines_before"] += before

        compressed = compress_skill(original, item)
        after = compressed.count("\n")
        stats["lines_after"] += after

        reduction = before - after
        if reduction > 0:
            stats["compressed"] += 1
            pct = (reduction / before * 100) if before > 0 else 0
            stats["details"].append((item, before, after, reduction, pct))
            if not dry_run:
                skill_md.write_text(compressed)

    return stats


def main():
    ap = argparse.ArgumentParser(description="Compress skills into practical reference cards")
    ap.add_argument("--skills-dir", default=str(SKILLS_DIR))
    ap.add_argument("--dry-run", action="store_true", help="Preview without writing")
    args = ap.parse_args()

    stats = compress_all(Path(args.skills_dir), dry_run=args.dry_run)
    mode = "DRY RUN" if args.dry_run else "APPLIED"
    reduction = stats["lines_before"] - stats["lines_after"]
    pct = (reduction / stats["lines_before"] * 100) if stats["lines_before"] > 0 else 0

    print(f"\n{'='*60}")
    print(f"  SKILL COMPRESSION — {mode}")
    print(f"{'='*60}")
    print(f"  Skills checked:    {stats['checked']}")
    print(f"  Skills compressed: {stats['compressed']}")
    print(f"  Lines before:      {stats['lines_before']:,}")
    print(f"  Lines after:       {stats['lines_after']:,}")
    print(f"  Lines removed:     {reduction:,} (-{pct:.0f}%)")
    print()

    # Top 20 reductions
    sorted_d = sorted(stats["details"], key=lambda x: x[3], reverse=True)
    for name, b, a, r, p in sorted_d[:20]:
        print(f"  {name}: {b} → {a} (-{r}, -{p:.0f}%)")
    if len(sorted_d) > 20:
        print(f"  ... and {len(sorted_d) - 20} more")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
