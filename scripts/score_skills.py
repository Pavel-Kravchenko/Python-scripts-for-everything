#!/usr/bin/env python3
"""Score all skills in Skills/ against a quality rubric (0-100).

Usage:
    python3 scripts/score_skills.py                    # default paths
    python3 scripts/score_skills.py --output report.csv
"""

import argparse
import csv
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

SKILLS_DIR = Path(__file__).resolve().parent.parent / "Skills"
GOLD_DIR = Path.home() / ".claude" / "skills"
REPORTS_DIR = Path(__file__).resolve().parent.parent / "reports"

# ── Frontmatter parsing ─────────────────────────────────────────────

def parse_frontmatter(content: str) -> dict:
    m = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if not m:
        return {"_has_frontmatter": False}
    fm = {"_has_frontmatter": True}
    for line in m.group(1).splitlines():
        if ":" in line:
            key, _, val = line.partition(":")
            fm[key.strip()] = val.strip().strip('"').strip("'")
    return fm

def body_after_frontmatter(content: str) -> str:
    m = re.match(r"^---\n.*?\n---\n?", content, re.DOTALL)
    return content[m.end():] if m else content

# ── Dataclass ────────────────────────────────────────────────────────

@dataclass
class SkillInfo:
    name: str
    path: str
    skill_type: str          # "directory" or "flat"
    frontmatter: dict = field(default_factory=dict)
    content: str = ""
    body: str = ""
    line_count: int = 0
    code_blocks: dict = field(default_factory=dict)
    h2_sections: list = field(default_factory=list)

# ── Discovery ────────────────────────────────────────────────────────

def discover_skills(skills_dir: Path) -> list[SkillInfo]:
    skills = []
    for item in sorted(os.listdir(skills_dir)):
        p = skills_dir / item
        # Directory-based
        skill_md = p / "SKILL.md"
        if p.is_dir() and skill_md.exists():
            content = skill_md.read_text(errors="replace")
            fm = parse_frontmatter(content)
            body = body_after_frontmatter(content)
            skills.append(SkillInfo(
                name=item, path=str(skill_md), skill_type="directory",
                frontmatter=fm, content=content, body=body,
                line_count=content.count("\n"),
                code_blocks=count_code_blocks(content),
                h2_sections=extract_h2_sections(body),
            ))
        # Flat .md
        elif p.is_file() and p.suffix == ".md" and item != "README.md":
            content = p.read_text(errors="replace")
            fm = parse_frontmatter(content)
            body = body_after_frontmatter(content)
            skills.append(SkillInfo(
                name=item.replace(".md", ""), path=str(p), skill_type="flat",
                frontmatter=fm, content=content, body=body,
                line_count=content.count("\n"),
                code_blocks=count_code_blocks(content),
                h2_sections=extract_h2_sections(body),
            ))
    return skills

def discover_gold_names(gold_dir: Path) -> set[str]:
    if not gold_dir.exists():
        return set()
    return {d for d in os.listdir(gold_dir) if (gold_dir / d).is_dir()}

# ── Helpers ──────────────────────────────────────────────────────────

def count_code_blocks(content: str) -> dict:
    counts = {"python": 0, "r": 0, "bash": 0, "other": 0, "untagged": 0}
    for m in re.finditer(r"```(\w*)", content):
        lang = m.group(1).lower()
        if lang in ("python", "py"):
            counts["python"] += 1
        elif lang in ("r",):
            counts["r"] += 1
        elif lang in ("bash", "sh", "shell", "zsh"):
            counts["bash"] += 1
        elif lang == "":
            counts["untagged"] += 1
        else:
            counts["other"] += 1
    return counts

def extract_h2_sections(body: str) -> list[tuple[str, int]]:
    """Return list of (section_name, content_line_count)."""
    sections = []
    lines = body.splitlines()
    current_name = None
    current_lines = 0
    for line in lines:
        if line.startswith("## "):
            if current_name is not None:
                sections.append((current_name, current_lines))
            current_name = line[3:].strip()
            current_lines = 0
        elif current_name is not None:
            if line.strip():
                current_lines += 1
    if current_name is not None:
        sections.append((current_name, current_lines))
    return sections

def total_code_block_count(cb: dict) -> int:
    return sum(cb.values())

# ── Scoring dimensions ───────────────────────────────────────────────

def score_metadata(s: SkillInfo) -> tuple[int, list[str]]:
    """20 points max."""
    score = 0
    issues = []
    fm = s.frontmatter

    if fm.get("_has_frontmatter"):
        score += 4
    else:
        issues.append("no_frontmatter")
        return score, issues

    if fm.get("name"):
        score += 3
    else:
        issues.append("no_name")

    desc = fm.get("description", "")
    if desc and len(desc) >= 50:
        score += 2
    elif desc:
        score += 1
        issues.append("short_description")
    else:
        issues.append("no_description")

    # Description quality: use-case language
    if re.search(r"Use (when|for)", desc, re.I):
        score += 4
    elif len(desc) >= 100 and not re.search(r"^(\*\*)?Tier\s+\d|^Split from|^By the end", desc):
        score += 2
    else:
        issues.append("bad_description")

    if fm.get("tool_type"):
        score += 2
    else:
        issues.append("no_tool_type")

    if fm.get("primary_tool"):
        score += 3
    else:
        issues.append("no_primary_tool")

    # Cap at 20 (extra from partial scoring)
    return min(score, 20), issues


def score_structure(s: SkillInfo) -> tuple[int, list[str]]:
    """20 points max."""
    score = 0
    issues = []
    body = s.body

    # Version Compatibility section (5 pts)
    if re.search(r"##\s*Version Compatibility", body):
        score += 5
    else:
        issues.append("no_version_compat")

    # Directory format (3 pts)
    if s.skill_type == "directory":
        score += 3
    else:
        issues.append("flat_format")

    # Duplicate H1 (2 pts — lose if duplicate exists)
    h1_matches = re.findall(r"^# .+", body, re.MULTILINE)
    if len(h1_matches) <= 1:
        score += 2
    else:
        issues.append("duplicate_h1")

    # Learning objectives (3 pts)
    if re.search(r"(Learning Objectives|you will be able to|When to [Uu]se)", body):
        score += 3
    else:
        issues.append("no_learning_objectives")

    # Prerequisites (2 pts)
    if re.search(r"[Pp]rerequisites?", body):
        score += 2
    else:
        issues.append("no_prerequisites")

    # Pitfalls / gotchas (3 pts)
    if re.search(r"(pitfall|gotcha|sticking point|common.*(mistake|error|issue))", body, re.I):
        score += 3
    else:
        issues.append("no_pitfalls")

    # Related / See Also (2 pts)
    if re.search(r"(Related|See [Aa]lso|Next [Ss]tep|Further)", body):
        score += 2

    return min(score, 20), issues


def score_code_quality(s: SkillInfo) -> tuple[int, list[str]]:
    """25 points max."""
    score = 0
    issues = []
    cb = s.code_blocks
    total_cb = total_code_block_count(cb)
    body = s.body

    # At least 3 code blocks (5 pts)
    if total_cb >= 5:
        score += 5
    elif total_cb >= 3:
        score += 3
    elif total_cb >= 1:
        score += 1
    else:
        issues.append("no_code_blocks")

    # Goal → Approach → Code pattern (8 pts)
    if re.search(r"\*\*Goal[:\*]", body) and re.search(r"\*\*Approach[:\*]", body):
        score += 8
    elif re.search(r"\*\*Goal[:\*]", body) or re.search(r"\*\*Approach[:\*]", body):
        score += 3
    else:
        issues.append("no_goal_approach")

    # Language specifiers on code blocks (3 pts)
    if cb.get("untagged", 0) == 0 and total_cb > 0:
        score += 3
    elif total_cb > 0 and cb.get("untagged", 0) < total_cb // 2:
        score += 1
    else:
        if total_cb > 0:
            issues.append("untagged_code_blocks")

    # R coverage for topics that should have it (4 pts)
    r_topics = r"(differential.expression|single.cell|statistic|biostatistic|regression|hypothesis|deseq|edger|seurat|ggplot|r.fundamental)"
    needs_r = bool(re.search(r_topics, s.name, re.I))
    if needs_r:
        if cb.get("r", 0) > 0:
            score += 4
        else:
            issues.append("needs_r_coverage")
    else:
        score += 4  # N/A — full marks

    # Docstrings or comments in code (3 pts)
    if re.search(r'"""', body) or re.search(r"'''", body):
        score += 3
    elif re.search(r"^#\s", body, re.MULTILINE):
        score += 2

    # At least one function def (2 pts)
    if re.search(r"def\s+\w+\(", body):
        score += 2
    else:
        issues.append("no_function_defs")

    return min(score, 25), issues


def score_content_depth(s: SkillInfo) -> tuple[int, list[str]]:
    """20 points max."""
    score = 0
    issues = []

    # Minimum length >= 100 lines (4 pts)
    if s.line_count >= 200:
        score += 4
    elif s.line_count >= 100:
        score += 3
    elif s.line_count >= 50:
        score += 1
    else:
        issues.append("too_short")

    # Not a stub (5 pts) — detect placeholder patterns
    stub_patterns = r"(> .*description|> .*TODO|> .*placeholder|> .*coming soon)"
    if re.search(stub_patterns, s.body, re.I):
        issues.append("is_stub")
    else:
        score += 5

    # Reference tables (3 pts)
    table_rows = len(re.findall(r"^\|.+\|", s.body, re.MULTILINE))
    if table_rows >= 5:
        score += 3
    elif table_rows >= 2:
        score += 1

    # Complexity / performance info for algo skills (3 pts)
    is_algo = s.name.startswith("algo-")
    if is_algo:
        if re.search(r"O\(", s.body):
            score += 3
    else:
        score += 3  # N/A

    # Multi-section depth (5 pts)
    substantive = [sec for sec in s.h2_sections if sec[1] >= 3]
    if len(substantive) >= 5:
        score += 5
    elif len(substantive) >= 3:
        score += 3
    elif len(substantive) >= 1:
        score += 1
    else:
        issues.append("shallow_sections")

    return min(score, 20), issues


def score_uniqueness(s: SkillInfo, all_skills: list[SkillInfo], gold_names: set[str]) -> tuple[int, list[str]]:
    """15 points max."""
    score = 0
    issues = []
    duplicate_of = ""

    # Format duplicate: old .md ↔ new directory (8 pts)
    if s.skill_type == "flat":
        # Check if a directory version exists
        for other in all_skills:
            if other.skill_type == "directory" and _names_overlap(s.name, other.name):
                issues.append("duplicate_of_directory")
                duplicate_of = other.name
                break
        else:
            score += 8
    else:
        score += 8

    # Near-duplicate of gold standard (4 pts)
    if s.name in gold_names:
        # Exact match is OK — it IS the installed version
        score += 4
    else:
        # Check for fuzzy match
        matched = False
        for gn in gold_names:
            if _names_overlap(s.name, gn) and s.name != gn:
                issues.append("near_gold_duplicate")
                duplicate_of = duplicate_of or f"gold:{gn}"
                matched = True
                break
        if not matched:
            score += 4

    # Fragment check (3 pts)
    desc = s.frontmatter.get("description", "")
    if "Split from" in desc or "split from" in desc:
        if s.line_count < 100:
            issues.append("is_fragment")
        else:
            score += 3
    else:
        score += 3

    return min(score, 15), issues, duplicate_of


def _names_overlap(a: str, b: str) -> bool:
    """Check if two skill names refer to the same topic."""
    stop = {"bio", "applied", "core", "python", "foundations", "algo", "ai", "science", "analysis", "for", "and", "the"}
    wa = {w for w in re.split(r"[-_]", a.lower()) if w and w not in stop and len(w) > 2}
    wb = {w for w in re.split(r"[-_]", b.lower()) if w and w not in stop and len(w) > 2}
    if not wa or not wb:
        return False
    overlap = len(wa & wb)
    return overlap >= 2 or (overlap == 1 and min(len(wa), len(wb)) == 1)

# ── Grading ──────────────────────────────────────────────────────────

def grade_and_recommend(total: int, issues: list[str]) -> tuple[str, str]:
    if "duplicate_of_directory" in issues:
        return ("D", "DELETE_DUPLICATE")

    if total >= 85:
        g = "A"
    elif total >= 70:
        g = "B"
    elif total >= 50:
        g = "C"
    elif total >= 30:
        g = "D"
    else:
        g = "F"

    if g in ("A", "B"):
        rec = "KEEP" if g == "A" else "FIX_MINOR"
    elif g == "C":
        rec = "FIX_MAJOR"
    elif "is_stub" in issues or "too_short" in issues:
        rec = "DELETE_LOW_QUALITY"
    else:
        rec = "REWRITE"

    return g, rec

# ── Main ─────────────────────────────────────────────────────────────

def score_all(skills_dir: Path, gold_dir: Path) -> list[dict]:
    skills = discover_skills(skills_dir)
    gold_names = discover_gold_names(gold_dir)

    rows = []
    for s in skills:
        meta_score, meta_issues = score_metadata(s)
        struct_score, struct_issues = score_structure(s)
        code_score, code_issues = score_code_quality(s)
        depth_score, depth_issues = score_content_depth(s)
        uniq_score, uniq_issues, dup_of = score_uniqueness(s, skills, gold_names)

        total = meta_score + struct_score + code_score + depth_score + uniq_score
        all_issues = meta_issues + struct_issues + code_issues + depth_issues + uniq_issues
        grade, rec = grade_and_recommend(total, all_issues)

        rows.append({
            "name": s.name,
            "type": s.skill_type,
            "path": s.path,
            "lines": s.line_count,
            "code_blocks": total_code_block_count(s.code_blocks),
            "metadata": meta_score,
            "structure": struct_score,
            "code_quality": code_score,
            "content_depth": depth_score,
            "uniqueness": uniq_score,
            "total": total,
            "grade": grade,
            "recommendation": rec,
            "duplicate_of": dup_of,
            "issues": "|".join(all_issues),
        })

    rows.sort(key=lambda r: r["total"])
    return rows


def write_csv(rows: list[dict], output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def print_summary(rows: list[dict]):
    grades = {}
    recs = {}
    for r in rows:
        grades[r["grade"]] = grades.get(r["grade"], 0) + 1
        recs[r["recommendation"]] = recs.get(r["recommendation"], 0) + 1

    total_score = sum(r["total"] for r in rows)
    avg = total_score / len(rows) if rows else 0

    print(f"\n{'='*60}")
    print(f"  SKILL QUALITY REPORT — {len(rows)} skills scored")
    print(f"{'='*60}")
    print(f"  Average score: {avg:.1f}/100")
    print()
    print("  Grade distribution:")
    for g in "ABCDF":
        n = grades.get(g, 0)
        bar = "#" * (n // 2)
        print(f"    {g}: {n:3d}  {bar}")
    print()
    print("  Recommendations:")
    for rec in ["KEEP", "FIX_MINOR", "FIX_MAJOR", "DELETE_DUPLICATE", "DELETE_LOW_QUALITY", "REWRITE"]:
        print(f"    {rec:<20s} {recs.get(rec, 0):3d}")
    print()
    print("  Bottom 10:")
    for r in rows[:10]:
        print(f"    {r['total']:3d}  {r['grade']}  {r['name'][:50]:<50s}  {r['recommendation']}")
    print(f"{'='*60}")


def main():
    ap = argparse.ArgumentParser(description="Score all skills against quality rubric")
    ap.add_argument("--skills-dir", default=str(SKILLS_DIR))
    ap.add_argument("--gold-dir", default=str(GOLD_DIR))
    ap.add_argument("--output", default=str(REPORTS_DIR / "skill_scores.csv"))
    args = ap.parse_args()

    rows = score_all(Path(args.skills_dir), Path(args.gold_dir))
    write_csv(rows, Path(args.output))
    print_summary(rows)
    print(f"\n  CSV written to: {args.output}")


if __name__ == "__main__":
    main()
