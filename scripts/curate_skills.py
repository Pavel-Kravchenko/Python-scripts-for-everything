#!/usr/bin/env python3
"""Orchestrate full skill curation: score → fix → dedup → re-score → report.

Usage:
    python3 scripts/curate_skills.py                     # full pipeline
    python3 scripts/curate_skills.py --threshold 50      # custom gate
    python3 scripts/curate_skills.py --no-delete          # skip deletion
"""

import argparse
import csv
import os
import shutil
from pathlib import Path

SKILLS_DIR = Path(__file__).resolve().parent.parent / "Skills"
GOLD_DIR = Path.home() / ".claude" / "skills"
REPORTS_DIR = Path(__file__).resolve().parent.parent / "reports"

# Import sibling modules
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_skills import score_all, write_csv, print_summary
from improve_skills import apply_fixes


def load_csv(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


def delete_duplicates_and_junk(skills_dir: Path, rows: list[dict], threshold: int) -> dict:
    """Delete skills marked DELETE_DUPLICATE or DELETE_LOW_QUALITY."""
    stats = {"deleted_duplicate": 0, "deleted_low_quality": 0, "migrated": 0}

    for r in rows:
        path = Path(r["path"])
        rec = r["recommendation"]

        if rec == "DELETE_DUPLICATE":
            if path.exists():
                if path.is_file():
                    os.remove(path)
                    print(f"  Deleted flat duplicate: {r['name']}")
                elif path.parent.is_dir() and path.name == "SKILL.md":
                    shutil.rmtree(path.parent)
                    print(f"  Deleted directory duplicate: {r['name']}")
                stats["deleted_duplicate"] += 1

        elif rec == "DELETE_LOW_QUALITY":
            total = int(r["total"])
            if total < threshold:
                if path.exists():
                    if path.is_file():
                        os.remove(path)
                    elif path.parent.is_dir() and path.name == "SKILL.md":
                        shutil.rmtree(path.parent)
                    print(f"  Deleted low-quality ({total}pts): {r['name']}")
                    stats["deleted_low_quality"] += 1

    return stats


def migrate_valuable_flat_files(skills_dir: Path, rows: list[dict], threshold: int) -> int:
    """Migrate flat .md files that scored above threshold to directory format."""
    migrated = 0
    for r in rows:
        if r["type"] != "flat":
            continue
        if r["recommendation"] in ("DELETE_DUPLICATE", "DELETE_LOW_QUALITY"):
            continue
        total = int(r["total"])
        if total < threshold:
            continue

        flat_path = Path(r["path"])
        if not flat_path.exists():
            continue

        dir_path = skills_dir / r["name"]
        if dir_path.exists():
            continue  # Already has directory version

        dir_path.mkdir(parents=True, exist_ok=True)
        shutil.move(str(flat_path), str(dir_path / "SKILL.md"))
        print(f"  Migrated: {r['name']}.md → {r['name']}/SKILL.md")
        migrated += 1

    return migrated


def compare_reports(before: list[dict], after: list[dict]):
    """Print before/after comparison."""
    before_map = {r["name"]: r for r in before}
    after_map = {r["name"]: r for r in after}

    improved = 0
    degraded = 0
    unchanged = 0
    removed = 0
    total_before = sum(int(r["total"]) for r in before)
    total_after = sum(int(r["total"]) for r in after)

    for name, b in before_map.items():
        if name not in after_map:
            removed += 1
            continue
        a = after_map[name]
        bt, at = int(b["total"]), int(a["total"])
        if at > bt:
            improved += 1
        elif at < bt:
            degraded += 1
        else:
            unchanged += 1

    avg_before = total_before / len(before) if before else 0
    avg_after = total_after / len(after) if after else 0

    # Grade distributions
    def grade_dist(rows):
        d = {}
        for r in rows:
            d[r["grade"]] = d.get(r["grade"], 0) + 1
        return d

    gb = grade_dist(before)
    ga = grade_dist(after)

    print(f"\n{'='*60}")
    print(f"  BEFORE vs AFTER COMPARISON")
    print(f"{'='*60}")
    print(f"  Skills:     {len(before):>4d}  →  {len(after):>4d}  ({len(after)-len(before):+d})")
    print(f"  Avg score:  {avg_before:>5.1f}  →  {avg_after:>5.1f}  ({avg_after-avg_before:+.1f})")
    print(f"  Improved:   {improved}")
    print(f"  Unchanged:  {unchanged}")
    print(f"  Degraded:   {degraded}")
    print(f"  Removed:    {removed}")
    print()
    print("  Grade distribution (before → after):")
    for g in "ABCDF":
        print(f"    {g}: {gb.get(g,0):3d} → {ga.get(g,0):3d}  ({ga.get(g,0)-gb.get(g,0):+d})")
    print(f"{'='*60}")


def main():
    ap = argparse.ArgumentParser(description="Full skill curation pipeline")
    ap.add_argument("--skills-dir", default=str(SKILLS_DIR))
    ap.add_argument("--gold-dir", default=str(GOLD_DIR))
    ap.add_argument("--threshold", type=int, default=50, help="Minimum score to keep")
    ap.add_argument("--no-delete", action="store_true", help="Skip deletion step")
    args = ap.parse_args()

    skills_dir = Path(args.skills_dir)
    gold_dir = Path(args.gold_dir)
    reports = REPORTS_DIR
    reports.mkdir(parents=True, exist_ok=True)

    # Phase 1: Baseline scoring
    print("\n▶ Phase 1: Baseline scoring...")
    before_rows = score_all(skills_dir, gold_dir)
    before_csv = reports / "skill_scores_before.csv"
    write_csv(before_rows, before_csv)
    print_summary(before_rows)

    # Phase 2: Auto-fix
    print("\n▶ Phase 2: Auto-fixing...")
    fixes = ["add_frontmatter", "version_compat", "primary_tool", "duplicate_h1"]
    stats, details = apply_fixes(skills_dir, fixes, dry_run=False)
    print(f"  Fixed {stats['modified']} skills ({stats['fixes_applied']} total fixes)")

    # Phase 3: Delete duplicates and junk
    if not args.no_delete:
        print("\n▶ Phase 3: Removing duplicates and low-quality skills...")
        del_stats = delete_duplicates_and_junk(skills_dir, before_rows, args.threshold)
        print(f"  Deleted {del_stats['deleted_duplicate']} duplicates, {del_stats['deleted_low_quality']} low-quality")

        # Phase 3b: Migrate remaining valuable flat files
        print("\n▶ Phase 3b: Migrating valuable flat files to directory format...")
        migrated = migrate_valuable_flat_files(skills_dir, before_rows, args.threshold)
        print(f"  Migrated {migrated} flat files to directory format")

    # Phase 4: Re-score
    print("\n▶ Phase 4: Re-scoring...")
    after_rows = score_all(skills_dir, gold_dir)
    after_csv = reports / "skill_scores_after.csv"
    write_csv(after_rows, after_csv)
    print_summary(after_rows)

    # Phase 5: Comparison
    compare_reports(before_rows, after_rows)

    print(f"\n  Reports: {before_csv}")
    print(f"           {after_csv}")


if __name__ == "__main__":
    main()
