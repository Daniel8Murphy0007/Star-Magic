#!/usr/bin/env python3
"""
Session 166 patch: Upgrade PAPER_400–421 SM Anchor headers to CVW v2.0.0 format.
Changes:
  1. Old header → new CVW v2.0.0 header
  2. Replace old footer tag with PAPER_642 cite line
"""

import os
import re

WHITEPAPER_DIR = os.path.join(os.path.dirname(__file__), "whitepapers")

OLD_HEADER = "## §SM Anchors — UQFF Predictions vs. Standard-Model Experiments"
NEW_HEADER = "## §SM Anchors — Standard Model Cross-Validation (G6 Gate, CVW v2.0.0)"

OLD_FOOTER = "*CVW Gate G6 — Session 164 patch*"
NEW_FOOTER = ("*CVW Gate G6 — Session 166 patch (CVW v2.0.0 upgrade)*\n\n"
              "*Cite PAPER_642 (`UQFFSMParameterBridgeMasterComparisonCalculator`) "
              "for full UQFF–SM bridge.*")

TARGET_RANGE = set(range(400, 422))  # 400–421 inclusive

patched = []
skipped = []
errors = []

for fname in sorted(os.listdir(WHITEPAPER_DIR)):
    if not fname.endswith(".md"):
        continue
    m = re.match(r"PAPER_(\d+)_", fname)
    if not m:
        continue
    num = int(m.group(1))
    if num not in TARGET_RANGE:
        continue

    fpath = os.path.join(WHITEPAPER_DIR, fname)
    try:
        with open(fpath, "r", encoding="utf-8") as f:
            content = f.read()
    except Exception as e:
        errors.append(f"{fname}: read error {e}")
        continue

    if OLD_HEADER not in content:
        skipped.append(f"{fname}: old header not found (may already be patched)")
        continue

    new_content = content.replace(OLD_HEADER, NEW_HEADER, 1)

    if OLD_FOOTER in new_content:
        new_content = new_content.replace(OLD_FOOTER, NEW_FOOTER, 1)

    try:
        with open(fpath, "w", encoding="utf-8") as f:
            f.write(new_content)
        patched.append(fname)
    except Exception as e:
        errors.append(f"{fname}: write error {e}")

print(f"\n=== patch_cvw_400_421.py complete ===")
print(f"  Patched : {len(patched)}")
for p in patched:
    print(f"    ✓ {p}")
if skipped:
    print(f"  Skipped : {len(skipped)}")
    for s in skipped:
        print(f"    ~ {s}")
if errors:
    print(f"  Errors  : {len(errors)}")
    for e in errors:
        print(f"    ✗ {e}")
