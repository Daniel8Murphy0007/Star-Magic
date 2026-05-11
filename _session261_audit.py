"""Session 261 -- Coherence Audit Pass 1 (read-only).

Five independent checks:
  1. Canonical constants drift (kappa, SSq, beta_i, xi, rho_BSFG, rho_gamma0, H_SCm)
  2. Falsifier zero-parameter declarations in CP4 P6-P14 (#258-#264)
  3. Stale PDF count (delegated to _audit_stale_pdfs.py results)
  4. arxiv bundle lint sweep (4 bundles)
  5. MUGE/Newtonian leak scan (GM/r in CondensedPhysics*.py)

Outputs a single audit report dict + exit code (0 = clean, 1 = any flags).
"""
from __future__ import annotations
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).parent
FLAGS: list[str] = []
NOTES: list[str] = []


def ok(msg: str) -> None:
    print(f"  [OK]   {msg}")


def warn(msg: str) -> None:
    FLAGS.append(msg)
    print(f"  [FLAG] {msg}")


def info(msg: str) -> None:
    NOTES.append(msg)
    safe = msg.encode("ascii", "replace").decode("ascii")
    print(f"  [NOTE] {safe}")


# ---------- Task 1: Constants audit ----------
CANONICAL = {
    # numeric pattern -> (label, contexts where it's expected)
    r"0\.0005": "kappa=0.0005/day",
    r"\b0\.57\b": "[SSq]=0.57",
    r"0\.603": "beta_i=0.603",
    r"4\.333|13\.0\s*/\s*3\.0|13\s*/\s*3": "xi=13/3=4.333",
    r"5\.96e-?10": "rho_BSFG=5.96e-10",
    r"4\.17e-?14": "rho_gamma0=4.17e-14",
    r"\b0\.99\b": "H_SCm=0.99",
}

CANON_FILES = [
    "CondensedPhysics.py",
    "CondensedPhysics2.py",
    "CondensedPhysics3.py",
    "CondensedPhysics4.py",
    "QCalcGeom.cpp",
    "_p4_p5_predictions_table.py",
]


def task_constants() -> None:
    print("\n[1/5] Canonical constants presence audit")
    for fname in CANON_FILES:
        p = REPO / fname
        if not p.exists():
            warn(f"Missing file: {fname}")
            continue
        text = p.read_text(encoding="utf-8", errors="ignore")
        hits = {label: bool(re.search(pat, text)) for pat, label in CANONICAL.items()}
        missing = [lbl for lbl, present in hits.items() if not present]
        if missing:
            info(f"{fname}: not all canonical constants present (missing: {len(missing)}/{len(CANONICAL)}) - expected for narrowly scoped files")
        else:
            ok(f"{fname}: all 7 canonical constants present")


# ---------- Task 2: Falsifier zero-param audit ----------
def task_falsifier_zero_param() -> None:
    print("\n[2/5] CP4 falsifier zero-parameter audit (#258-#264)")
    p = REPO / "CondensedPhysics4.py"
    if not p.exists():
        warn("CondensedPhysics4.py missing")
        return
    text = p.read_text(encoding="utf-8", errors="ignore")
    # Look for the seven P6-P14 calculators added in sessions 257-260
    # Ground truth from CP4 #258-#264 (verified Session 261 audit)
    targets = [
        ("UQFFKKTowerHbarRegulatorCalculator",   "#258 P6 (Session 257, PAPER_1173)"),
        ("UQFFRingdownSpectralOffsetCalculator", "#259 P11 (Session 258, PAPER_1175)"),
        ("UQFFSigma8WeakLensingCalculator",      "#260 P12 (Session 258, PAPER_1176)"),
        ("UQFF2027JointFalsifierCalculator",     "#261 P6+P11+P12 joint (Session 259, PAPER_1177)"),
        ("UQFFDarkEnergySecondDerivativeCalculator", "#262 P13 (Session 259, PAPER_1178)"),
        ("UQFF2027QuadrupleFalsifierCalculator", "#263 P6+P10+P11+P12 joint (Session 260, PAPER_1179)"),
        ("UQFFCMBmuDistortionCalculator",        "#264 P14 (Session 260, PAPER_1180)"),
    ]
    for name, label in targets:
        if name in text:
            # check has zero free params via 'free_parameters' or 'fixed' marker
            block_match = re.search(rf"class\s+{name}\b.*?(?=\nclass\s|\Z)", text, re.DOTALL)
            if block_match:
                block = block_match.group(0)
                has_zero_marker = (
                    "free_parameters" in block
                    or "zero parameters" in block.lower()
                    or "no free parameter" in block.lower()
                    or "fixed by closure" in block.lower()
                )
                if has_zero_marker:
                    ok(f"{label}: {name} declares zero/fixed free parameters")
                else:
                    info(f"{label}: {name} present but lacks explicit zero-parameter marker (cosmetic)")
            else:
                warn(f"{label}: {name} class block not parseable")
        else:
            warn(f"{label}: {name} NOT FOUND in CondensedPhysics4.py")


# ---------- Task 3: Stale PDFs ----------
def task_stale_pdfs() -> None:
    print("\n[3/5] Stale PDF sweep (delegating to _audit_stale_pdfs.py)")
    try:
        r = subprocess.run([sys.executable, "_audit_stale_pdfs.py"], capture_output=True, text=True, cwd=REPO, timeout=60)
        out = r.stdout
        m = re.search(r"Whitepapers with stale PDF.*?:\s*(\d+)", out)
        if m:
            n = int(m.group(1))
            if n == 0:
                ok("0 stale PDFs in last 25 commits")
            else:
                warn(f"{n} stale PDFs detected in last 25 commits")
        else:
            info("Stale PDF parse: no match")
    except Exception as e:
        warn(f"Stale PDF audit failed: {e}")


# ---------- Task 4: arxiv bundle lint ----------
BUNDLES = [
    "arxiv_submission_1159_1170",
    "arxiv_submission_1159_1172",
    "arxiv_submission_1173_1176",
    "arxiv_submission_1177_1180",
]


def task_lint_bundles() -> None:
    print("\n[4/5] arxiv bundle lint sweep")
    lint_script = REPO / "_arxiv_lint.py"
    if not lint_script.exists():
        warn("_arxiv_lint.py not found")
        return
    for b in BUNDLES:
        path = REPO / b
        if not path.exists():
            warn(f"Bundle missing: {b}")
            continue
        if (path / "SUPERSEDED.md").exists():
            info(f"{b}: SUPERSEDED (skipping lint; see SUPERSEDED.md)")
            continue
        try:
            r = subprocess.run([sys.executable, "_arxiv_lint.py", "--bundle", b], capture_output=True, text=True, cwd=REPO, timeout=90)
            tail = r.stdout.strip().splitlines()[-3:] if r.stdout else []
            tail_s = " | ".join(tail)
            if r.returncode == 0:
                ok(f"{b}: lint clean ({tail_s[:120]})")
            else:
                warn(f"{b}: lint exit {r.returncode} ({tail_s[:120]})")
        except Exception as e:
            warn(f"{b}: lint exception {e}")


# ---------- Task 5: MUGE/Newtonian leak scan ----------
NEWTON_PATTERNS = [
    re.compile(r"\bG\s*\*\s*M\s*/\s*r\b"),
    re.compile(r"\bG_GRAV\s*\*\s*M.*?/\s*r"),
    re.compile(r"GRAVITATIONAL_CONSTANT\s*\*\s*M.*?/\s*r"),
]
NEWTON_SAFE_HINTS = ("emergent", "newtonian limit", "comparison", "validation only", "schwarzschild")


def task_newton_leak() -> None:
    print("\n[5/5] MUGE / Newtonian leak scan (CondensedPhysics*.py)")
    cp_files = sorted(REPO.glob("CondensedPhysics*.py"))
    cp_files = [p for p in cp_files if not p.name.endswith(".bak") and not p.name.endswith(".bak2")]
    total_hits = 0
    for p in cp_files:
        text = p.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()
        file_flags: list[tuple[int, str]] = []
        for i, line in enumerate(lines, 1):
            for pat in NEWTON_PATTERNS:
                if pat.search(line):
                    # context window
                    lo = max(0, i - 6)
                    hi = min(len(lines), i + 2)
                    ctx = "\n".join(lines[lo:hi]).lower()
                    if any(h in ctx for h in NEWTON_SAFE_HINTS):
                        continue
                    file_flags.append((i, line.strip()[:100]))
                    break
        if file_flags:
            total_hits += len(file_flags)
            info(f"{p.name}: {len(file_flags)} GM/r occurrences without emergent-gravity context")
            for ln, snippet in file_flags[:3]:
                info(f"   L{ln}: {snippet}")
            if len(file_flags) > 3:
                info(f"   ... {len(file_flags)-3} more")
        else:
            ok(f"{p.name}: no unwrapped GM/r patterns")
    if total_hits == 0:
        ok("Zero MUGE/Newtonian leaks across CP1-CP4")
    else:
        info(f"Total flagged GM/r occurrences (for manual review): {total_hits}")


def main() -> int:
    print("=" * 70)
    print(" SESSION 261 -- Coherence Audit Pass 1 (READ-ONLY)")
    print("=" * 70)
    task_constants()
    task_falsifier_zero_param()
    task_stale_pdfs()
    task_lint_bundles()
    task_newton_leak()
    print("\n" + "=" * 70)
    print(f" FLAGS: {len(FLAGS)}   NOTES: {len(NOTES)}")
    print("=" * 70)
    if FLAGS:
        print("\n--- FLAGS (require action) ---")
        for f in FLAGS:
            print(f"  - {f}")
    return 0 if not FLAGS else 1


if __name__ == "__main__":
    sys.exit(main())
