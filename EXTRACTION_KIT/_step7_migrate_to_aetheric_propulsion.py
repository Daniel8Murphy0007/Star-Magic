"""
_step7_migrate_to_aetheric_propulsion.py
Phase Step 7 migration — copies the Aetheric-Propulsion package contents from
EXTRACTION_KIT/new_repo_layout/ into a target Aetheric-Propulsion repo checkout.

Usage:
    python EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py /path/to/Aetheric-Propulsion

The script also refreshes EXTRACTION_KIT/new_repo_layout/ from Star-Magic's current
state (so future `git push` of EXTRACTION_KIT/ stays in sync with the calculator + dispatch).

Idempotent in both directions.
"""
import shutil
import sys
from pathlib import Path

STAR_MAGIC_ROOT = Path(__file__).resolve().parent.parent
KIT_ROOT = STAR_MAGIC_ROOT / "EXTRACTION_KIT"
PACKAGE_ROOT = KIT_ROOT / "new_repo_layout" / "aetheric_propulsion"

# Files copied from Star-Magic ROOT into the new package
PACKAGE_FILES = [
    ("uqff_pure_calculator.py",      "calculator.py"),
    ("assimilation_dispatch.py",     "assimilation_dispatch.py"),
    ("qcalcgeom_solver.py",          "qcalcgeom_solver.py"),
    ("provenance_recorder.py",       "provenance_recorder.py"),
    ("_build_overdetermination_views.py", "_build_overdetermination_views.py"),
]
PACKAGE_DIRS = [
    ("geometry_backends",            "geometry_backends"),
    ("numeric_backends",             "numeric_backends"),
]

# Top-level repo files (placed at NEW_REPO_ROOT/, not inside the package)
TOP_LEVEL_FILES = [
    "pyproject.toml",
    "README.md",
    "LICENSE-AGPL-3.0.txt",
    "COMMERCIAL.md",
    ".github/workflows/release.yml",
    "tests/test_smoke.py",
    "aetheric_propulsion/__init__.py",
]

# Data files: docs + maps that ship with the package
DATA_FILES = [
    ("OVERDETERMINATION_MAP.csv",       "data/OVERDETERMINATION_MAP.csv"),
    ("OVERDETERMINATION_WIDE.csv",      "data/OVERDETERMINATION_WIDE.csv"),
    ("OVERDETERMINATION_MAP.md",        "data/OVERDETERMINATION_MAP.md"),
    ("ASSIMILATION_GEOMETRY_ATLAS.md",  "data/ASSIMILATION_GEOMETRY_ATLAS.md"),
]


def refresh_kit_from_star_magic():
    """Copy current Star-Magic files into EXTRACTION_KIT/new_repo_layout/aetheric_propulsion/."""
    PACKAGE_ROOT.mkdir(parents=True, exist_ok=True)
    n = 0
    for src_name, dst_name in PACKAGE_FILES:
        src = STAR_MAGIC_ROOT / src_name
        if src.exists():
            dst = PACKAGE_ROOT / dst_name
            shutil.copyfile(src, dst)
            n += 1
    for src_dir, dst_dir in PACKAGE_DIRS:
        src = STAR_MAGIC_ROOT / src_dir
        if src.exists():
            dst = PACKAGE_ROOT / dst_dir
            # FUSE-mount-friendly: copy file-by-file rather than rmtree+copytree
            dst.mkdir(parents=True, exist_ok=True)
            for f in src.iterdir():
                if f.is_file() and "__pycache__" not in f.parts:
                    shutil.copyfile(f, dst / f.name)
            n += 1
    # Data files
    data_dst = KIT_ROOT / "new_repo_layout" / "data"
    data_dst.mkdir(parents=True, exist_ok=True)
    for src_name, dst_rel in DATA_FILES:
        src = STAR_MAGIC_ROOT / src_name
        if src.exists():
            dst = KIT_ROOT / "new_repo_layout" / dst_rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            n += 1
    return n


def migrate_to_target(target_root):
    """Copy EXTRACTION_KIT/new_repo_layout/ contents into target Aetheric-Propulsion checkout."""
    target = Path(target_root).resolve()
    if not target.exists():
        print(f"ERROR: target directory {target} does not exist")
        return 1
    src_root = KIT_ROOT / "new_repo_layout"
    if not src_root.exists():
        print(f"ERROR: {src_root} not found; run refresh first")
        return 1
    n = 0
    for item in src_root.rglob("*"):
        if item.is_file():
            rel = item.relative_to(src_root)
            dst = target / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(item, dst)
            n += 1
    print(f"Migrated {n} files to {target}")
    return 0


def main():
    refresh_count = refresh_kit_from_star_magic()
    print(f"Refreshed EXTRACTION_KIT/new_repo_layout/ from Star-Magic ({refresh_count} files+dirs)")
    if len(sys.argv) > 1:
        return migrate_to_target(sys.argv[1])
    print("(no target dir provided; only refreshed the kit)")
    print("To migrate: python EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py /path/to/Aetheric-Propulsion")
    return 0


if __name__ == "__main__":
    sys.exit(main())
