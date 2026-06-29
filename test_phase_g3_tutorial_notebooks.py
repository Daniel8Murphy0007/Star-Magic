"""
test_phase_g3_tutorial_notebooks.py
Phase G3 regression — verify the 10 per-domain tutorial notebooks parse cleanly and
their code cells execute without raising. Lightweight test (no Jupyter kernel needed
in CI): we extract code cells via JSON and exec() them in a single namespace per notebook.
"""
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NB_DIR = ROOT / "notebooks"

DOMAINS = ["SI", "SM", "LCDM", "astro", "GR", "chem", "CM", "bio", "geo", "KK"]


def execute_notebook(path):
    """Load + execute all code cells in a notebook; return (ok, error_msg, n_code_cells)."""
    try:
        nb = json.loads(path.read_text(encoding='utf-8'))
    except Exception as e:
        return False, f"JSON parse failed: {e}", 0

    if nb.get("nbformat") != 4:
        return False, f"nbformat is {nb.get('nbformat')} (expected 4)", 0

    cells = nb.get("cells", [])
    code_cells = [c for c in cells if c.get("cell_type") == "code"]
    if not code_cells:
        return False, "no code cells", 0

    ns = {}
    for i, cell in enumerate(code_cells):
        src = cell.get("source", [])
        if isinstance(src, list):
            src = "".join(src) if src else ""
        try:
            exec(compile(src, f"{path.name}#cell{i}", "exec"), ns)
        except Exception as e:
            return False, f"cell {i} raised {type(e).__name__}: {e}", len(code_cells)

    return True, None, len(code_cells)


def main():
    print("=" * 72)
    print("PHASE G3 — Per-domain tutorial notebook regression")
    print("=" * 72)

    all_pass = True
    sys.path.insert(0, str(ROOT))

    for dom in DOMAINS:
        nbs = list(NB_DIR.glob(f"1?_assimilation_{dom}.ipynb"))
        if not nbs:
            print(f"FAIL: no notebook for domain {dom}")
            all_pass = False
            continue
        nb = nbs[0]
        ok, err, n_cells = execute_notebook(nb)
        if ok:
            print(f"  PASS  {nb.name:<40s}  {n_cells} code cells executed")
        else:
            print(f"  FAIL  {nb.name:<40s}  {err}")
            all_pass = False

    print()
    if all_pass:
        print("PHASE G3 SUCCESS CRITERION MET. 10 / 10 notebooks executable.")
        return 0
    print("PHASE G3 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
