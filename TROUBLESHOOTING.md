# UQFF Star-Magic — Troubleshooting Guide

Common issues and their fixes. If you don't find your issue here, open one
at https://github.com/Daniel8Murphy0007/Star-Magic/issues with the exact
command you ran and the full error output.

---

## Installation

### `pip install uqff` fails with "Could not find a version that satisfies"

UQFF requires **Python 3.10 or newer**. Check your Python version:

```bash
python --version
```

If you're on 3.9 or older, install a newer Python (3.10+ recommended).
The package is published as a pure-Python wheel that works on any OS.

### `pip install uqff` succeeds but `import uqff_pure_calculator` fails

After install, the importable name is `uqff_pure_calculator`, NOT `uqff`.
The PyPI package name (`uqff`) and the importable module name
(`uqff_pure_calculator`) are different — this is a setuptools convention.

```python
# Correct:
import uqff_pure_calculator as u

# Incorrect — fails:
import uqff
```

### After install, the `uqff` CLI command isn't recognized

The `uqff` console script should be installed when you `pip install uqff`.
If your shell can't find it:

- **Windows**: the script is at `C:\Python\Python312\Scripts\uqff.exe`.
  Add `Scripts/` to your PATH, or use the full path.
- **macOS/Linux**: the script is at `~/.local/bin/uqff` if you used
  `pip install --user`. Add `~/.local/bin` to your PATH:
  ```bash
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
  source ~/.bashrc
  ```
- **Any platform**: as a fallback, run the module directly:
  ```bash
  python -m uqff_cli version
  ```

### Optional C++ extension `uqff_core` not building

The C++ pybind11 extension (`uqff_core`) is **not required** for the pure-
Python calculator. The `setup.py` auto-skips building it if `pybind11` is
not installed.

If you DO want the C++ acceleration:

```bash
pip install pybind11
pip install --no-binary :all: uqff   # force compilation
```

Build requires a C++17 compiler (gcc 7+, clang 5+, MSVC 2017+).

To skip the C++ build even when pybind11 is present:

```bash
UQFF_SKIP_CPP=1 pip install uqff
```

---

## Dispatch & API

### `calculate_paradox` returns `{'value': None}` for a known paradox name

**Dispatch keys are LOWERCASE.** The dispatcher normalizes input via
`.lower().strip()`, so `"hubble_tension"` works but `"Hubble_Tension"` does
not (silently returns None).

```python
# Correct:
u.calculate_paradox({"paradox": "hubble_tension"})

# Returns {'value': None} silently — wrong case:
u.calculate_paradox({"paradox": "Hubble_Tension"})
```

To find the correct case-sensitive name:

```bash
uqff list --filter hubble
```

Or in Python:

```python
import uqff_pure_calculator as u
matches = [k for k in u.PARADOX_TO_CLOSURE if "hubble" in k.lower()]
```

### `AttributeError: module 'uqff_pure_calculator' has no attribute 'calculate_xyz'`

The public surface is exactly 34 `calculate_*` functions. List them:

```bash
uqff surfaces
```

Or in Python:

```python
import uqff_pure_calculator as u
[n for n in dir(u) if n.startswith("calculate_")]
```

Anything not in that list is internal and may move or change without
notice. Use `calculate_paradox` for closure lookups, not internal helpers.

### Closure returns `{'value': X}` but X is a dict not a number

Most `calculate_*` surfaces return structured dicts (with `observables`,
`magic_numbers`, `ker_chain`, etc.) — only a few primitive-computation
surfaces return scalar values. See `INPUT_DOMAINS.md` for the full per-
surface return contract.

---

## Encoding & line endings (Windows-specific)

### `SyntaxError: source code string cannot contain null bytes`

Some Windows edit-save cycles append null padding to `.py` files. Strip
them with:

```bash
python scripts/ci_strip_nulls.py uqff_fidelity_tests.py
```

Or manually:

```python
import pathlib
p = pathlib.Path("uqff_fidelity_tests.py")
p.write_bytes(p.read_bytes().replace(b"\x00", b""))
```

### `UnicodeEncodeError` when running fidelity gate on Windows

If your console doesn't support UTF-8 output:

```bash
set PYTHONIOENCODING=utf-8
python uqff_fidelity_tests.py
```

Or in PowerShell:

```powershell
$env:PYTHONIOENCODING = "utf-8"
python uqff_fidelity_tests.py
```

### Git Bash on Windows: `cd C:\Users\...` says "No such file or directory"

Git Bash interprets `\` as escape characters. Use forward slashes:

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
```

Or quote the path:

```bash
cd "C:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic"
```

---

## Fidelity gate failures

### `KK regulator` or other "EXACT regression" tests fail with `err < 1e-15`

This is a known floating-point precision issue. The values agree to 15+
digits but the test tolerance was set tighter than float64 can represent.
Already fixed in v5.27.0. If you see it on an older version:

```bash
pip install --upgrade uqff
```

### Test fails on Python 3.12 but passes on 3.10

Python 3.12 ships a different libm than 3.10. If a test fails by ~1e-16
relative error, it's a float64 rounding artifact, not a physics issue. We
fixed 5 such tests in v5.27.0; if you find another, file an issue.

### `ModuleNotFoundError: No module named 'uqff_pure_calculator'` when running gate

The gate file must live in the same directory as `uqff_pure_calculator.py`
OR be installed via `pip install -e .` from the repo. The CLI `uqff gate`
command handles this lookup automatically.

---

## CI / GitHub Actions failures

### Smoke job fails in ~20 seconds

Almost always one of:
- `uqff_fidelity_tests.py` or `uqff_pure_calculator.py` not committed
- Test file has CRLF line endings (need LF for CI)
- Hard-coded numeric assertion changed (e.g., closure count drifted)

Read the "Diagnostics" step output in the failed run — it lists which
files are present and their sizes.

### PyPI publish fails with `invalid-publisher`

The Trusted Publisher entry on PyPI doesn't match the OIDC claims from
GitHub. Verify exactly:
- Owner: `Daniel8Murphy0007` (case-sensitive)
- Repository name: `Star-Magic` (capital S, M, hyphen not underscore)
- Workflow name: `release.yml` (lowercase, `.yml` not `.yaml`)
- Environment name: `pypi` (lowercase)

PyPI's add-publisher form sometimes silently fails to save. Verify by
scrolling DOWN past the form to the "Pending publishers" section — your
entry must appear there in a table, otherwise the save didn't go through.

---

## Performance

### `import uqff_pure_calculator` takes a few seconds

The single-file module is 2.66 MB / 48,000 lines and ~600 unique closure
functions, all defined at import time. Cold import takes 2-5 seconds
depending on disk speed; warm import (cached .pyc) is < 1 second.

A modular refactor (Tier-3 H1 in `PRODUCTION_ROADMAP.md`) would amortize
import cost across submodules but is a multi-week effort.

### Calculator memory footprint

Approximately 25 MB resident after `import uqff_pure_calculator`. The
PARADOX_TO_CLOSURE dispatch table dominates (~12 MB). The 794 closures
share 616 unique function objects (178 aliases save memory).

---

## Reporting issues

Please include:

1. UQFF version: `uqff version`
2. Python version: `python --version`
3. OS + version (e.g., Windows 11 24H2, macOS 14.5, Ubuntu 22.04)
4. The exact command you ran
5. The full output (stderr included)
6. Whether the issue is reproducible

File at: https://github.com/Daniel8Murphy0007/Star-Magic/issues

Use the bug-report template — it has placeholders for all of the above.
