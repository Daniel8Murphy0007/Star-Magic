# UQFF Star-Magic — Installation Guide

Three install paths depending on your use case.

---

## Path 1: From PyPI (Recommended for users)

The simplest. Works on Windows, macOS, Linux. Python 3.10+ required.

```bash
pip install uqff
```

That's it. Verify:

```bash
uqff version
```

Expected output:

```
uqff 5.27.1
  closures: 794
  truly_independent_primitives: 9
  derivative_primitives: 2
  python: 3.x.y
```

Use it in Python:

```python
import uqff_pure_calculator as u

# Fetch any closure
result = u.calculate_paradox({"paradox": "hubble_tension"})
print(result["value"])

# All 7 nuclear magic numbers, exactly
magic = u.calculate_nuclear_magic({})["value"]["magic_numbers"]
print(magic)  # {'shell_1_He': 2, 'shell_2_O': 8, ...}

# Status report
print(u.calculate_status_report({})["value"]["summary"])
```

---

## Path 2: From source (Recommended for contributors)

If you want to run the fidelity gate, add closures, or modify the
calculator, install from source instead.

```bash
git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
cd Star-Magic
pip install -e .
```

The `-e` flag (editable install) means changes to the source are picked
up without reinstalling.

For documentation building:

```bash
pip install -e ".[docs]"
cd docs
make html              # Linux/macOS
make.bat html          # Windows
# Open docs/_build/html/index.html in a browser
```

For running the fidelity gate:

```bash
python uqff_fidelity_tests.py
# Or via CLI:
uqff gate
```

Expected: `TOTAL: 857 passed, 0 failed`.

---

## Path 3: Docker container (Recommended for isolation)

If you want UQFF in a sandboxed container without affecting your Python
environment:

```bash
# Pull and run latest image
docker run --rm uqff/uqff:5.27.1 version

# Run a specific closure prediction
docker run --rm uqff/uqff:5.27.1 predict hubble_tension --json

# Interactive Python shell with uqff loaded
docker run --rm -it --entrypoint python uqff/uqff:5.27.1
```

Build the image locally from source:

```bash
git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
cd Star-Magic
docker build -t uqff:5.27.1 .
docker run --rm uqff:5.27.1 version
```

---

## Optional: C++ acceleration (`uqff_core`)

The `uqff_core` pybind11 extension provides C++-accelerated versions of
the hot computational kernels. Performance gain: ~5-10× on compute-heavy
workloads. **Not required** for any standard use.

```bash
pip install pybind11
git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
cd Star-Magic
pip install -e .
python -c "import uqff_core; print(uqff_core.__doc__)"
```

To explicitly skip C++ even when pybind11 is installed:

```bash
UQFF_SKIP_CPP=1 pip install -e .
```

---

## Verifying your install

After install, run this comprehensive sanity check:

```bash
# Version + headline numbers
uqff version

# Fetch a known closure
uqff predict hubble_tension --json

# List available closures
uqff list --filter holmlid

# Production status summary
uqff status
```

If `uqff` command is not found, see TROUBLESHOOTING.md "After install, the
`uqff` CLI command isn't recognized."

---

## Upgrading

```bash
pip install --upgrade uqff
```

Or to a specific version:

```bash
pip install uqff==5.27.1
```

---

## Uninstalling

```bash
pip uninstall uqff
```

This removes the package and the `uqff` console script. No system files are
touched.

---

## Requirements

| Component | Version |
|---|---|
| Python | >= 3.10 (tested on 3.10, 3.11, 3.12, 3.13) |
| pip | >= 21.0 (most modern pip versions work) |
| OS | Any (Windows, macOS, Linux; uqff is pure Python) |
| Runtime dependencies | NONE (uses only Python standard library) |
| Optional (C++ extension) | pybind11 + C++17 compiler |
| Optional (docs build) | sphinx, sphinx-rtd-theme |
| Optional (test coverage) | coverage >= 7.0 |

---

## License

UQFF is dual-licensed: AGPL-3.0-or-later (free for academic / research /
non-commercial use) OR commercial license (proprietary deployments).

See `LICENSE`, `LICENSE-AGPL-3.0.txt`, `COMMERCIAL.md`, and the install
output footer for full terms.

---

## Need help?

- TROUBLESHOOTING.md — common install + runtime issues
- FAQ.md — answers to recurring questions
- GitHub Issues: https://github.com/Daniel8Murphy0007/Star-Magic/issues
- Email (author): `daniel.murphy00@enrgyone.com`
