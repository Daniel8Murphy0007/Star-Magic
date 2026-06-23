---
name: Bug report
about: Report something broken or unexpected
title: "[BUG] "
labels: bug
assignees: ''
---

## Describe the bug

A clear, concise description of what's broken.

## To reproduce

Steps to reproduce the behavior:
1. `pip install uqff==<version>`
2. Run `python -c "..."` (paste the exact code)
3. See error: `<paste the error here>`

## Expected behavior

What did you expect to happen?

## Actual behavior

What happened instead?

## Environment

- **UQFF version**: (run `uqff version` and paste output)
- **Python version**: (run `python --version`)
- **Operating system**: (e.g., Windows 11 24H2, macOS 14.5, Ubuntu 22.04)
- **Install method**: (e.g., `pip install uqff`, source clone, conda, etc.)

## Full error output

```
Paste the complete stderr + stdout here, including any tracebacks.
Strip secrets/PII if necessary.
```

## Additional context

- Have you tried `pip install --upgrade uqff` to ensure you're on the latest?
- Does the issue reproduce on a fresh Python environment?
- Is the fidelity gate passing locally? (`uqff gate`)

## Reading list before reporting

- [ ] I checked `TROUBLESHOOTING.md` for my issue
- [ ] I searched existing issues for duplicates
- [ ] I'm on the latest UQFF version
