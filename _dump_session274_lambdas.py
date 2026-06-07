"""Session 274 raw derivative-function dump.

For every entry in _PA_S274_CLOSURE_REGISTRY whose form() is a real callable
(not None), print:
  - S_id
  - paper_tag
  - label
  - the verbatim lambda source line as it appears in uqff_pure_calculator.py
  - the actual numeric value returned by calling form()

No percentages. No comparisons. Just the raw computation each closed-form
derivative function produces from the locked vacuum-ledger primitives.
"""
from __future__ import annotations

import re
from pathlib import Path

import uqff_pure_calculator as upc


SRC = Path(upc.__file__).read_text(encoding="utf-8").splitlines()


def find_source_line(s_id: str) -> str:
    needle = f'"{s_id}":'
    for line in SRC:
        if needle in line:
            return line.strip()
    return "<source line not located>"


def extract_form_text(src_line: str) -> str:
    m = re.search(r'\("[^"]+",\s*"[^"]+",\s*(.+?),\s*(?:-?[\d\.eE+]|None)', src_line)
    if not m:
        return "<form not parseable>"
    return m.group(1).strip()


def fmt_value(v) -> str:
    if v is None:
        return "None"
    if isinstance(v, bool):
        return repr(v)
    if isinstance(v, int):
        return f"{v}"
    if isinstance(v, float):
        if v == 0.0:
            return "0.0"
        a = abs(v)
        if a >= 1e6 or a < 1e-3:
            return f"{v:.10e}"
        return f"{v:.10g}"
    return repr(v)


def main() -> None:
    reg = upc._PA_S274_CLOSURE_REGISTRY
    transcribed = [(sid, t) for sid, t in reg.items() if t[2] is not None]
    print(f"Session 274 transcribed lambdas: {len(transcribed)} / {len(reg)} entries")
    print("Locked primitive snapshot:")
    primitives = {
        "BETA_I":        upc.BETA_I,
        "BETA0_DPM":     upc.BETA0_DPM,
        "S_26":          upc.S_26,
        "S26_DPM":       upc.S26_DPM,
        "S26_3":         upc.S26_3,
        "PHI_RESONANCE": upc.PHI_RESONANCE,
        "SSQ":           upc.SSQ,
        "D_CRIT":        upc.D_CRIT,
        "D_BSFG":        upc.D_BSFG,
        "TRZ":           upc.TRZ,
        "A_26":          upc.A_26,
        "RHO_SCM":       upc.RHO_SCM,
        "RHO_UA":        upc.RHO_UA,
        "OMEGA_SCM":     upc.OMEGA_SCM,
        "F_THZ":         upc.F_THZ,
        "PLANCK_H":      upc.PLANCK_H,
        "C_LIGHT":       upc.C_LIGHT,
        "EV_J":          upc.EV_J,
        "K_B":           upc.K_B,
        "G_NEWTON":      upc.G_NEWTON,
    }
    for k, v in primitives.items():
        print(f"    {k:<14} = {fmt_value(v)}")
    print()
    print("=" * 110)

    current_paper = None
    for sid, (tag, label, form, _pred, _obs, unit) in transcribed:
        if tag != current_paper:
            current_paper = tag
            print()
            print(f"---- PAPER_{tag} ----")
        src_line = find_source_line(sid)
        form_text = extract_form_text(src_line)
        try:
            value = form()
        except Exception as exc:
            value = f"<RAISED: {exc}>"
        print(f"{sid}  {label}  [{unit}]")
        print(f"    form  = {form_text}")
        print(f"    value = {fmt_value(value)}")


if __name__ == "__main__":
    main()
