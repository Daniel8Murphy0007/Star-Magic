#!/usr/bin/env python3
"""Append CP4 #259 (UQFFComparedToGW170817Calculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #259 — UQFFComparedToGW170817Calculator\n# PAPER_675 | Session 173 | appended by _append_cp4_259.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFComparedToGW170817Calculator:\n    """\n    CP4 #259 — PAPER_675 | UQFF vs GW170817 NS-NS merger — GRB delay, tidal deformability, inspiral strain\n    CPP_MODULE: UQFFComparedToGW170817\n    """\n    ENTRY   = 259\n    PAPER   = "PAPER_675"\n    CPP_MODULE = "UQFFComparedToGW170817"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'m1_kg\', \'float\', 2.704e+30, \'NS1 mass kg (1.36 Msun)\'),\n        (\'m2_kg\', \'float\', 2.327e+30, \'NS2 mass kg (1.17 Msun)\'),\n        (\'distance_m\', \'float\', 1.234e+24, \'40 Mpc in m\'),\n        (\'frequency_hz\', \'float\', 1500.0, \'GW frequency Hz\'),\n    ]\n    PRIMARY_OUTPUT = "h_uqff"\n    PRIMARY_INPUT  = "frequency_hz"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'m1_kg\': 2.704e+30, \'m2_kg\': 2.327e+30, \'distance_m\': 1.234e+24, \'frequency_hz\': 1500.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        m1_kg = params.get(\'m1_kg\', 2.704e+30)\n        m2_kg = params.get(\'m2_kg\', 2.327e+30)\n        distance_m = params.get(\'distance_m\', 1.234e+24)\n        frequency_hz = params.get(\'frequency_hz\', 1500.0)\n        def chirp_mass(m1, m2):\n            return (m1*m2)**0.6 / (m1+m2)**0.2\n        def h_GR(Mc, d, f):\n            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2\n        M  = m1_kg + m2_kg; Mc = chirp_mass(m1_kg, m2_kg)\n        r  = r_s(M); TH = T_H(M)\n        S1 = S_SCm(M, TH, r)\n        S2 = S_Um_f(M, r, 1e8, frequency_hz)\n        h_gr   = h_GR(Mc, distance_m, frequency_hz)\n        h_uqff = (1.0-F_TRZ)*S1*S2*h_gr\n        grb_delay = 1.7*(1.0 + F_TRZ*(RHO_UA/RHO_SCM))\n        kappa_tidal_ratio = 1.0 - F_TRZ*RHO_SCM/RHO_UA\n        result = {\'h_gr\': h_gr, \'h_uqff\': h_uqff, \'grb_delay_UQFF_s\': grb_delay,\n                   \'kappa_tidal_ratio\': kappa_tidal_ratio}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFComparedToGW170817Calculator" in existing:
        print(f"SKIP: UQFFComparedToGW170817Calculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #259 UQFFComparedToGW170817Calculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFComparedToGW170817Calculator" in tail or "UQFFComparedToGW170817Calculator" in existing:
        print(f"VERIFIED: UQFFComparedToGW170817Calculator present")
    else:
        print(f"WARNING: Could not verify UQFFComparedToGW170817Calculator — check file")

if __name__ == "__main__":
    main()
