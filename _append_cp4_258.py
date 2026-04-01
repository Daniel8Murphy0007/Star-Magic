#!/usr/bin/env python3
"""Append CP4 #258 (UQFFComparedToLIGODataCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #258 — UQFFComparedToLIGODataCalculator\n# PAPER_674 | Session 173 | appended by _append_cp4_258.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFComparedToLIGODataCalculator:\n    """\n    CP4 #258 — PAPER_674 | UQFF vs General LIGO Dataset — chirp mass, h_UQFF(f), phase shift, frequency sweep\n    CPP_MODULE: UQFFComparedToLIGOData\n    """\n    ENTRY   = 258\n    PAPER   = "PAPER_674"\n    CPP_MODULE = "UQFFComparedToLIGOData"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'m1_kg\', \'float\', 7.16e+31, \'primary BH mass kg\'),\n        (\'m2_kg\', \'float\', 5.77e+31, \'secondary BH mass kg\'),\n        (\'distance_m\', \'float\', 1.27e+25, \'luminosity distance m\'),\n        (\'frequency_hz\', \'float\', 150.0, \'GW frequency Hz\'),\n    ]\n    PRIMARY_OUTPUT = "h_uqff"\n    PRIMARY_INPUT  = "frequency_hz"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'m1_kg\': 7.16e+31, \'m2_kg\': 5.77e+31, \'distance_m\': 1.27e+25, \'frequency_hz\': 150.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        m1_kg = params.get(\'m1_kg\', 7.16e+31)\n        m2_kg = params.get(\'m2_kg\', 5.77e+31)\n        distance_m = params.get(\'distance_m\', 1.27e+25)\n        frequency_hz = params.get(\'frequency_hz\', 150.0)\n        def chirp_mass(m1, m2):\n            return (m1*m2)**0.6 / (m1+m2)**0.2\n        def h_GR(Mc, d, f):\n            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2\n\n        M  = m1_kg + m2_kg\n        Mc = chirp_mass(m1_kg, m2_kg)\n        r  = r_s(M); TH = T_H(M)\n        S1 = S_SCm(M, TH, r)\n        S2 = S_Um_f(M, r, 1e8, frequency_hz)\n        h_gr   = h_GR(Mc, distance_m, frequency_hz)\n        h_uqff = (1.0 - F_TRZ) * S1 * S2 * h_gr\n        delta_phi = KAPPA * F_TRZ * 1e8\n        result = {\'h_gr\': h_gr, \'h_uqff\': h_uqff, \'S_SCm\': S1, \'S_Um\': S2,\n                   \'chirp_mass_kg\': Mc, \'delta_phi_rad\': delta_phi}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFComparedToLIGODataCalculator" in existing:
        print(f"SKIP: UQFFComparedToLIGODataCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #258 UQFFComparedToLIGODataCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFComparedToLIGODataCalculator" in tail or "UQFFComparedToLIGODataCalculator" in existing:
        print(f"VERIFIED: UQFFComparedToLIGODataCalculator present")
    else:
        print(f"WARNING: Could not verify UQFFComparedToLIGODataCalculator — check file")

if __name__ == "__main__":
    main()
