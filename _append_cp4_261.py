#!/usr/bin/env python3
"""Append CP4 #261 (UQFFPredictionsForLISACalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #261 — UQFFPredictionsForLISACalculator\n# PAPER_677 | Session 173 | appended by _append_cp4_261.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFPredictionsForLISACalculator:\n    """\n    CP4 #261 — PAPER_677 | UQFF predictions for LISA space GW observatory — h_LISA_UQFF, Omega_GW, EMRI rate\n    CPP_MODULE: UQFFPredictionsForLISA\n    """\n    ENTRY   = 261\n    PAPER   = "PAPER_677"\n    CPP_MODULE = "UQFFPredictionsForLISA"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'chirp_mass_kg\', \'float\', 1.989e+38, \'SMBH chirp mass (1e8 Msun)\'),\n        (\'distance_m\', \'float\', 1e+27, \'source distance m\'),\n        (\'frequency_hz\', \'float\', 0.001, \'GW frequency in LISA band Hz\'),\n        (\'T_eff_K\', \'float\', 2.73, \'effective temperature K\'),\n    ]\n    PRIMARY_OUTPUT = "h_uqff_lisa"\n    PRIMARY_INPUT  = "frequency_hz"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'chirp_mass_kg\': 1.989e+38, \'distance_m\': 1e+27, \'frequency_hz\': 0.001, \'T_eff_K\': 2.73}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        chirp_mass_kg = params.get(\'chirp_mass_kg\', 1.989e+38)\n        distance_m = params.get(\'distance_m\', 1e+27)\n        frequency_hz = params.get(\'frequency_hz\', 0.001)\n        T_eff_K = params.get(\'T_eff_K\', 2.73)\n        def h_GR_SMBH(Mc, d, f):\n            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2\n        L_LISA = 2.5e9\n        S_UA   = max(0.0, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))\n        M_tot  = chirp_mass_kg * 4.0**0.2\n        r  = r_s(M_tot); TH = T_H(M_tot)\n        S1 = S_SCm(M_tot, TH, r)\n        h_gr   = h_GR_SMBH(chirp_mass_kg, distance_m, frequency_hz)\n        h_uqff_lisa = (1.0-F_TRZ)*S_UA*S1*h_gr\n        omega_GW_ratio = (RHO_UA/9.47e-27)**F_TRZ\n        EMRI_boost = 1.0 + F_TRZ*(RHO_UA/RHO_SCM)\n        result = {\'h_uqff_lisa\': h_uqff_lisa, \'S_UA_LISA\': S_UA,\n                   \'omega_GW_UQFF_ratio\': omega_GW_ratio, \'EMRI_rate_boost\': EMRI_boost}\n        def chirp_mass(m1, m2): return (m1*m2)**0.6/(m1+m2)**0.2\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFPredictionsForLISACalculator" in existing:
        print(f"SKIP: UQFFPredictionsForLISACalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #261 UQFFPredictionsForLISACalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFPredictionsForLISACalculator" in tail or "UQFFPredictionsForLISACalculator" in existing:
        print(f"VERIFIED: UQFFPredictionsForLISACalculator present")
    else:
        print(f"WARNING: Could not verify UQFFPredictionsForLISACalculator — check file")

if __name__ == "__main__":
    main()
