#!/usr/bin/env python3
"""Append CP4 #263 (AetherSuperfluidDynamicsCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #263 — AetherSuperfluidDynamicsCalculator\n# PAPER_679 | Session 173 | appended by _append_cp4_263.py\n# ========================================================================\n'
CLASS_SRC = 'class AetherSuperfluidDynamicsCalculator:\n    """\n    CP4 #263 — PAPER_679 | UQFF aether superfluid — m_UA, healing length, sound speed, vortex circulation, g_eff\n    CPP_MODULE: AetherSuperfluidDynamics\n    """\n    ENTRY   = 263\n    PAPER   = "PAPER_679"\n    CPP_MODULE = "AetherSuperfluidDynamics"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'radius_m\', \'float\', 10000000000.0, \'radial distance from BH m\'),\n        (\'mass_kg\', \'float\', 8.548e+36, \'BH mass kg (Sgr A*)\'),\n        (\'n_UA\', \'float\', 2.95e+31, \'aether number density m^-3\'),\n    ]\n    PRIMARY_OUTPUT = "g_eff"\n    PRIMARY_INPUT  = "radius_m"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'radius_m\': 10000000000.0, \'mass_kg\': 8.548e+36, \'n_UA\': 2.95e+31}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        radius_m = params.get(\'radius_m\', 10000000000.0)\n        mass_kg = params.get(\'mass_kg\', 8.548e+36)\n        n_UA = params.get(\'n_UA\', 2.95e+31)\n        H0_SI = 2.184e-18\n        M_UA  = HBAR*H0_SI/C**2\n        G_UA  = 1.0e-10\n        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)\n        xi_UA = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)\n        c_UA  = math.sqrt(G_UA*N_UA/M_UA)\n        kv    = 2.0*PI*HBAR/M_UA\n        g_Newton = G*mass_kg/radius_m**2\n        boost = 1.0+(c_UA**2/C**2)*F_TRZ*(RHO_UA/RHO_SCM)\n        g_eff = g_Newton*boost\n        result = {\'healing_length_m\': xi_UA, \'sound_speed_m_s\': c_UA,\n                   \'vortex_circulation\': kv, \'g_eff\': g_eff,\n                   \'g_Newton\': g_Newton, \'m_UA_kg\': M_UA}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class AetherSuperfluidDynamicsCalculator" in existing:
        print(f"SKIP: AetherSuperfluidDynamicsCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #263 AetherSuperfluidDynamicsCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class AetherSuperfluidDynamicsCalculator" in tail or "AetherSuperfluidDynamicsCalculator" in existing:
        print(f"VERIFIED: AetherSuperfluidDynamicsCalculator present")
    else:
        print(f"WARNING: Could not verify AetherSuperfluidDynamicsCalculator — check file")

if __name__ == "__main__":
    main()
