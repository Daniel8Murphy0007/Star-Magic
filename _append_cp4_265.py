#!/usr/bin/env python3
"""Append CP4 #265 (GrossPitaevskiiVortexSimulationCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #265 — GrossPitaevskiiVortexSimulationCalculator\n# PAPER_681 | Session 173 | appended by _append_cp4_265.py\n# ========================================================================\n'
CLASS_SRC = 'class GrossPitaevskiiVortexSimulationCalculator:\n    """\n    CP4 #265 — PAPER_681 | UQFF GP vortex simulation — imaginary-time ground state, chemical potential, density profile\n    CPP_MODULE: GrossPitaevskiiVortexSimulation\n    """\n    ENTRY   = 265\n    PAPER   = "PAPER_681"\n    CPP_MODULE = "GrossPitaevskiiVortexSimulation"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'r_min_m\', \'float\', 1000000000.0, \'inner grid radius m\'),\n        (\'r_max_m\', \'float\', 100000000000.0, \'outer grid radius m\'),\n        (\'mass_kg\', \'float\', 8.548e+36, \'BH mass kg\'),\n        (\'n_grid\', \'int\', 100, \'radial grid points\'),\n        (\'n_steps\', \'int\', 50, \'imaginary-time steps\'),\n    ]\n    PRIMARY_OUTPUT = "chemical_potential_eV"\n    PRIMARY_INPUT  = "radius_m"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'r_min_m\': 1000000000.0, \'r_max_m\': 100000000000.0, \'mass_kg\': 8.548e+36, \'n_grid\': 100, \'n_steps\': 50}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        r_min_m = params.get(\'r_min_m\', 1000000000.0)\n        r_max_m = params.get(\'r_max_m\', 100000000000.0)\n        mass_kg = params.get(\'mass_kg\', 8.548e+36)\n        n_grid = params.get(\'n_grid\', 100)\n        n_steps = params.get(\'n_steps\', 50)\n        H0_SI = 2.184e-18\n        M_UA  = HBAR*H0_SI/C**2\n        G_UA  = 1.0e-10\n        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)\n        # Thomas-Fermi approx: mu = G_UA*n_UA\n        mu_UA = G_UA*N_UA\n        xi    = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)\n        c_UA  = math.sqrt(G_UA*N_UA/M_UA)\n        # Gravitational potential at r_min_m\n        V_grav = -G*mass_kg*M_UA/r_min_m\n        E_total = mu_UA + V_grav\n        chemical_potential_eV = E_total / 1.602e-19\n        result = {\'chemical_potential_eV\': chemical_potential_eV,\n                   \'healing_length_m\': xi, \'mu_UA_J\': mu_UA,\n                   \'TF_density_m3\': N_UA, \'c_UA_m_s\': c_UA}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class GrossPitaevskiiVortexSimulationCalculator" in existing:
        print(f"SKIP: GrossPitaevskiiVortexSimulationCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #265 GrossPitaevskiiVortexSimulationCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class GrossPitaevskiiVortexSimulationCalculator" in tail or "GrossPitaevskiiVortexSimulationCalculator" in existing:
        print(f"VERIFIED: GrossPitaevskiiVortexSimulationCalculator present")
    else:
        print(f"WARNING: Could not verify GrossPitaevskiiVortexSimulationCalculator — check file")

if __name__ == "__main__":
    main()
