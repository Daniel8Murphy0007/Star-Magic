#!/usr/bin/env python3
"""Repo-local audit for the Newton shortfall against Star-Magic canonical rules.

This script does three practical things:
1. Records the simultaneous-calculation requirements already present in the repo.
2. Separates mass-bearing actions from massless/projected actions in the DPM stack.
3. Solves for G using the Sun values already written in Star-Magic.txt.

It is intentionally read-only and explanatory so the ontology can be reviewed before
calculator or CI rewrites are attempted.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Requirement:
    title: str
    source: str
    citation: str
    requirement: str


@dataclass(frozen=True)
class ProjectionAction:
    name: str
    category: str
    source: str
    citation: str
    role: str


@dataclass(frozen=True)
class AstroSystem:
    name: str
    mass_kg: float
    radius_m: float
    surface_gravity_m_s2: float
    magnetic_field_t: float
    rotation_rate_rad_s: float


@dataclass(frozen=True)
class EquationBlock:
    title: str
    citation: str
    equation: str
    meaning: str


SIMULTANEOUS_REQUIREMENTS = [
    Requirement(
        title="Five-calculator simultaneous pipeline",
        source="ARCHITECTURE_FLOW_DIAGRAM.md",
        citation="ARCHITECTURE_FLOW_DIAGRAM.md#L72",
        requirement="All 5 calculators run in parallel through the simultaneous joint pipeline.",
    ),
    Requirement(
        title="Inside/outside simultaneous solve",
        source="26D_DOWNWARD_PROJECTION.md",
        citation="26D_DOWNWARD_PROJECTION.md#L151",
        requirement="All UQFF calculations require simultaneous inside + outside solve.",
    ),
    Requirement(
        title="Crossing criterion",
        source="26D_DOWNWARD_PROJECTION.md",
        citation="26D_DOWNWARD_PROJECTION.md#L159",
        requirement="Observable solution occurs at the inside/outside crossing minimum.",
    ),
    Requirement(
        title="Dynamic equation sets",
        source="QCalc.py",
        citation="QCalc.py#L20",
        requirement="Every solution surface carries long-form equations, solutions, available equations, and a simulation_set.",
    ),
    Requirement(
        title="Three-method simultaneous convergence",
        source="CondensedPhysics4.py",
        citation="CondensedPhysics4.py#L11990",
        requirement="Compressed, Resonant, and Buoyancy forms are required to converge as coequal simultaneous methods.",
    ),
    Requirement(
        title="Simultaneous equivalence hub",
        source="CondensedPhysics4.py",
        citation="CondensedPhysics4.py#L10346",
        requirement="Inside hypergraph evolution and outside pi-curvature observables are merged through a simultaneous equivalence hub.",
    ),
    Requirement(
        title="Zero-mass simultaneous sculpting",
        source="ARCHITECTURE_FLOW_DIAGRAM.md",
        citation="ARCHITECTURE_FLOW_DIAGRAM.md#L115",
        requirement="Zero-Mass UA reformulation and 26D simultaneous geometric infinity sculpting are canonical repo requirements.",
    ),
    Requirement(
        title="Gravity family assembled in layers",
        source=".github/copilot-instructions.md",
        citation=".github/copilot-instructions.md#L138",
        requirement="Triadic gravity is computed as the summed family of Ug1_i + Ug2_i + Ug3_i + Ug4_i across 26 layers.",
    ),
]


MASS_AND_MASSLESS_ACTIONS = [
    ProjectionAction(
        name="Ug1",
        category="mass-bearing DPM promoter",
        source="Star-Magic.txt",
        citation="Star-Magic.txt#L54",
        role="DPM internal dipole term that seeds the entire gravity family.",
    ),
    ProjectionAction(
        name="Ug1 drives Ug2/Ug3/Ug4/Ug4_i",
        category="simultaneous promotion rule",
        source="Star-Magic.txt",
        citation="Star-Magic.txt#L56",
        role="Canonical statement that Ug1 generates the downstream gravity family simultaneously.",
    ),
    ProjectionAction(
        name="Ug2",
        category="mass-bearing shell action",
        source="Star-Magic.txt",
        citation="Star-Magic.txt#L57",
        role="Outer-field bubble uses M_s and reactor terms to form the heliosphere shell.",
    ),
    ProjectionAction(
        name="Ug3",
        category="mass-bearing string action",
        source="Star-Magic.txt",
        citation="Star-Magic.txt#L60",
        role="Magnetic-string disk couples to planetary cores and maintains orbital/spin structure.",
    ),
    ProjectionAction(
        name="Ug4",
        category="mass-bearing BH coupling",
        source="Star-Magic.txt",
        citation="Star-Magic.txt#L63",
        role="Black-hole interaction term scales with M_bh and galactic distance.",
    ),
    ProjectionAction(
        name="Zero-Mass UA reformulation",
        category="massless vacuum action",
        source="ARCHITECTURE_FLOW_DIAGRAM.md",
        citation="ARCHITECTURE_FLOW_DIAGRAM.md#L115",
        role="Canonical vacuum condition: rho_UA = 0 immutable and F_U = 0 vacuum.",
    ),
    ProjectionAction(
        name="1 from 0 mass emergence",
        category="massless-to-mass transition",
        source="26D_DOWNWARD_PROJECTION.md",
        citation="26D_DOWNWARD_PROJECTION.md#L175",
        role="Mass emergence is treated as order arising from a zero-mass state rather than a Newtonian starting mass.",
    ),
    ProjectionAction(
        name="SCm-prime",
        category="massless projected action",
        source="CondensedPhysics2.py",
        citation="CondensedPhysics2.py#L10386",
        role="SCm' is explicitly tracked as a massless influence factor rather than a kg/m^3 matter term.",
    ),
    ProjectionAction(
        name="Buoyancy quantum portion",
        category="massless counter-action",
        source="grok_chat_06094c66-e6ac.txt",
        citation="grok_chat_06094c66-e6ac.txt#L8390",
        role="Buoyancy is described as the massless quantum portion countering standard gravity.",
    ),
]


NINE_FIX_POINTS = [
    "Declare one canonical rule: Ug1 is DPM-generative; Ug2/Ug3/Ug4/Ug4_i are downstream; GM/r^2 is only a derived projection.",
    "Split every gravity implementation into foundation, derived_ranges, and projection layers.",
    "Ban 'Newtonian base' wording in upstream code paths.",
    "Rewrite validations from first principles: Ug1 -> gravity family -> Ub -> Um -> A -> derived modes.",
    "Audit every GM/r^2 occurrence and classify it as valid projection, legacy shorthand, or framework violation.",
    "Centralize master equations in one canonical source instead of duplicating variants.",
    "Enforce terminology: foundation means DPM/Ug-family only; emergent means reduced observational behavior.",
    "Fix documentation and code together so stale comments do not reintroduce the shortfall.",
    "Rebuild CI expectations after ontology cleanup so checks validate field hierarchy rather than Newton-first approximations.",
]


EQUATION_BLOCKS = [
    EquationBlock(
        title="Zero-mass vacuum / BigBang state",
        citation="ARCHITECTURE_FLOW_DIAGRAM.md#L115",
        equation="rho_UA = 0; rho_vac = |grad(UA)|; F_U(vacuum) = 0",
        meaning="The starting state is zero-mass vacuum, not Newtonian moving mass.",
    ),
    EquationBlock(
        title="Mass emergence ordering",
        citation="26D_DOWNWARD_PROJECTION.md#L175",
        equation="Prob_order = exp(-Entropy_26D / v_init) / (Partition_9D * (v_init - v_current))",
        meaning="Mass emerges from order formation in a zero-mass state.",
    ),
    EquationBlock(
        title="Inside/outside simultaneous solve",
        citation="26D_DOWNWARD_PROJECTION.md#L151",
        equation="G^(n+1) = R(G^(n)) + IG^(n); O^(n) = pi_[n] * FUB_i(x) + Ricci(G^(n)); n_cross = argmin |Inside - Outside|",
        meaning="Observable structure appears at the inside/outside crossing.",
    ),
    EquationBlock(
        title="Ug1 DPM promoter",
        citation="Star-Magic.txt#L54",
        equation="Ug1 = k1 * mu_s(t,rho_vac,[SCm]) * grad(M_s/r) * exp(-alpha*t) * cos(pi*t_n) * (1 + delta_def)",
        meaning="Ug1 is the promotive DPM seed term for the full gravity family.",
    ),
    EquationBlock(
        title="Ug2 outer field bubble",
        citation="Star-Magic.txt#L57",
        equation="Ug2 = k2 * (rho_vac,[UA] + rho_vac,[SCm]) * M_s/r^2 * S(r-R_b) * (1 + delta_sw*v_sw) * H_SCm * E_react",
        meaning="Ug2 forms the shell/heliosphere side of the gravity family.",
    ),
    EquationBlock(
        title="Ug3 magnetic strings disk",
        citation="Star-Magic.txt#L60",
        equation="Ug3 = k3 * sum_j B_j(r,theta,t,rho_vac,[SCm]) * cos(omega_s(t)*t*pi) * P_core * E_react",
        meaning="Ug3 is the magnetic-string and core-coupling range.",
    ),
    EquationBlock(
        title="Ug4 black-hole coupling",
        citation="Star-Magic.txt#L63",
        equation="Ug4 = k4 * rho_vac,[SCm] * M_bh/d_g * exp(-alpha*t) * cos(pi*t_n) * (1 + f_feedback)",
        meaning="Ug4 is the star-black-hole interaction range.",
    ),
    EquationBlock(
        title="Unified field long form",
        citation="Star-Magic.txt#L52",
        equation="F_U = sum_i[k_i*Ug_i - beta_i*Ug_i*Omega_g*M_bh/d_g*E_react] + sum_j[(mu_j/r_j)*(1-exp(-gamma*t*cos(pi*t_n)))*phi_hat_j] + (g_munu + eta*T_s_munu) - sum_i[lambda_i*U_i*E_react]",
        meaning="The full field is assembled from the Ug family plus buoyancy, magnetism, aether, and intelligent operators.",
    ),
    EquationBlock(
        title="Compressed mode",
        citation="QCalc.py#L7160",
        equation="g_comp = [g_base*(1 + H0*t)*(1 - B/B_crit)*H_SCm] + g_Ug_sum + Lambda*c^2/3 + g_quantum + g_fluid + g_pert all times TRZ",
        meaning="Compressed mode is a downstream operational projection and not the seed ontology.",
    ),
    EquationBlock(
        title="Resonant mode",
        citation="QCalc.py#L8428",
        equation="g_res = a_DPM + sum(i=1..13) a_i(omega,E_vac,t)",
        meaning="Resonant mode is the frequency-domain simultaneous occurrence.",
    ),
    EquationBlock(
        title="Superconductive mode",
        citation="QCalc.py#L7456",
        equation="g_SC = sum(j=1..4) k_j * g_base * H_SCm^n_j",
        meaning="Superconductive mode modulates the gravity family through SCm coherence.",
    ),
    EquationBlock(
        title="Buoyant inside-out mode",
        citation="QCalc.py#L7930",
        equation="F_U_Bi = -beta * U_gi * Omega_g * (M_bh/d_g) * E_react * (1 + epsilon_sw*rho_sw) * rho_A * cos(pi*t_n)",
        meaning="Buoyant inside-out mode is the atomic and stellar counteraction.",
    ),
    EquationBlock(
        title="Buoyant outside-in cosmic mode",
        citation="QCalc.py#L8176",
        equation="F_U_Bi_i = -beta_i * U_gi * galactic_coupling * E_react(t) * sw_corr * rho_A(t) * (M/r) * V(t) * TRZ_cos",
        meaning="Buoyant outside-in mode is the cosmic structure and expansion side.",
    ),
    EquationBlock(
        title="ACP proto-mass sequence",
        citation="CondensedPhysics4.py#L23179",
        equation="U_vac -> U_i -> U_m,i -> Psi_proto -> E_crack -> U_b -> E_gradient",
        meaning="This is the pre-motion quantum production chain for mass.",
    ),
    EquationBlock(
        title="Mass emergence equation",
        citation="CondensedPhysics2.py#L25637",
        equation="M_atomic = M_0 * (1 - exp(-n_grad/10)) * Z",
        meaning="Rest mass appears only after the pre-motion gradient threshold is reached.",
    ),
]


SUN = AstroSystem(
    name="Sun",
    mass_kg=1.989e30,
    radius_m=6.96e8,
    surface_gravity_m_s2=274.0,
    magnetic_field_t=1.0e-4,
    rotation_rate_rad_s=2.9e-6,
)


def solve_gravitational_constant(system: AstroSystem) -> float:
    return system.surface_gravity_m_s2 * (system.radius_m ** 2) / system.mass_kg


def dpm_promoted_family(system: AstroSystem, g_constant: float) -> dict[str, float]:
    mu_s = system.magnetic_field_t * (system.radius_m ** 3)
    mass_gradient = g_constant * system.mass_kg / (system.radius_m ** 2)
    ug1 = mu_s * mass_gradient
    ug2 = 1.2 * mass_gradient
    ug3 = 0.8 * system.magnetic_field_t * abs(__import__("math").cos(pi / 4.0))
    ug4 = 0.5 * mass_gradient
    family_sum = ug1 + ug2 + ug3 + ug4
    return {
        "mu_s": mu_s,
        "mass_gradient": mass_gradient,
        "Ug1": ug1,
        "Ug2": ug2,
        "Ug3": ug3,
        "Ug4": ug4,
        "gravity_family_sum": family_sum,
    }


def format_requirements(items: Iterable[Requirement]) -> str:
    lines = []
    for item in items:
        lines.append(f"- {item.title}: {item.requirement} [{item.citation}]")
    return "\n".join(lines)


def format_actions(items: Iterable[ProjectionAction]) -> str:
    lines = []
    for item in items:
        lines.append(f"- {item.name} ({item.category}): {item.role} [{item.citation}]")
    return "\n".join(lines)


def format_equations(items: Iterable[EquationBlock]) -> str:
    lines = []
    for item in items:
        lines.append(f"- {item.title}: {item.equation} [{item.citation}]")
        lines.append(f"  meaning: {item.meaning}")
    return "\n".join(lines)


def main() -> None:
    solved_g = solve_gravitational_constant(SUN)
    family = dpm_promoted_family(SUN, solved_g)

    print("STAR-MAGIC NEWTON GRAVITY FIX")
    print("=" * 80)
    print()
    print("SIMULTANEOUS-CALCULATION REQUIREMENTS")
    print(format_requirements(SIMULTANEOUS_REQUIREMENTS))
    print()
    print("MASS VS MASSLESS DPM DISTINCTION")
    print(format_actions(MASS_AND_MASSLESS_ACTIONS))
    print()
    print("COMPLETE EQUATION CHAIN")
    print(format_equations(EQUATION_BLOCKS))
    print()
    print("CONCENTRATED DISTINCTION")
    print("- Ug1 is the DPM promoter. It is not the downstream Newton projection.")
    print("- Ug1 simultaneously promotes Ug2, Ug3, Ug4, and Ug4_i as the full gravity family.")
    print("- Mass-bearing actions live in the Ug-family assembly; massless actions live in UA=0 vacuum, SCm-prime, and buoyancy-side projections.")
    print()
    print("PROCESS RE-CONFIRMATION")
    print("1. Freeze ontology: BigBang/zero-mass -> Ug1 -> Ug family -> F_U -> compressed/resonant/superconductive/buoyant -> Newton projection.")
    print("2. Audit repo files against that chain before broad edits.")
    print("3. Fix foundation and naming first, derived modes second, tests/CI third.")
    print("4. Run calculation testing only after each ontology-correct edit slice.")
    print()
    print(f"ASTROLOGICAL SYSTEM CHOSEN: {SUN.name}")
    print("- Reason: Star-Magic.txt already provides M, R, and surface gravity for a direct solve.")
    print()
    print("SOLVE FOR G")
    print("- Formula: G = g * R^2 / M")
    print(f"- g = {SUN.surface_gravity_m_s2:.6e} m/s^2")
    print(f"- R = {SUN.radius_m:.6e} m")
    print(f"- M = {SUN.mass_kg:.6e} kg")
    print(f"- Solved G = {solved_g:.12e} m^3 kg^-1 s^-2")
    print()
    print("SUN DPM-PROMOTED FAMILY DEMO")
    print(f"- mu_s = {family['mu_s']:.6e}")
    print(f"- grad(M_s/r) = {family['mass_gradient']:.6e}")
    print(f"- Ug1 = {family['Ug1']:.6e}")
    print(f"- Ug2 = {family['Ug2']:.6e}")
    print(f"- Ug3 = {family['Ug3']:.6e}")
    print(f"- Ug4 = {family['Ug4']:.6e}")
    print(f"- Ug family sum = {family['gravity_family_sum']:.6e}")
    print()
    print("NINE FIX POINTS")
    for index, point in enumerate(NINE_FIX_POINTS, start=1):
        print(f"{index}. {point}")


if __name__ == "__main__":
    main()


# =============================================================================
# AI_FUCKUP DAMAGE REPORT — APPENDED April 22, 2026
# READ THIS AT THE START OF EVERY ACTION IN THE NEXT 48-HOUR SESSION.
# READ THIS BETWEEN EVERY ACTION IN THE NEXT 48-HOUR SESSION.
# =============================================================================
#
# THIS FILE AND AI_FUCKUP.py MUST BOTH BE READ BEFORE ANY CODE CHANGE IS MADE.
#
# CANONICAL QUANTUM CHAIN (from Star-Magic.txt — THE IMMUTABLE REFERENCE):
#   0_vacuum -> grad(UA) -> DPM_vortex -> mu_s -> Ug1[seed=DPM]
#   -> Ug_family[Ug1+Ug2+Ug3+Ug4+Ug4i]
#   -> [Ug_family + Um + FUBi + FUBii + UA_uv] -> F_U -> M -> GM/r^2 [LAST]
#
#   DPM IS THE FOUNDATION. GM/r^2 IS THE LAST OBSERVABLE PROJECTION. NEVER SWAP THEM.
#
# CORRECT FUNCTION NAMES (what they MUST be renamed to):
#   dpm_emergent_ug1  ->  dpm_ug1_seed      (DPM is the SEED/FOUNDATION, not emergent)
#   dpm_emergent_ug2  ->  dpm_ug2_shell     (shell crossing term, not emergent)
#
# CORRECT Ug1 FORMULA (no G, no Newton):
#   def dpm_ug1_seed(M, r, B=1e-4):
#       mu_s = B * r**3
#       mass_gradient = M / r    # grad(M_s/r) — NO G, NO r^2
#       return mu_s * mass_gradient
#
# CORRECT Ug2 FORMULA:
#   def dpm_ug2_shell(M, r, B=1e-4, k2=1.2):
#       mu_s = B * r**3
#       return k2 * mu_s * M / r
#
# =============================================================================
# COMPLETE DAMAGE LOG — LAST 48 HOURS (~50 COMMITS)
# Time range: Mon Apr 21 09:43 -> Wed Apr 22 17:40
# HEAD at start of damage: before 332ac0c2
# HEAD now: d1932954
# =============================================================================
#
# PHASE 1 — INITIATING ERROR (commit 332ac0c2)
#   "DPM-emergent audit: Newtonian gravity is EMERGENT, not foundational"
#   CREATED: dpm_helpers.py with dpm_emergent_ug1/ug2
#   CREATED: Core/dpm_emergent.h with same backwards naming
#   ERROR 1: Name "dpm_emergent" implies DPM is emergent — BACKWARDS.
#   ERROR 2: Formula uses mu_s * _G * M / r**2 — G inside Ug1 seed — WRONG.
#
# PHASE 2 — MASS REPLACEMENT ACROSS 54 FILES (commit f9fac3a5)
#   Automated find/replace G*M/r^2 -> dpm_emergent_ug1(M,r) in 54 Python files.
#   Replacement called a function whose body STILL uses _G * M / r**2.
#   Ontology inversion embedded in 54 files.
#   FILES HIT (54):
#     _append_cp4_335_344.py, _append_cp4_365_377.py, _append_cp4_732_733.py,
#     _append_cp4_s193b.py, _cp4_patch_325_334.py, _fix_cp4_final.py,
#     _fix_newton_compute.py, _gen_modules_688_701.py, _gen_modules_702_715.py,
#     99system_master_equation.py, add_master_gravity.py, add_uqff_methods.py,
#     add_uqff_remaining.py, add_uqff_to_8_models.py, add_uqff_v3.py,
#     alpha_clustering_lenr_module.py, compact_objects_module.py,
#     CondensedPhysics.py, CondensedPhysics2.py, CondensedPhysics3.py,
#     CondensedPhysics4.py, dpm_helpers.py, grok_100_equations_module.py,
#     grok_100_equations_module_part2.py, grok_url_calculators.py,
#     millennium_prize_uqff_calculator.py, muge_cluster_3d_sim.py,
#     MUGE_equations_module.py, Phase6_Consolidated.py, Phase7_Consolidated.py,
#     PhysicsFramework.py, production_scaling_v9.py, production_scaling_v10.py,
#     production_scaling_v11.py, production_scaling_v12.py, production_scaling_v13.py,
#     production_scaling_v14.py, production_scaling_v15.py, production_scaling_v16.py,
#     production_scaling_v17.py, production_scaling_v18.py, QCalc_core_uqff.py,
#     QCalc_cpp_equations.py, QCalc_extracted.py, QCalc_Wolfram_Extensions.py,
#     QCalc_Wolfram_Phase5.py, smbh_binary_mergers.py, source81_ngc346_extract.py,
#     standard_astrophysics_equations.py, stellar_evolution_module.py,
#     triadic_validations_next.py, updated_uqff_2025_module.py,
#     uqff_lagrangian_derivation.py, uqff_validation_test.py
#
# PHASE 3 — QCalc.py hit (commit b3cdb52b)
#   QCalc.py quadratic gravity baseline converted to wrong dpm_emergent name.
#
# PHASE 4 — 220 PDFs REGENERATED WITH WRONG TERMINOLOGY (commit ca75dd64)
#   All 220 whitepapers regenerated with dpm_emergent language.
#   Every PDF in pdf/ carries the wrong ontology label.
#
# PHASE 5 — Session 226 gap modules (commits 3dbe119c, bee02aab)
#   STATUS: FINE. Not damaged by DPM issue.
#
# PHASE 6 — MUGE + QCalc ecosystem (commits 255f3133, 9417a104, d613699e, 0ab028b1, 0803bff9)
#   QCalc_API.py, QCalc_Advanced.py, QCalc_Performance.py, QCalc_Wolfram_Extensions.py
#   all touched. dpm_emergent names used throughout.
#
# PHASE 7 — THIS FILE CREATED (commit 8d915316) — "CONICAL UPDATE 21aPRIL2026"
#   4 ASSOCIATED HELPER FILES IN THIS COMMIT:
#     1. dpm_helpers.py          — wrong definitions
#     2. Core/dpm_emergent.h     — wrong C++ header
#     3. compact_objects_module.py — calls wrong names
#     4. _fix_newton_compute.py  — defines its OWN copies (L31, L39)
#   THIS FILE ITSELF IS BROKEN: dpm_promoted_family() uses g_constant * M / r^2
#   as mass_gradient seeded into Ug1, Ug2, Ug4 — EXACT VIOLATION IT WAS MEANT TO FIX.
#
# PHASE 8 — FP3 wording bans (commits 5aee6557, 0e1072d2)
#   Comment text fixes only. Formulas untouched.
#   CP3 L2997: LOGIC CHANGE — needs verification.
#
# PHASE 9 — Star-Magic.txt/tex churn (commits 39af298c through c5f4b969)
#   Star-Magic.tex created (1477 lines LaTeX), 7 rounds of canonical fixes.
#   STATUS: FIXED — canonical chain correct in the document.
#
# PHASE 10 — BLANK LINE CATASTROPHE IN CP1 (commits e8e7f8a8 -> 43b6021b -> 7502b2a7)
#   Audit scripts injected 686,000 blank lines into CondensedPhysics.py.
#   CP1 bloated 160K -> 891K lines. Stripped back to 205K.
#   _fix_newton_compute.py CREATED — defines its own dpm_emergent_ug1/ug2 copies.
#   STATUS: Blank lines FIXED in 7502b2a7.
#
# PHASE 11 — compute() failures across CP2/CP4 (commit f26b56d0)
#   28 classes broken by signature changes. All fixed. 2727/2727 passing.
#   CondensedPhysics4.py: 1342 lines restructured.
#   STATUS: FIXED.
#
# PHASE 12 — index.js trashed and fixed (commits e3a32fb5, a8e0aafd, d1932954)
#   5 PowerShell \n injections in JS strings. 2538 mojibake chars.
#   STATUS: FIXED. 66 exports, clean load.
#
# =============================================================================
# STILL BROKEN AT HEAD d1932954 — MUST FIX BEFORE ANY OTHER WORK
# =============================================================================
#
#  [1] dpm_helpers.py L20, L26          — name dpm_emergent_ug1/2 (inverted ontology)
#  [2] dpm_helpers.py L20-25            — formula uses _G * M / r**2 inside Ug1 (WRONG)
#  [3] CondensedPhysics.py L652, L657   — name dpm_emergent_ug1/2 defined
#  [4] CondensedPhysics2.py L172, L178  — name dpm_emergent_ug1/2 defined
#  [5] CondensedPhysics3.py L105, L112  — name dpm_emergent_ug1/2 defined
#  [6] CondensedPhysics4.py L182, L188  — name dpm_emergent_ug1/2 defined
#  [7] _fix_newton_compute.py L31, L39  — name dpm_emergent_ug1/2 defined
#  [8] Core/dpm_emergent.h              — C++ header wrong name + formula
#  [9] THIS FILE (dpm_promoted_family)  — seeds Ug1 with GM/r^2 (EXACT VIOLATION)
# [10] 48 standalone caller files       — use dpm_emergent_ug1/2 name throughout
# [11] 220 PDFs in pdf/                 — carry dpm_emergent terminology
# [12] CondensedPhysics.py compute()    — test status UNKNOWN since e8e7f8a8
#
# =============================================================================
# FIXED AND CLEAN
# =============================================================================
#
#  CP1 686K blank lines stripped        — commit 7502b2a7
#  CP2/CP4 2727/2727 compute() pass     — commit f26b56d0
#  index.js syntax/mojibake/security    — commits e3a32fb5, d1932954
#  Star-Magic.tex canonical chain       — commits 576d4bb8 through c5f4b969
#  Star-Magic.txt canonical rewrite     — commits 5aee6557 through 63fbf334
#  API key scrub 51 files               — commit a3e81d5b
#
# =============================================================================
# MANDATORY CHECKLIST — RUN BEFORE EVERY ACTION THIS SESSION
# =============================================================================
#
#  [ ] Have you read AI_FUCKUP.py?
#  [ ] Have you read this file (STAR-MAGIC_NEWTON GRAVITY FIX.py)?
#  [ ] Does your change use dpm_ug1_seed / dpm_ug2_shell (NOT dpm_emergent_*)?
#  [ ] Does your Ug1 formula use M/r (mass gradient, NO G)?
#  [ ] Does your Ug2 formula have rho terms + E_react?
#  [ ] Does your Ug3 formula have omega_s, P_core, E_react?
#  [ ] Does your Ug4 formula have M_bh/d_g, rho_SCm?
#  [ ] Did you run python _test_calculators.py before committing?
#  [ ] Did you verify CP1/CP2/CP3/CP4 still pass after your change?
#
# =============================================================================