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