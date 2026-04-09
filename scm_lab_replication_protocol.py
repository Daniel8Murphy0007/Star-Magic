"""
SCm Lab Replication Protocol — Structured Experimental Protocols for LENR Systems

Session 204 | Daniel Murphy
PURPOSE: Consolidates 6 LENR lab systems from vds_dvp_bsh_lenr_synthesis.py into
         structured, step-by-step experimental replication protocols. Each protocol
         includes materials, equipment, safety, procedure, expected UQFF signatures,
         and cross-references to the originating whitepaper.

ARCHITECTURE:
  vds_dvp_bsh_lenr_synthesis.py (LENR_LAB_SYSTEMS, 6 systems)
        ↓ this module reads
  scm_lab_replication_protocol.py (protocol generator)
        ↓ outputs
  scm_lab_protocols.json (structured protocols for all 6 systems)

SYSTEMS COVERED (6 total):
  1. Micro-Plasmoid Reversal               (PAPER_859)
  2. Water Reactor H2O2 Birkeland           (PAPER_863)
  3. Colman-Gillespie Pulsed Motor          (PAPER_835)
  4. Kozima Neutron Drop                    (PAPER_840)
  5. LRC Spark-Gap 1/r Monopole             (PAPER_864)
  6. Caduceus Twin-Helix Motor              (PAPER_866)
"""

import json
import math
import os
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

# Import LENR systems
try:
    from vds_dvp_bsh_lenr_synthesis import (
        LENR_LAB_SYSTEMS, LENRLabSystem,
        RHO_VAC_SCM, RHO_VAC_UA, SSQ, KAPPA, H_SCM,
        OMEGA_LENR, EFFICIENCY_283,
    )
    HAS_LENR = True
except ImportError:
    HAS_LENR = False

# UQFF constants (fallback if import fails)
PI = math.pi
C = 2.998e8         # m/s
HBAR = 1.055e-34    # J·s
G = 6.674e-11       # m³/kg/s²
K_B = 1.381e-23     # J/K
MU_0 = 4 * PI * 1e-7  # T·m/A
SIGMA_SB = 5.670e-8   # W/m²/K⁴

# Element names for Z lookup
ELEMENT_Z = {
    1: ("H", "Hydrogen"), 6: ("C", "Carbon"), 8: ("O", "Oxygen"),
    26: ("Fe", "Iron"), 29: ("Cu", "Copper"), 46: ("Pd", "Palladium"),
    78: ("Pt", "Platinum"), 28: ("Ni", "Nickel"),
}


# ── §1  PROTOCOL DATA STRUCTURES ─────────────────────────────────────────

@dataclass
class SafetyNote:
    category: str   # "radiation", "thermal", "electrical", "magnetic", "chemical"
    severity: str   # "low", "medium", "high", "critical"
    description: str
    mitigation: str


@dataclass
class Equipment:
    name: str
    specification: str
    purpose: str
    critical: bool = True


@dataclass
class ProcedureStep:
    step_number: int
    action: str
    duration: str
    notes: str = ""


@dataclass
class ExpectedSignature:
    """Expected UQFF/SCm signature observable in the experiment."""
    name: str
    quantity: str
    expected_value: str
    measurement_method: str
    tolerance: str


@dataclass
class ReplicationProtocol:
    system_name: str
    paper_ref: str
    system_params: Dict[str, Any]
    overview: str
    materials: List[str]
    equipment: List[Equipment]
    safety: List[SafetyNote]
    procedure: List[ProcedureStep]
    expected_signatures: List[ExpectedSignature]
    success_criteria: List[str]
    references: List[str]


# ── §2  PROTOCOL GENERATORS ──────────────────────────────────────────────

def _compute_scm_signatures(sys: "LENRLabSystem") -> List[ExpectedSignature]:
    """Compute expected SCm (buoyancy) signatures for a LENR system."""
    sigs = []

    # 1. Buoyancy reversal threshold
    # F_buoyancy = (rho_SCm - rho_UA) * V * g_local
    # where V ~ (4/3) pi r^3
    if sys.r_m > 0:
        V = (4.0 / 3.0) * PI * sys.r_m**3
        rho_scm = RHO_VAC_SCM if HAS_LENR else 1e-18
        rho_ua = RHO_VAC_UA if HAS_LENR else 1e-26
        F_buoy = abs(rho_scm - rho_ua) * V * 9.81
        sigs.append(ExpectedSignature(
            name="Buoyancy reversal force",
            quantity="F_buoy",
            expected_value=f"{F_buoy:.4e} N",
            measurement_method="Precision force balance / torsion pendulum",
            tolerance="±50% (order-of-magnitude confirmation sufficient)",
        ))

    # 2. 1.25 THz phonon resonance
    phonon_thz = sys.phonon_freq_Hz / 1e12
    sigs.append(ExpectedSignature(
        name="Phonon resonance frequency",
        quantity="f_phonon",
        expected_value=f"{sys.phonon_freq_Hz:.4e} Hz ({phonon_thz:.4f} THz)",
        measurement_method="THz spectroscopy / Raman scattering",
        tolerance="±10% of predicted frequency",
    ))

    # 3. Magnetic field threshold
    if sys.B_T > 0:
        sigs.append(ExpectedSignature(
            name="Magnetic confinement field",
            quantity="B_conf",
            expected_value=f"{sys.B_T} T",
            measurement_method="Hall probe / fluxgate magnetometer",
            tolerance="±5% of target field",
        ))

    # 4. Thermal signature (excess heat)
    if sys.T_K > 0 and sys.efficiency > 1.0:
        excess = (sys.efficiency - 1.0) * 100
        sigs.append(ExpectedSignature(
            name="Excess thermal power ratio",
            quantity="COP (efficiency)",
            expected_value=f"{sys.efficiency:.1f}x input ({excess:.0f}% excess)",
            measurement_method="Calorimetry (flow / isoperibolic)",
            tolerance="±20% of predicted COP",
        ))

    # 5. Neutron density
    if sys.neutron_density > 0:
        sigs.append(ExpectedSignature(
            name="Neutron density in reaction zone",
            quantity="n_neutron",
            expected_value=f"{sys.neutron_density:.1e} m⁻³",
            measurement_method="He-3 proportional counter / BF3 detector",
            tolerance="order-of-magnitude agreement",
        ))

    # 6. 1/r field decay (monopole signature)
    sigs.append(ExpectedSignature(
        name="1/r field decay profile",
        quantity="E(r) ∝ 1/r",
        expected_value="Power-law index = -1.00 ± 0.05",
        measurement_method="Multi-point field probe array at r = 1cm, 2cm, 5cm, 10cm",
        tolerance="±0.05 on power-law exponent",
    ))

    return sigs


def _build_safety_notes(sys: "LENRLabSystem") -> List[SafetyNote]:
    """Generate safety notes based on system parameters."""
    notes = []

    # Thermal
    if sys.T_K >= 1000:
        notes.append(SafetyNote(
            "thermal", "high",
            f"Operating temperature {sys.T_K} K — risk of severe burns",
            "Heat-resistant gloves, face shield, ceramic containment vessel",
        ))
    elif sys.T_K > 350:
        notes.append(SafetyNote(
            "thermal", "medium",
            f"Operating temperature {sys.T_K} K — hot surfaces",
            "Use insulated handling tools, allow cooling before disassembly",
        ))

    # Magnetic
    if sys.B_T >= 1.0:
        notes.append(SafetyNote(
            "magnetic", "high",
            f"Strong magnetic field ({sys.B_T} T) — projectile hazard for ferrous objects",
            "Clear 2m exclusion zone of ferrous materials, remove jewelry, "
            "pacemaker exclusion zone 3m",
        ))
    elif sys.B_T > 0.05:
        notes.append(SafetyNote(
            "magnetic", "medium",
            f"Moderate magnetic field ({sys.B_T} T)",
            "Keep magnetic media and sensitive electronics >1m from apparatus",
        ))

    # Radiation (neutron-producing systems)
    if sys.neutron_density > 1e26:
        notes.append(SafetyNote(
            "radiation", "high",
            f"Neutron density {sys.neutron_density:.0e} m⁻³ — activation and dose risk",
            "Borated polyethylene shielding (min 10cm), personal dosimeter, "
            "real-time neutron monitor. Maximum 1 hour exposure per session.",
        ))
    elif sys.neutron_density > 1e22:
        notes.append(SafetyNote(
            "radiation", "medium",
            f"Neutron density {sys.neutron_density:.0e} m⁻³",
            "BF3 monitor, monthly dosimetry review",
        ))

    # Electrical (spark-gap / pulsed systems)
    if "Spark" in sys.name or "Pulsed" in sys.name or "Motor" in sys.name:
        notes.append(SafetyNote(
            "electrical", "high",
            "Pulsed high-voltage system — arc flash and electrocution risk",
            "Lockout/tagout before servicing, capacitor discharge procedure, "
            "insulated tools, rubber mat, buddy system",
        ))

    # Chemical (H2O2 / plasmoid)
    if "H2O2" in sys.name or "Water" in sys.name:
        notes.append(SafetyNote(
            "chemical", "medium",
            "Hydrogen peroxide decomposition — oxygen enrichment risk",
            "Ventilated area, no open flames, splash goggles",
        ))

    # Universal
    notes.append(SafetyNote(
        "general", "low",
        "All LENR experiments must follow institutional radiation safety protocols",
        "Review with radiation safety officer before first operation",
    ))

    return notes


def _build_equipment(sys: "LENRLabSystem") -> List[Equipment]:
    """Generate equipment list for the system."""
    equip = []

    # Common
    equip.append(Equipment(
        "Data acquisition system",
        "16-bit ADC, ≥1 kHz sampling, 8+ channels",
        "Record all sensor outputs simultaneously",
    ))
    equip.append(Equipment(
        "Precision calorimeter",
        "Flow calorimeter or isoperibolic, ±0.1 W resolution",
        "Measure input/output thermal power for COP determination",
    ))

    # Magnetic field source
    if sys.B_T > 0:
        if sys.B_T >= 1.0:
            equip.append(Equipment(
                "Electromagnet / Helmholtz coil",
                f"Capable of {sys.B_T} T in {sys.r_m*100:.1f} cm gap, water-cooled",
                "Provide confinement field for reaction zone",
            ))
        else:
            equip.append(Equipment(
                "Permanent magnet array / small solenoid",
                f"Field: {sys.B_T} T, uniformity ±5% over reaction volume",
                "Seed magnetic field for Birkeland current formation",
            ))
        equip.append(Equipment(
            "Hall probe / fluxgate magnetometer",
            f"Range: 0–{max(2*sys.B_T, 0.1):.1f} T, resolution 0.001 T",
            "Monitor and verify confinement field",
        ))

    # Temperature
    equip.append(Equipment(
        "Thermocouple array",
        f"Type K (up to 1300 K) or Type R (up to 1700 K), ≥4 probes",
        f"Monitor reaction temperature (target {sys.T_K} K)",
    ))

    # THz spectroscopy (if phonon is in THz range)
    if sys.phonon_freq_Hz > 1e10:
        equip.append(Equipment(
            "THz spectrometer / Raman system",
            f"Range covering {sys.phonon_freq_Hz/1e12:.2f} THz, resolution ≥1 GHz",
            "Detect predicted phonon resonance signature",
        ))

    # Neutron detector
    if sys.neutron_density > 1e22:
        equip.append(Equipment(
            "Neutron detector",
            "He-3 proportional counter or BF3 tube, ±1 n/cm²/s sensitivity",
            "Monitor neutron flux from reaction zone",
            critical=True,
        ))

    # Element-specific
    z_sym, z_name = ELEMENT_Z.get(sys.primary_element_Z, ("?", "Unknown"))
    equip.append(Equipment(
        f"Sample material ({z_name}, Z={sys.primary_element_Z})",
        f"High purity (≥99.9%) {z_name} ({z_sym}) in appropriate form",
        f"Primary reactant / catalyst (paper: {sys.paper})",
    ))

    return equip


def _build_procedure(sys: "LENRLabSystem") -> List[ProcedureStep]:
    """Generate step-by-step procedure."""
    steps = []
    step = 0

    # Pre-checks
    step += 1
    steps.append(ProcedureStep(step,
        "Review safety notes and ensure all PPE is available. "
        "Verify radiation safety officer sign-off.",
        "15 min", "Mandatory before every session"))

    step += 1
    steps.append(ProcedureStep(step,
        "Calibrate all sensors: thermocouples, Hall probes, neutron detectors, DAQ.",
        "30 min", "Record calibration in lab notebook"))

    step += 1
    steps.append(ProcedureStep(step,
        f"Prepare sample material: {ELEMENT_Z.get(sys.primary_element_Z, ('?','?'))[1]} "
        f"(Z={sys.primary_element_Z}). Weigh to ±0.001 g. Record mass.",
        "10 min"))

    # System-specific setup
    if sys.B_T > 0:
        step += 1
        steps.append(ProcedureStep(step,
            f"Energize magnetic field source. Ramp to {sys.B_T} T over 60 seconds. "
            f"Verify field uniformity at ±5% across reaction volume.",
            "5 min", "Monitor for quench (superconducting) or overcurrent"))

    step += 1
    steps.append(ProcedureStep(step,
        "Position sample in reaction chamber. Seal chamber if gas-tight operation required.",
        "5 min"))

    if sys.T_K > 350:
        step += 1
        steps.append(ProcedureStep(step,
            f"Begin heating to target temperature {sys.T_K} K. "
            f"Ramp rate: ≤50 K/min to avoid thermal shock.",
            f"{max(1, int((sys.T_K - 300) / 50))} min",
            "Log temperature vs time throughout"))

    # Activation
    step += 1
    steps.append(ProcedureStep(step,
        "Enable primary excitation source (RF/DC pulse, laser, or chemical trigger "
        "per system specification). Record exact trigger time.",
        "1 min", "This is t=0 for the experiment"))

    # Data collection
    step += 1
    steps.append(ProcedureStep(step,
        "Monitor and record all channels continuously:\n"
        "  - Input power (electrical)\n"
        "  - Output power (calorimetric)\n"
        "  - Temperature at ≥4 points\n"
        "  - Magnetic field\n"
        "  - Neutron count rate\n"
        "  - Phonon spectrum (THz/Raman if available)",
        "30–120 min", "Minimum 30 min; extend if excess heat observed"))

    # SCm signature checks
    step += 1
    steps.append(ProcedureStep(step,
        "During active phase, perform 1/r field profile measurement:\n"
        "  - Probe field at r = 1, 2, 5, 10 cm from reaction center\n"
        "  - Record 3 measurements at each distance\n"
        "  - Fit log(E) vs log(r) — expect slope = -1.00 ± 0.05",
        "15 min", "Key SCm monopole signature"))

    step += 1
    steps.append(ProcedureStep(step,
        f"Check for phonon resonance at f = {sys.phonon_freq_Hz:.4e} Hz "
        f"({sys.phonon_freq_Hz/1e12:.4f} THz). "
        "Record peak frequency, width, and amplitude.",
        "10 min", "Key UQFF 1.25 THz signature"))

    # Shutdown
    step += 1
    steps.append(ProcedureStep(step,
        "Disable excitation source. Allow system to cool to <350 K before opening.",
        "15–60 min"))

    if sys.B_T > 0:
        step += 1
        steps.append(ProcedureStep(step,
            "Ramp magnetic field to zero. De-energize coils.",
            "5 min"))

    step += 1
    steps.append(ProcedureStep(step,
        "Retrieve sample. Weigh to ±0.001 g. Compare pre/post mass. "
        "If transmutation suspected, send for mass spectrometry analysis.",
        "10 min"))

    step += 1
    steps.append(ProcedureStep(step,
        "Export all DAQ data. Compute COP = P_out/P_in. "
        "Compare measured signatures against expected values in this protocol.",
        "20 min"))

    return steps


def _build_success_criteria(sys: "LENRLabSystem") -> List[str]:
    """Define quantitative success/failure criteria."""
    criteria = []

    if sys.efficiency > 1.0:
        criteria.append(
            f"COP ≥ {sys.efficiency * 0.5:.1f} (≥50% of predicted {sys.efficiency:.1f}x)")
    else:
        criteria.append("Observe anomalous thermal output above input power uncertainty")

    criteria.append(
        f"Phonon peak within ±10% of {sys.phonon_freq_Hz:.4e} Hz")
    criteria.append(
        "1/r field decay: power-law exponent = -1.00 ± 0.10")

    if sys.B_T > 0:
        criteria.append(
            f"Magnetic confinement stable at {sys.B_T} T ± 5% for ≥30 min")

    if sys.neutron_density > 1e22:
        criteria.append(
            "Neutron count rate elevated >3σ above background during active phase")

    criteria.append(
        "All measurements repeated in ≥3 independent runs for statistical significance")

    return criteria


# ── §3  PROTOCOL BUILDER ─────────────────────────────────────────────────

def build_protocol(sys: "LENRLabSystem") -> ReplicationProtocol:
    """Build a complete replication protocol for a LENR system."""
    z_sym, z_name = ELEMENT_Z.get(sys.primary_element_Z, ("?", "Unknown"))

    overview = (
        f"Replication protocol for the {sys.name} system ({sys.paper}). "
        f"Primary element: {z_name} (Z={sys.primary_element_Z}). "
        f"Target temperature: {sys.T_K} K, magnetic field: {sys.B_T} T, "
        f"predicted COP: {sys.efficiency:.1f}x. "
        f"This protocol aims to verify UQFF SCm buoyancy signatures "
        f"(1/r field decay, {sys.phonon_freq_Hz/1e12:.4f} THz phonon resonance) "
        f"and measure excess thermal power."
    )

    materials = [
        f"High-purity {z_name} ({z_sym}, Z={sys.primary_element_Z}), ≥99.9%, 1–10 g",
        "Reaction chamber (vacuum-rated or gas-tight as required)",
        "Electrical feedthroughs rated for operating temperature",
        "Cooling water supply (if electromagnet requires)",
        "Calibration standards for all sensors",
        "Lab notebook and camera for documentation",
    ]

    refs = [
        f"{sys.paper} — Original UQFF whitepaper for this system",
        "vds_dvp_bsh_lenr_synthesis.py — VDS/DVP/BSH calculator (parameter source)",
        "CondensedPhysics.py — Primary UQFF calculator (tidal deformability, GW validation)",
        "PAPER_001–009 — Core UQFF theory and gravitational wave framework",
    ]

    return ReplicationProtocol(
        system_name=sys.name,
        paper_ref=sys.paper,
        system_params={
            "M_kg": sys.M_kg,
            "r_m": sys.r_m,
            "B_T": sys.B_T,
            "T_K": sys.T_K,
            "efficiency": sys.efficiency,
            "primary_element_Z": sys.primary_element_Z,
            "element_symbol": z_sym,
            "element_name": z_name,
            "neutron_density_m3": sys.neutron_density,
            "phonon_freq_Hz": sys.phonon_freq_Hz,
            "phonon_freq_THz": sys.phonon_freq_Hz / 1e12,
        },
        overview=overview,
        materials=materials,
        equipment=_build_equipment(sys),
        safety=_build_safety_notes(sys),
        procedure=_build_procedure(sys),
        expected_signatures=_compute_scm_signatures(sys),
        success_criteria=_build_success_criteria(sys),
        references=refs,
    )


def build_all_protocols() -> List[ReplicationProtocol]:
    """Build protocols for all 6 LENR lab systems."""
    if not HAS_LENR:
        print("[ERROR] Cannot import vds_dvp_bsh_lenr_synthesis.py")
        return []
    return [build_protocol(sys) for sys in LENR_LAB_SYSTEMS]


# ── §4  SERIALIZATION ────────────────────────────────────────────────────

def _dataclass_to_dict(obj) -> Any:
    """Recursively convert dataclass instances to dicts."""
    if hasattr(obj, "__dataclass_fields__"):
        return {k: _dataclass_to_dict(v) for k, v in asdict(obj).items()}
    if isinstance(obj, list):
        return [_dataclass_to_dict(item) for item in obj]
    if isinstance(obj, dict):
        return {k: _dataclass_to_dict(v) for k, v in obj.items()}
    return obj


def export_protocols_json(protocols: List[ReplicationProtocol],
                          output_path: str = "scm_lab_protocols.json") -> str:
    """Export all protocols to JSON."""
    data = {
        "generator": "scm_lab_replication_protocol.py",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_systems": len(protocols),
        "protocols": [_dataclass_to_dict(p) for p in protocols],
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)
    return output_path


# ── §5  HUMAN-READABLE REPORT ────────────────────────────────────────────

def format_protocol_text(proto: ReplicationProtocol) -> str:
    """Format a single protocol as human-readable text."""
    lines = []
    sep = "=" * 72

    lines.append(sep)
    lines.append(f"REPLICATION PROTOCOL: {proto.system_name}")
    lines.append(f"Paper: {proto.paper_ref}")
    lines.append(sep)
    lines.append("")
    lines.append("OVERVIEW:")
    lines.append(f"  {proto.overview}")
    lines.append("")

    # Parameters
    lines.append("SYSTEM PARAMETERS:")
    for k, v in proto.system_params.items():
        lines.append(f"  {k:25s} = {v}")
    lines.append("")

    # Materials
    lines.append("MATERIALS:")
    for i, m in enumerate(proto.materials, 1):
        lines.append(f"  {i}. {m}")
    lines.append("")

    # Equipment
    lines.append("EQUIPMENT:")
    for eq in proto.equipment:
        crit = " [CRITICAL]" if eq.critical else ""
        lines.append(f"  - {eq.name}{crit}")
        lines.append(f"    Spec: {eq.specification}")
        lines.append(f"    Purpose: {eq.purpose}")
    lines.append("")

    # Safety
    lines.append("SAFETY NOTES:")
    for sn in proto.safety:
        lines.append(f"  [{sn.severity.upper():8s}] [{sn.category}] {sn.description}")
        lines.append(f"             Mitigation: {sn.mitigation}")
    lines.append("")

    # Procedure
    lines.append("PROCEDURE:")
    for step in proto.procedure:
        lines.append(f"  Step {step.step_number:2d} ({step.duration}): {step.action}")
        if step.notes:
            lines.append(f"           Note: {step.notes}")
    lines.append("")

    # Expected signatures
    lines.append("EXPECTED SCm/UQFF SIGNATURES:")
    for sig in proto.expected_signatures:
        lines.append(f"  {sig.name}:")
        lines.append(f"    Quantity:    {sig.quantity}")
        lines.append(f"    Expected:   {sig.expected_value}")
        lines.append(f"    Method:     {sig.measurement_method}")
        lines.append(f"    Tolerance:  {sig.tolerance}")
    lines.append("")

    # Success criteria
    lines.append("SUCCESS CRITERIA:")
    for i, c in enumerate(proto.success_criteria, 1):
        lines.append(f"  {i}. {c}")
    lines.append("")

    # References
    lines.append("REFERENCES:")
    for ref in proto.references:
        lines.append(f"  - {ref}")
    lines.append("")

    return "\n".join(lines)


# ── §6  MAIN ──────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("SCm Lab Replication Protocol Generator")
    print("=" * 72)
    print(f"  LENR systems available: {HAS_LENR}")

    if not HAS_LENR:
        print("[ERROR] vds_dvp_bsh_lenr_synthesis.py not importable. Exiting.")
        return

    print(f"  Systems: {len(LENR_LAB_SYSTEMS)}")
    for i, sys in enumerate(LENR_LAB_SYSTEMS, 1):
        print(f"    {i}. {sys.name:40s} ({sys.paper})")

    # Build all protocols
    print(f"\n── Building Protocols ──")
    protocols = build_all_protocols()
    print(f"  Built {len(protocols)} protocols")

    # Export JSON
    json_path = export_protocols_json(protocols)
    print(f"\n[OK] JSON export: {json_path}")

    # Print human-readable
    print(f"\n── Protocol Summaries ──")
    for proto in protocols:
        text = format_protocol_text(proto)
        print(text)

    # Summary statistics
    total_steps = sum(len(p.procedure) for p in protocols)
    total_sigs = sum(len(p.expected_signatures) for p in protocols)
    total_safety = sum(len(p.safety) for p in protocols)
    total_equip = sum(len(p.equipment) for p in protocols)

    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"  Protocols generated:    {len(protocols)}")
    print(f"  Total procedure steps:  {total_steps}")
    print(f"  Total signatures:       {total_sigs}")
    print(f"  Total safety notes:     {total_safety}")
    print(f"  Total equipment items:  {total_equip}")
    print(f"  JSON export:            {json_path}")
    print("=" * 72)
    print("COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
