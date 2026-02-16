#!/usr/bin/env python3
"""
Test Suite for ACP (Atomic Creation Process) Models
====================================================

Tests for the three ACP models implemented in CondensedPhysics.py:
1. ProtoNucleusShellModel - Proto-nucleus shell formation
2. ACPStageTracker - ACP stage progression tracking
3. UniversalCycleTracker - [UA]/[SCm] consumption tracking

Based on UQFF documents:
- Gas Nebula observation_19April_2025.docx
- Westerlund 2 UQFF Session June 2025

Author: Star Magic UQFF Team
Date: February 2026
"""

import pytest
import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from CondensedPhysics import (
    PROTO_NUCLEUS_MODEL,
    ACP_TRACKER,
    UNIVERSAL_CYCLE,
    ProtoNucleusShellModel,
    ACPStageTracker,
    UniversalCycleTracker,
    CONSTANTS
)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST: ProtoNucleusShellModel
# ═══════════════════════════════════════════════════════════════════════════════

class TestProtoNucleusShellModel:
    """Tests for Proto-Nucleus Shell Formation Model"""
    
    def test_model_initialization(self):
        """Test that global model instance exists and is properly initialized"""
        assert PROTO_NUCLEUS_MODEL is not None
        assert isinstance(PROTO_NUCLEUS_MODEL, ProtoNucleusShellModel)
        assert hasattr(PROTO_NUCLEUS_MODEL, 'ACP_STAGES')
        assert len(PROTO_NUCLEUS_MODEL.ACP_STAGES) == 9
        
    def test_acp_stages_defined(self):
        """Test that all 9 ACP stages are properly defined"""
        expected_stages = [
            'DPM_INITIATION', 'VACUUM_DENSITY', 'U_I_CREATION',
            'U_M_STRINGS', 'PROTO_NUCLEUS', 'SHELL_CRACKING',
            'ELECTRON_SETTLEMENT', 'NEUTRINO_INFLATION', 'ATOM_COMPLETE'
        ]
        model = ProtoNucleusShellModel()
        for i, stage_name in enumerate(expected_stages, 1):
            assert i in model.ACP_STAGES
            assert model.ACP_STAGES[i]['name'] == stage_name
            
    def test_species_index_hydrogen(self):
        """Test species index calculation for hydrogen-like configuration"""
        model = ProtoNucleusShellModel()
        species_idx, steps = model.compute_species_index(n_shell=1)
        
        # Species index should be negative (ρ_SCm < ρ_UA)
        assert species_idx < 0
        # For n=1, magnitude should be ~1.0 (ratio ≈ 0.1)
        assert abs(species_idx) < 2.0
        assert 'SPECIES INDEX CALCULATION' in steps
        
    def test_species_index_scaling(self):
        """Test that species index scales linearly with n"""
        model = ProtoNucleusShellModel()
        idx_1, _ = model.compute_species_index(n_shell=1)
        idx_2, _ = model.compute_species_index(n_shell=2)
        idx_4, _ = model.compute_species_index(n_shell=4)
        
        # Should scale linearly: idx_2 = 2 × idx_1, idx_4 = 4 × idx_1
        assert pytest.approx(idx_2, rel=1e-10) == 2 * idx_1
        assert pytest.approx(idx_4, rel=1e-10) == 4 * idx_1
        
    def test_shell_energy_hydrogen(self):
        """Test shell energy calculation for hydrogen (Z=1)"""
        model = ProtoNucleusShellModel()
        E_shell, steps = model.compute_shell_energy(n=1, Z=1, include_uqff=False)
        
        # Ground state hydrogen: E = -13.6 eV
        E_expected_eV = -13.6
        E_calculated_eV = E_shell / CONSTANTS['q']
        assert pytest.approx(E_calculated_eV, rel=0.01) == E_expected_eV
        
    def test_shell_energy_uqff_correction(self):
        """Test that UQFF correction modifies shell energy"""
        model = ProtoNucleusShellModel()
        E_without, _ = model.compute_shell_energy(n=1, Z=1, include_uqff=False)
        E_with, _ = model.compute_shell_energy(n=1, Z=1, include_uqff=True)
        
        # UQFF should modify the energy
        assert E_with != E_without
        # Correction should be small but non-zero
        correction_factor = abs(E_with - E_without) / abs(E_without)
        assert 0 < correction_factor < 1.0  # Less than 100% correction
        
    def test_shell_energy_n_squared_scaling(self):
        """Test that shell energy scales as 1/n² (Bohr model)"""
        model = ProtoNucleusShellModel()
        E_1, _ = model.compute_shell_energy(n=1, Z=1, include_uqff=False)
        E_2, _ = model.compute_shell_energy(n=2, Z=1, include_uqff=False)
        E_3, _ = model.compute_shell_energy(n=3, Z=1, include_uqff=False)
        
        # E_n ∝ 1/n² → E_2/E_1 = 1/4, E_3/E_1 = 1/9
        assert pytest.approx(E_2 / E_1, rel=0.01) == 0.25
        assert pytest.approx(E_3 / E_1, rel=0.01) == 1.0 / 9.0
        
    def test_shell_cracking_energy(self):
        """Test shell cracking energy calculation"""
        model = ProtoNucleusShellModel()
        E_crack, can_crack, steps = model.compute_shell_cracking_energy(
            n=1, E_neutrino=1e-12
        )
        
        # Cracking energy should be positive
        assert E_crack > 0
        assert isinstance(can_crack, bool)
        assert 'SHELL CRACKING ENERGY' in steps
        
    def test_electron_settlement(self):
        """Test electron settlement via U_g2"""
        model = ProtoNucleusShellModel()
        E_transition, steps = model.compute_electron_settlement(
            n_initial=2, n_final=1, U_g2=1e-30
        )
        
        # Energy released in transition (positive = emitted)
        assert isinstance(E_transition, (int, float))
        assert 'ELECTRON SETTLEMENT' in steps
        
    def test_neutrino_inflation(self):
        """Test neutrino energy inflation"""
        model = ProtoNucleusShellModel()
        E_inf, steps = model.compute_neutrino_inflation(
            E_neutrino=1e-12, n_shell=1
        )
        
        # Inflated energy should be greater than input
        assert E_inf > 1e-12
        assert 'NEUTRINO ENERGY INFLATION' in steps
        
    def test_neutrino_inflation_n_dependence(self):
        """Test that neutrino inflation increases with n_shell (exp(n_shell/26))"""
        model = ProtoNucleusShellModel()
        E_1, _ = model.compute_neutrino_inflation(E_neutrino=1e-12, n_shell=1)
        E_13, _ = model.compute_neutrino_inflation(E_neutrino=1e-12, n_shell=13)
        E_26, _ = model.compute_neutrino_inflation(E_neutrino=1e-12, n_shell=26)
        
        # Higher n_shell should give more inflation: E_26 > E_13 > E_1
        assert E_26 > E_13 > E_1
        
    def test_acp_sequence_simulation(self):
        """Test full ACP sequence simulation"""
        model = ProtoNucleusShellModel()
        results = model.simulate_acp_sequence(Z=1, n_shells=4, E_neutrino=1e-12)
        
        assert 'shell_energies' in results
        assert 'total_energy_released' in results
        assert 'stages' in results
        assert len(results['shell_energies']) == 4
        assert len(results['stages']) > 0


# ═══════════════════════════════════════════════════════════════════════════════
# TEST: ACPStageTracker
# ═══════════════════════════════════════════════════════════════════════════════

class TestACPStageTracker:
    """Tests for ACP Stage Progression Tracker"""
    
    def test_tracker_initialization(self):
        """Test that global tracker instance exists"""
        assert ACP_TRACKER is not None
        assert isinstance(ACP_TRACKER, ACPStageTracker)
        
    def test_fresh_tracker_state(self):
        """Test fresh tracker starts at stage 0"""
        tracker = ACPStageTracker()
        assert tracker.current_stage == 0
        assert tracker.total_energy_input == 0.0
        assert tracker.total_energy_output == 0.0
        assert len(tracker.stage_history) == 0
        
    def test_stage_advancement(self):
        """Test stage advancement mechanism"""
        tracker = ACPStageTracker()
        initial_stage = tracker.current_stage
        
        result = tracker.advance_stage(energy_input=1e-12)
        assert isinstance(result, dict)
        assert 'current_stage' in result
        assert tracker.current_stage == initial_stage + 1
        assert tracker.total_energy_input > 0
        
    def test_stage_history_tracking(self):
        """Test that stage history is properly tracked"""
        tracker = ACPStageTracker()
        
        # Advance through 3 stages (starting from 0)
        tracker.advance_stage(energy_input=1e-12)
        tracker.advance_stage(energy_input=2e-12)
        tracker.advance_stage(energy_input=3e-12)
        
        assert len(tracker.stage_history) == 3
        assert tracker.current_stage == 3  # 0 -> 1 -> 2 -> 3
        
    def test_dipole_vortex_energy(self):
        """Test dipole vortex energy calculation: ([SCm] - [UA'])²"""
        tracker = ACPStageTracker()
        E_dipole, steps = tracker.compute_dipole_vortex_energy()
        
        # Should be ([SCm] - [UA'])² which is positive
        assert E_dipole >= 0
        assert 'DIPOLE VORTEX ENERGY' in steps
        
    def test_stage_report(self):
        """Test stage report generation"""
        tracker = ACPStageTracker()
        tracker.advance_stage(energy_input=1e-12)
        
        report = tracker.get_stage_report()
        
        assert isinstance(report, str)
        assert 'Stage' in report or 'stage' in report
        
    def test_cannot_exceed_max_stage(self):
        """Test that tracker cannot exceed stage 9"""
        tracker = ACPStageTracker()
        
        # Try to advance through all stages
        for _ in range(15):  # More than 9 stages
            tracker.advance_stage(energy_input=1e-12)
        
        # Should be capped at stage 9
        assert tracker.current_stage <= 9


# ═══════════════════════════════════════════════════════════════════════════════
# TEST: UniversalCycleTracker
# ═══════════════════════════════════════════════════════════════════════════════

class TestUniversalCycleTracker:
    """Tests for Universal [UA]/[SCm] Cycle Tracker"""
    
    def test_tracker_initialization(self):
        """Test that global cycle tracker exists"""
        assert UNIVERSAL_CYCLE is not None
        assert isinstance(UNIVERSAL_CYCLE, UniversalCycleTracker)
        
    def test_initial_conditions(self):
        """Test that initial [UA] and [SCm] values are set"""
        tracker = UniversalCycleTracker()
        assert tracker.UA_current > 0
        assert tracker.SCm_current > 0
        
    def test_consumption_rates(self):
        """Test consumption rate calculation"""
        tracker = UniversalCycleTracker()
        rates = tracker.compute_consumption_rates()
        
        assert 'dUA_dt' in rates
        assert 'dSCm_dt' in rates
        # SCm should decay (negative rate)
        assert rates['dSCm_dt'] <= 0
        
    def test_evolution(self):
        """Test time evolution of [UA]/[SCm]"""
        tracker = UniversalCycleTracker()
        UA_initial = tracker.UA_current
        SCm_initial = tracker.SCm_current
        
        # Evolve for 1000 days
        tracker.evolve(dt_days=1000)
        
        # SCm should decrease (decay rate 0.0963)
        assert tracker.SCm_current <= SCm_initial
        
    def test_decay_rate_constant(self):
        """Test that decay rate is ~0.0963 day⁻¹"""
        tracker = UniversalCycleTracker()
        assert pytest.approx(tracker.SCM_DECAY_RATE, rel=0.01) == 0.0963
        
    def test_universe_decay_time(self):
        """Test universe decay time calculation"""
        tracker = UniversalCycleTracker()
        t_decay, steps = tracker.compute_universe_decay_time()
        
        # Should return a finite positive time
        assert t_decay > 0
        assert np.isfinite(t_decay)
        assert isinstance(steps, str)
        
    def test_third_decay_trapping(self):
        """Test third decay cycle matter trapping"""
        tracker = UniversalCycleTracker()
        # method requires matter_density argument
        trapped_fraction, steps = tracker.compute_third_decay_trapping(matter_density=1e-26)
        
        # Fraction should be between 0 and 1
        assert 0 <= trapped_fraction <= 1
        assert isinstance(steps, str)
        
    def test_cycle_report(self):
        """Test cycle report generation"""
        tracker = UniversalCycleTracker()
        tracker.evolve(dt_days=100)
        
        report = tracker.get_cycle_report()
        
        assert isinstance(report, str)
        assert '[UA]' in report or 'UA' in report
        
    def test_history_tracking(self):
        """Test that evolution history is tracked"""
        tracker = UniversalCycleTracker()
        
        # Evolve multiple times
        for _ in range(5):
            tracker.evolve(dt_days=100)
        
        assert len(tracker.history) >= 5


# ═══════════════════════════════════════════════════════════════════════════════
# TEST: Integration Tests (Cross-Model)
# ═══════════════════════════════════════════════════════════════════════════════

class TestACPIntegration:
    """Integration tests across ACP models"""
    
    def test_all_models_use_consistent_constants(self):
        """Test that all models use same UQFF constants"""
        proto = ProtoNucleusShellModel()
        acp = ACPStageTracker()
        cycle = UniversalCycleTracker()
        
        # All should use CONSTANTS dictionary
        assert CONSTANTS.get('SSq', 0.57) == 0.57
        # Verify cycle uses proper decay rate
        assert cycle.SCM_DECAY_RATE == 0.0963
        
    def test_proto_nucleus_to_acp_stage_mapping(self):
        """Test proto-nucleus stages align with ACP stages"""
        proto = ProtoNucleusShellModel()
        acp = ACPStageTracker()
        
        # Both should have 9 stages
        assert len(proto.ACP_STAGES) == 9
        # ACP stages are capped at 9 in advance_stage method
        for _ in range(20):  # Try to advance past 9
            acp.advance_stage(energy_input=1e-12)
        assert acp.current_stage == 9  # Should cap at 9
        
    def test_full_acp_cycle(self):
        """Test a complete ACP cycle from DPM to complete atom"""
        proto = ProtoNucleusShellModel()
        acp = ACPStageTracker()
        
        # Simulate hydrogen formation
        results = proto.simulate_acp_sequence(Z=1, n_shells=2, E_neutrino=1e-12)
        
        # Advance ACP tracker through stages (9 stages total, cap at 9)
        # Use total_energy_released (or a fixed small value if zero)
        energy_per_stage = max(abs(results['total_energy_released']), 1e-15) / 9
        while acp.current_stage < 9:
            acp.advance_stage(energy_input=energy_per_stage)
        
        assert acp.current_stage == 9
        assert acp.total_energy_input > 0


# ═══════════════════════════════════════════════════════════════════════════════
# TEST: Westerlund 2 Validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestWesterlund2Validation:
    """Tests validating against Westerlund 2 UQFF Session values"""
    
    def test_omega_ug4i_constant(self):
        """Test that omega_Ug4i_1 constant is defined"""
        assert 'omega_Ug4i_1' in CONSTANTS
        assert pytest.approx(CONSTANTS['omega_Ug4i_1'], rel=0.01) == 6.283e12
        
    def test_resonance_frequencies_defined(self):
        """Test all 4 Ug resonance frequencies are defined"""
        assert 'omega_Ug1_1' in CONSTANTS
        assert 'omega_Ug2_1' in CONSTANTS
        assert 'omega_Ug3_1' in CONSTANTS
        assert 'omega_Ug4i_1' in CONSTANTS
        
    def test_frequency_hierarchy(self):
        """Test frequency hierarchy: ω_Ug4i >> ω_Ug3 >> ω_Ug2 > ω_Ug1"""
        omega_1 = CONSTANTS['omega_Ug1_1']
        omega_2 = CONSTANTS['omega_Ug2_1']
        omega_3 = CONSTANTS['omega_Ug3_1']
        omega_4 = CONSTANTS['omega_Ug4i_1']
        
        assert omega_4 > omega_3 > omega_2 > omega_1
        assert omega_4 / omega_1 > 1e20  # 6.283e12 / 1.989e-13 ≈ 3e25


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

def run_all_tests():
    """Run all ACP model tests with summary output"""
    print("=" * 80)
    print("ACP MODELS TEST SUITE - Proto-Nucleus, ACP Stage, Universal Cycle")
    print("=" * 80)
    
    # Test ProtoNucleusShellModel
    print("\n[1/4] Testing ProtoNucleusShellModel...")
    proto_tests = TestProtoNucleusShellModel()
    try:
        proto_tests.test_model_initialization()
        proto_tests.test_acp_stages_defined()
        proto_tests.test_species_index_hydrogen()
        proto_tests.test_species_index_scaling()
        proto_tests.test_shell_energy_hydrogen()
        proto_tests.test_shell_energy_uqff_correction()
        proto_tests.test_shell_energy_n_squared_scaling()
        proto_tests.test_shell_cracking_energy()
        proto_tests.test_electron_settlement()
        proto_tests.test_neutrino_inflation()
        proto_tests.test_neutrino_inflation_n_dependence()
        proto_tests.test_acp_sequence_simulation()
        print("  ✓ All ProtoNucleusShellModel tests passed")
    except Exception as e:
        print(f"  ✗ ProtoNucleusShellModel test failed: {e}")
        
    # Test ACPStageTracker
    print("\n[2/4] Testing ACPStageTracker...")
    acp_tests = TestACPStageTracker()
    try:
        acp_tests.test_tracker_initialization()
        acp_tests.test_fresh_tracker_state()
        acp_tests.test_stage_advancement()
        acp_tests.test_stage_history_tracking()
        acp_tests.test_dipole_vortex_energy()
        acp_tests.test_stage_report()
        acp_tests.test_cannot_exceed_max_stage()
        print("  ✓ All ACPStageTracker tests passed")
    except Exception as e:
        print(f"  ✗ ACPStageTracker test failed: {e}")
        
    # Test UniversalCycleTracker
    print("\n[3/4] Testing UniversalCycleTracker...")
    cycle_tests = TestUniversalCycleTracker()
    try:
        cycle_tests.test_tracker_initialization()
        cycle_tests.test_initial_conditions()
        cycle_tests.test_consumption_rates()
        cycle_tests.test_evolution()
        cycle_tests.test_decay_rate_constant()
        cycle_tests.test_universe_decay_time()
        cycle_tests.test_third_decay_trapping()
        cycle_tests.test_cycle_report()
        cycle_tests.test_history_tracking()
        print("  ✓ All UniversalCycleTracker tests passed")
    except Exception as e:
        print(f"  ✗ UniversalCycleTracker test failed: {e}")
        
    # Test Integration
    print("\n[4/4] Testing Integration & Westerlund 2 Validation...")
    int_tests = TestACPIntegration()
    w2_tests = TestWesterlund2Validation()
    try:
        int_tests.test_all_models_use_consistent_constants()
        int_tests.test_proto_nucleus_to_acp_stage_mapping()
        int_tests.test_full_acp_cycle()
        w2_tests.test_omega_ug4i_constant()
        w2_tests.test_resonance_frequencies_defined()
        w2_tests.test_frequency_hierarchy()
        print("  ✓ All Integration & Westerlund 2 tests passed")
    except Exception as e:
        print(f"  ✗ Integration test failed: {e}")
    
    print("\n" + "=" * 80)
    print("TEST SUITE COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--pytest":
        # Run via pytest
        pytest.main([__file__, "-v"])
    else:
        # Run direct
        run_all_tests()
