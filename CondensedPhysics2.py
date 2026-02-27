#!/usr/bin/env python3
"""
CondensedPhysics2.py - UQFF Extended Calculators Module
========================================================

Pipeline extension of CondensedPhysics.py for scalable calculator additions.
This module imports base infrastructure from CondensedPhysics.py and provides
new Calculator classes starting with UFT Orb Analysis_10/11.

ARCHITECTURE:
    CondensedPhysics.py  → Foundation (1011 base classes, frozen)
    CondensedPhysics2.py → Extension 1 (Orb Analysis + new systems)
    CondensedPhysics3.py → Extension 2 (future overflow)
    
    __init__.py aggregates all exports into unified API

CAPACITY: ~500-600 Calculator classes (~80-100K lines)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: February 27, 2026
"""

import math
import numpy as np
from typing import Dict, List, Optional, Any

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORT BASE INFRASTRUCTURE FROM CONDENSEDHYSICS.PY
# ═══════════════════════════════════════════════════════════════════════════════

# Import shared constants and base classes
from CondensedPhysics import (
    SCIPY_AVAILABLE,
    SYMPY_AVAILABLE,
    PHYSICS_FRAMEWORK_AVAILABLE,
)

# NumPy 2.0 compatibility
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid


# ═══════════════════════════════════════════════════════════════════════════════
# CONDENSEDPHYSICS2 REGISTRY TRACKING
# ═══════════════════════════════════════════════════════════════════════════════

CP2_VERSION = "2.0.0"
CP2_CLASS_COUNT = 0  # Updated dynamically


# ============================================================================
# UFT ORB ANALYSIS_10 CALCULATORS (8 Calculator Classes)
# Source: Grok UFT Orb Analysis_10 (March 4, 2025)
# Photos: #35-#37 of 496 Red Dwarf Reactor Plasma Orb infrared images
# Physics: 36 frames (~1.08 s), cyclical convection pattern:
#          lower right (#35) → upper left (#36) → upper right (#37)
#          ~1.62 J total energy, non-local "spooky action" via [Ug3] and [Um]
# ============================================================================

# UFT Orb Analysis_10 Parameters
ORB_ANALYSIS_10_PARAMS = {
    'r_reactor': 0.0889,           # m (3.5 in diameter / 2)
    'M_s': 0.5e-3,                 # kg (plasma orb mass 0.5 g)
    'omega_s': 2 * math.pi * 6000, # rad/s (6000 Hz field resonance)
    'T_base': 366,                 # K (bulb base temperature)
    'T_top': 288,                  # K (ambient top temperature)
    'B_s': 1e-3,                   # T (magnetic field from H2 bubbles)
    'SCm': 1e15,                   # kg/m³ (hypothetical density at base)
    'UA': 1e-11,                   # C (trapped Aether charge)
    'dt_frame': 0.03,              # s (frame interval)
    'n_frames_orb10': 36,          # total frames (#1-#34 minus #14, +#35-#37)
    'total_time_orb10': 1.08,      # s (36 frames × 0.03 s)
    'E_total_orb10': 1.62,         # J (total energy ~45 mJ/frame × 36)
    'E_react': 1e15,               # W/m³ (reactivity with thermal decay)
    'gamma_decay': 0.001,          # decay constant
    'v_plasmoid': 0.5,             # m/s (plasmoid motion speed)
    'n_H2_bubbles': 15,            # average hydrogen bubbles (12-18)
}


class ThirtySixFrameSequenceCalculator:
    """
    UFT Orb Analysis_10: 36-frame sequence energy calculator.
    
    Extends the video analysis to 36 frames (Photos #1-#34 minus #14, plus #35-#37),
    covering ~1.08 s at 33.3 fps. Calculates cumulative energy and temporal evolution
    of plasmoid dynamics through the extended sequence.
    
    Physics:
        E_total = Σᵢ[E_frame_i] = 36 × 45 mJ ≈ 1.62 J
        t_total = n_frames × dt = 36 × 0.03 s = 1.08 s
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_10_PARAMS['n_frames_orb10'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_10_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.045)  # J (~45 mJ/frame)
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_10_PARAMS['gamma_decay'])
        
        # Total sequence time
        t_total = n_frames * dt
        
        # Cumulative energy with decay correction
        E_cumulative = 0.0
        frame_energies = []
        for i in range(n_frames):
            t_n = i * dt
            E_frame = E_per_frame * math.exp(-gamma * t_n)
            E_cumulative += E_frame
            frame_energies.append({
                'frame': i + 1,
                't': round(t_n, 4),
                'E': round(E_frame, 6)
            })
        
        # Average power
        P_avg = E_cumulative / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_cumulative': round(E_cumulative, 6),
            'P_avg': round(P_avg, 6),
            'E_per_frame_avg': round(E_cumulative / n_frames, 6) if n_frames > 0 else 0.0,
            'frame_energies': frame_energies[:5],  # First 5 for preview
            'equation': 'E_total = Σᵢ[E_frame_i × e^(-γtᵢ)]',
            'source': 'Grok UFT Orb Analysis_10 36-Frame Sequence (March 4, 2025)'
        }


class CyclicalConvectionPatternCalculator:
    """
    UFT Orb Analysis_10: Cyclical convection pattern tracking calculator.
    
    Tracks the cyclical concentration shifts: lower right (#35) → upper left (#36)
    → upper right (#37). Models the spatial redistribution as a convection cycle
    driven by thermal gradients and non-local effects.
    
    Physics:
        Cycle_path: LR → UL → UR (3-point convection loop)
        θ_rotation = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]
        Cycle_period ≈ 3 frames × 0.03 s = 0.09 s
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Cyclical quadrant positions (normalized 0-1 coordinates)
        cycle_sequence = dataset.get('cycle_sequence', [
            {'frame': 35, 'quadrant': 'lower_right', 'x_c': 0.75, 'y_c': 0.25},
            {'frame': 36, 'quadrant': 'upper_left', 'x_c': 0.25, 'y_c': 0.75},
            {'frame': 37, 'quadrant': 'upper_right', 'x_c': 0.75, 'y_c': 0.75},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_10_PARAMS['dt_frame'])
        
        # Track cycle transitions
        transitions = []
        total_angle = 0.0
        total_distance = 0.0
        
        for i in range(1, len(cycle_sequence)):
            prev = cycle_sequence[i - 1]
            curr = cycle_sequence[i]
            
            dx = curr['x_c'] - prev['x_c']
            dy = curr['y_c'] - prev['y_c']
            distance = math.sqrt(dx**2 + dy**2)
            theta = math.atan2(dy, dx) * 180 / math.pi
            
            total_angle += theta
            total_distance += distance
            
            transitions.append({
                'from_frame': prev['frame'],
                'to_frame': curr['frame'],
                'transition': f"{prev['quadrant']} → {curr['quadrant']}",
                'dx': round(dx, 4),
                'dy': round(dy, 4),
                'distance': round(distance, 4),
                'theta_deg': round(theta, 2)
            })
        
        # Cycle characteristics
        n_transitions = len(transitions)
        cycle_period = n_transitions * dt
        avg_angular_velocity = total_angle / cycle_period if cycle_period > 0 else 0.0
        
        return {
            'cycle_sequence': cycle_sequence,
            'transitions': transitions,
            'total_distance': round(total_distance, 4),
            'total_angle_deg': round(total_angle, 2),
            'cycle_period': round(cycle_period, 4),
            'avg_angular_velocity': round(avg_angular_velocity, 2),
            'pattern': 'LR → UL → UR (cyclical convection)',
            'equation': 'θ_cycle = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]',
            'source': 'Grok UFT Orb Analysis_10 Cyclical Convection (March 4, 2025)'
        }


class Orb10RefinedFUCalculator:
    """
    UFT Orb Analysis_10: Refined Unified Field F_U calculator.
    
    Computes the complete unified field equation refined with Photos #35-#37 data,
    extending the temporal coverage to 1.08 s with 36 frames.
    
    Physics:
        F_U = Σᵢ[kᵢ·Ugᵢ(r,t,Mₛ,ωₛ,Tₛ,Bₛ,SCm,UA,tₙ) - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react]
              + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ]
              + (gμν + η·Tₛμν(UA,SCm,ρA))
              
    Parameters updated for 36-frame sequence:
        r = 0.0889 m, Mₛ = 0.5 g, ωₛ = 2π×6000 rad/s
        Tₛ = 366→288 K, Bₛ = 10⁻³ T
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r = dataset.get('r_reactor', ORB_ANALYSIS_10_PARAMS['r_reactor'])
        M_s = dataset.get('M_s', ORB_ANALYSIS_10_PARAMS['M_s'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_10_PARAMS['omega_s'])
        T_base = dataset.get('T_base', ORB_ANALYSIS_10_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_10_PARAMS['T_top'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_10_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_10_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_10_PARAMS['UA'])
        t = dataset.get('t', 0.54)  # midpoint of 1.08 s sequence
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_10_PARAMS['gamma_decay'])
        E_react = dataset.get('E_react', ORB_ANALYSIS_10_PARAMS['E_react'])
        
        t_n = t  # normalized time
        
        # Ug_1: Internal dipole
        k_1 = 1.5e-4
        Ug1 = k_1 * (M_s / r) * math.exp(-gamma * t) * math.cos(math.pi * t_n) * (1 + 0.01 * math.sin(gamma * t))
        
        # Ug_2: Outer field bubble
        k_2 = 1.2
        Ug2 = k_2 * (UA + UA) * M_s / (r**2) * SCm * math.exp(-gamma * t)
        
        # Ug_3: Magnetic strings
        k_3 = 1.8
        Ug3 = k_3 * B_s * math.cos(omega_s * t * math.pi) * SCm * math.exp(-gamma * t)
        
        # Ub_i: Universal buoyancy
        beta_i = 0.8
        Omega_g = 7.3e-16  # rad/s
        M_bh = 8.15e36     # kg (reference black hole mass)
        d_g = 2.55e20      # m (galactic distance)
        Ubi = -beta_i * (Ug1 + Ug2 + Ug3) * Omega_g * (M_bh / d_g) * E_react * math.cos(math.pi * t_n)
        
        # Um: Universal magnetism
        mu_j = 1e-4
        Um = (mu_j / r) * (1 - math.exp(-gamma * t * math.cos(math.pi * t_n))) * SCm * math.exp(-gamma * t)
        
        # A_μν: Cosmic Aether tensor contribution
        eta = 1e-22
        rho_A = 1e-23
        A_munu = 1.0 + eta * (UA * SCm * rho_A) * t_n
        
        # Total F_U
        F_U = Ug1 + Ug2 + Ug3 + Ubi + Um + A_munu
        
        return {
            'Ug1': Ug1,
            'Ug2': Ug2,
            'Ug3': Ug3,
            'Ubi': Ubi,
            'Um': Um,
            'A_munu': A_munu,
            'F_U_total': F_U,
            'parameters': {
                'r': r, 'M_s': M_s, 'omega_s': omega_s,
                'T_gradient': f'{T_base}→{T_top} K',
                'B_s': B_s, 'SCm': SCm, 'UA': UA, 't': t
            },
            'equation': 'F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] + (gμν + η·Tₛμν)',
            'source': 'Grok UFT Orb Analysis_10 Refined F_U (March 4, 2025)'
        }


class SpookyActionNonLocalTransferCalculator:
    """
    UFT Orb Analysis_10: Spooky action non-local energy transfer calculator.
    
    Models the non-local "spooky action" observed in cyclical shifts, where plasmoid
    concentrations jump across quadrants faster than classical diffusion would predict.
    Driven by [Ug3] (magnetic strings) and [Um] (universal magnetism).
    
    Physics:
        v_apparent = Δr / Δt (exceeds local diffusion limit)
        E_nonlocal = ε·[Ug3 + Um]·cos(ωt)
        Correlation_length = λ·(SCm/UA)^(1/2)
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        dt = dataset.get('dt_frame', ORB_ANALYSIS_10_PARAMS['dt_frame'])
        v_plasmoid = dataset.get('v_plasmoid', ORB_ANALYSIS_10_PARAMS['v_plasmoid'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_10_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_10_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_10_PARAMS['UA'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_10_PARAMS['omega_s'])
        t = dataset.get('t', 0.54)
        
        # Classical diffusion limit (for comparison)
        D_classical = 1e-9  # m²/s (typical for viscous oil)
        l_diffusion = math.sqrt(2 * D_classical * dt)  # diffusion length in one frame
        
        # Actual displacement observed
        delta_r_observed = v_plasmoid * dt
        
        # Ratio (how much faster than diffusion)
        speedup_ratio = delta_r_observed / l_diffusion if l_diffusion > 0 else float('inf')
        
        # Ug3 contribution (magnetic strings)
        k_3 = 1.8
        gamma = ORB_ANALYSIS_10_PARAMS['gamma_decay']
        Ug3 = k_3 * B_s * math.cos(omega_s * t * math.pi) * SCm * math.exp(-gamma * t)
        
        # Um contribution (universal magnetism)
        mu_j = 1e-4
        r = ORB_ANALYSIS_10_PARAMS['r_reactor']
        Um = (mu_j / r) * (1 - math.exp(-gamma * t * math.cos(math.pi * t))) * SCm * math.exp(-gamma * t)
        
        # Non-local energy
        epsilon = 1e-10
        E_nonlocal = epsilon * (Ug3 + Um) * math.cos(omega_s * t)
        
        # Correlation length
        lambda_corr = 1e-6  # m (characteristic length scale)
        correlation_length = lambda_corr * math.sqrt(SCm / UA) if UA > 0 else 0.0
        
        return {
            'v_plasmoid': v_plasmoid,
            'delta_r_per_frame': round(delta_r_observed, 6),
            'l_diffusion_classical': round(l_diffusion, 9),
            'speedup_vs_diffusion': round(speedup_ratio, 2),
            'Ug3': Ug3,
            'Um': Um,
            'E_nonlocal': E_nonlocal,
            'correlation_length': correlation_length,
            'mechanism': 'Ug3 (magnetic strings) + Um (universal magnetism)',
            'equation': 'E_nonlocal = ε·[Ug3 + Um]·cos(ωt)',
            'source': 'Grok UFT Orb Analysis_10 Spooky Action (March 4, 2025)'
        }


class ThermalGradientDrivenDynamicsCalculator:
    """
    UFT Orb Analysis_10: Thermal gradient-driven dynamics calculator.
    
    Models plasmoid dynamics driven by the thermal gradient (366K base → 288K top).
    The 78 K temperature difference creates convective forces that drive the cyclical
    concentration patterns observed in Photos #35-#37.
    
    Physics:
        ΔT = T_base - T_top = 366 - 288 = 78 K
        F_thermal = ρ·β·ΔT·V·g
        v_thermal ∝ √(β·ΔT·g·L)
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        T_base = dataset.get('T_base', ORB_ANALYSIS_10_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_10_PARAMS['T_top'])
        L = dataset.get('L_reactor', 0.254)  # m (10 in height)
        
        # Oil properties
        rho_oil = dataset.get('rho_oil', 850)  # kg/m³
        beta_thermal = dataset.get('beta_thermal', 7e-4)  # K⁻¹
        g = 9.81  # m/s²
        
        Delta_T = T_base - T_top
        
        # Thermal-driven velocity estimate
        v_thermal = math.sqrt(beta_thermal * Delta_T * g * L)
        
        # Grashof number (ratio of buoyancy to viscous forces)
        nu_oil = dataset.get('nu_oil', 1e-5)  # m²/s
        Gr = (g * beta_thermal * Delta_T * L**3) / (nu_oil**2)
        
        # Thermal force on plasmoid
        V_plasmoid = (4/3) * math.pi * (0.001)**3  # 1 mm radius sphere
        F_thermal = rho_oil * beta_thermal * Delta_T * V_plasmoid * g
        
        # Temperature profile (linear approximation)
        z_positions = [0.0, 0.127, 0.254]  # bottom, middle, top
        T_profile = [T_base - (T_base - T_top) * z / L for z in z_positions]
        
        return {
            'Delta_T': Delta_T,
            'v_thermal': round(v_thermal, 4),
            'Grashof_number': Gr,
            'F_thermal': F_thermal,
            'T_profile': [{'z': z, 'T': round(T, 1)} for z, T in zip(z_positions, T_profile)],
            'gradient': round(Delta_T / L, 2),
            'equation': 'v_thermal ∝ √(β·ΔT·g·L)',
            'source': 'Grok UFT Orb Analysis_10 Thermal Dynamics (March 4, 2025)'
        }


class QuadrantTransitionTrackerCalculator:
    """
    UFT Orb Analysis_10: Quadrant transition tracker calculator.
    
    Tracks all quadrant transitions across the 36-frame sequence, identifying
    patterns in plasmoid concentration redistribution. Maps transitions to
    identify dominant flow paths and cycle frequencies.
    
    Physics:
        Transition_matrix[i,j] = P(quadrant_j | quadrant_i)
        Dominant_path = argmax(transition_counts)
        Cycle_frequency = n_complete_cycles / t_total
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Recent quadrant sequence from photos
        quadrant_sequence = dataset.get('quadrant_sequence', [
            # Orb_8: Photos #30-#31
            {'frame': 30, 'quadrant': 'upper_right'},
            {'frame': 31, 'quadrant': 'lower_left'},
            # Orb_9: Photos #32-#34
            {'frame': 32, 'quadrant': 'upper_left'},
            {'frame': 33, 'quadrant': 'upper_right'},
            {'frame': 34, 'quadrant': 'upper_right'},
            # Orb_10: Photos #35-#37
            {'frame': 35, 'quadrant': 'lower_right'},
            {'frame': 36, 'quadrant': 'upper_left'},
            {'frame': 37, 'quadrant': 'upper_right'},
        ])
        
        # Count transitions
        transition_counts = {}
        for i in range(1, len(quadrant_sequence)):
            prev_q = quadrant_sequence[i-1]['quadrant']
            curr_q = quadrant_sequence[i]['quadrant']
            key = f"{prev_q} → {curr_q}"
            transition_counts[key] = transition_counts.get(key, 0) + 1
        
        # Find dominant transition
        dominant = max(transition_counts.items(), key=lambda x: x[1]) if transition_counts else ('none', 0)
        
        # Total transitions
        total_transitions = sum(transition_counts.values())
        
        # Quadrant visit counts
        quadrant_visits = {}
        for entry in quadrant_sequence:
            q = entry['quadrant']
            quadrant_visits[q] = quadrant_visits.get(q, 0) + 1
        
        return {
            'quadrant_sequence': quadrant_sequence,
            'transition_counts': transition_counts,
            'dominant_transition': dominant[0],
            'dominant_count': dominant[1],
            'total_transitions': total_transitions,
            'quadrant_visits': quadrant_visits,
            'most_visited': max(quadrant_visits.items(), key=lambda x: x[1])[0] if quadrant_visits else 'none',
            'equation': 'P(Q_j|Q_i) = count(Q_i→Q_j) / count(Q_i)',
            'source': 'Grok UFT Orb Analysis_10 Quadrant Tracker (March 4, 2025)'
        }


class ACEDCEModulatedEnergyCalculator:
    """
    UFT Orb Analysis_10: ACE/DCE modulated energy calculator.
    
    Calculates energy contributions from [ACE/DCE] (Aether Creation/Destruction 
    Concentration Events) modulated by [SCm] reactivity and 65W bulb input.
    Brightness fluctuations (~1 mJ/spot) are tied to these events.
    
    Physics:
        E_ACE = k_ACE × [SCm] × [UA] × (1 - e^(-γt))
        E_DCE = k_DCE × [SCm] × [UA] × e^(-γt)
        E_total_events = Σᵢ[E_ACE_i + E_DCE_i]
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_10_PARAMS['n_frames_orb10'])
        E_per_spot = dataset.get('E_per_spot', 1e-3)  # J (1 mJ per spot)
        n_spots_avg = dataset.get('n_spots_avg', 45)  # average spots per frame
        SCm = dataset.get('SCm', ORB_ANALYSIS_10_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_10_PARAMS['UA'])
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_10_PARAMS['gamma_decay'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_10_PARAMS['dt_frame'])
        P_bulb = dataset.get('P_bulb', 65)  # W
        
        k_ACE = 1e-30  # coupling constant
        k_DCE = 0.8e-30
        
        E_ACE_total = 0.0
        E_DCE_total = 0.0
        event_log = []
        
        for i in range(n_frames):
            t = i * dt
            
            # ACE (creation events - increase with time)
            E_ACE = k_ACE * SCm * UA * (1 - math.exp(-gamma * t))
            E_ACE_total += E_ACE
            
            # DCE (destruction events - decrease with time)
            E_DCE = k_DCE * SCm * UA * math.exp(-gamma * t)
            E_DCE_total += E_DCE
            
            if i < 5:  # Log first 5 frames
                event_log.append({
                    'frame': i + 1,
                    'E_ACE': E_ACE,
                    'E_DCE': E_DCE,
                    'ratio': round(E_ACE / E_DCE, 4) if E_DCE > 0 else float('inf')
                })
        
        # Observable energy from spots
        E_spots_observable = n_frames * n_spots_avg * E_per_spot
        
        return {
            'n_frames': n_frames,
            'E_ACE_total': E_ACE_total,
            'E_DCE_total': E_DCE_total,
            'E_events_total': E_ACE_total + E_DCE_total,
            'E_spots_observable': E_spots_observable,
            'n_spots_avg': n_spots_avg,
            'event_log': event_log,
            'P_bulb_input': P_bulb,
            'equation': 'E_total = Σᵢ[E_ACE + E_DCE] with E_ACE ∝ (1-e^(-γt)), E_DCE ∝ e^(-γt)',
            'source': 'Grok UFT Orb Analysis_10 ACE/DCE Energy (March 4, 2025)'
        }


class MagneticBubbleConfinementCalculator:
    """
    UFT Orb Analysis_10: Magnetic bubble confinement calculator.
    
    Models how the 12-18 hydrogen bubbles provide magnetic confinement (~10⁻³ T)
    that anchors plasmoid paths and enhances [SCm]/[UA] reactivity. The even
    spacing of bubbles creates field uniformity.
    
    Physics:
        B_total = Σⱼ[Bⱼ(r)] (superposition of bubble fields)
        E_confinement = B²/(2μ₀) × V
        τ_confinement = r²/(D × β_confinement)
        
    Source: Grok UFT Orb Analysis_10 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_bubbles = dataset.get('n_H2_bubbles', ORB_ANALYSIS_10_PARAMS['n_H2_bubbles'])
        B_per_bubble = dataset.get('B_s', ORB_ANALYSIS_10_PARAMS['B_s'])
        r_reactor = dataset.get('r_reactor', ORB_ANALYSIS_10_PARAMS['r_reactor'])
        L_reactor = dataset.get('L_reactor', 0.254)  # m
        
        mu_0 = 4 * math.pi * 1e-7  # H/m
        
        # Total field (simplified superposition at center)
        B_total = B_per_bubble * math.sqrt(n_bubbles)  # statistical enhancement
        
        # Confinement volume (cylindrical reactor core)
        V_confinement = math.pi * r_reactor**2 * L_reactor
        
        # Magnetic energy density
        u_magnetic = B_total**2 / (2 * mu_0)
        
        # Total magnetic energy
        E_confinement = u_magnetic * V_confinement
        
        # Confinement time estimate
        D_plasmoid = 1e-9  # m²/s (effective diffusion)
        beta_confinement = 0.1  # field effectiveness factor
        tau_confinement = r_reactor**2 / (D_plasmoid * beta_confinement) if beta_confinement > 0 else float('inf')
        
        # Bubble spacing (assumed uniform)
        bubble_spacing = L_reactor / n_bubbles
        
        return {
            'n_bubbles': n_bubbles,
            'B_per_bubble': B_per_bubble,
            'B_total_center': round(B_total, 6),
            'V_confinement': round(V_confinement, 8),
            'u_magnetic': round(u_magnetic, 6),
            'E_confinement': round(E_confinement, 9),
            'tau_confinement': round(tau_confinement, 2),
            'bubble_spacing': round(bubble_spacing, 4),
            'field_uniformity': 'enhanced by even bubble spacing',
            'equation': 'E_conf = B²/(2μ₀)×V, τ = r²/(D×β)',
            'source': 'Grok UFT Orb Analysis_10 Bubble Confinement (March 4, 2025)'
        }


# UFT Orb Analysis_10 registry dict
ORB_ANALYSIS_10_CALCULATORS = {
    'ThirtySixFrameSequenceCalculator': ThirtySixFrameSequenceCalculator(),
    'CyclicalConvectionPatternCalculator': CyclicalConvectionPatternCalculator(),
    'Orb10RefinedFUCalculator': Orb10RefinedFUCalculator(),
    'SpookyActionNonLocalTransferCalculator': SpookyActionNonLocalTransferCalculator(),
    'ThermalGradientDrivenDynamicsCalculator': ThermalGradientDrivenDynamicsCalculator(),
    'QuadrantTransitionTrackerCalculator': QuadrantTransitionTrackerCalculator(),
    'ACEDCEModulatedEnergyCalculator': ACEDCEModulatedEnergyCalculator(),
    'MagneticBubbleConfinementCalculator': MagneticBubbleConfinementCalculator(),
}


# ============================================================================
# UFT ORB ANALYSIS_11 CALCULATORS (9 Calculator Classes)
# Source: Grok UFT Orb Analysis_11 (March 4, 2025)
# Photos: #38-#40 of 496 Red Dwarf Reactor Plasma Orb infrared images
# Physics: 39 frames (~1.17 s), counter-clockwise diagonal cycle:
#          upper left (#38) → upper right (#39) → lower left (#40)
#          ~1.755 J total energy, intelligent plasmoid behavior
# ============================================================================

# UFT Orb Analysis_11 Parameters
ORB_ANALYSIS_11_PARAMS = {
    'r_reactor': 0.0889,           # m (3.5 in diameter / 2)
    'M_s': 0.5e-3,                 # kg (plasma orb mass 0.5 g)
    'omega_s': 2 * math.pi * 6000, # rad/s (6000 Hz field resonance)
    'T_base': 366,                 # K (bulb base temperature)
    'T_top': 288,                  # K (ambient top temperature)
    'B_s': 1e-3,                   # T (magnetic field from H2 bubbles)
    'SCm': 1e15,                   # kg/m³ (hypothetical density at base)
    'UA': 1e-11,                   # C (trapped Aether charge)
    'dt_frame': 0.03,              # s (frame interval)
    'n_frames_orb11': 39,          # total frames (#1-#37 minus #14, +#38-#40)
    'total_time_orb11': 1.17,      # s (39 frames × 0.03 s)
    'E_total_orb11': 1.755,        # J (total energy ~45 mJ/frame × 39)
    'E_react': 1e15,               # W/m³ (reactivity with thermal decay)
    'gamma_decay': 0.001,          # decay constant
    'v_plasmoid': 0.5,             # m/s (plasmoid motion speed)
    'n_H2_bubbles': 15,            # average hydrogen bubbles (12-18)
}


class ThirtyNineFrameSequenceCalculator:
    """
    UFT Orb Analysis_11: 39-frame sequence energy calculator.
    
    Extends the video analysis to 39 frames (Photos #1-#37 minus #14, plus #38-#40),
    covering ~1.17 s at 33.3 fps. Calculates cumulative energy and temporal evolution
    of plasmoid dynamics through the extended sequence.
    
    Physics:
        E_total = Σᵢ[E_frame_i] = 39 × 45 mJ ≈ 1.755 J
        t_total = n_frames × dt = 39 × 0.03 s = 1.17 s
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_11_PARAMS['n_frames_orb11'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_11_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.045)  # J (~45 mJ/frame)
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_11_PARAMS['gamma_decay'])
        
        # Total sequence time
        t_total = n_frames * dt
        
        # Cumulative energy with decay correction
        E_cumulative = 0.0
        frame_energies = []
        for i in range(n_frames):
            t_n = i * dt
            E_frame = E_per_frame * math.exp(-gamma * t_n)
            E_cumulative += E_frame
            frame_energies.append({
                'frame': i + 1,
                't': round(t_n, 4),
                'E': round(E_frame, 6)
            })
        
        # Average power
        P_avg = E_cumulative / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_cumulative': round(E_cumulative, 6),
            'P_avg': round(P_avg, 6),
            'E_per_frame_avg': round(E_cumulative / n_frames, 6) if n_frames > 0 else 0.0,
            'frame_energies': frame_energies[:5],  # First 5 for preview
            'equation': 'E_total = Σᵢ[E_frame_i × e^(-γtᵢ)]',
            'source': 'Grok UFT Orb Analysis_11 39-Frame Sequence (March 4, 2025)'
        }


class CounterClockwiseDiagonalCycleCalculator:
    """
    UFT Orb Analysis_11: Counter-clockwise diagonal cycle calculator.
    
    Tracks the counter-clockwise diagonal cycle: upper left (#38) → upper right (#39)
    → lower left (#40). This pattern suggests a rotational convection mode with
    diagonal transitions across quadrants.
    
    Physics:
        Cycle_path: UL → UR → LL (counter-clockwise diagonal)
        θ_total = Σᵢ[θᵢ] (accumulated rotation angle)
        ω_cycle = θ_total / (n_transitions × dt)
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Diagonal cycle positions (normalized 0-1 coordinates)
        cycle_sequence = dataset.get('cycle_sequence', [
            {'frame': 38, 'quadrant': 'upper_left', 'x_c': 0.25, 'y_c': 0.75},
            {'frame': 39, 'quadrant': 'upper_right', 'x_c': 0.75, 'y_c': 0.75},
            {'frame': 40, 'quadrant': 'lower_left', 'x_c': 0.25, 'y_c': 0.25},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_11_PARAMS['dt_frame'])
        
        # Track diagonal transitions
        transitions = []
        total_angle = 0.0
        total_distance = 0.0
        
        for i in range(1, len(cycle_sequence)):
            prev = cycle_sequence[i - 1]
            curr = cycle_sequence[i]
            
            dx = curr['x_c'] - prev['x_c']
            dy = curr['y_c'] - prev['y_c']
            distance = math.sqrt(dx**2 + dy**2)
            theta = math.atan2(dy, dx) * 180 / math.pi
            
            total_angle += theta
            total_distance += distance
            
            # Classify transition type
            if abs(dx) > 0.3 and abs(dy) > 0.3:
                trans_type = 'diagonal'
            elif abs(dx) > 0.3:
                trans_type = 'horizontal'
            else:
                trans_type = 'vertical'
            
            transitions.append({
                'from_frame': prev['frame'],
                'to_frame': curr['frame'],
                'transition': f"{prev['quadrant']} → {curr['quadrant']}",
                'type': trans_type,
                'dx': round(dx, 4),
                'dy': round(dy, 4),
                'distance': round(distance, 4),
                'theta_deg': round(theta, 2)
            })
        
        # Cycle angular velocity
        n_transitions = len(transitions)
        cycle_period = n_transitions * dt
        omega_cycle = total_angle / cycle_period if cycle_period > 0 else 0.0
        
        return {
            'cycle_sequence': cycle_sequence,
            'transitions': transitions,
            'total_distance': round(total_distance, 4),
            'total_angle_deg': round(total_angle, 2),
            'cycle_period': round(cycle_period, 4),
            'omega_cycle': round(omega_cycle, 2),
            'rotation_sense': 'counter-clockwise',
            'pattern': 'UL → UR → LL (diagonal counter-clockwise)',
            'equation': 'ω_cycle = θ_total / (n × dt)',
            'source': 'Grok UFT Orb Analysis_11 Diagonal Cycle (March 4, 2025)'
        }


class Orb11RefinedFUCalculator:
    """
    UFT Orb Analysis_11: Refined Unified Field F_U calculator.
    
    Computes the complete unified field equation refined with Photos #38-#40 data,
    extending the temporal coverage to 1.17 s with 39 frames.
    
    Physics:
        F_U = Σᵢ[kᵢ·Ugᵢ(r,t,Mₛ,ωₛ,Tₛ,Bₛ,SCm,UA,tₙ) - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react]
              + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ]
              + (gμν + η·Tₛμν(UA,SCm,ρA))
              
    Parameters updated for 39-frame sequence:
        r = 0.0889 m, Mₛ = 0.5 g, ωₛ = 2π×6000 rad/s
        Tₛ = 366→288 K, Bₛ = 10⁻³ T, t_total = 1.17 s
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r = dataset.get('r_reactor', ORB_ANALYSIS_11_PARAMS['r_reactor'])
        M_s = dataset.get('M_s', ORB_ANALYSIS_11_PARAMS['M_s'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_11_PARAMS['omega_s'])
        T_base = dataset.get('T_base', ORB_ANALYSIS_11_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_11_PARAMS['T_top'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_11_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_11_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_11_PARAMS['UA'])
        t = dataset.get('t', 0.585)  # midpoint of 1.17 s sequence
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_11_PARAMS['gamma_decay'])
        E_react = dataset.get('E_react', ORB_ANALYSIS_11_PARAMS['E_react'])
        
        t_n = t  # normalized time
        
        # Ug_1: Internal dipole
        k_1 = 1.5e-4
        Ug1 = k_1 * (M_s / r) * math.exp(-gamma * t) * math.cos(math.pi * t_n) * (1 + 0.01 * math.sin(gamma * t))
        
        # Ug_2: Outer field bubble
        k_2 = 1.2
        Ug2 = k_2 * (UA + UA) * M_s / (r**2) * SCm * math.exp(-gamma * t)
        
        # Ug_3: Magnetic strings
        k_3 = 1.8
        Ug3 = k_3 * B_s * math.cos(omega_s * t * math.pi) * SCm * math.exp(-gamma * t)
        
        # Ub_i: Universal buoyancy
        beta_i = 0.8
        Omega_g = 7.3e-16  # rad/s
        M_bh = 8.15e36     # kg (reference black hole mass)
        d_g = 2.55e20      # m (galactic distance)
        Ubi = -beta_i * (Ug1 + Ug2 + Ug3) * Omega_g * (M_bh / d_g) * E_react * math.cos(math.pi * t_n)
        
        # Um: Universal magnetism
        mu_j = 1e-4
        Um = (mu_j / r) * (1 - math.exp(-gamma * t * math.cos(math.pi * t_n))) * SCm * math.exp(-gamma * t)
        
        # A_μν: Cosmic Aether tensor contribution
        eta = 1e-22
        rho_A = 1e-23
        A_munu = 1.0 + eta * (UA * SCm * rho_A) * t_n
        
        # Total F_U
        F_U = Ug1 + Ug2 + Ug3 + Ubi + Um + A_munu
        
        return {
            'Ug1': Ug1,
            'Ug2': Ug2,
            'Ug3': Ug3,
            'Ubi': Ubi,
            'Um': Um,
            'A_munu': A_munu,
            'F_U_total': F_U,
            'parameters': {
                'r': r, 'M_s': M_s, 'omega_s': omega_s,
                'T_gradient': f'{T_base}→{T_top} K',
                'B_s': B_s, 'SCm': SCm, 'UA': UA, 't': t
            },
            'equation': 'F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] + (gμν + η·Tₛμν)',
            'source': 'Grok UFT Orb Analysis_11 Refined F_U (March 4, 2025)'
        }


class ExtendedCyclePatternAnalyzerCalculator:
    """
    UFT Orb Analysis_11: Extended cycle pattern analyzer calculator.
    
    Analyzes the extended cycle patterns across multiple Orb analysis sessions,
    identifying recurring convection modes and transition frequencies.
    Combines data from Orb_8 through Orb_11.
    
    Physics:
        Cycle_frequency = n_complete_cycles / t_total
        Dominant_mode = mode(transition_types)
        Pattern_correlation = Σᵢ[θᵢ·θᵢ₊₁] / n
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Extended quadrant history from Orb_8 through Orb_11
        full_sequence = dataset.get('full_sequence', [
            # Orb_8: Photos #30-#31
            {'frame': 30, 'quadrant': 'upper_right', 'orb': 8},
            {'frame': 31, 'quadrant': 'lower_left', 'orb': 8},
            # Orb_9: Photos #32-#34
            {'frame': 32, 'quadrant': 'upper_left', 'orb': 9},
            {'frame': 33, 'quadrant': 'upper_right', 'orb': 9},
            {'frame': 34, 'quadrant': 'upper_right', 'orb': 9},
            # Orb_10: Photos #35-#37
            {'frame': 35, 'quadrant': 'lower_right', 'orb': 10},
            {'frame': 36, 'quadrant': 'upper_left', 'orb': 10},
            {'frame': 37, 'quadrant': 'upper_right', 'orb': 10},
            # Orb_11: Photos #38-#40
            {'frame': 38, 'quadrant': 'upper_left', 'orb': 11},
            {'frame': 39, 'quadrant': 'upper_right', 'orb': 11},
            {'frame': 40, 'quadrant': 'lower_left', 'orb': 11},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_11_PARAMS['dt_frame'])
        
        # Count transition types
        transition_types = {}
        quadrant_visits = {}
        
        for i in range(1, len(full_sequence)):
            prev = full_sequence[i-1]['quadrant']
            curr = full_sequence[i]['quadrant']
            
            # Count transitions
            trans = f"{prev} → {curr}"
            transition_types[trans] = transition_types.get(trans, 0) + 1
            
            # Count visits
            quadrant_visits[curr] = quadrant_visits.get(curr, 0) + 1
        
        # First quadrant
        quadrant_visits[full_sequence[0]['quadrant']] = quadrant_visits.get(
            full_sequence[0]['quadrant'], 0) + 1
        
        # Find dominant patterns
        dominant_transition = max(transition_types.items(), key=lambda x: x[1]) if transition_types else ('none', 0)
        most_visited = max(quadrant_visits.items(), key=lambda x: x[1]) if quadrant_visits else ('none', 0)
        
        # Total frames and time
        n_frames = len(full_sequence)
        t_total = n_frames * dt
        
        # Transition frequency
        n_transitions = len(full_sequence) - 1
        trans_frequency = n_transitions / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames_analyzed': n_frames,
            't_total': round(t_total, 4),
            'n_transitions': n_transitions,
            'transition_frequency': round(trans_frequency, 2),
            'transition_types': transition_types,
            'dominant_transition': dominant_transition[0],
            'dominant_count': dominant_transition[1],
            'quadrant_visits': quadrant_visits,
            'most_visited_quadrant': most_visited[0],
            'orb_sessions_covered': [8, 9, 10, 11],
            'equation': 'f_trans = n_transitions / t_total',
            'source': 'Grok UFT Orb Analysis_11 Extended Pattern (March 4, 2025)'
        }


class IntelligentPlasmoidBehaviorCalculator:
    """
    UFT Orb Analysis_11: Intelligent plasmoid behavior quantifier calculator.
    
    Quantifies the "intelligent" behavior traits observed in plasmoids:
    - Independence in pass-throughs (non-interfering)
    - Rare rotational energy transfers
    - Multi-axial rotations and spin drift
    - Shape-shifting (elongation, contraction)
    
    Physics:
        I_independence = 1 - P(collision)
        I_transfer = n_transfers / n_encounters
        I_intelligence = f(I_independence, I_transfer, complexity)
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Observed behavior parameters
        n_plasmoids = dataset.get('n_plasmoids', 15)  # average count
        n_encounters = dataset.get('n_encounters', 50)  # total close encounters
        n_collisions = dataset.get('n_collisions', 2)  # actual collisions
        n_energy_transfers = dataset.get('n_energy_transfers', 3)  # rotation transfers
        
        # Shape-shifting observations
        n_elongations = dataset.get('n_elongations', 20)
        n_contractions = dataset.get('n_contractions', 15)
        n_shapeshifts = n_elongations + n_contractions
        
        # Multi-axial rotation observations
        n_rotation_modes = dataset.get('n_rotation_modes', 3)  # axes observed
        n_reversals = dataset.get('n_reversals', 8)  # instant reversals
        
        # Calculate independence index
        I_independence = 1 - (n_collisions / n_encounters) if n_encounters > 0 else 1.0
        
        # Calculate transfer rarity index
        I_transfer_rarity = 1 - (n_energy_transfers / n_encounters) if n_encounters > 0 else 1.0
        
        # Calculate complexity index (normalized)
        complexity_raw = n_shapeshifts + n_reversals + n_rotation_modes * 5
        complexity_max = 100
        I_complexity = min(complexity_raw / complexity_max, 1.0)
        
        # Aggregate intelligence index
        I_intelligence = (I_independence * 0.4 + I_transfer_rarity * 0.3 + I_complexity * 0.3)
        
        # Behavior classification
        if I_intelligence > 0.9:
            classification = 'highly_intelligent'
        elif I_intelligence > 0.7:
            classification = 'intelligent'
        elif I_intelligence > 0.5:
            classification = 'semi_intelligent'
        else:
            classification = 'reactive'
        
        return {
            'n_plasmoids': n_plasmoids,
            'n_encounters': n_encounters,
            'n_collisions': n_collisions,
            'I_independence': round(I_independence, 4),
            'n_energy_transfers': n_energy_transfers,
            'I_transfer_rarity': round(I_transfer_rarity, 4),
            'n_shapeshifts': n_shapeshifts,
            'n_reversals': n_reversals,
            'n_rotation_modes': n_rotation_modes,
            'I_complexity': round(I_complexity, 4),
            'I_intelligence': round(I_intelligence, 4),
            'classification': classification,
            'equation': 'I = 0.4×I_indep + 0.3×I_transfer + 0.3×I_complex',
            'source': 'Grok UFT Orb Analysis_11 Intelligent Behavior (March 4, 2025)'
        }


class BulbDrivenPlasmaEnergeticsCalculator:
    """
    UFT Orb Analysis_11: 65W bulb-driven plasma energetics calculator.
    
    Models the energy input from the 65W incandescent bulb at the reactor bottom
    and its influence on plasmoid dynamics. Calculates thermal/electromagnetic
    coupling and energy transfer efficiency.
    
    Physics:
        P_input = 65 W (bulb power)
        E_thermal = P × t × η_thermal
        E_EM = P × t × η_EM × (λ_IR / λ_total)
        η_plasma = E_plasma / E_input
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        P_bulb = dataset.get('P_bulb', 65)  # W
        t_total = dataset.get('t_total', ORB_ANALYSIS_11_PARAMS['total_time_orb11'])
        eta_thermal = dataset.get('eta_thermal', 0.90)  # 90% to heat
        eta_EM = dataset.get('eta_EM', 0.10)  # 10% to light
        eta_IR = dataset.get('eta_IR', 0.60)  # 60% of light is IR
        
        # Total energy input
        E_input = P_bulb * t_total
        
        # Thermal energy (heats oil)
        E_thermal = E_input * eta_thermal
        
        # Electromagnetic energy (light output)
        E_light = E_input * eta_EM
        
        # IR component (most relevant for plasma)
        E_IR = E_light * eta_IR
        
        # Observed plasma energy
        E_plasma_observed = ORB_ANALYSIS_11_PARAMS['E_total_orb11']
        
        # Overall plasma coupling efficiency
        eta_plasma = E_plasma_observed / E_input if E_input > 0 else 0.0
        
        # Power density at plasma region
        V_plasma = math.pi * (ORB_ANALYSIS_11_PARAMS['r_reactor'])**2 * 0.05  # 5 cm plasma height
        P_density = P_bulb / V_plasma if V_plasma > 0 else 0.0
        
        return {
            'P_bulb': P_bulb,
            't_total': round(t_total, 4),
            'E_input_total': round(E_input, 4),
            'E_thermal': round(E_thermal, 4),
            'E_light': round(E_light, 4),
            'E_IR': round(E_IR, 4),
            'E_plasma_observed': E_plasma_observed,
            'eta_plasma_coupling': round(eta_plasma, 6),
            'P_density_plasma': round(P_density, 2),
            'efficiency_breakdown': {
                'thermal': eta_thermal,
                'EM': eta_EM,
                'IR_fraction': eta_IR
            },
            'equation': 'η_plasma = E_plasma / (P × t)',
            'source': 'Grok UFT Orb Analysis_11 Bulb Energetics (March 4, 2025)'
        }


class WaxCapCoolingDynamicsCalculator:
    """
    UFT Orb Analysis_11: Paraffin wax cap cooling dynamics calculator.
    
    Models the cooling effect of the 3-inch paraffin wax cap emulsified with
    red mercury. The cap maintains viscosity at 180-215°F and provides the
    thermal gradient that drives convection.
    
    Physics:
        Q_cooling = m_wax × c_wax × ΔT
        R_thermal = L / (k × A)
        T(z) = T_base - (T_base - T_cap) × (z/L)
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Wax cap properties
        h_cap = dataset.get('h_cap', 0.0762)  # m (3 in)
        r_reactor = dataset.get('r_reactor', ORB_ANALYSIS_11_PARAMS['r_reactor'])
        T_wax = dataset.get('T_wax', 370)  # K (~180-215°F average)
        
        # Thermal properties of paraffin wax
        rho_wax = dataset.get('rho_wax', 900)  # kg/m³
        c_wax = dataset.get('c_wax', 2900)  # J/(kg·K)
        k_wax = dataset.get('k_wax', 0.25)  # W/(m·K)
        
        # Reactor temperatures
        T_base = ORB_ANALYSIS_11_PARAMS['T_base']
        T_top = ORB_ANALYSIS_11_PARAMS['T_top']
        
        # Cap volume and mass
        V_cap = math.pi * r_reactor**2 * h_cap
        m_cap = rho_wax * V_cap
        
        # Cross-sectional area
        A_cross = math.pi * r_reactor**2
        
        # Thermal resistance of cap
        R_thermal = h_cap / (k_wax * A_cross) if (k_wax * A_cross) > 0 else float('inf')
        
        # Heat flow through cap
        Delta_T = T_wax - T_top
        Q_dot = Delta_T / R_thermal if R_thermal > 0 else 0.0
        
        # Cooling capacity (energy stored in cap)
        Q_capacity = m_cap * c_wax * (T_wax - T_top)
        
        # Temperature profile (linear approximation)
        L_reactor = 0.254  # m (10 in)
        z_positions = [0.0, 0.127, 0.254]
        T_profile = [T_base - (T_base - T_top) * z / L_reactor for z in z_positions]
        
        return {
            'h_cap': h_cap,
            'm_cap': round(m_cap, 6),
            'V_cap': round(V_cap, 8),
            'T_wax': T_wax,
            'R_thermal': round(R_thermal, 4),
            'Q_dot_heat_flow': round(Q_dot, 4),
            'Q_capacity': round(Q_capacity, 2),
            'Delta_T_across_cap': Delta_T,
            'T_profile': [{'z': z, 'T': round(T, 1)} for z, T in zip(z_positions, T_profile)],
            'equation': 'Q̇ = ΔT / R_thermal, R = L/(k×A)',
            'source': 'Grok UFT Orb Analysis_11 Wax Cap Cooling (March 4, 2025)'
        }


class FieldGeneratorResonanceCouplingCalculator:
    """
    UFT Orb Analysis_11: Field generator resonance coupling calculator.
    
    Models the 6000 Hz field resonance and its coupling to plasmoid dynamics.
    The resonance drives the "flapping" behavior observed in the Field Generator
    Experiment and influences [Ug3] and [Um] components.
    
    Physics:
        f_resonance = 6000 Hz
        ω = 2π × f
        E_resonance = ½ × L_eff × I² × sin²(ωt)
        Coupling = k_coupling × [Ug3 + Um] × cos(ωt)
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        f_resonance = dataset.get('f_resonance', 6000)  # Hz
        omega = 2 * math.pi * f_resonance
        t = dataset.get('t', 0.5)  # sample time
        
        # Effective inductance (from hydrogen bubble array)
        L_eff = dataset.get('L_eff', 1e-6)  # H (approximate)
        I_peak = dataset.get('I_peak', 0.1)  # A (approximate current)
        
        # Resonance energy
        E_resonance = 0.5 * L_eff * I_peak**2 * (math.sin(omega * t))**2
        
        # Field coupling coefficients
        k_coupling = 1e-10
        B_s = ORB_ANALYSIS_11_PARAMS['B_s']
        SCm = ORB_ANALYSIS_11_PARAMS['SCm']
        gamma = ORB_ANALYSIS_11_PARAMS['gamma_decay']
        
        # Ug3 + Um contribution
        k_3 = 1.8
        Ug3 = k_3 * B_s * math.cos(omega * t * math.pi) * SCm * math.exp(-gamma * t)
        
        mu_j = 1e-4
        r = ORB_ANALYSIS_11_PARAMS['r_reactor']
        Um = (mu_j / r) * (1 - math.exp(-gamma * t)) * SCm * math.exp(-gamma * t)
        
        # Coupling strength
        Coupling = k_coupling * (Ug3 + Um) * math.cos(omega * t)
        
        # Flapping amplitude estimate
        A_flap = abs(Coupling) * 1e3  # normalized
        
        return {
            'f_resonance': f_resonance,
            'omega': round(omega, 2),
            'E_resonance': E_resonance,
            'L_eff': L_eff,
            'I_peak': I_peak,
            'Ug3': Ug3,
            'Um': Um,
            'Coupling_strength': Coupling,
            'A_flap_normalized': round(A_flap, 6),
            'equation': 'Coupling = k×[Ug3+Um]×cos(ωt)',
            'source': 'Grok UFT Orb Analysis_11 Resonance Coupling (March 4, 2025)'
        }


class TotalEnergyBudgetCalculator:
    """
    UFT Orb Analysis_11: Total energy budget calculator.
    
    Computes the complete energy budget across all 39 frames, accounting for
    bulb input, thermal losses, plasma energy, and [ACE/DCE] events.
    
    Physics:
        E_in = P_bulb × t
        E_out = E_thermal_loss + E_radiation + E_plasma
        E_balance = E_in - E_out
        η_overall = E_plasma / E_in
        
    Source: Grok UFT Orb Analysis_11 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        P_bulb = dataset.get('P_bulb', 65)  # W
        t_total = dataset.get('t_total', ORB_ANALYSIS_11_PARAMS['total_time_orb11'])
        
        # Energy input
        E_in = P_bulb * t_total
        
        # Thermal losses (estimated)
        eta_thermal_loss = dataset.get('eta_thermal_loss', 0.3)  # 30% lost to environment
        E_thermal_loss = E_in * eta_thermal_loss
        
        # Radiation losses (visible + IR leaving system)
        eta_radiation = dataset.get('eta_radiation', 0.05)  # 5% radiated
        E_radiation = E_in * eta_radiation
        
        # Plasma energy (observed)
        E_plasma = ORB_ANALYSIS_11_PARAMS['E_total_orb11']
        
        # Energy retained in oil/wax
        E_stored = E_in - E_thermal_loss - E_radiation - E_plasma
        
        # Overall efficiency
        eta_overall = E_plasma / E_in if E_in > 0 else 0.0
        
        return {
            'P_bulb': P_bulb,
            't_total': round(t_total, 4),
            'E_input': round(E_in, 4),
            'E_thermal_loss': round(E_thermal_loss, 4),
            'E_radiation': round(E_radiation, 4),
            'E_plasma': E_plasma,
            'E_stored': round(E_stored, 4),
            'eta_overall': round(eta_overall, 6),
            'energy_breakdown': {
                'input': E_in,
                'thermal_loss': E_thermal_loss,
                'radiation': E_radiation,
                'plasma': E_plasma,
                'stored': E_stored
            },
            'equation': 'E_balance = E_in - E_loss - E_plasma',
            'source': 'Grok UFT Orb Analysis_11 Energy Budget (March 4, 2025)'
        }


# UFT Orb Analysis_11 registry dict
ORB_ANALYSIS_11_CALCULATORS = {
    'ThirtyNineFrameSequenceCalculator': ThirtyNineFrameSequenceCalculator(),
    'CounterClockwiseDiagonalCycleCalculator': CounterClockwiseDiagonalCycleCalculator(),
    'Orb11RefinedFUCalculator': Orb11RefinedFUCalculator(),
    'ExtendedCyclePatternAnalyzerCalculator': ExtendedCyclePatternAnalyzerCalculator(),
    'IntelligentPlasmoidBehaviorCalculator': IntelligentPlasmoidBehaviorCalculator(),
    'BulbDrivenPlasmaEnergeticsCalculator': BulbDrivenPlasmaEnergeticsCalculator(),
    'WaxCapCoolingDynamicsCalculator': WaxCapCoolingDynamicsCalculator(),
    'FieldGeneratorResonanceCouplingCalculator': FieldGeneratorResonanceCouplingCalculator(),
    'TotalEnergyBudgetCalculator': TotalEnergyBudgetCalculator(),
}


# ============================================================================
# UFT ORB ANALYSIS_12 CALCULATORS (7 Calculator Classes)
# Source: Grok UFT Orb Analysis_12 (March 4, 2025)
# Photos: #41-#43 of 496 Red Dwarf Reactor Plasma Orb infrared images
# Physics: 42 frames (~1.26 s), cyclical convection pattern:
#          lower right (#41) → upper left (#42) → upper right (#43)
#          ~1.89 J total energy, hydrogen bubble path anchoring
# ============================================================================

# UFT Orb Analysis_12 Parameters
ORB_ANALYSIS_12_PARAMS = {
    'r_reactor': 0.0889,           # m (3.5 in diameter / 2)
    'M_s': 0.5e-3,                 # kg (plasma orb mass 0.5 g)
    'omega_s': 2 * math.pi * 6000, # rad/s (6000 Hz field resonance)
    'T_base': 366,                 # K (bulb base temperature)
    'T_top': 288,                  # K (ambient top temperature)
    'B_s': 1e-3,                   # T (magnetic field from H2 bubbles)
    'SCm': 1e15,                   # kg/m³ (hypothetical density at base)
    'UA': 1e-11,                   # C (trapped Aether charge)
    'dt_frame': 0.03,              # s (frame interval)
    'n_frames_orb12': 42,          # total frames (#1-#43 minus #14)
    'total_time_orb12': 1.26,      # s (42 frames × 0.03 s)
    'E_total_orb12': 1.89,         # J (total energy ~45 mJ/frame × 42)
    'E_react': 1e15,               # W/m³ (reactivity with thermal decay)
    'gamma_decay': 0.001,          # decay constant
    'v_plasmoid': 0.5,             # m/s (plasmoid motion speed)
    'n_H2_bubbles': 15,            # average hydrogen bubbles (12-18)
}


class FortyTwoFrameSequenceCalculator:
    """
    UFT Orb Analysis_12: 42-frame sequence energy calculator.
    
    Extends the video analysis to 42 frames (Photos #1-#43 minus #14),
    covering ~1.26 s at 33.3 fps. Calculates cumulative energy and temporal
    evolution of plasmoid dynamics through the extended sequence.
    
    Physics:
        E_total = Σᵢ[E_frame_i] = 42 × 45 mJ ≈ 1.89 J
        t_total = n_frames × dt = 42 × 0.03 s = 1.26 s
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_12_PARAMS['n_frames_orb12'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_12_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.045)  # J (~45 mJ/frame)
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_12_PARAMS['gamma_decay'])
        
        t_total = n_frames * dt
        
        E_cumulative = 0.0
        frame_energies = []
        for i in range(n_frames):
            t_n = i * dt
            E_frame = E_per_frame * math.exp(-gamma * t_n)
            E_cumulative += E_frame
            frame_energies.append({
                'frame': i + 1,
                't': round(t_n, 4),
                'E': round(E_frame, 6)
            })
        
        P_avg = E_cumulative / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_cumulative': round(E_cumulative, 6),
            'P_avg': round(P_avg, 6),
            'E_per_frame_avg': round(E_cumulative / n_frames, 6) if n_frames > 0 else 0.0,
            'frame_energies': frame_energies[:5],
            'equation': 'E_total = Σᵢ[E_frame_i × e^(-γtᵢ)]',
            'source': 'Grok UFT Orb Analysis_12 42-Frame Sequence (March 4, 2025)'
        }


class CyclicalConvectionOrb12Calculator:
    """
    UFT Orb Analysis_12: Cyclical convection pattern calculator for Photos #41-#43.
    
    Tracks the cyclical concentration shifts: lower right (#41) → upper left (#42)
    → upper right (#43). This pattern reflects ongoing convection driven by
    thermal gradients and [Ug3]/[Um] fields.
    
    Physics:
        Cycle_path: LR → UL → UR (3-point convection loop)
        θ_rotation = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]
        Cycle_period ≈ 3 frames × 0.03 s = 0.09 s
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        cycle_sequence = dataset.get('cycle_sequence', [
            {'frame': 41, 'quadrant': 'lower_right', 'x_c': 0.75, 'y_c': 0.25},
            {'frame': 42, 'quadrant': 'upper_left', 'x_c': 0.25, 'y_c': 0.75},
            {'frame': 43, 'quadrant': 'upper_right', 'x_c': 0.75, 'y_c': 0.75},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_12_PARAMS['dt_frame'])
        
        transitions = []
        total_angle = 0.0
        total_distance = 0.0
        
        for i in range(1, len(cycle_sequence)):
            prev = cycle_sequence[i - 1]
            curr = cycle_sequence[i]
            
            dx = curr['x_c'] - prev['x_c']
            dy = curr['y_c'] - prev['y_c']
            distance = math.sqrt(dx**2 + dy**2)
            theta = math.atan2(dy, dx) * 180 / math.pi
            
            total_angle += theta
            total_distance += distance
            
            transitions.append({
                'from_frame': prev['frame'],
                'to_frame': curr['frame'],
                'transition': f"{prev['quadrant']} → {curr['quadrant']}",
                'dx': round(dx, 4),
                'dy': round(dy, 4),
                'distance': round(distance, 4),
                'theta_deg': round(theta, 2)
            })
        
        n_transitions = len(transitions)
        cycle_period = n_transitions * dt
        avg_angular_velocity = total_angle / cycle_period if cycle_period > 0 else 0.0
        
        return {
            'cycle_sequence': cycle_sequence,
            'transitions': transitions,
            'total_distance': round(total_distance, 4),
            'total_angle_deg': round(total_angle, 2),
            'cycle_period': round(cycle_period, 4),
            'avg_angular_velocity': round(avg_angular_velocity, 2),
            'pattern': 'LR → UL → UR (cyclical convection)',
            'equation': 'θ_cycle = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]',
            'source': 'Grok UFT Orb Analysis_12 Cyclical Convection (March 4, 2025)'
        }


class Orb12RefinedFUCalculator:
    """
    UFT Orb Analysis_12: Refined Unified Field F_U calculator.
    
    Computes the complete unified field equation refined with Photos #41-#43 data,
    extending the temporal coverage to 1.26 s with 42 frames.
    
    Physics:
        F_U = Σᵢ[kᵢ·Ugᵢ(r,t,Mₛ,ωₛ,Tₛ,Bₛ,SCm,UA,tₙ) - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react]
              + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ]
              + (gμν + η·Tₛμν(UA,SCm,ρA))
              
    Parameters updated for 42-frame sequence:
        t_total = 1.26 s, E_total = 1.89 J
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r = dataset.get('r_reactor', ORB_ANALYSIS_12_PARAMS['r_reactor'])
        M_s = dataset.get('M_s', ORB_ANALYSIS_12_PARAMS['M_s'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_12_PARAMS['omega_s'])
        T_base = dataset.get('T_base', ORB_ANALYSIS_12_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_12_PARAMS['T_top'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_12_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_12_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_12_PARAMS['UA'])
        t = dataset.get('t', 0.63)  # midpoint of 1.26 s sequence
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_12_PARAMS['gamma_decay'])
        E_react = dataset.get('E_react', ORB_ANALYSIS_12_PARAMS['E_react'])
        
        t_n = t
        
        # Ug_1: Internal dipole
        k_1 = 1.5e-4
        Ug1 = k_1 * (M_s / r) * math.exp(-gamma * t) * math.cos(math.pi * t_n) * (1 + 0.01 * math.sin(gamma * t))
        
        # Ug_2: Outer field bubble
        k_2 = 1.2
        Ug2 = k_2 * (UA + UA) * M_s / (r**2) * SCm * math.exp(-gamma * t)
        
        # Ug_3: Magnetic strings
        k_3 = 1.8
        Ug3 = k_3 * B_s * math.cos(omega_s * t * math.pi) * SCm * math.exp(-gamma * t)
        
        # Ub_i: Universal buoyancy
        beta_i = 0.8
        Omega_g = 7.3e-16
        M_bh = 8.15e36
        d_g = 2.55e20
        Ubi = -beta_i * (Ug1 + Ug2 + Ug3) * Omega_g * (M_bh / d_g) * E_react * math.cos(math.pi * t_n)
        
        # Um: Universal magnetism
        mu_j = 1e-4
        Um = (mu_j / r) * (1 - math.exp(-gamma * t * math.cos(math.pi * t_n))) * SCm * math.exp(-gamma * t)
        
        # A_μν: Cosmic Aether tensor
        eta = 1e-22
        rho_A = 1e-23
        A_munu = 1.0 + eta * (UA * SCm * rho_A) * t_n
        
        F_U = Ug1 + Ug2 + Ug3 + Ubi + Um + A_munu
        
        return {
            'Ug1': Ug1,
            'Ug2': Ug2,
            'Ug3': Ug3,
            'Ubi': Ubi,
            'Um': Um,
            'A_munu': A_munu,
            'F_U_total': F_U,
            'parameters': {
                'r': r, 'M_s': M_s, 'omega_s': omega_s,
                'T_gradient': f'{T_base}→{T_top} K',
                'B_s': B_s, 'SCm': SCm, 'UA': UA, 't': t
            },
            'equation': 'F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] + (gμν + η·Tₛμν)',
            'source': 'Grok UFT Orb Analysis_12 Refined F_U (March 4, 2025)'
        }


class FullSequenceQuadrantTrackerCalculator:
    """
    UFT Orb Analysis_12: Full sequence quadrant tracker calculator.
    
    Tracks quadrant concentration shifts across the entire sequence
    (Photos #1-#43, minus #14), identifying dominant patterns and
    transition frequencies.
    
    Physics:
        Transition_matrix[i,j] = count(Q_i → Q_j)
        Dominant_mode = argmax(quadrant_visits)
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Extended quadrant history from Orb_8 through Orb_12
        full_sequence = dataset.get('full_sequence', [
            # Orb_8: Photos #30-#31
            {'frame': 30, 'quadrant': 'upper_right', 'orb': 8},
            {'frame': 31, 'quadrant': 'lower_left', 'orb': 8},
            # Orb_9: Photos #32-#34
            {'frame': 32, 'quadrant': 'upper_left', 'orb': 9},
            {'frame': 33, 'quadrant': 'upper_right', 'orb': 9},
            {'frame': 34, 'quadrant': 'upper_right', 'orb': 9},
            # Orb_10: Photos #35-#37
            {'frame': 35, 'quadrant': 'lower_right', 'orb': 10},
            {'frame': 36, 'quadrant': 'upper_left', 'orb': 10},
            {'frame': 37, 'quadrant': 'upper_right', 'orb': 10},
            # Orb_11: Photos #38-#40
            {'frame': 38, 'quadrant': 'upper_left', 'orb': 11},
            {'frame': 39, 'quadrant': 'upper_right', 'orb': 11},
            {'frame': 40, 'quadrant': 'lower_left', 'orb': 11},
            # Orb_12: Photos #41-#43
            {'frame': 41, 'quadrant': 'lower_right', 'orb': 12},
            {'frame': 42, 'quadrant': 'upper_left', 'orb': 12},
            {'frame': 43, 'quadrant': 'upper_right', 'orb': 12},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_12_PARAMS['dt_frame'])
        
        # Count transitions and visits
        transition_types = {}
        quadrant_visits = {}
        
        for i in range(1, len(full_sequence)):
            prev = full_sequence[i-1]['quadrant']
            curr = full_sequence[i]['quadrant']
            
            trans = f"{prev} → {curr}"
            transition_types[trans] = transition_types.get(trans, 0) + 1
            quadrant_visits[curr] = quadrant_visits.get(curr, 0) + 1
        
        quadrant_visits[full_sequence[0]['quadrant']] = quadrant_visits.get(
            full_sequence[0]['quadrant'], 0) + 1
        
        dominant_transition = max(transition_types.items(), key=lambda x: x[1]) if transition_types else ('none', 0)
        most_visited = max(quadrant_visits.items(), key=lambda x: x[1]) if quadrant_visits else ('none', 0)
        
        n_frames = len(full_sequence)
        t_total = n_frames * dt
        n_transitions = n_frames - 1
        trans_frequency = n_transitions / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames_analyzed': n_frames,
            't_total': round(t_total, 4),
            'n_transitions': n_transitions,
            'transition_frequency': round(trans_frequency, 2),
            'transition_types': transition_types,
            'dominant_transition': dominant_transition[0],
            'quadrant_visits': quadrant_visits,
            'most_visited_quadrant': most_visited[0],
            'orb_sessions_covered': [8, 9, 10, 11, 12],
            'equation': 'f_trans = n_transitions / t_total',
            'source': 'Grok UFT Orb Analysis_12 Full Sequence Tracker (March 4, 2025)'
        }


class HydrogenBubblePathAnchorCalculator:
    """
    UFT Orb Analysis_12: Hydrogen bubble path anchoring calculator.
    
    Models how the 12-18 evenly-spaced hydrogen bubbles anchor plasmoid paths,
    enhancing [SCm] and [UA] reactivity via magnetic fields (~10⁻³ T).
    
    Physics:
        B_anchor = Σⱼ[B_bubble_j × exp(-r/λ_anchor)]
        Path_stability = τ_anchor / τ_diffusion
        Field_uniformity = 1 - σ(B)/⟨B⟩
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_bubbles = dataset.get('n_H2_bubbles', ORB_ANALYSIS_12_PARAMS['n_H2_bubbles'])
        B_per_bubble = dataset.get('B_s', ORB_ANALYSIS_12_PARAMS['B_s'])
        r_reactor = dataset.get('r_reactor', ORB_ANALYSIS_12_PARAMS['r_reactor'])
        L_reactor = dataset.get('L_reactor', 0.254)  # m (10 in)
        
        # Bubble spacing (evenly distributed)
        bubble_spacing = L_reactor / n_bubbles
        
        # Characteristic anchoring length
        lambda_anchor = bubble_spacing / 2
        
        # Total anchor field at center
        B_anchor_total = 0.0
        bubble_positions = []
        for j in range(n_bubbles):
            z_bubble = bubble_spacing * (j + 0.5)
            z_center = L_reactor / 2
            r_from_center = abs(z_bubble - z_center)
            B_contribution = B_per_bubble * math.exp(-r_from_center / lambda_anchor)
            B_anchor_total += B_contribution
            bubble_positions.append({
                'bubble': j + 1,
                'z_position': round(z_bubble, 4),
                'B_contribution': round(B_contribution, 8)
            })
        
        # Diffusion timescale
        D_plasmoid = 1e-9  # m²/s
        tau_diffusion = r_reactor**2 / (4 * D_plasmoid)
        
        # Anchor timescale (magnetic confinement)
        mu_0 = 4 * math.pi * 1e-7
        m_plasmoid = ORB_ANALYSIS_12_PARAMS['M_s']
        tau_anchor = m_plasmoid * r_reactor / (B_anchor_total**2 / mu_0) if B_anchor_total > 0 else 0.0
        
        # Path stability ratio
        path_stability = tau_anchor / tau_diffusion if tau_diffusion > 0 else 0.0
        
        # Field uniformity estimate (1.0 = perfectly uniform)
        B_values = [bp['B_contribution'] for bp in bubble_positions]
        B_mean = np.mean(B_values)
        B_std = np.std(B_values)
        field_uniformity = 1 - (B_std / B_mean) if B_mean > 0 else 0.0
        
        return {
            'n_bubbles': n_bubbles,
            'bubble_spacing': round(bubble_spacing, 4),
            'lambda_anchor': round(lambda_anchor, 4),
            'B_anchor_total': round(B_anchor_total, 8),
            'tau_diffusion': round(tau_diffusion, 2),
            'tau_anchor': tau_anchor,
            'path_stability': path_stability,
            'field_uniformity': round(field_uniformity, 4),
            'bubble_positions': bubble_positions[:5],
            'equation': 'B_anchor = Σⱼ[Bⱼ × e^(-r/λ)]',
            'source': 'Grok UFT Orb Analysis_12 H2 Bubble Path Anchor (March 4, 2025)'
        }


class ThermalConvectionCycleCalculator:
    """
    UFT Orb Analysis_12: Thermal convection cycle calculator.
    
    Models the thermal-driven cyclical convection patterns observed across
    the sequence, driven by the 366K→288K gradient and [Ub] effects.
    
    Physics:
        Ra = (g × β × ΔT × L³) / (ν × α) (Rayleigh number)
        n_cycles = t_total / τ_cycle
        ω_convection = 2π / τ_cycle
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        T_base = dataset.get('T_base', ORB_ANALYSIS_12_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_12_PARAMS['T_top'])
        L = dataset.get('L_reactor', 0.254)  # m
        t_total = dataset.get('t_total', ORB_ANALYSIS_12_PARAMS['total_time_orb12'])
        
        # Oil properties
        rho_oil = dataset.get('rho_oil', 850)  # kg/m³
        beta_thermal = dataset.get('beta_thermal', 7e-4)  # K⁻¹
        nu_oil = dataset.get('nu_oil', 1e-5)  # m²/s (kinematic viscosity)
        alpha_oil = dataset.get('alpha_oil', 1e-7)  # m²/s (thermal diffusivity)
        g = 9.81  # m/s²
        
        Delta_T = T_base - T_top
        
        # Rayleigh number
        Ra = (g * beta_thermal * Delta_T * L**3) / (nu_oil * alpha_oil)
        
        # Convection regime
        Ra_critical = 1708  # onset of convection
        if Ra > 1e6:
            regime = 'turbulent'
        elif Ra > Ra_critical:
            regime = 'laminar_convective'
        else:
            regime = 'conductive'
        
        # Convection velocity estimate
        v_convection = math.sqrt(g * beta_thermal * Delta_T * L) if Ra > Ra_critical else 0.0
        
        # Cycle period estimate (single convection roll)
        tau_cycle = 2 * L / v_convection if v_convection > 0 else float('inf')
        
        # Number of complete cycles
        n_cycles = t_total / tau_cycle if tau_cycle > 0 else 0.0
        
        # Angular velocity of convection
        omega_convection = 2 * math.pi / tau_cycle if tau_cycle > 0 else 0.0
        
        return {
            'Delta_T': Delta_T,
            'Rayleigh_number': Ra,
            'Ra_critical': Ra_critical,
            'regime': regime,
            'v_convection': round(v_convection, 4),
            'tau_cycle': round(tau_cycle, 4) if tau_cycle != float('inf') else 'inf',
            'n_cycles_in_sequence': round(n_cycles, 2) if n_cycles != 0.0 else 0.0,
            'omega_convection': round(omega_convection, 4) if omega_convection != 0.0 else 0.0,
            'equation': 'Ra = (gβΔTL³)/(να)',
            'source': 'Grok UFT Orb Analysis_12 Thermal Convection (March 4, 2025)'
        }


class CumulativeEnergyAnalyzerCalculator:
    """
    UFT Orb Analysis_12: Cumulative energy analyzer calculator.
    
    Tracks total energy accumulation across the full 42-frame sequence,
    comparing observed vs theoretical predictions from F_U.
    
    Physics:
        E_observed = n_frames × E_per_frame = 42 × 45 mJ ≈ 1.89 J
        E_theoretical = ∫F_U dt over 1.26 s
        η_model = E_observed / E_theoretical
        
    Source: Grok UFT Orb Analysis_12 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_12_PARAMS['n_frames_orb12'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_12_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.045)  # J
        P_bulb = dataset.get('P_bulb', 65)  # W
        
        t_total = n_frames * dt
        
        # Observed energy
        E_observed = n_frames * E_per_frame
        
        # Input energy from bulb
        E_input = P_bulb * t_total
        
        # Theoretical energy (simplified F_U integral)
        # Using average F_U contribution scaled to reactor
        k_FU = dataset.get('k_FU', 1e-5)  # coupling constant
        E_theoretical = k_FU * E_input
        
        # Model efficiency
        eta_model = E_observed / E_theoretical if E_theoretical > 0 else 0.0
        
        # Energy breakdown by Orb session
        energy_by_orb = {
            'Orb_10': {'frames': 36, 'E': 1.62},
            'Orb_11': {'frames': 39, 'E': 1.755},
            'Orb_12': {'frames': 42, 'E': 1.89},
        }
        
        # Delta from previous session
        delta_E = energy_by_orb['Orb_12']['E'] - energy_by_orb['Orb_11']['E']
        delta_frames = energy_by_orb['Orb_12']['frames'] - energy_by_orb['Orb_11']['frames']
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_observed': round(E_observed, 4),
            'E_input_bulb': round(E_input, 4),
            'E_theoretical': round(E_theoretical, 6),
            'eta_model': round(eta_model, 4),
            'energy_by_orb': energy_by_orb,
            'delta_E_from_Orb11': round(delta_E, 4),
            'delta_frames_from_Orb11': delta_frames,
            'E_per_frame_avg': round(E_observed / n_frames, 6) if n_frames > 0 else 0.0,
            'equation': 'E_obs = n × E_frame, η = E_obs / E_theory',
            'source': 'Grok UFT Orb Analysis_12 Energy Analyzer (March 4, 2025)'
        }


# UFT Orb Analysis_12 registry dict
ORB_ANALYSIS_12_CALCULATORS = {
    'FortyTwoFrameSequenceCalculator': FortyTwoFrameSequenceCalculator(),
    'CyclicalConvectionOrb12Calculator': CyclicalConvectionOrb12Calculator(),
    'Orb12RefinedFUCalculator': Orb12RefinedFUCalculator(),
    'FullSequenceQuadrantTrackerCalculator': FullSequenceQuadrantTrackerCalculator(),
    'HydrogenBubblePathAnchorCalculator': HydrogenBubblePathAnchorCalculator(),
    'ThermalConvectionCycleCalculator': ThermalConvectionCycleCalculator(),
    'CumulativeEnergyAnalyzerCalculator': CumulativeEnergyAnalyzerCalculator(),
}


# ============================================================================
# UFT ORB ANALYSIS_13 CALCULATORS (6 Calculator Classes)
# Source: Grok UFT Orb Analysis_13 (March 4, 2025)
# Photos: #44-#45 of 496 Red Dwarf Reactor Plasma Orb infrared images
# Physics: 44 frames (~1.32 s at 33.3 fps), ~2.02 J total energy
#          Diagonal shift: lower left (#44) → upper right (#45)
#          Intelligent quantum plasmoid dynamics, celestial-like behavior
# ============================================================================

# UFT Orb Analysis_13 Parameters
ORB_ANALYSIS_13_PARAMS = {
    'r_reactor': 0.0889,           # m (3.5 in diameter / 2)
    'M_s': 0.5e-3,                 # kg (plasma orb mass 0.5 g)
    'omega_s': 2 * math.pi * 6000, # rad/s (6000 Hz field resonance)
    'T_base': 366,                 # K (bulb base temperature)
    'T_top': 288,                  # K (ambient top temperature)
    'B_s': 1e-3,                   # T (magnetic field from H2 bubbles)
    'SCm': 1e15,                   # kg/m³ (hypothetical density at base)
    'UA': 1e-11,                   # C (trapped Aether charge)
    'dt_frame': 0.03,              # s (frame interval)
    'n_frames_orb13': 44,          # total frames (#1-#45 minus #14)
    'total_time_orb13': 1.32,      # s (44 frames × 0.03 s)
    'E_total_orb13': 2.02,         # J (total energy ~45 mJ/frame × 44)
    'E_react': 1e15,               # W/m³ (reactivity with thermal decay)
    'gamma_decay': 0.001,          # decay constant
    'v_plasmoid': 0.5,             # m/s (plasmoid motion speed)
    'n_H2_bubbles': 15,            # average hydrogen bubbles (12-18)
}


class FortyFourFrameSequenceCalculator:
    """
    UFT Orb Analysis_13: 44-frame sequence energy calculator.
    
    Extends the video analysis to 44 frames (Photos #1-#45 minus #14),
    covering ~1.32 s at 33.3 fps. Calculates cumulative energy and temporal
    evolution of plasmoid dynamics through the extended sequence.
    
    Physics:
        E_total = Σᵢ[E_frame_i] = 44 × 45 mJ ≈ 2.02 J
        t_total = n_frames × dt = 44 × 0.03 s = 1.32 s
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_13_PARAMS['n_frames_orb13'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_13_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.046)  # J (~46 mJ/frame for Orb_13)
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_13_PARAMS['gamma_decay'])
        
        t_total = n_frames * dt
        
        E_cumulative = 0.0
        frame_energies = []
        for i in range(n_frames):
            t_n = i * dt
            E_frame = E_per_frame * math.exp(-gamma * t_n)
            E_cumulative += E_frame
            frame_energies.append({
                'frame': i + 1,
                't': round(t_n, 4),
                'E': round(E_frame, 6)
            })
        
        P_avg = E_cumulative / t_total if t_total > 0 else 0.0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_cumulative': round(E_cumulative, 6),
            'P_avg': round(P_avg, 6),
            'E_per_frame_avg': round(E_cumulative / n_frames, 6) if n_frames > 0 else 0.0,
            'frame_energies': frame_energies[:5],
            'equation': 'E_total = Σᵢ[E_frame_i × e^(-γtᵢ)]',
            'source': 'Grok UFT Orb Analysis_13 44-Frame Sequence (March 4, 2025)'
        }


class DiagonalShiftOrb13Calculator:
    """
    UFT Orb Analysis_13: Diagonal shift pattern calculator for Photos #44-#45.
    
    Tracks the diagonal concentration shift: lower left (#44) → upper right (#45).
    This pattern reflects convection driven by thermal gradients and [Ug3]/[Um] fields.
    
    Physics:
        Diagonal_path: LL → UR (2-point diagonal shift)
        θ_diagonal = arctan2(Δy, Δx) ≈ 45° (ideal diagonal)
        Shift_distance = √(Δx² + Δy²)
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        shift_sequence = dataset.get('shift_sequence', [
            {'frame': 44, 'quadrant': 'lower_left', 'x_c': 0.25, 'y_c': 0.25},
            {'frame': 45, 'quadrant': 'upper_right', 'x_c': 0.75, 'y_c': 0.75},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_13_PARAMS['dt_frame'])
        
        if len(shift_sequence) < 2:
            return {'error': 'Need at least 2 frames for shift analysis'}
        
        prev = shift_sequence[0]
        curr = shift_sequence[1]
        
        dx = curr['x_c'] - prev['x_c']
        dy = curr['y_c'] - prev['y_c']
        distance = math.sqrt(dx**2 + dy**2)
        theta = math.atan2(dy, dx) * 180 / math.pi
        
        # Velocity estimate (normalized distance per frame time)
        t_shift = dt  # single frame transition
        v_shift = distance / t_shift if t_shift > 0 else 0.0
        
        # Check if it's a true diagonal (within 10° of 45° or 135°)
        is_diagonal = abs(abs(theta) - 45) < 10 or abs(abs(theta) - 135) < 10
        
        return {
            'shift_sequence': shift_sequence,
            'from_quadrant': prev['quadrant'],
            'to_quadrant': curr['quadrant'],
            'dx': round(dx, 4),
            'dy': round(dy, 4),
            'distance': round(distance, 4),
            'theta_deg': round(theta, 2),
            'is_diagonal': is_diagonal,
            'shift_velocity': round(v_shift, 4),
            'pattern': 'LL → UR (diagonal shift)',
            'equation': 'θ_shift = arctan2(Δy, Δx)',
            'source': 'Grok UFT Orb Analysis_13 Diagonal Shift (March 4, 2025)'
        }


class Orb13RefinedFUCalculator:
    """
    UFT Orb Analysis_13: Refined Unified Field F_U calculator.
    
    Computes the complete unified field equation refined with Photos #44-#45 data,
    extending the temporal coverage to 1.32 s with 44 frames.
    
    Physics:
        F_U = Σᵢ[kᵢ·Ugᵢ(r,t,Mₛ,ωₛ,Tₛ,Bₛ,SCm,UA,tₙ) - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react]
              + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ]
              + (gμν + η·Tₛμν(UA,SCm,ρA))
              
    Parameters updated for 44-frame sequence:
        t_total = 1.32 s, E_total = 2.02 J
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r = dataset.get('r_reactor', ORB_ANALYSIS_13_PARAMS['r_reactor'])
        M_s = dataset.get('M_s', ORB_ANALYSIS_13_PARAMS['M_s'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_13_PARAMS['omega_s'])
        T_base = dataset.get('T_base', ORB_ANALYSIS_13_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_13_PARAMS['T_top'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_13_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_13_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_13_PARAMS['UA'])
        t = dataset.get('t', 0.66)  # midpoint of 1.32 s sequence
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_13_PARAMS['gamma_decay'])
        E_react = dataset.get('E_react', ORB_ANALYSIS_13_PARAMS['E_react'])
        
        t_n = t
        
        # Ug_1: Internal dipole
        k_1 = 1.5e-4
        Ug1 = k_1 * (M_s / r) * math.exp(-gamma * t) * math.cos(math.pi * t_n) * (1 + 0.01 * math.sin(gamma * t))
        
        # Ug_2: Outer field bubble
        k_2 = 1.2
        Ug2 = k_2 * (UA + UA) * M_s / (r**2) * SCm * math.exp(-gamma * t)
        
        # Ug_3: Magnetic strings
        k_3 = 1.8
        Ug3 = k_3 * B_s * math.cos(omega_s * t * math.pi) * SCm * math.exp(-gamma * t)
        
        # Ub_i: Universal buoyancy
        beta_i = 0.8
        Omega_g = 7.3e-16
        M_bh = 8.15e36
        d_g = 2.55e20
        Ubi = -beta_i * (Ug1 + Ug2 + Ug3) * Omega_g * (M_bh / d_g) * E_react * math.cos(math.pi * t_n)
        
        # Um: Universal magnetism
        mu_j = 1e-4
        Um = (mu_j / r) * (1 - math.exp(-gamma * t * math.cos(math.pi * t_n))) * SCm * math.exp(-gamma * t)
        
        # A_μν: Cosmic Aether tensor
        eta = 1e-22
        rho_A = 1e-23
        A_munu = 1.0 + eta * (UA * SCm * rho_A) * t_n
        
        F_U = Ug1 + Ug2 + Ug3 + Ubi + Um + A_munu
        
        return {
            'Ug1': Ug1,
            'Ug2': Ug2,
            'Ug3': Ug3,
            'Ubi': Ubi,
            'Um': Um,
            'A_munu': A_munu,
            'F_U_total': F_U,
            'parameters': {
                'r': r, 'M_s': M_s, 'omega_s': omega_s,
                'T_gradient': f'{T_base}→{T_top} K',
                'B_s': B_s, 'SCm': SCm, 'UA': UA, 't': t
            },
            'equation': 'F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] + (gμν + η·Tₛμν)',
            'source': 'Grok UFT Orb Analysis_13 Refined F_U (March 4, 2025)'
        }


class IntelligentQuantumPlasmoidCalculator:
    """
    UFT Orb Analysis_13: Intelligent quantum plasmoid behavior calculator.
    
    Models the observed intelligent behavior of quantum plasmoids:
    - Unique spins and shape-shifting (elongation)
    - Multi-axial rotations and spin drift
    - Instant reverse rotations
    - Non-interfering pass-throughs
    - Rare rotational energy transfers
    - Celestial-like dynamics
    
    Physics:
        Intelligence_metric = (spin_uniqueness × shape_shift_factor × independence_score)
        ω_plasmoid = Σᵢ[ωᵢ × (1 - δᵢⱼ)]  (superposition of independent spins)
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Plasmoid characteristics
        n_plasmoids = dataset.get('n_plasmoids', 15)  # average observed
        spin_rate_range = dataset.get('spin_rate_range', (1, 10))  # rotations/sec
        shape_shift_factor = dataset.get('shape_shift_factor', 0.8)  # elongation degree
        independence_score = dataset.get('independence_score', 0.95)  # non-interference rate
        energy_transfer_rate = dataset.get('energy_transfer_rate', 0.05)  # rare transfers
        
        # Generate mock plasmoid ensemble
        plasmoids = []
        total_spin = 0.0
        for i in range(n_plasmoids):
            omega_i = spin_rate_range[0] + (spin_rate_range[1] - spin_rate_range[0]) * (i / n_plasmoids)
            # Multi-axial rotation: 3 axes
            axes = [
                {'axis': 'x', 'omega': omega_i * (1 + 0.1 * math.sin(i))},
                {'axis': 'y', 'omega': omega_i * (1 - 0.1 * math.cos(i))},
                {'axis': 'z', 'omega': omega_i * 0.5},
            ]
            total_spin += sum([ax['omega'] for ax in axes])
            plasmoids.append({
                'id': i + 1,
                'primary_omega': round(omega_i, 2),
                'shape': 'elongated' if i % 3 == 0 else 'spherical',
                'brightness': 'sun-like' if i % 4 == 0 else 'normal',
            })
        
        # Intelligence metric
        spin_uniqueness = 1 - (sum([p['primary_omega'] for p in plasmoids]) / (n_plasmoids * spin_rate_range[1]))
        intelligence_metric = spin_uniqueness * shape_shift_factor * independence_score
        
        # Celestial behavior indicators
        celestial_features = {
            'multi_axial_rotation': True,
            'spin_drift': True,
            'instant_reverse': True,
            'non_interfering_passthrough': True,
            'rare_energy_transfer': True,
        }
        celestial_score = sum(celestial_features.values()) / len(celestial_features)
        
        return {
            'n_plasmoids': n_plasmoids,
            'total_spin_sum': round(total_spin, 2),
            'spin_uniqueness': round(spin_uniqueness, 4),
            'shape_shift_factor': shape_shift_factor,
            'independence_score': independence_score,
            'energy_transfer_rate': energy_transfer_rate,
            'intelligence_metric': round(intelligence_metric, 4),
            'celestial_features': celestial_features,
            'celestial_score': celestial_score,
            'sample_plasmoids': plasmoids[:5],
            'equation': 'I_metric = spin_uniqueness × shape_shift × independence',
            'source': 'Grok UFT Orb Analysis_13 Intelligent Plasmoid (March 4, 2025)'
        }


class Orb13EnergyProgressionCalculator:
    """
    UFT Orb Analysis_13: Energy progression tracker across Orb sessions.
    
    Tracks the cumulative energy growth from Orb_10 through Orb_13,
    showing the linear progression as more frames are analyzed.
    
    Physics:
        E_orb(n) = n × E_per_frame ≈ n × 45~46 mJ
        dE/dn ≈ 45-46 mJ/frame (nearly constant)
        
    Progression:
        Orb_10: 36 frames → 1.62 J
        Orb_11: 39 frames → 1.755 J
        Orb_12: 42 frames → 1.89 J
        Orb_13: 44 frames → 2.02 J
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        orb_progression = [
            {'orb': 10, 'n_frames': 36, 'E_total': 1.62, 'E_per_frame': 0.045},
            {'orb': 11, 'n_frames': 39, 'E_total': 1.755, 'E_per_frame': 0.045},
            {'orb': 12, 'n_frames': 42, 'E_total': 1.89, 'E_per_frame': 0.045},
            {'orb': 13, 'n_frames': 44, 'E_total': 2.02, 'E_per_frame': 0.046},
        ]
        
        # Calculate deltas
        for i in range(1, len(orb_progression)):
            orb_progression[i]['delta_frames'] = orb_progression[i]['n_frames'] - orb_progression[i-1]['n_frames']
            orb_progression[i]['delta_E'] = round(orb_progression[i]['E_total'] - orb_progression[i-1]['E_total'], 4)
        orb_progression[0]['delta_frames'] = 0
        orb_progression[0]['delta_E'] = 0.0
        
        # Linear fit estimate
        n_values = [o['n_frames'] for o in orb_progression]
        E_values = [o['E_total'] for o in orb_progression]
        n_mean = np.mean(n_values)
        E_mean = np.mean(E_values)
        
        # Slope: dE/dn
        numerator = sum([(n - n_mean) * (E - E_mean) for n, E in zip(n_values, E_values)])
        denominator = sum([(n - n_mean)**2 for n in n_values])
        slope = numerator / denominator if denominator > 0 else 0.0
        intercept = E_mean - slope * n_mean
        
        # Predicted energy at n=496 (full video)
        E_predicted_496 = slope * 496 + intercept
        
        return {
            'orb_progression': orb_progression,
            'slope_dE_dn': round(slope, 6),
            'intercept': round(intercept, 4),
            'E_predicted_496_frames': round(E_predicted_496, 2),
            'current_orb': 13,
            'current_frames': 44,
            'current_E_total': 2.02,
            'equation': 'E(n) = slope × n + intercept',
            'source': 'Grok UFT Orb Analysis_13 Energy Progression (March 4, 2025)'
        }


class WaxCapCoolingSCMCalculator:
    """
    UFT Orb Analysis_13: Wax cap cooling and [SCm] stability calculator.
    
    Models how the paraffin wax/red mercury cap (180-215°F) cools the
    plasma and stabilizes [SCm] reactivity as plasmoids rise.
    
    Physics:
        Q_cooling = m_wax × c_wax × ΔT
        τ_cooling = m_wax × c_wax / (h × A_surface)
        SCm_stability = SCm_base × exp(-z/λ_stability)
        
    Source: Grok UFT Orb Analysis_13 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        # Wax cap parameters
        T_wax_hot = dataset.get('T_wax_hot', 373)  # K (215°F ≈ 101°C)
        T_wax_cool = dataset.get('T_wax_cool', 355)  # K (180°F ≈ 82°C)
        T_ambient = dataset.get('T_ambient', ORB_ANALYSIS_13_PARAMS['T_top'])  # 288 K
        
        # Wax properties
        m_wax = dataset.get('m_wax', 0.5)  # kg (estimate)
        c_wax = dataset.get('c_wax', 2500)  # J/(kg·K) paraffin specific heat
        h_conv = dataset.get('h_conv', 10)  # W/(m²·K) natural convection
        
        # Geometry
        r_reactor = dataset.get('r_reactor', ORB_ANALYSIS_13_PARAMS['r_reactor'])
        cap_thickness = dataset.get('cap_thickness', 0.0762)  # m (3 in)
        A_surface = math.pi * r_reactor**2  # top surface area
        L_reactor = dataset.get('L_reactor', 0.254)  # m (10 in)
        
        # Cooling calculations
        Delta_T_wax = (T_wax_hot + T_wax_cool) / 2 - T_ambient  # average wax to ambient
        Q_cooling = m_wax * c_wax * Delta_T_wax  # total cooling capacity
        tau_cooling = m_wax * c_wax / (h_conv * A_surface) if h_conv * A_surface > 0 else 0.0
        
        # SCm stability profile
        SCm_base = dataset.get('SCm', ORB_ANALYSIS_13_PARAMS['SCm'])
        lambda_stability = L_reactor / 3  # characteristic decay length
        z_positions = [0, L_reactor/4, L_reactor/2, 3*L_reactor/4, L_reactor]
        SCm_profile = []
        for z in z_positions:
            SCm_z = SCm_base * math.exp(-z / lambda_stability)
            SCm_profile.append({'z': round(z, 4), 'SCm': SCm_z})
        
        # Plasmoid stabilization factor
        stabilization_factor = 1 - math.exp(-tau_cooling / 10)  # asymptotic to 1
        
        return {
            'T_wax_range': f'{T_wax_cool}→{T_wax_hot} K',
            'Delta_T_wax': round(Delta_T_wax, 2),
            'Q_cooling_capacity': round(Q_cooling, 2),
            'tau_cooling': round(tau_cooling, 2),
            'SCm_base': SCm_base,
            'lambda_stability': round(lambda_stability, 4),
            'SCm_profile': SCm_profile,
            'stabilization_factor': round(stabilization_factor, 4),
            'equation': 'SCm(z) = SCm_base × e^(-z/λ)',
            'source': 'Grok UFT Orb Analysis_13 Wax Cap Cooling (March 4, 2025)'
        }


# UFT Orb Analysis_13 registry dict
ORB_ANALYSIS_13_CALCULATORS = {
    'FortyFourFrameSequenceCalculator': FortyFourFrameSequenceCalculator(),
    'DiagonalShiftOrb13Calculator': DiagonalShiftOrb13Calculator(),
    'Orb13RefinedFUCalculator': Orb13RefinedFUCalculator(),
    'IntelligentQuantumPlasmoidCalculator': IntelligentQuantumPlasmoidCalculator(),
    'Orb13EnergyProgressionCalculator': Orb13EnergyProgressionCalculator(),
    'WaxCapCoolingSCMCalculator': WaxCapCoolingSCMCalculator(),
}


# ═══════════════════════════════════════════════════════════════════════════════
# CONDENSEDPHYSICS2 AGGREGATED REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

# Combine all CP2 calculators into master registry
CP2_CALCULATORS = {
    **ORB_ANALYSIS_10_CALCULATORS,
    **ORB_ANALYSIS_11_CALCULATORS,
    **ORB_ANALYSIS_12_CALCULATORS,
    **ORB_ANALYSIS_13_CALCULATORS,
}

# Update class count
CP2_CLASS_COUNT = len(CP2_CALCULATORS)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE EXPORTS
# ═══════════════════════════════════════════════════════════════════════════════

__all__ = [
    # Version & metadata
    'CP2_VERSION',
    'CP2_CLASS_COUNT',
    
    # Orb Analysis_10 (8 classes)
    'ORB_ANALYSIS_10_PARAMS',
    'ThirtySixFrameSequenceCalculator',
    'CyclicalConvectionPatternCalculator',
    'Orb10RefinedFUCalculator',
    'SpookyActionNonLocalTransferCalculator',
    'ThermalGradientDrivenDynamicsCalculator',
    'QuadrantTransitionTrackerCalculator',
    'ACEDCEModulatedEnergyCalculator',
    'MagneticBubbleConfinementCalculator',
    'ORB_ANALYSIS_10_CALCULATORS',
    
    # Orb Analysis_11 (9 classes)
    'ORB_ANALYSIS_11_PARAMS',
    'ThirtyNineFrameSequenceCalculator',
    'CounterClockwiseDiagonalCycleCalculator',
    'Orb11RefinedFUCalculator',
    'ExtendedCyclePatternAnalyzerCalculator',
    'IntelligentPlasmoidBehaviorCalculator',
    'BulbDrivenPlasmaEnergeticsCalculator',
    'WaxCapCoolingDynamicsCalculator',
    'FieldGeneratorResonanceCouplingCalculator',
    'TotalEnergyBudgetCalculator',
    'ORB_ANALYSIS_11_CALCULATORS',
    
    # Orb Analysis_12 (7 classes)
    'ORB_ANALYSIS_12_PARAMS',
    'FortyTwoFrameSequenceCalculator',
    'CyclicalConvectionOrb12Calculator',
    'Orb12RefinedFUCalculator',
    'FullSequenceQuadrantTrackerCalculator',
    'HydrogenBubblePathAnchorCalculator',
    'ThermalConvectionCycleCalculator',
    'CumulativeEnergyAnalyzerCalculator',
    'ORB_ANALYSIS_12_CALCULATORS',
    
    # Orb Analysis_13 (6 classes)
    'ORB_ANALYSIS_13_PARAMS',
    'FortyFourFrameSequenceCalculator',
    'DiagonalShiftOrb13Calculator',
    'Orb13RefinedFUCalculator',
    'IntelligentQuantumPlasmoidCalculator',
    'Orb13EnergyProgressionCalculator',
    'WaxCapCoolingSCMCalculator',
    'ORB_ANALYSIS_13_CALCULATORS',
    
    # Aggregated registry
    'CP2_CALCULATORS',
]
