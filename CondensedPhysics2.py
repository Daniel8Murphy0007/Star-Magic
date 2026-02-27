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


# ============================================================================
# UFT ORB ANALYSIS_14 CALCULATORS (6 Calculator Classes)
# Source: Grok UFT Orb Analysis_14 (March 4, 2025)
# Photos: #46-#48 of 496 Red Dwarf Reactor Plasma Orb infrared images
# Physics: 47 frames (~1.41 s at 33.3 fps), ~2.15 J total energy
#          Cyclical convection: UL (#46) → LR (#47) → UR (#48)
#          Reinforced cyclical pattern with ~0.66 s half-cycle
# ============================================================================

# UFT Orb Analysis_14 Parameters
ORB_ANALYSIS_14_PARAMS = {
    'r_reactor': 0.0889,           # m (3.5 in diameter / 2)
    'M_s': 0.5e-3,                 # kg (plasma orb mass 0.5 g)
    'omega_s': 2 * math.pi * 6000, # rad/s (6000 Hz field resonance)
    'T_base': 366,                 # K (bulb base temperature)
    'T_top': 288,                  # K (ambient top temperature)
    'B_s': 1e-3,                   # T (magnetic field from H2 bubbles)
    'SCm': 1e15,                   # kg/m³ (hypothetical density at base)
    'UA': 1e-11,                   # C (trapped Aether charge)
    'dt_frame': 0.03,              # s (frame interval)
    'n_frames_orb14': 47,          # total frames (#1-#48 minus #14)
    'total_time_orb14': 1.41,      # s (47 frames × 0.03 s)
    'E_total_orb14': 2.15,         # J (total energy ~45.7 mJ/frame × 47)
    'E_react': 1e15,               # W/m³ (reactivity with thermal decay)
    'gamma_decay': 0.001,          # decay constant
    'v_plasmoid': 0.5,             # m/s (plasmoid motion speed)
    'n_H2_bubbles': 15,            # average hydrogen bubbles (12-18)
    'half_cycle_period': 0.66,     # s (cyclical oscillation half-cycle)
}


class FortySevenFrameSequenceCalculator:
    """
    UFT Orb Analysis_14: 47-frame sequence energy calculator.
    
    Extends the video analysis to 47 frames (Photos #1-#48 minus #14),
    covering ~1.41 s at 33.3 fps. Calculates cumulative energy and temporal
    evolution of plasmoid dynamics through the extended sequence.
    
    Physics:
        E_total = Σᵢ[E_frame_i] = 47 × 45.7 mJ ≈ 2.15 J
        t_total = n_frames × dt = 47 × 0.03 s = 1.41 s
        P_output ≈ 1.53 W (scaled from 65W input)
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        n_frames = dataset.get('n_frames', ORB_ANALYSIS_14_PARAMS['n_frames_orb14'])
        dt = dataset.get('dt_frame', ORB_ANALYSIS_14_PARAMS['dt_frame'])
        E_per_frame = dataset.get('E_per_frame', 0.0457)  # J (~45.7 mJ/frame)
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_14_PARAMS['gamma_decay'])
        
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
        
        # Efficiency estimate
        P_input = 65  # W bulb power
        efficiency = P_avg / P_input if P_input > 0 else 0.0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 4),
            'E_cumulative': round(E_cumulative, 6),
            'P_avg': round(P_avg, 6),
            'P_input': P_input,
            'efficiency': round(efficiency, 6),
            'E_per_frame_avg': round(E_cumulative / n_frames, 6) if n_frames > 0 else 0.0,
            'frame_energies': frame_energies[:5],
            'equation': 'E_total = Σᵢ[E_frame_i × e^(-γtᵢ)]',
            'source': 'Grok UFT Orb Analysis_14 47-Frame Sequence (March 4, 2025)'
        }


class CyclicalConvectionOrb14Calculator:
    """
    UFT Orb Analysis_14: Cyclical convection pattern calculator for Photos #46-#48.
    
    Tracks the cyclical concentration shifts: UL (#46) → LR (#47) → UR (#48).
    This pattern reinforces the ~0.66 s half-cycle oscillation driven by
    thermal gradients and [Ub] effects, with [Ug3] and [Um] modulating paths.
    
    Physics:
        Cycle_path: UL → LR → UR (3-point convection)
        Half_cycle ≈ 0.66 s
        θ_rotation = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        cycle_sequence = dataset.get('cycle_sequence', [
            {'frame': 46, 'quadrant': 'upper_left', 'x_c': 0.25, 'y_c': 0.75},
            {'frame': 47, 'quadrant': 'lower_right', 'x_c': 0.75, 'y_c': 0.25},
            {'frame': 48, 'quadrant': 'upper_right', 'x_c': 0.75, 'y_c': 0.75},
        ])
        
        dt = dataset.get('dt_frame', ORB_ANALYSIS_14_PARAMS['dt_frame'])
        
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
        
        # Half-cycle estimate from overall sequence
        half_cycle = ORB_ANALYSIS_14_PARAMS['half_cycle_period']
        
        return {
            'cycle_sequence': cycle_sequence,
            'transitions': transitions,
            'total_distance': round(total_distance, 4),
            'total_angle_deg': round(total_angle, 2),
            'cycle_period_3frame': round(cycle_period, 4),
            'estimated_half_cycle': half_cycle,
            'pattern': 'UL → LR → UR (cyclical convection)',
            'equation': 'θ_cycle = Σᵢ[arctan2(Δyᵢ, Δxᵢ)]',
            'source': 'Grok UFT Orb Analysis_14 Cyclical Convection (March 4, 2025)'
        }


class Orb14RefinedFUCalculator:
    """
    UFT Orb Analysis_14: Refined Unified Field F_U calculator.
    
    Computes the complete unified field equation refined with Photos #46-#48 data,
    extending the temporal coverage to 1.41 s with 47 frames.
    
    Physics:
        F_U = Σᵢ[kᵢ·Ugᵢ(r,t,Mₛ,ωₛ,Tₛ,Bₛ,SCm,UA,tₙ) - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react]
              + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ]
              + (gμν + η·Tₛμν(UA,SCm,ρA))
              
    Parameters updated for 47-frame sequence:
        t_total = 1.41 s, E_total = 2.15 J
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r = dataset.get('r_reactor', ORB_ANALYSIS_14_PARAMS['r_reactor'])
        M_s = dataset.get('M_s', ORB_ANALYSIS_14_PARAMS['M_s'])
        omega_s = dataset.get('omega_s', ORB_ANALYSIS_14_PARAMS['omega_s'])
        T_base = dataset.get('T_base', ORB_ANALYSIS_14_PARAMS['T_base'])
        T_top = dataset.get('T_top', ORB_ANALYSIS_14_PARAMS['T_top'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_14_PARAMS['B_s'])
        SCm = dataset.get('SCm', ORB_ANALYSIS_14_PARAMS['SCm'])
        UA = dataset.get('UA', ORB_ANALYSIS_14_PARAMS['UA'])
        t = dataset.get('t', 0.705)  # midpoint of 1.41 s sequence
        gamma = dataset.get('gamma_decay', ORB_ANALYSIS_14_PARAMS['gamma_decay'])
        E_react = dataset.get('E_react', ORB_ANALYSIS_14_PARAMS['E_react'])
        
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
            'source': 'Grok UFT Orb Analysis_14 Refined F_U (March 4, 2025)'
        }


class HalfCycleOscillationCalculator:
    """
    UFT Orb Analysis_14: Half-cycle oscillation calculator.
    
    Models the ~0.66 s half-cycle oscillation observed in the cyclical
    convection pattern, driven by [Ub] and thermal gradients.
    
    Physics:
        τ_half = L / v_convection ≈ 0.66 s
        f_oscillation = 1 / (2 × τ_half)
        ω_cycle = 2π × f_oscillation
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        L_reactor = dataset.get('L_reactor', 0.254)  # m (10 in)
        v_plasmoid = dataset.get('v_plasmoid', ORB_ANALYSIS_14_PARAMS['v_plasmoid'])
        t_total = dataset.get('t_total', ORB_ANALYSIS_14_PARAMS['total_time_orb14'])
        
        # Half-cycle period
        tau_half = L_reactor / v_plasmoid if v_plasmoid > 0 else 0.0
        
        # Full cycle period
        tau_full = 2 * tau_half
        
        # Oscillation frequency
        f_osc = 1 / tau_full if tau_full > 0 else 0.0
        
        # Angular frequency
        omega_cycle = 2 * math.pi * f_osc
        
        # Number of cycles in sequence
        n_cycles = t_total / tau_full if tau_full > 0 else 0.0
        
        # Phase at end of sequence
        phi_end = 2 * math.pi * n_cycles
        
        return {
            'L_reactor': L_reactor,
            'v_plasmoid': v_plasmoid,
            'tau_half': round(tau_half, 4),
            'tau_full': round(tau_full, 4),
            'f_oscillation': round(f_osc, 4),
            'omega_cycle': round(omega_cycle, 4),
            'n_cycles_in_sequence': round(n_cycles, 2),
            'phi_end_rad': round(phi_end, 4),
            'phi_end_deg': round(phi_end * 180 / math.pi, 2),
            'equation': 'τ_half = L/v, f = 1/(2τ)',
            'source': 'Grok UFT Orb Analysis_14 Half-Cycle Oscillation (March 4, 2025)'
        }


class CelestialDynamicsComparisonCalculator:
    """
    UFT Orb Analysis_14: Celestial dynamics comparison calculator.
    
    Compares the observed plasmoid dynamics to celestial phenomena:
    - Red dwarf stability via [R_orbit]
    - Black hole jet analogs via [USm]_jet
    - Stellar spin and multi-axial rotation
    
    Physics:
        Similarity_metric = f(spin, shape_shift, brightness, independence)
        R_orbit_analog = r_reactor × (T_plasmoid / T_star)^(1/4)
        USm_jet = v_plasmoid × B_s / μ₀
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        r_reactor = dataset.get('r_reactor', ORB_ANALYSIS_14_PARAMS['r_reactor'])
        v_plasmoid = dataset.get('v_plasmoid', ORB_ANALYSIS_14_PARAMS['v_plasmoid'])
        B_s = dataset.get('B_s', ORB_ANALYSIS_14_PARAMS['B_s'])
        T_plasmoid = dataset.get('T_plasmoid', 3000)  # K (approx plasma temp)
        T_red_dwarf = dataset.get('T_red_dwarf', 3500)  # K (typical red dwarf)
        
        mu_0 = 4 * math.pi * 1e-7  # H/m
        
        # R_orbit analog (scaling relation)
        R_orbit_analog = r_reactor * (T_plasmoid / T_red_dwarf)**(1/4)
        
        # USm_jet analog (magnetic jet velocity-field product)
        USm_jet = v_plasmoid * B_s / mu_0
        
        # Celestial comparison features
        celestial_features = {
            'multi_axial_rotation': {'observed': True, 'celestial_analog': 'stellar rotation'},
            'spin_drift': {'observed': True, 'celestial_analog': 'precession'},
            'instant_reverse': {'observed': True, 'celestial_analog': 'magnetic reconnection'},
            'non_interfering_passthrough': {'observed': True, 'celestial_analog': 'collisionless plasma'},
            'shape_shifting': {'observed': True, 'celestial_analog': 'stellar pulsation'},
            'brightness_variation': {'observed': True, 'celestial_analog': 'stellar variability'},
        }
        
        # Similarity score
        similarity_score = sum([1 for f in celestial_features.values() if f['observed']]) / len(celestial_features)
        
        return {
            'R_orbit_analog': round(R_orbit_analog, 6),
            'USm_jet': round(USm_jet, 4),
            'T_plasmoid': T_plasmoid,
            'T_red_dwarf': T_red_dwarf,
            'celestial_features': celestial_features,
            'similarity_score': similarity_score,
            'cosmic_applications': [
                'Red dwarf stability modeling',
                'Black hole jet analogs',
                'Stellar plasma dynamics'
            ],
            'equation': 'R_orbit = r × (T_p/T_s)^(1/4), USm = vB/μ₀',
            'source': 'Grok UFT Orb Analysis_14 Celestial Dynamics (March 4, 2025)'
        }


class Orb14EnergyEfficiencyCalculator:
    """
    UFT Orb Analysis_14: Energy efficiency and power scaling calculator.
    
    Analyzes the energy output efficiency of the reactor system:
    - 65W input from incandescent bulb
    - ~2.15 J output over 1.41 s (~1.53 W average)
    - Efficiency and scaling factors
    
    Physics:
        η = P_out / P_in
        P_out = E_total / t_total ≈ 1.53 W
        
    Source: Grok UFT Orb Analysis_14 (March 4, 2025)
    """
    
    def compute(self, dataset: dict) -> dict:
        E_total = dataset.get('E_total', ORB_ANALYSIS_14_PARAMS['E_total_orb14'])
        t_total = dataset.get('t_total', ORB_ANALYSIS_14_PARAMS['total_time_orb14'])
        P_input = dataset.get('P_input', 65)  # W (bulb power)
        
        # Average output power
        P_output = E_total / t_total if t_total > 0 else 0.0
        
        # Efficiency
        efficiency = P_output / P_input if P_input > 0 else 0.0
        
        # Energy input during sequence
        E_input = P_input * t_total
        
        # Energy ratio (captured fraction)
        energy_ratio = E_total / E_input if E_input > 0 else 0.0
        
        # Projected full video energy (496 frames)
        n_full = 496 - 1  # minus #14
        E_full_projected = (E_total / 47) * n_full
        t_full_projected = n_full * ORB_ANALYSIS_14_PARAMS['dt_frame']
        
        return {
            'E_total': round(E_total, 4),
            't_total': round(t_total, 4),
            'P_input': P_input,
            'P_output': round(P_output, 4),
            'efficiency': round(efficiency, 6),
            'E_input': round(E_input, 4),
            'energy_ratio': round(energy_ratio, 6),
            'n_full_frames': n_full,
            'E_full_projected': round(E_full_projected, 2),
            't_full_projected': round(t_full_projected, 2),
            'equation': 'η = P_out/P_in, P_out = E/t',
            'source': 'Grok UFT Orb Analysis_14 Energy Efficiency (March 4, 2025)'
        }


# UFT Orb Analysis_14 registry dict
ORB_ANALYSIS_14_CALCULATORS = {
    'FortySevenFrameSequenceCalculator': FortySevenFrameSequenceCalculator(),
    'CyclicalConvectionOrb14Calculator': CyclicalConvectionOrb14Calculator(),
    'Orb14RefinedFUCalculator': Orb14RefinedFUCalculator(),
    'HalfCycleOscillationCalculator': HalfCycleOscillationCalculator(),
    'CelestialDynamicsComparisonCalculator': CelestialDynamicsComparisonCalculator(),
    'Orb14EnergyEfficiencyCalculator': Orb14EnergyEfficiencyCalculator(),
}


# ═══════════════════════════════════════════════════════════════════════════════
# UFT ORB ANALYSIS_15: DATASET COMPLETION & ERROR REDUCTION
# Photos #49-#50, Photo #14 restored, 500 total images, 2 min 45 s
# Plasmoid spin rate: 0.15 rot/s (~9 RPM), 6.6 s cycle
# Error reduction: >50% across all F_U components
# Exp_2 preview: 4,943 images, 1 hr 15 min 7 sec
# Source: https://grok.com/share/bGVnYWN5LWNvcHk_ee173e65-63cd-41de-9a62-f4153aba3858
# ═══════════════════════════════════════════════════════════════════════════════

# UFT Orb Analysis_15 consolidated parameters
ORB_ANALYSIS_15_PARAMS = {
    # Dataset completion (500 images total)
    'n_batches': 50,                   # 50 batches of 10 images each
    'images_per_batch': 10,            # 10 images per batch
    'n_frames_total': 500,             # Total frames in completed dataset
    'n_frames_exp1': 500,              # Exp_1 total
    't_total_exp1': 165.0,             # 2 min 45 sec in seconds
    'dt_frame_exp1': 0.33,             # seconds per frame (Exp_1)
    
    # Photos #49-#50 quadrant sequence
    'photo_49_quadrant': 'UL',         # Upper Left concentration
    'photo_50_quadrant': 'LR',         # Lower Right concentration
    'photo_14_quadrant': 'LL',         # Lower Left (restored)
    
    # Plasmoid spin dynamics
    'spin_rate': 0.15,                 # rotations/second
    'spin_rpm': 9.0,                   # ~9 rotations per minute
    'cycle_time': 6.6,                 # seconds per full cycle
    'n_cycles_exp1': 25,               # ~25 cycles in Exp_1 (165 s / 6.6 s)
    'v_spot': 0.5,                     # m/s plasmoid velocity
    
    # Error reduction achievements
    'error_before_Ug1': 0.12,          # ±12% before
    'error_after_Ug1': 0.05,           # ±5% after (>50% reduction)
    'error_before_Ug2': 0.10,          # ±10% before
    'error_after_Ug2': 0.04,           # ±4% after
    'error_before_Ug3': 0.15,          # ±15% before
    'error_after_Ug3': 0.05,           # ±5% after
    'error_before_Ubi': 0.13,          # ±13% before
    'error_after_Ubi': 0.06,           # ±6% after
    'error_before_Um': 0.14,           # ±14% before
    'error_after_Um': 0.05,            # ±5% after
    'error_before_Amuv': 0.12,         # ±12% before
    'error_after_Amuv': 0.04,          # ±4% after
    
    # Energy metrics (finalized)
    'E_total_360_frames': 7.10,        # Joules over 360 frames
    'P_output': 0.43,                  # Watts output
    'efficiency': 0.109,               # ~10.9% of 65W input
    'P_input': 65.0,                   # 65W bulb
    
    # Exp_2 preview parameters
    'n_frames_exp2': 4943,             # Exp_2 total images
    'dt_frame_exp2': 0.033,            # 33 ms per frame
    't_total_exp2': 4507.0,            # 1 hr 15 min 7 sec in seconds
    'n_cycles_exp2': 682,              # ~682 cycles in Exp_2
    
    # Reactor constants (consolidated)
    'r_reactor': 0.0889,               # m
    'M_s': 0.5e-3,                     # kg (0.5 g plasma orb)
    'omega_s': 2 * 3.14159 * 6000,     # rad/s (6000 Hz)
    'T_base': 366.0,                   # K
    'T_top': 288.0,                    # K
    'B_s': 1e-3,                       # T
    'SCm': 1e15,                       # kg/m³
    'UA': 1e-11,                       # C
    'gamma': 0.001,                    # decay constant
    'eta_tensor': 1e-22,               # Aether tensor coupling
    'rho_A': 1e-23,                    # Aether density
}


class FiveHundredFrameDatasetCalculator:
    """
    Complete 500-frame dataset energy calculator.
    
    500 images × 0.33 s/frame = 165 seconds (2 min 45 s)
    Photo #14 (10 images) restored to complete sequence.
    
    Physics: E_total = n_frames × E_per_frame
    Energy budget: ~7.10 J over 360 frames scaled to 500 frames
    """
    
    def compute(self, n_frames: int = 500, E_per_frame: float = None) -> dict:
        """
        Compute total energy for complete 500-frame dataset.
        
        Args:
            n_frames: Total frame count (default 500)
            E_per_frame: Energy per frame (J), default from params
        
        Returns:
            Complete dataset energy analysis
        """
        if E_per_frame is None:
            # Derived from 7.10 J / 360 frames
            E_per_frame = ORB_ANALYSIS_15_PARAMS['E_total_360_frames'] / 360
        
        dt = ORB_ANALYSIS_15_PARAMS['dt_frame_exp1']
        t_total = n_frames * dt
        E_total = n_frames * E_per_frame
        
        # Power output
        P_output = E_total / t_total if t_total > 0 else 0
        P_input = ORB_ANALYSIS_15_PARAMS['P_input']
        efficiency = P_output / P_input
        
        # Cycle count
        T_cycle = ORB_ANALYSIS_15_PARAMS['cycle_time']
        n_cycles = t_total / T_cycle
        
        return {
            'n_frames': n_frames,
            'E_per_frame': round(E_per_frame, 6),
            't_total': round(t_total, 2),
            't_formatted': f'{int(t_total // 60)} min {t_total % 60:.1f} s',
            'E_total': round(E_total, 3),
            'P_output': round(P_output, 4),
            'P_input': P_input,
            'efficiency': round(efficiency, 4),
            'efficiency_percent': round(efficiency * 100, 2),
            'T_cycle': T_cycle,
            'n_cycles': round(n_cycles, 1),
            'equation': 'E_total = n × E_frame, P = E/t, η = P_out/P_in',
            'source': 'Grok UFT Orb Analysis_15 Dataset Completion (March 4, 2025)'
        }


class PlasmoidSpinRateCalculator:
    """
    Plasmoid spin rate and rotational dynamics calculator.
    
    Observed: ~0.15 rotations/second (~9 RPM)
    Cycle time: 6.6 seconds per full convection cycle
    
    Physics: ω_spin = 2π / T_cycle
             f_spin = 1 / T_cycle
             RPM = 60 × f_spin
    """
    
    def compute(self, T_cycle: float = None, t_experiment: float = None) -> dict:
        """
        Compute plasmoid spin rates and rotation counts.
        
        Args:
            T_cycle: Cycle period (s), default 6.6 s
            t_experiment: Total experiment time (s), default 165 s
        
        Returns:
            Spin dynamics analysis
        """
        if T_cycle is None:
            T_cycle = ORB_ANALYSIS_15_PARAMS['cycle_time']
        if t_experiment is None:
            t_experiment = ORB_ANALYSIS_15_PARAMS['t_total_exp1']
        
        import math
        
        # Spin frequency and rate
        f_spin = 1.0 / T_cycle
        omega_spin = 2 * math.pi * f_spin  # rad/s
        rpm = 60 * f_spin
        
        # Total rotations in experiment
        n_rotations = t_experiment / T_cycle
        
        # Velocity at reactor radius
        r = ORB_ANALYSIS_15_PARAMS['r_reactor']
        v_tangential = omega_spin * r
        v_observed = ORB_ANALYSIS_15_PARAMS['v_spot']
        
        return {
            'T_cycle': T_cycle,
            'f_spin': round(f_spin, 4),
            'omega_spin': round(omega_spin, 4),
            'rpm': round(rpm, 2),
            't_experiment': t_experiment,
            'n_rotations': round(n_rotations, 1),
            'r_reactor': r,
            'v_tangential': round(v_tangential, 4),
            'v_observed': v_observed,
            'v_ratio': round(v_tangential / v_observed, 3) if v_observed > 0 else 0,
            'equation': 'f = 1/T, ω = 2πf, RPM = 60f, v = ωr',
            'source': 'Grok UFT Orb Analysis_15 Spin Dynamics (March 4, 2025)'
        }


class ErrorReductionValidatorCalculator:
    """
    Error reduction validation calculator for F_U components.
    
    Achieves >50% error reduction across all components:
    - Ug_1: 12% → 5%
    - Ug_2: 10% → 4%
    - Ug_3: 15% → 5%
    - Ub_i: 13% → 6%
    - Um: 14% → 5%
    - A_μν: 12% → 4%
    
    Physics: Error reduction ratio = (ε_before - ε_after) / ε_before
    """
    
    def compute(self, component: str = None) -> dict:
        """
        Compute error reduction metrics for F_U components.
        
        Args:
            component: Specific component (Ug1, Ug2, Ug3, Ubi, Um, Amuv) or None for all
        
        Returns:
            Error reduction analysis
        """
        p = ORB_ANALYSIS_15_PARAMS
        
        components = {
            'Ug1': {'before': p['error_before_Ug1'], 'after': p['error_after_Ug1']},
            'Ug2': {'before': p['error_before_Ug2'], 'after': p['error_after_Ug2']},
            'Ug3': {'before': p['error_before_Ug3'], 'after': p['error_after_Ug3']},
            'Ubi': {'before': p['error_before_Ubi'], 'after': p['error_after_Ubi']},
            'Um': {'before': p['error_before_Um'], 'after': p['error_after_Um']},
            'Amuv': {'before': p['error_before_Amuv'], 'after': p['error_after_Amuv']},
        }
        
        results = {}
        for name, errors in components.items():
            if component is not None and name != component:
                continue
                
            before = errors['before']
            after = errors['after']
            reduction = (before - after) / before if before > 0 else 0
            reduction_pct = reduction * 100
            
            results[name] = {
                'error_before': round(before * 100, 1),
                'error_after': round(after * 100, 1),
                'reduction_ratio': round(reduction, 3),
                'reduction_percent': round(reduction_pct, 1),
                'meets_50pct_target': reduction >= 0.50,
            }
        
        # Aggregate stats
        all_reductions = [r['reduction_ratio'] for r in results.values()]
        avg_reduction = sum(all_reductions) / len(all_reductions) if all_reductions else 0
        all_meet_target = all(r['meets_50pct_target'] for r in results.values())
        
        return {
            'components': results,
            'avg_reduction_percent': round(avg_reduction * 100, 1),
            'all_meet_50pct_target': all_meet_target,
            'n_components': len(results),
            'equation': 'ε_reduction = (ε_before - ε_after) / ε_before',
            'source': 'Grok UFT Orb Analysis_15 Error Validation (March 4, 2025)'
        }


class Exp2PreviewCalculator:
    """
    Red Dwarf Reactor Exp_2 preview parameter calculator.
    
    Dataset: 4,943 images at 33 ms intervals
    Duration: 1 hr 15 min 7 sec (4,507 seconds)
    Cycles: ~682 convection cycles
    
    Physics: Scales Exp_1 parameters to 10× larger dataset
    """
    
    def compute(self, n_frames: int = None, dt_frame: float = None) -> dict:
        """
        Compute Exp_2 preview parameters and scaling factors.
        
        Args:
            n_frames: Exp_2 frame count (default 4943)
            dt_frame: Frame interval (default 0.033 s)
        
        Returns:
            Exp_2 parameter preview
        """
        p = ORB_ANALYSIS_15_PARAMS
        
        if n_frames is None:
            n_frames = p['n_frames_exp2']
        if dt_frame is None:
            dt_frame = p['dt_frame_exp2']
        
        t_total = n_frames * dt_frame
        T_cycle = p['cycle_time']
        n_cycles = t_total / T_cycle
        
        # Scale from Exp_1
        n_exp1 = p['n_frames_exp1']
        t_exp1 = p['t_total_exp1']
        E_exp1 = p['E_total_360_frames'] * (n_exp1 / 360)
        
        scale_factor = n_frames / n_exp1
        E_projected = E_exp1 * scale_factor
        
        # Frames per cycle
        frames_per_cycle = T_cycle / dt_frame
        
        return {
            'n_frames': n_frames,
            'dt_frame': dt_frame,
            't_total': round(t_total, 2),
            't_formatted': f'{int(t_total // 3600)} hr {int((t_total % 3600) // 60)} min {t_total % 60:.0f} s',
            'T_cycle': T_cycle,
            'n_cycles': round(n_cycles, 0),
            'frames_per_cycle': round(frames_per_cycle, 1),
            'scale_factor': round(scale_factor, 2),
            'E_exp1': round(E_exp1, 2),
            'E_projected': round(E_projected, 1),
            'n_batches': int(n_frames // 10),
            'expected_error_reduction': '≤3% (from ≤5% with 10× data)',
            'equation': 't = n × dt, n_cycles = t/T_cycle, scale = n_exp2/n_exp1',
            'source': 'Grok UFT Orb Analysis_15 Exp_2 Preview (March 4, 2025)'
        }


class QuadrantSequenceOrb15Calculator:
    """
    Photo #49-#50 quadrant sequence completion calculator.
    
    Quadrant pattern:
    - #48: UR (Upper Right) → #49: UL (Upper Left) → #50: LR (Lower Right)
    - #14 restored: LL (Lower Left)
    - Full sequence: 500 images completing the cycle
    
    Physics: θ = arctan2(Δy, Δx) for transition angle
             v = Δr / Δt for transition velocity
    """
    
    QUADRANT_CENTERS = {
        'UL': (0.25, 0.75),  # Upper Left (x, y normalized)
        'UR': (0.75, 0.75),  # Upper Right
        'LL': (0.25, 0.25),  # Lower Left
        'LR': (0.75, 0.25),  # Lower Right
    }
    
    def compute(self, sequence: list = None) -> dict:
        """
        Compute quadrant transition dynamics for Photos #48-#50.
        
        Args:
            sequence: List of quadrant codes, default ['UR', 'UL', 'LR']
        
        Returns:
            Quadrant sequence analysis
        """
        import math
        
        if sequence is None:
            sequence = ['UR', 'UL', 'LR']  # #48 → #49 → #50
        
        p = ORB_ANALYSIS_15_PARAMS
        r = p['r_reactor']
        dt = p['dt_frame_exp1']
        
        transitions = []
        for i in range(len(sequence) - 1):
            q1, q2 = sequence[i], sequence[i+1]
            c1 = self.QUADRANT_CENTERS.get(q1, (0.5, 0.5))
            c2 = self.QUADRANT_CENTERS.get(q2, (0.5, 0.5))
            
            # Distance in normalized units, scale to reactor
            dx = (c2[0] - c1[0]) * 2 * r
            dy = (c2[1] - c1[1]) * 2 * r
            dist = math.sqrt(dx**2 + dy**2)
            
            # Angle
            theta = math.degrees(math.atan2(dy, dx))
            
            # Velocity (10 frames per batch × dt)
            dt_batch = 10 * dt
            v = dist / dt_batch if dt_batch > 0 else 0
            
            transitions.append({
                'from': q1,
                'to': q2,
                'dx': round(dx, 4),
                'dy': round(dy, 4),
                'distance': round(dist, 4),
                'theta_deg': round(theta, 1),
                'dt': round(dt_batch, 3),
                'velocity': round(v, 4),
            })
        
        # Total path length
        total_dist = sum(t['distance'] for t in transitions)
        avg_velocity = sum(t['velocity'] for t in transitions) / len(transitions) if transitions else 0
        
        return {
            'sequence': sequence,
            'transitions': transitions,
            'total_distance': round(total_dist, 4),
            'avg_velocity': round(avg_velocity, 4),
            'v_observed': p['v_spot'],
            'r_reactor': r,
            'equation': 'd = √(Δx² + Δy²), θ = atan2(Δy, Δx), v = d/Δt',
            'source': 'Grok UFT Orb Analysis_15 Quadrant Sequence (March 4, 2025)'
        }


class FinalizedFURefinementCalculator:
    """
    Finalized F_U refinement calculator with consolidated parameters.
    
    F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] 
        + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] 
        + (gμν + η·Tₛμν)
    
    All parameters locked to ≤5% error after 500-image validation.
    """
    
    def compute(self, t: float = 0.0, t_n: float = 0.0) -> dict:
        """
        Compute finalized F_U field value with validated parameters.
        
        Args:
            t: Time (s) from experiment start
            t_n: Normalized time increment (t - t_0)
        
        Returns:
            Finalized F_U computation with all components
        """
        import math
        
        p = ORB_ANALYSIS_15_PARAMS
        
        # Validated parameters (≤5% error)
        r = p['r_reactor']
        M_s = p['M_s']
        omega_s = p['omega_s']
        T_base = p['T_base']
        T_top = p['T_top']
        B_s = p['B_s']
        SCm = p['SCm']
        UA = p['UA']
        gamma = p['gamma']
        eta = p['eta_tensor']
        rho_A = p['rho_A']
        
        # Coupling constants (validated)
        k_1 = 1.5e-4
        k_2 = 1.2
        k_3 = 1.8
        beta_i = 0.8
        mu_j = 1e-4
        
        # E_react with exponential decay
        E_react_0 = 1e15
        E_react = E_react_0 * math.exp(-0.001 * t)
        
        # Cosmological factors (scaled for reactor)
        Omega_g = 7.3e-16
        M_bh = 8.15e36
        d_g = 2.55e20
        
        # Time-dependent terms
        cos_pi_tn = math.cos(math.pi * t_n) if t_n > 0 else 1.0
        exp_decay = math.exp(-gamma * t)
        sin_mod = math.sin(0.001 * t)
        
        # Ug_1: Internal dipole (±5% error)
        grad_Ms_r = M_s / r
        Ug_1 = k_1 * grad_Ms_r * exp_decay * cos_pi_tn * (1 + 0.01 * sin_mod)
        
        # Ug_2: Outer field bubble (±4% error)
        S_boundary = 1.0  # At boundary r=R
        field_factor = (UA + UA) * M_s / (r**2)
        pert = 1 + 0.01 * 5e5
        Ug_2 = k_2 * field_factor * S_boundary * pert * SCm * exp_decay
        
        # Ug_3: Magnetic strings (±5% error)
        cos_omega = math.cos(omega_s * t * math.pi)
        Ug_3 = k_3 * B_s * cos_omega * SCm * exp_decay
        
        # Ub_i: Universal buoyancy (±6% error)
        buoyancy_factor = Omega_g * (M_bh / d_g)
        quantum_corr = 1 + 0.001 * 1e-21
        Ub_i = -beta_i * (Ug_1 + Ug_2 + Ug_3) * buoyancy_factor * quantum_corr * SCm * cos_pi_tn
        
        # Um: Universal magnetism (±5% error)
        Um = (mu_j / r) * (1 - exp_decay * cos_pi_tn) * SCm * exp_decay
        
        # A_μν: Aether tensor metric (±4% error)
        g_munu = 1.0  # Flat metric baseline
        T_s_munu = UA * SCm * rho_A * t_n if t_n > 0 else 0
        A_munu = g_munu + eta * T_s_munu
        
        # Total F_U
        sum_Ug = Ug_1 + Ug_2 + Ug_3
        sum_correction = beta_i * sum_Ug * Omega_g * (M_bh / d_g) * E_react
        sum_field = sum_Ug - sum_correction
        
        F_U = sum_field + Um + Ub_i + A_munu
        
        return {
            't': t,
            't_n': t_n,
            'Ug_1': Ug_1,
            'Ug_2': Ug_2,
            'Ug_3': Ug_3,
            'Ub_i': Ub_i,
            'Um': Um,
            'A_munu': A_munu,
            'E_react': E_react,
            'F_U': F_U,
            'error_bound': '≤5% across all components',
            'parameters': {
                'r': r, 'M_s': M_s, 'omega_s': omega_s,
                'T_grad': f'{T_base}K → {T_top}K',
                'B_s': B_s, 'SCm': SCm, 'UA': UA,
            },
            'equation': 'F_U = Σᵢ[kᵢ·Ugᵢ - βᵢ·Ugᵢ·Ωg·(Mbh/dg)·E_react] + Σⱼ[μⱼ/rⱼ·(1-e^(-γt·cos(πtₙ)))·φ̂ⱼ] + (gμν + η·Tₛμν)',
            'source': 'Grok UFT Orb Analysis_15 Finalized F_U (March 4, 2025)'
        }


# UFT Orb Analysis_15 registry dict
ORB_ANALYSIS_15_CALCULATORS = {
    'FiveHundredFrameDatasetCalculator': FiveHundredFrameDatasetCalculator(),
    'PlasmoidSpinRateCalculator': PlasmoidSpinRateCalculator(),
    'ErrorReductionValidatorCalculator': ErrorReductionValidatorCalculator(),
    'Exp2PreviewCalculator': Exp2PreviewCalculator(),
    'QuadrantSequenceOrb15Calculator': QuadrantSequenceOrb15Calculator(),
    'FinalizedFURefinementCalculator': FinalizedFURefinementCalculator(),
}


# ═══════════════════════════════════════════════════════════════════════════════
# UFT ORB ANALYSIS_16: EXP_2 BATCH 1 (Photos #1-#15)
# Exp_2: 4,943 images at 33 ms (0.03 s/frame), 33.3 fps
# First 15 frames analyzed (~0.45 s duration)
# Plasmoid species classification: Standard, ACE/DCE, Non-Local
# Comparison to Maxwell's, Einstein's, QFT equations
# Source: https://grok.com/share/bGVnYWN5LWNvcHk_8c320a8f-d1a0-4d5c-8421-4023097835f3
# ═══════════════════════════════════════════════════════════════════════════════

# UFT Orb Analysis_16 / Exp_2 Batch 1 parameters
ORB_ANALYSIS_16_PARAMS = {
    # Exp_2 Batch 1 dataset
    'exp_id': 'UFT_Exp_2_1_03Mar2025',
    'n_frames_batch1': 15,             # Photos #1-#15
    'fps_exp2': 33.3,                  # Frames per second
    'dt_frame_exp2': 0.03,             # seconds per frame
    't_batch1': 0.45,                  # 15 × 0.03 s
    
    # Plasmoid counts and distribution
    'n_spots_per_frame': 45,           # ~45 plasmoids/frame
    'spot_count_error': 0.05,          # ±5% counting error
    'spot_intensity_range': (0.1, 1.0),  # mJ per spot
    'spot_size_range': (0.5, 2.0),     # mm diameter
    
    # Circulation dynamics
    'v_plasmoid': 0.5,                 # m/s
    'T_cycle_33fps': 3.3,              # seconds per cycle (at 33 fps)
    'spin_rate': 0.15,                 # rotations/second
    
    # Energy metrics (Batch 1)
    'E_per_10_frames': 0.19,           # Joules
    'E_batch1_total': 0.57,            # ~0.57 J over 30 frames (15 × 2 data points)
    'efficiency_batch1': 0.0029,       # 0.29% of 65W
    
    # Plasmoid species fractions
    'frac_standard': 0.80,             # ~80% standard plasmoids
    'frac_ace_dce': 0.15,              # ~15% ACE/DCE events
    'frac_non_local': 0.05,            # ~5% non-local entities
    
    # Standard physics comparison errors
    'maxwell_initial_error': 0.15,     # ±15% before
    'maxwell_refined_error': 0.05,     # ±5% after
    'einstein_initial_error': 0.10,    # ±10% before
    'einstein_refined_error': 0.05,    # ±5% after
    'qft_initial_error': 0.15,         # ±15% before
    'qft_refined_error': 0.05,         # ±5% after
    
    # Reactor constants (inherited)
    'r_reactor': 0.0889,               # m
    'M_s': 0.5e-3,                     # kg
    'omega_s': 2 * 3.14159 * 6000,     # rad/s
    'T_base': 366.0,                   # K
    'T_top': 288.0,                    # K
    'B_s': 1e-3,                       # T
    'SCm': 1e15,                       # kg/m³
    'UA': 1e-11,                       # C
}


class PlasmoidSpeciesClassifierCalculator:
    """
    Plasmoid species classification calculator.
    
    Three species identified in Exp_2:
    1. Standard Plasmoids: ~80%, small uniform (~1 mJ, 1 mm), steady spins
    2. ACE/DCE Events: ~15%, larger brighter (~2 mJ, 2 mm), plasma excitation
    3. Non-Local Entities: ~5%, sudden positional jumps, ghost-like behavior
    
    Physics: Species fraction = N_type / N_total
    """
    
    def compute(self, n_total: int = None, fracs: dict = None) -> dict:
        """
        Classify plasmoid species distribution.
        
        Args:
            n_total: Total plasmoid count per frame (default 45)
            fracs: Custom species fractions dict
        
        Returns:
            Species classification with counts and characteristics
        """
        p = ORB_ANALYSIS_16_PARAMS
        
        if n_total is None:
            n_total = p['n_spots_per_frame']
        if fracs is None:
            fracs = {
                'standard': p['frac_standard'],
                'ace_dce': p['frac_ace_dce'],
                'non_local': p['frac_non_local'],
            }
        
        species = {
            'standard': {
                'fraction': fracs['standard'],
                'count': int(n_total * fracs['standard']),
                'intensity_mJ': 1.0,
                'size_mm': 1.0,
                'spin_rate': p['spin_rate'],
                'driver': '[Um] magnetism + 6000 Hz resonance',
                'description': 'Small, uniform spots with steady spins',
            },
            'ace_dce': {
                'fraction': fracs['ace_dce'],
                'count': int(n_total * fracs['ace_dce']),
                'intensity_mJ': 2.0,
                'size_mm': 2.0,
                'spin_rate': p['spin_rate'] * 1.2,  # Slightly faster
                'driver': '[Ug_1] dipole + [SCm] dust + H2 bubbles',
                'description': 'Larger, brighter spots indicating plasma excitation',
            },
            'non_local': {
                'fraction': fracs['non_local'],
                'count': int(n_total * fracs['non_local']),
                'intensity_mJ': 1.5,
                'size_mm': 1.5,
                'spin_rate': None,  # Discontinuous motion
                'driver': '[UA] Aether non-locality',
                'description': 'Sudden positional jumps, ghost-like behavior',
            },
        }
        
        # Validate fractions sum
        total_frac = sum(fracs.values())
        
        return {
            'n_total': n_total,
            'species': species,
            'total_fraction': round(total_frac, 2),
            'fraction_valid': abs(total_frac - 1.0) < 0.01,
            'dominant_species': 'standard',
            'equation': 'f_i = N_i / N_total',
            'source': 'Grok UFT Orb Analysis_16 Plasmoid Species (March 4, 2025)'
        }


class CirculationPatternExp2Calculator:
    """
    Exp_2 plasmoid circulation pattern analyzer.
    
    Cyclic convection: ~3.3 s cycle at 0.5 m/s
    Quadrant sequence: upper → lower → upper (Photos #1-#15)
    
    Physics: Driven by [Ub] buoyancy, thermal gradient, [Ug_3]/[Um]
             Non-local jumps via [UA]
    """
    
    def compute(self, n_frames: int = None, photo_sequence: list = None) -> dict:
        """
        Analyze circulation patterns over frame sequence.
        
        Args:
            n_frames: Number of frames analyzed (default 15)
            photo_sequence: List of quadrant concentrations per batch
        
        Returns:
            Circulation pattern analysis
        """
        p = ORB_ANALYSIS_16_PARAMS
        
        if n_frames is None:
            n_frames = p['n_frames_batch1']
        
        if photo_sequence is None:
            # Default sequence from analysis: upper → lower → upper → lower_right
            photo_sequence = [
                {'photos': '#1-#3', 'quadrant': 'upper', 'trend': 'initial distribution'},
                {'photos': '#4-#8', 'quadrant': 'upper_left → lower_left', 'trend': 'lateral shift'},
                {'photos': '#9-#12', 'quadrant': 'upper → lower → upper', 'trend': 'cyclic + non-local'},
                {'photos': '#13-#15', 'quadrant': 'lower_right', 'trend': 'concentration'},
            ]
        
        dt = p['dt_frame_exp2']
        t_total = n_frames * dt
        T_cycle = p['T_cycle_33fps']
        v = p['v_plasmoid']
        
        # Cycles completed
        n_cycles = t_total / T_cycle
        
        # Distance traveled
        d_total = v * t_total
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 3),
            'dt_frame': dt,
            'T_cycle': T_cycle,
            'n_cycles': round(n_cycles, 3),
            'v_plasmoid': v,
            'd_total': round(d_total, 3),
            'pattern_sequence': photo_sequence,
            'drivers': ['[Ub] buoyancy', 'thermal gradient 366→288 K', '[Ug_3] magnetic strings', '[Um] magnetism', '[UA] non-local jumps'],
            'equation': 'd = v × t, n_cycles = t / T_cycle',
            'source': 'Grok UFT Orb Analysis_16 Circulation (March 4, 2025)'
        }


class StandardPhysicsComparisonCalculator:
    """
    Comparison to standard physics equations calculator.
    
    Compares F_U model to:
    1. Maxwell's equations (electromagnetism)
    2. Einstein's field equations (general relativity)
    3. QFT Lagrangian (quantum field theory)
    
    Shows error reduction from initial to refined values.
    """
    
    def compute(self, component: str = None) -> dict:
        """
        Compare F_U to standard physics equations.
        
        Args:
            component: Specific comparison (maxwell, einstein, qft) or None for all
        
        Returns:
            Standard physics comparison analysis
        """
        p = ORB_ANALYSIS_16_PARAMS
        
        comparisons = {
            'maxwell': {
                'equation': '∇·E = ρ/ε₀, ∇×B = μ₀J + μ₀ε₀∂E/∂t',
                'F_U_analog': '[Um] magnetism + [Ug_3] magnetic strings',
                'deviation_cause': '[SCm] and [UA] add unverified terms',
                'error_initial': p['maxwell_initial_error'],
                'error_refined': p['maxwell_refined_error'],
                'reduction': (p['maxwell_initial_error'] - p['maxwell_refined_error']) / p['maxwell_initial_error'],
            },
            'einstein': {
                'equation': 'Gμν = (8πG/c⁴)Tμν',
                'F_U_analog': '[Ug_1] dipole + [Ug_2] outer field + [Ub_i] buoyancy',
                'deviation_cause': '[SCm] density and [UA] charge not in GR',
                'error_initial': p['einstein_initial_error'],
                'error_refined': p['einstein_refined_error'],
                'reduction': (p['einstein_initial_error'] - p['einstein_refined_error']) / p['einstein_initial_error'],
            },
            'qft': {
                'equation': 'L_QCD = Σ_q q̄(iγ^μD_μ - m_q)q - ¼G^a_μνG_a^μν',
                'F_U_analog': '[UA] Aether + non-local plasmoid behavior',
                'deviation_cause': '[UA]/[SCm] lack QFT experimental basis',
                'error_initial': p['qft_initial_error'],
                'error_refined': p['qft_refined_error'],
                'reduction': (p['qft_initial_error'] - p['qft_refined_error']) / p['qft_initial_error'],
            },
        }
        
        if component is not None and component in comparisons:
            comparisons = {component: comparisons[component]}
        
        # Summary stats
        all_reductions = [c['reduction'] for c in comparisons.values()]
        avg_reduction = sum(all_reductions) / len(all_reductions) if all_reductions else 0
        
        return {
            'comparisons': comparisons,
            'avg_error_reduction': round(avg_reduction * 100, 1),
            'all_refined_to_5pct': all(c['error_refined'] <= 0.05 for c in comparisons.values()),
            'note': 'F_U empirically fits observations but uses speculative [SCm]/[UA] entities',
            'source': 'Grok UFT Orb Analysis_16 Standard Physics Comparison (March 4, 2025)'
        }


class Exp2Batch1EnergyCalculator:
    """
    Exp_2 Batch 1 (Photos #1-#15) energy budget calculator.
    
    Energy: ~0.19 J per 10 frames, ~0.57 J over 30 frames
    Efficiency: 0.29% of 65W input (1.95% × measurement uncertainty)
    
    Physics: E = n_frames × E_per_frame
             η = E / (P_in × t)
    """
    
    def compute(self, n_frames: int = None) -> dict:
        """
        Compute energy budget for Exp_2 Batch 1.
        
        Args:
            n_frames: Number of frames (default 15)
        
        Returns:
            Energy budget analysis
        """
        p = ORB_ANALYSIS_16_PARAMS
        
        if n_frames is None:
            n_frames = p['n_frames_batch1']
        
        dt = p['dt_frame_exp2']
        t_total = n_frames * dt
        
        # Energy per 10-frame batch
        E_per_10 = p['E_per_10_frames']
        E_per_frame = E_per_10 / 10
        
        # Total energy for batch
        E_total = n_frames * E_per_frame
        
        # Efficiency
        P_in = 65.0  # Watts
        E_input = P_in * t_total
        efficiency = E_total / E_input if E_input > 0 else 0
        
        # Per-spot energy
        n_spots = p['n_spots_per_frame']
        E_per_spot = E_per_frame / n_spots if n_spots > 0 else 0
        
        # Compare to classical plasma efficiency
        classical_efficiency = 0.002  # ~0.1-0.2%
        efficiency_ratio = efficiency / classical_efficiency if classical_efficiency > 0 else 0
        
        return {
            'n_frames': n_frames,
            't_total': round(t_total, 3),
            'E_per_frame': round(E_per_frame, 4),
            'E_per_10_frames': E_per_10,
            'E_total': round(E_total, 4),
            'n_spots_per_frame': n_spots,
            'E_per_spot_mJ': round(E_per_spot * 1000, 3),
            'P_input': P_in,
            'E_input': round(E_input, 3),
            'efficiency': round(efficiency, 5),
            'efficiency_pct': round(efficiency * 100, 3),
            'classical_efficiency': classical_efficiency,
            'exceeds_classical_by': round(efficiency_ratio, 2),
            'equation': 'E = n × E_frame, η = E_out / (P_in × t)',
            'source': 'Grok UFT Orb Analysis_16 Batch 1 Energy (March 4, 2025)'
        }


class NavierStokesPlasmaFlowCalculator:
    """
    Navier-Stokes based plasma flow dynamics calculator.
    
    Standard fluid equation:
    ρ(∂v/∂t + v·∇v) = -∇P + μ∇²v + ρg
    
    Parameters: ρ ≈ 10³ kg/m³ (oil), μ ≈ 0.01 Pa·s, g ≈ 9.8 m/s²
    Observed: v ≈ 0.5 m/s with 5-7% thermal gradient error
    """
    
    def compute(self, v: float = None, T_gradient: tuple = None) -> dict:
        """
        Compute plasma flow dynamics using Navier-Stokes framework.
        
        Args:
            v: Flow velocity (m/s), default 0.5
            T_gradient: (T_base, T_top) in K
        
        Returns:
            Flow dynamics analysis
        """
        import math
        
        p = ORB_ANALYSIS_16_PARAMS
        
        if v is None:
            v = p['v_plasmoid']
        if T_gradient is None:
            T_gradient = (p['T_base'], p['T_top'])
        
        # Oil properties
        rho = 1e3  # kg/m³
        mu = 0.01  # Pa·s (dynamic viscosity)
        g = 9.8    # m/s²
        
        # Reynolds number
        L = p['r_reactor'] * 2  # Characteristic length
        Re = rho * v * L / mu
        
        # Thermal buoyancy
        T_base, T_top = T_gradient
        delta_T = T_base - T_top
        alpha = 7e-4  # Thermal expansion coefficient (1/K)
        
        # Buoyancy-driven velocity estimate
        Ra = (g * alpha * delta_T * L**3) / (mu / rho * 1.4e-7)  # Rayleigh number (approx)
        v_buoyancy = math.sqrt(g * alpha * delta_T * L) if delta_T > 0 else 0
        
        # Error analysis
        thermal_error = 0.06  # 5-7% error from thermal gradients
        v_error_range = (v * (1 - thermal_error), v * (1 + thermal_error))
        
        return {
            'v_observed': v,
            'v_buoyancy_estimate': round(v_buoyancy, 3),
            'v_match': abs(v - v_buoyancy) / v < 0.5 if v > 0 else False,
            'rho': rho,
            'mu': mu,
            'g': g,
            'L': round(L, 4),
            'Re': round(Re, 1),
            'flow_regime': 'laminar' if Re < 2300 else 'turbulent',
            'T_gradient': f'{T_base}K → {T_top}K',
            'delta_T': delta_T,
            'thermal_error_pct': round(thermal_error * 100, 1),
            'v_error_range': (round(v_error_range[0], 3), round(v_error_range[1], 3)),
            'equation': 'ρ(∂v/∂t + v·∇v) = -∇P + μ∇²v + ρg',
            'source': 'Grok UFT Orb Analysis_16 Navier-Stokes (March 4, 2025)'
        }


class PlanckBlackbodyValidatorCalculator:
    """
    Planck blackbody radiation validator for plasma background.
    
    B(λ,T) = (2hc²/λ⁵) × 1/(e^(hc/λkT) - 1)
    
    Observed: Reddish-orange at 2500-4000 K, λ ≈ 0.7-10 µm
    Deviation: 2-3% from oil/wax scattering
    """
    
    def compute(self, T: float = None, wavelength_range: tuple = None) -> dict:
        """
        Validate plasma background against Planck blackbody spectrum.
        
        Args:
            T: Temperature (K), default 3000 K (midrange)
            wavelength_range: (λ_min, λ_max) in µm
        
        Returns:
            Blackbody validation analysis
        """
        import math
        
        if T is None:
            T = 3000.0  # K (midrange)
        if wavelength_range is None:
            wavelength_range = (0.7, 10.0)  # µm
        
        # Constants
        h = 6.626e-34  # J·s
        c = 3e8        # m/s
        k = 1.38e-23   # J/K
        
        # Wien's displacement law: λ_peak = b/T
        b = 2.898e-3  # m·K
        lambda_peak_um = (b / T) * 1e6
        
        # Calculate spectral radiance at peak
        lambda_peak = b / T
        exp_term = math.exp(h * c / (lambda_peak * k * T)) - 1
        B_peak = (2 * h * c**2 / lambda_peak**5) / exp_term
        
        # Color interpretation
        if 0.6 < lambda_peak_um < 0.8:
            color = 'red'
        elif 0.8 < lambda_peak_um < 1.2:
            color = 'near-infrared'
        elif lambda_peak_um < 0.6:
            color = 'orange-yellow'
        else:
            color = 'infrared'
        
        # Deviation from oil/wax scattering
        scattering_deviation = 0.025  # 2-3%
        
        return {
            'T': T,
            'T_range': '2500-4000 K',
            'wavelength_range_um': wavelength_range,
            'lambda_peak_um': round(lambda_peak_um, 3),
            'expected_color': color,
            'observed_color': 'reddish-orange',
            'color_match': color in ['red', 'near-infrared'],
            'B_peak_W_m2_sr_m': B_peak,
            'scattering_deviation': round(scattering_deviation * 100, 1),
            'deviation_source': 'oil and wax scattering',
            'equation': 'B(λ,T) = (2hc²/λ⁵) / (e^(hc/λkT) - 1)',
            'source': 'Grok UFT Orb Analysis_16 Planck Blackbody (March 4, 2025)'
        }


# UFT Orb Analysis_16 registry dict
ORB_ANALYSIS_16_CALCULATORS = {
    'PlasmoidSpeciesClassifierCalculator': PlasmoidSpeciesClassifierCalculator(),
    'CirculationPatternExp2Calculator': CirculationPatternExp2Calculator(),
    'StandardPhysicsComparisonCalculator': StandardPhysicsComparisonCalculator(),
    'Exp2Batch1EnergyCalculator': Exp2Batch1EnergyCalculator(),
    'NavierStokesPlasmaFlowCalculator': NavierStokesPlasmaFlowCalculator(),
    'PlanckBlackbodyValidatorCalculator': PlanckBlackbodyValidatorCalculator(),
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
    **ORB_ANALYSIS_14_CALCULATORS,
    **ORB_ANALYSIS_15_CALCULATORS,
    **ORB_ANALYSIS_16_CALCULATORS,
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
    
    # Orb Analysis_14 (6 classes)
    'ORB_ANALYSIS_14_PARAMS',
    'FortySevenFrameSequenceCalculator',
    'CyclicalConvectionOrb14Calculator',
    'Orb14RefinedFUCalculator',
    'HalfCycleOscillationCalculator',
    'CelestialDynamicsComparisonCalculator',
    'Orb14EnergyEfficiencyCalculator',
    'ORB_ANALYSIS_14_CALCULATORS',
    
    # Orb Analysis_15 (6 classes)
    'ORB_ANALYSIS_15_PARAMS',
    'FiveHundredFrameDatasetCalculator',
    'PlasmoidSpinRateCalculator',
    'ErrorReductionValidatorCalculator',
    'Exp2PreviewCalculator',
    'QuadrantSequenceOrb15Calculator',
    'FinalizedFURefinementCalculator',
    'ORB_ANALYSIS_15_CALCULATORS',
    
    # Orb Analysis_16 / Exp_2 Batch 1 (6 classes)
    'ORB_ANALYSIS_16_PARAMS',
    'PlasmoidSpeciesClassifierCalculator',
    'CirculationPatternExp2Calculator',
    'StandardPhysicsComparisonCalculator',
    'Exp2Batch1EnergyCalculator',
    'NavierStokesPlasmaFlowCalculator',
    'PlanckBlackbodyValidatorCalculator',
    'ORB_ANALYSIS_16_CALCULATORS',
    
    # Aggregated registry
    'CP2_CALCULATORS',
]
