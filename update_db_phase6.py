#!/usr/bin/env python3
"""
update_db_phase6.py - Phase 6 Database Integration
==================================================

Adds Phase 6 equations (SOURCE70-71, 80) to symbolic database.

Phase 6 Scope:
- SOURCE70: M51 Whirlpool Galaxy (11 functions)
- SOURCE71: NGC1316 Fornax A (11 functions)
- SOURCE80: SMBH Binary (9 functions)

Total: 31 functions across 3 galactic systems

Author: Daniel T. Murphy
Date: February 14, 2026
"""

import sys
from SymbolicDB import SymbolicDatabase, EquationMetadata

def add_phase6_equations():
    """Add all Phase 6 equations to the database"""
    
    db = SymbolicDatabase('uqff_equations.db')
    
    # Phase 6 equations using EquationMetadata
    equations = []
    
    # ═══════════════════════════════════════════════════════════════
    # SOURCE70: M51 Whirlpool Galaxy (11 equations)
    # ═══════════════════════════════════════════════════════════════
    
    equations.append(EquationMetadata(
        id='source70_m51_gravity',
        sympy_expr='g_M51(r,t,z)',
        latex=r'g_{M51} = g_{base}(1+H)(1-B/B_c)(1+F_{env}) + \sum U_g + \Lambda c^2/3 + U_i + Q + F + DM',
        category='astrophysics.galaxy_dynamics',
        subcategory='spiral_galaxies',
        parameters=['M', 'r', 'z', 't', 'SFR', 'M_NGC5195', 'B'],
        units='m/s²',
        source_file='source70.cpp',
        source_function='M51UQFFModule::computeG',
        description='M51 Whirlpool Galaxy total gravity with NGC5195 tidal interaction',
        refs='source70.cpp M51UQFFModule',
        self_expand=True,
        version='2.0-Enhanced'
    ))
    
    for func_name, desc in [
        ('hubble', 'Time-dependent Hubble parameter'),
        ('fenv', 'Environmental forces (tidal + star formation)'),
        ('ui', 'Vacuum concentration energy')
        {
            'equation_id': 'source70_m51_fenv',
            'latex': r'F_{env} = F_{tidal} + F_{SF} = \frac{G M_{NGC5195}}{d^2} + k_{SF} \cdot SFR',
            'sympy_expr': 'F_env(t)',
            'description': 'Environmental forces on M51: tidal from NGC5195 + star formation pressure',
            'domain': 'Galaxy Dynamics',
            'complexity': 0.70,
            'source_file': 'source70.cpp',
            'tags': 'environmental,tidal,star_formation,M51,NGC5195'
        },
        {
            'equation_id': 'source70_m51_ui',
            'latex': r'U_i = \lambda_I \frac{\rho_{SCm}}{\rho_{UA}} \omega_i \cos(\pi t_n)(1 + F_{RZ})',
            'sympy_expr': 'U_i(t)',
            'description': 'Vacuum concentration energy for M51',
            'domain': 'Quantum Field Theory',
            'complexity': 0.75,
            'source_file': 'source70.cpp',
            'tags': 'vacuum,energy,quantum,M51'
        },
        {
            'equation_id': 'source70_m51_psi_integral',
            'latex': r'\int |\psi_{spiral}|^2 dV = \int |A e^{-r^2/(2\sigma^2)} e^{i(m\theta - \omega t)}|^2',
            'sympy_expr': 'psi_integral(r,theta,t)',
            'description': 'Spiral arm wavefunction integral for M51',
            'domain': 'Galactic Structure',
            'complexity': 0.80,
            'source_file': 'source70.cpp',
            'tags': 'spiral,wavefunction,M51,density_wave'
        },
        {
            'equation_id': 'source70_m51_quantum',
            'latex': r'Q = \frac{\hbar}{\Delta_{unc}} \int |\psi|^2 \frac{2\pi}{t_{Hubble}}',
            'sympy_expr': 'quantum_term(t_H, r)',
            'description': 'Quantum correction term for M51',
            'domain': 'Quantum Field Theory',
            'complexity': 0.85,
            'source_file': 'source70.cpp',
            'tags': 'quantum,uncertainty,M51'
        },
        {
            'equation_id': 'source70_m51_fluid',
            'latex': r'F_{fluid} = \rho_{fluid} V g_{base}',
            'sympy_expr': 'fluid_term(g)',
            'description': 'Fluid coupling term for M51 ISM',
            'domain': 'Fluid Dynamics',
            'complexity': 0.50,
            'source_file': 'source70.cpp',
            'tags': 'fluid,ISM,M51'
        },
        {
            'equation_id': 'source70_m51_dm',
            'latex': r'DM = (M_{vis} + M_{DM})(\frac{\delta\rho}{\rho} + \frac{3GM}{r^3})',
            'sympy_expr': 'dm_term(r)',
            'description': 'Dark matter halo contribution to M51',
            'domain': 'Dark Matter',
            'complexity': 0.70,
            'source_file': 'source70.cpp',
            'tags': 'dark_matter,halo,M51'
        },
        {
            'equation_id': 'source70_m51_ugsum',
            'latex': r'\sum U_g = U_{g1} + U_{g2} + U_{g3}^\prime + U_{g4}',
            'sympy_expr': 'Ug_sum(r)',
            'description': 'Sum of magnetic dipole, superconductor, external, and reactive terms for M51',
            'domain': 'Unified Field',
            'complexity': 0.80,
            'source_file': 'source70.cpp',
            'tags': 'Ug,magnetic,superconductor,M51'
        },
        {
            'equation_id': 'source70_m51_msf',
            'latex': r'M(t) = M_0 (1 + \frac{SFR \cdot t}{M_0})',
            'sympy_expr': 'M_sf(t)',
            'description': 'Star formation mass growth for M51',
            'domain': 'Stellar Astrophysics',
            'complexity': 0.40,
            'source_file': 'source70.cpp',
            'tags': 'star_formation,mass,M51'
        },
        {
            'equation_id': 'source70_m51_rt',
            'latex': r'r(t) = r_0 + v_r t',
            'sympy_expr': 'r_t(t)',
            'description': 'Radial expansion for M51',
            'domain': 'Kinematics',
            'complexity': 0.20,
            'source_file': 'source70.cpp',
            'tags': 'expansion,kinematics,M51'
        },
        
        # ═══════════════════════════════════════════════════════
        # SOURCE71: NGC1316 Fornax A (11 equations)
        # ═══════════════════════════════════════════════════════
        {
            'equation_id': 'source71_ngc1316_gravity',
            'latex': r'g_{NGC1316} = g_{base}(1+H)(1-B/B_c)(1+F_{env}) + \sum U_g + \Lambda c^2/3 + U_i + Q + F_{dust} + DM',
            'sympy_expr': 'g_NGC1316(r,t,z)',
            'description': 'NGC1316 Fornax A total gravity with post-merger dynamics',
            'domain': 'Galaxy Dynamics',
            'complexity': 0.95,
            'source_file': 'source71.cpp',
            'tags': 'NGC1316,Fornax_A,merger,radio_galaxy,AGN'
        },
        {
            'equation_id': 'source71_ngc1316_hubble',
            'latex': r'H(t,z) = H_0 \sqrt{\Omega_m (1+z)^3 + \Omega_\Lambda}',
            'sympy_expr': 'H_tz(z)',
            'description': 'Time-dependent Hubble parameter for NGC1316 at z=0.005',
            'domain': 'Cosmology',
            'complexity': 0.60,
            'source_file': 'source71.cpp',
            'tags': 'Hubble,expansion,cosmology,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_fenv',
            'latex': r'F_{env} = F_{tidal} + F_{cluster} = \frac{G M_{spiral}}{d^2} e^{-t/t_{merge}} + \frac{G M_{cluster}}{r^2}',
            'sympy_expr': 'F_env(t)',
            'description': 'Environmental forces on NGC1316: merger tidal + cluster disruption',
            'domain': 'Galaxy Dynamics',
            'complexity': 0.75,
            'source_file': 'source71.cpp',
            'tags': 'environmental,tidal,merger,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_mmerge',
            'latex': r'M_{merge}(t) = M_{spiral} (1 - e^{-t/t_{merge}})',
            'sympy_expr': 'M_merge(t)',
            'description': 'Merger mass accretion evolution for NGC1316',
            'domain': 'Galaxy Dynamics',
            'complexity': 0.60,
            'source_file': 'source71.cpp',
            'tags': 'merger,accretion,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_ui',
            'latex': r'U_i = \lambda_I \frac{\rho_{SCm}}{\rho_{UA}} \omega_i \cos(\pi t_n)(1 + F_{RZ})',
            'sympy_expr': 'U_i(t)',
            'description': 'Vacuum concentration energy for NGC1316',
            'domain': 'Quantum Field Theory',
            'complexity': 0.75,
            'source_file': 'source71.cpp',
            'tags': 'vacuum,energy,quantum,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_psi_integral',
            'latex': r'\int |\psi_{dust}|^2 dV',
            'sympy_expr': 'psi_integral(r,t)',
            'description': 'Dust lane wavefunction integral for NGC1316',
            'domain': 'Galactic Structure',
            'complexity': 0.75,
            'source_file': 'source71.cpp',
            'tags': 'dust,wavefunction,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_quantum',
            'latex': r'Q = \frac{\hbar}{\Delta_{unc}} \int |\psi|^2 \frac{2\pi}{t_{Hubble}}',
            'sympy_expr': 'quantum_term(t_H, r)',
            'description': 'Quantum correction term for NGC1316',
            'domain': 'Quantum Field Theory',
            'complexity': 0.85,
            'source_file': 'source71.cpp',
            'tags': 'quantum,uncertainty,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_fluid',
            'latex': r'F_{fluid} = \rho_{dust} V g_{base}',
            'sympy_expr': 'fluid_term(g)',
            'description': 'Dust coupling term for NGC1316',
            'domain': 'Fluid Dynamics',
            'complexity': 0.50,
            'source_file': 'source71.cpp',
            'tags': 'dust,fluid,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_dm',
            'latex': r'DM = (M_{vis} + M_{DM})(\frac{\delta\rho}{\rho} + \frac{3GM}{r^3})',
            'sympy_expr': 'dm_term(r)',
            'description': 'Dark matter halo contribution to NGC1316',
            'domain': 'Dark Matter',
            'complexity': 0.70,
            'source_file': 'source71.cpp',
            'tags': 'dark_matter,halo,NGC1316'
        },
        {
            'equation_id': 'source71_ngc1316_ugsum',
            'latex': r'\sum U_g = U_{g1} + U_{g2} + U_{g3}^{ext} + U_{g4}',
            'sympy_expr': 'Ug_sum(r)',
            'description': 'Sum of magnetic dipole, superconductor, external, and reactive terms for NGC1316',
            'domain': 'Unified Field',
            'complexity': 0.80,
            'source_file': 'source71.cpp',
            'tags': 'Ug,magnetic,superconductor,NGC1316,AGN'
        },
        {
            'equation_id': 'source71_ngc1316_rt',
            'latex': r'r(t) = r_0 + v_r t',
            'sympy_expr': 'r_t(t)',
            'description': 'Radial expansion for NGC1316',
            'domain': 'Kinematics',
            'complexity': 0.20,
            'source_file': 'source71.cpp',
            'tags': 'expansion,kinematics,NGC1316'
        },
        
        # ═══════════════════════════════════════════════════════
        # SOURCE80: SMBH Binary Coalescence (9 equations)
        # ═══════════════════════════════════════════════════════
        {
            'equation_id': 'source80_smbh_binary_gravity',
            'latex': r'g = \sum_i f_i \cdot \frac{\lambda_P}{2\pi} = (f_{super} + f_{fluid} + f_{quantum} + f_{Aether} + f_{react} + f_{res} + f_{DPM} + f_{THz} + f_{Ug4i}) \frac{\lambda_P}{2\pi}',
            'sympy_expr': 'g_SMBH(t,r)',
            'description': 'SMBH binary frequency-based total gravity (revolutionary approach)',
            'domain': 'Gravitational Waves',
            'complexity': 0.98,
            'source_file': 'source80.cpp',
            'tags': 'SMBH,binary,frequency,coalescence,LISA,Aether'
        },
        {
            'equation_id': 'source80_freq_super',
            'latex': r'f_{super}(t) = f_{super,0} e^{-t/t_{coal}}',
            'sympy_expr': 'f_super(t)',
            'description': 'Super frequency exponential decay during SMBH coalescence',
            'domain': 'Gravitational Waves',
            'complexity': 0.50,
            'source_file': 'source80.cpp',
            'tags': 'frequency,coalescence,SMBH,decay'
        },
        {
            'equation_id': 'source80_freq_fluid',
            'latex': r'f_{fluid}(\rho) = f_{fluid,0} \frac{\rho}{\rho_0}',
            'sympy_expr': 'f_fluid(rho)',
            'description': 'Fluid frequency density modulation for SMBH binary',
            'domain': 'Fluid Dynamics',
            'complexity': 0.40,
            'source_file': 'source80.cpp',
            'tags': 'frequency,fluid,density,SMBH'
        },
        {
            'equation_id': 'source80_freq_quantum',
            'latex': r'f_{quantum}(\Delta) = \frac{f_{quantum,0}}{\Delta_{unc}}',
            'sympy_expr': 'f_quantum(unc)',
            'description': 'Quantum frequency with uncertainty principle',
            'domain': 'Quantum Field Theory',
            'complexity': 0.60,
            'source_file': 'source80.cpp',
            'tags': 'frequency,quantum,uncertainty,SMBH'
        },
        {
            'equation_id': 'source80_freq_aether',
            'latex': r'f_{Aether} = 1.576 \times 10^{-35} \text{ Hz}',
            'sympy_expr': 'f_Aether',
            'description': 'Aether frequency constant (replaces dark energy/Lambda)',
            'domain': 'Cosmology',
            'complexity': 0.30,
            'source_file': 'source80.cpp',
            'tags': 'Aether,dark_energy,cosmology,SMBH,revolutionary'
        },
        {
            'equation_id': 'source80_freq_react',
            'latex': r'f_{react}(t) = f_{react,0} \cos(\omega t)',
            'sympy_expr': 'f_react(t)',
            'description': 'Reactive frequency oscillation (Ug4i)',
            'domain': 'Unified Field',
            'complexity': 0.50,
            'source_file': 'source80.cpp',
            'tags': 'frequency,reactive,Ug4i,SMBH'
        },
        {
            'equation_id': 'source80_freq_resonance',
            'latex': r'f_{res} = 2\pi f_{super} |\psi|^2',
            'sympy_expr': 'f_res(t)',
            'description': 'Resonance frequency from wavefunction',
            'domain': 'Quantum Field Theory',
            'complexity': 0.80,
            'source_file': 'source80.cpp',
            'tags': 'frequency,resonance,wavefunction,SMBH'
        },
        {
            'equation_id': 'source80_freq_dpm',
            'latex': r'f_{DPM} = f_{DPM,0} \frac{\rho_{vac,plasm}}{c}',
            'sympy_expr': 'f_DPM',
            'description': 'Dipole Phase Modulation frequency',
            'domain': 'Unified Field',
            'complexity': 0.70,
            'source_file': 'source80.cpp',
            'tags': 'frequency,DPM,dipole,SMBH'
        },
        {
            'equation_id': 'source80_freq_thz',
            'latex': r'f_{THz}(t) = f_{THz,0} \sin(\omega t)',
            'sympy_expr': 'f_THz(t)',
            'description': 'THz hole pipeline frequency oscillation',
            'domain': 'Unified Field',
            'complexity': 0.50,
            'source_file': 'source80.cpp',
            'tags': 'frequency,THz,oscillation,SMBH'
        },
    ]
    
    # Add all equations
    print("\n" + "="*80)
    print("PHASE 6 DATABASE UPDATE - Galaxy-Scale Extraction")
    print("="*80)
    print(f"\nAdding {len(equations)} equations from Phase 6 (SOURCE70-71, 80)...")
    print("\nBreakdown:")
    print("  SOURCE70 (M51 Whirlpool): 11 equations")
    print("  SOURCE71 (NGC1316 Fornax A): 11 equations")
    print("  SOURCE80 (SMBH Binary): 9 equations")
    print(f"  Total: 31 equations\n")
    
    success, failure = db.batch_add_equations(equations)
    
    # Get updated count
    total = db.count_equations()
    
    print("\n" + "="*80)
    print("PHASE 6 UPDATE COMPLETE")
    print("="*80)
    print(f"✓ Successfully added: {success} equations")
    if failure > 0:
        print(f"✗ Failed to add: {failure} equations")
    print(f"✓ Total database equations: {total}")
    print(f"✓ Expected: 137 (Phase 1-5) + 31 (Phase 6) = 168")
    print("="*80)
    
    return success, failure, total

if __name__ == '__main__':
    try:
        success_count, fail_count, total_count = add_phase6_equations()
        
        if fail_count > 0:
            print(f"\n⚠ Warning: {fail_count} equations failed to add")
            sys.exit(1)
        else:
            print(f"\n✓ Phase 6 database update successful!")
            print(f"✓ Database now contains {total_count} total equations")
            sys.exit(0)
            
    except Exception as e:
        print(f"\n✗ Error during Phase 6 database update: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
