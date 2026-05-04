#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_calc_tests.py  —  Realistic calculator tests across CP1-CP4
Runs representative calculators from each module with physically motivated inputs.
"""
import sys, traceback
sys.stdout.reconfigure(errors='replace')

import warnings
warnings.filterwarnings('ignore')

import dpm_vacuum_manifold as DPM
import CondensedPhysics   as CP1
import CondensedPhysics2  as CP2
import CondensedPhysics3  as CP3
import CondensedPhysics4  as CP4

def run(label, fn):
    try:
        result = fn()
        if isinstance(result, dict):
            # Print first few numeric values
            nums = {k: v for k, v in result.items() if isinstance(v, (int, float))}
            preview = ', '.join(f'{k}={v:.4g}' for k, v in list(nums.items())[:4])
            print(f'  OK  {label}: {preview}')
        elif isinstance(result, (int, float)):
            print(f'  OK  {label}: {result:.6g}')
        elif isinstance(result, (list, tuple)) and len(result) > 0:
            print(f'  OK  {label}: [{result[0]:.4g} ... len={len(result)}]')
        else:
            print(f'  OK  {label}: {str(result)[:80]}')
        return True
    except Exception as e:
        print(f'  ERR {label}: {e}')
        for l in traceback.format_exc().splitlines()[-4:]:
            print(f'      {l}')
        return False

print('=' * 72)
print('DPM VACUUM MANIFOLD — Quantum Chain Constants')
print('=' * 72)

rho_scm, _ = DPM.derive_from_quantum_chain(n_levels=26, f_SCm=0.57)
rho_ua,  _ = DPM.derive_from_quantum_chain(n_levels=26, f_SCm=5.7)
print(f'  RHO_VAC_SCM = {rho_scm:.6f} J/m³  (canonical: 633333.333)')
print(f'  RHO_VAC_UA  = {rho_ua:.6f} J/m³  (canonical: 6333333.333)')
print(f'  DPM ratio   = {rho_ua/rho_scm:.6f}  (canonical: 10.0)')
print(f'  E_CRACK     = {DPM.E_CRACK:.6e} J  (canonical: 9.9862e+22)')
print(f'  S26_3       = {DPM.S26_3:.6e}  (canonical: 1.4531e26)')
print(f'  KER_SCm     = {DPM.KER_SCm:.6e} J')

# VDS and F_U_Bi_i numerics
vds = DPM.vds_numerical(terms=26)
print(f'  VDS(26)     = {vds:.6e}')
fubi = DPM.compute_F_U_Bi_i_numerical(M_bh=1.989e30, r=6.96e8)
print(f'  F_U_Bi_i    = {fubi:.6e} (Solar params)')

print()
print('=' * 72)
print('CP1 — CondensedPhysics  (representative calculators)')
print('=' * 72)

# Solar system dataset
solar = {
    'M': 1.989e30, 'R': 6.96e8, 'r': 1.496e11,
    'B0': 1e-4, 'omega0': 2.87e-6, 'v': 4.4e5,
    'z': 0.0, 'T': 5778.0, 'L': 3.828e26,
    'rho': 1.4e3, 'SFR': 0.0, 'age': 4.6e9 * 3.156e7,
    't': 0.0, 't_n': 0.0,
}

# Find and run representative CP1 calculators
cp1_tried = 0
cp1_ok = 0
for name in dir(CP1):
    obj = getattr(CP1, name, None)
    if obj is None: continue
    if not (isinstance(obj, type) and 'Calculator' in name and 'Base' not in name):
        continue
    if cp1_tried >= 8: break
    cp1_tried += 1
    try:
        inst = obj()
        if hasattr(inst, 'compute'):
            res = inst.compute(solar)
            if isinstance(res, dict):
                nums = {k: v for k, v in res.items() if isinstance(v, (int, float))}
                preview = ', '.join(f'{k}={v:.4g}' for k, v in list(nums.items())[:3])
                print(f'  OK  CP1.{name}: {preview}')
            else:
                print(f'  OK  CP1.{name}: {str(res)[:60]}')
            cp1_ok += 1
        elif hasattr(inst, 'calculate'):
            res = inst.calculate(solar)
            print(f'  OK  CP1.{name}: {str(res)[:60]}')
            cp1_ok += 1
    except Exception as e:
        print(f'  ERR CP1.{name}: {e}')
print(f'  CP1 summary: {cp1_ok}/{cp1_tried} calculators passed')

print()
print('=' * 72)
print('CP2 — CondensedPhysics2  (representative calculators)')
print('=' * 72)

# Sagittarius A* dataset
sgrA = {
    'M': 4.1e6 * 1.989e30, 'R': 1.2e10, 'r': 2.5e20,
    'B0': 8e-3, 'omega0': 1.5e-4, 'v': 1e7,
    'z': 0.0, 'T': 1e7, 'L': 1e31,
    'rho': 1e-16, 'SFR': 0.0, 'age': 1.3e10 * 3.156e7,
    't': 0.0, 't_n': 0.0,
}

cp2_tried = 0
cp2_ok = 0
for name in dir(CP2):
    obj = getattr(CP2, name, None)
    if obj is None: continue
    if not (isinstance(obj, type) and 'Calculator' in name and 'Base' not in name):
        continue
    if cp2_tried >= 6: break
    cp2_tried += 1
    try:
        inst = obj()
        if hasattr(inst, 'compute'):
            res = inst.compute(sgrA)
            if isinstance(res, dict):
                nums = {k: v for k, v in res.items() if isinstance(v, (int, float))}
                preview = ', '.join(f'{k}={v:.4g}' for k, v in list(nums.items())[:3])
                print(f'  OK  CP2.{name}: {preview}')
            else:
                print(f'  OK  CP2.{name}: {str(res)[:60]}')
            cp2_ok += 1
    except Exception as e:
        print(f'  ERR CP2.{name}: {e}')
print(f'  CP2 summary: {cp2_ok}/{cp2_tried} calculators passed')

print()
print('=' * 72)
print('CP3 — CondensedPhysics3')
print('=' * 72)

cp3_tried = 0
cp3_ok = 0
for name in dir(CP3):
    obj = getattr(CP3, name, None)
    if obj is None: continue
    if not (isinstance(obj, type) and 'Calculator' in name and 'Base' not in name):
        continue
    if cp3_tried >= 4: break
    cp3_tried += 1
    try:
        inst = obj()
        if hasattr(inst, 'compute'):
            res = inst.compute(solar)
            if isinstance(res, dict):
                nums = {k: v for k, v in res.items() if isinstance(v, (int, float))}
                preview = ', '.join(f'{k}={v:.4g}' for k, v in list(nums.items())[:3])
                print(f'  OK  CP3.{name}: {preview}')
            else:
                print(f'  OK  CP3.{name}: {str(res)[:60]}')
            cp3_ok += 1
    except Exception as e:
        print(f'  ERR CP3.{name}: {e}')
print(f'  CP3 summary: {cp3_ok}/{cp3_tried} calculators passed')

print()
print('=' * 72)
print('CP4 — CondensedPhysics4')
print('=' * 72)

cp4_tried = 0
cp4_ok = 0
for name in dir(CP4):
    obj = getattr(CP4, name, None)
    if obj is None: continue
    if not (isinstance(obj, type) and 'Calculator' in name and 'Base' not in name):
        continue
    if cp4_tried >= 4: break
    cp4_tried += 1
    try:
        inst = obj()
        if hasattr(inst, 'compute'):
            res = inst.compute(solar)
            if isinstance(res, dict):
                nums = {k: v for k, v in res.items() if isinstance(v, (int, float))}
                preview = ', '.join(f'{k}={v:.4g}' for k, v in list(nums.items())[:3])
                print(f'  OK  CP4.{name}: {preview}')
            else:
                print(f'  OK  CP4.{name}: {str(res)[:60]}')
            cp4_ok += 1
    except Exception as e:
        print(f'  ERR CP4.{name}: {e}')
print(f'  CP4 summary: {cp4_ok}/{cp4_tried} calculators passed')

print()
print('=' * 72)
print('LENR / SCm direct functions')
print('=' * 72)

run('parkhomov_excess_heat(Ni-H default)',
    lambda: DPM.parkhomov_excess_heat())
run('pons_fleischmann_excess_heat(Pd-D default)',
    lambda: DPM.pons_fleischmann_excess_heat())
run('mckubre_lenr(PdD=0.9, V=1e-6)',
    lambda: DPM.mckubre_lenr(PdD_loading=0.9, volume=1e-6, t_n=-100.0))
run('coleman_guillespie_scm',
    lambda: DPM.coleman_guillespie_scm())
run('qgp_energy_density_scm(T=1e11)',
    lambda: DPM.qgp_energy_density_scm(T_plasma=1e11))
run('s26_3_from_vds()',
    lambda: DPM.s26_3_from_vds())

print()
print('=' * 72)
print('ALL TESTS COMPLETE')
print('=' * 72)
