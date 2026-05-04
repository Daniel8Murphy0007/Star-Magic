#!/usr/bin/env python3
import sys, time, importlib, traceback
sys.stdout.reconfigure(errors='replace')

mods = ['dpm_vacuum_manifold','QCalc','CondensedPhysics','CondensedPhysics2','CondensedPhysics3','CondensedPhysics4']
for mod in mods:
    t0 = time.time()
    try:
        m = importlib.import_module(mod)
        print(f'OK  {mod} ({time.time()-t0:.1f}s)')
    except Exception as e:
        print(f'ERR {mod}: {e}')
        lines = traceback.format_exc().splitlines()
        for l in lines[-8:]:
            print('   ', l)
