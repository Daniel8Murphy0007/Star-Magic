#!/usr/bin/env python3
"""
Validation script for Drawings 1-31 Physics Models.
Tests FRB, Whittaker, Big Bang, Plasma Shield, BH Phases, and THz Holes.
"""

from CondensedPhysics import (
    FRB_MODEL, WHITTAKER_MODEL, BIG_BANG_MODEL, 
    PLASMA_SHIELD_MODEL, BH_PHASES_MODEL, THZ_HOLES_MODEL
)

def main():
    print('=' * 80)
    print('DRAWINGS 1-31 PHYSICS MODELS VALIDATION')
    print('=' * 80)

    results = []

    # 1. FRB Model
    print('\n[1/6] Fast Radio Burst Model (Drawing 1)...')
    passed, tests, summary = FRB_MODEL.validate_FRB_model()
    results.append(('FRB Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    # 2. Whittaker Decomposition Model
    print('\n[2/6] Whittaker Decomposition Model (Drawing 30)...')
    passed, tests, summary = WHITTAKER_MODEL.validate_Whittaker_model()
    results.append(('Whittaker Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    # 3. Big Bang Origin Model
    print('\n[3/6] Big Bang Origin Model (Drawings 14, 20)...')
    passed, tests, summary = BIG_BANG_MODEL.validate_BigBang_model()
    results.append(('Big Bang Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    # 4. Plasma Shield-Capture Model
    print('\n[4/6] Plasma Shield-Capture Model (Drawings 21, 28, 29)...')
    passed, tests, summary = PLASMA_SHIELD_MODEL.validate_plasma_model()
    results.append(('Plasma Shield Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    # 5. Black Hole Phases Model
    print('\n[5/6] Black Hole Phases Model (Drawings 5-9)...')
    passed, tests, summary = BH_PHASES_MODEL.validate_BH_phases_model()
    results.append(('BH Phases Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    # 6. Terahertz Holes Model
    print('\n[6/6] Terahertz Holes Model (Drawing 24)...')
    passed, tests, summary = THZ_HOLES_MODEL.validate_THz_model()
    results.append(('THz Holes Model', passed, len(tests), sum(t['passed'] for t in tests)))
    status = 'PASS' if passed else 'FAIL'
    print(f'    Result: {status} ({results[-1][3]}/{results[-1][2]} tests)')

    print('\n' + '=' * 80)
    print('FINAL SUMMARY')
    print('=' * 80)
    
    total_models = len(results)
    passed_models = sum(1 for r in results if r[1])
    total_tests = sum(r[2] for r in results)
    passed_tests = sum(r[3] for r in results)

    for name, passed, count, p_count in results:
        status = 'PASS' if passed else 'FAIL'
        mark = '[X]' if passed else '[ ]'
        print(f'  {mark} {name}: {status} ({p_count}/{count})')

    print(f'\n  Models: {passed_models}/{total_models} PASSED')
    print(f'  Tests:  {passed_tests}/{total_tests} PASSED')
    print('=' * 80)

    if passed_models == total_models:
        print('ALL MODELS VALIDATED SUCCESSFULLY!')
    else:
        print(f'WARNING: {total_models - passed_models} model(s) failed')
    print('=' * 80)


if __name__ == '__main__':
    main()
