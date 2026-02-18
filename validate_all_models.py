#!/usr/bin/env python3
"""Final validation of all 10 May 2025 Document UQFF Models."""

from CondensedPhysics import (
    NGC2264Model,
    UGC10214Model,
    NGC4676Model,
    RedSpiderNebulaModel,
    NGC3372Model,
    AGCarinaeModel,
    M42Model,
    TarantulaNebulaModel,
    NGC2841Model,
    MysticMountainModel
)

def main():
    models = [
        NGC2264Model,
        UGC10214Model,
        NGC4676Model,
        RedSpiderNebulaModel,
        NGC3372Model,
        AGCarinaeModel,
        M42Model,
        TarantulaNebulaModel,
        NGC2841Model,
        MysticMountainModel
    ]

    print("=" * 60)
    print("FINAL VALIDATION: 10 May 2025 Document UQFF Models")
    print("=" * 60)

    total_passed = 0
    total_tests = 0
    all_pass = True

    for model_cls in models:
        r = model_cls.run_tests()
        passed = sum(1 for t in r["tests"] if t["passed"])
        total = len(r["tests"])
        total_passed += passed
        total_tests += total
        status = "PASS" if r["all_passed"] else "FAIL"
        if not r["all_passed"]:
            all_pass = False
        print(f'  {r["class"]}: {passed}/{total} [{status}]')
        print()

    print("=" * 60)
    print(f"TOTAL: {total_passed}/{total_tests} tests passed")
    print(f'ALL 10 MODELS: {"COMPLETE" if all_pass else "INCOMPLETE"}')
    print("=" * 60)

if __name__ == "__main__":
    main()
