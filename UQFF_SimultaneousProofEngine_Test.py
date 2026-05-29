#!/usr/bin/env python3
"""
UQFF_SimultaneousProofEngine_Test.py

THIS FILE RECORDS THE OUTPUT.

Per user's explicit correction:
- UQFF_SimultaneousProofEngine.py contains the calculator
  (all calculate_* functions, run_live_calculations_with_timestamps, and all live math logic)
- UQFF_SimultaneousProofEngine_Test.py RECORDS the output produced by that calculator.

This file does NOT contain calculator logic.
It does NOT implement any formulas.
Its sole job is to invoke the calculator that lives in the main engine file
and record the results with precise timestamps.

Every time this file is executed, it produces official recorded output
of the live calculations at the exact moment of recording.
"""

import datetime
import pprint
from UQFF_SimultaneousProofEngine import run_live_calculations_with_timestamps


def record_output_from_calculator(num_records: int = 3):
    """
    Records the output produced by the calculator in UQFF_SimultaneousProofEngine.py.

    Each record contains:
    - Exact wall time the recording was made
    - The complete payload returned by the calculator
    - Clear markers so the recorded block is unambiguous
    """
    print("=" * 90)
    print("UQFF_SimultaneousProofEngine_Test.py")
    print("ROLE: RECORDER OF OUTPUT")
    print("CALCULATOR SOURCE: UQFF_SimultaneousProofEngine.py (contains all live math)")
    print("This file records the results. It does not calculate.")
    print("=" * 90)

    for i in range(1, num_records + 1):
        record_time = datetime.datetime.now(datetime.timezone.utc).isoformat() + 'Z'

        print(f"\n{'='*90}")
        print(f"RECORDED OUTPUT #{i:02d}")
        print(f"RECORDED AT (UTC): {record_time}")
        print(f"SOURCE CALCULATOR: UQFF_SimultaneousProofEngine.py")
        print(f"{'='*90}")

        # Invoke the calculator that lives in the main file and capture what it returns
        calculator_payload = run_live_calculations_with_timestamps()

        print("\n--- RECORDED PAYLOAD (exact output from the calculator at the moment above) ---")
        pprint.pprint(calculator_payload, width=120, sort_dicts=False)
        print("--- END OF RECORDED PAYLOAD ---")

        print(f"\n{'-'*90}")
        print(f"END OF RECORDED OUTPUT #{i:02d}")
        print(f"{'-'*90}")

    print("\n" + "=" * 90)
    print("RECORDING SESSION COMPLETE")
    print("All output above was produced by the calculator in UQFF_SimultaneousProofEngine.py")
    print("and recorded by this file at the timestamps shown.")
    print("No calculation logic exists in this recorder file.")
    print("=" * 90)


if __name__ == '__main__':
    # Default behavior: record 3 separate timestamped outputs from the calculator
    record_output_from_calculator(num_records=3)
