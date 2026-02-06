"""
121-System Full UQFF Validation Runner
Automated test of all 100 systems with new unified field equation
Generates comprehensive validation report
"""

import subprocess
import sys
import time
from pathlib import Path

def run_validation():
    """Run MAIN_1_CoAnQi.exe with Option 2 (parallel calculation)"""
    
    exe_path = Path("build_msvc/Release/MAIN_1_CoAnQi.exe")
    
    if not exe_path.exists():
        print(f"❌ ERROR: Executable not found at {exe_path}")
        return False
    
    print("=" * 80)
    print("121-SYSTEM FULL UQFF VALIDATION")
    print("=" * 80)
    print(f"Executable: {exe_path}")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    try:
        # Run executable with Option 2 piped in
        print("Launching calculator with Option 2 (Calculate ALL systems parallel)...")
        
        # Use Popen for interactive process
        process = subprocess.Popen(
            [str(exe_path)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Wait for menu to appear (look for "Enter choice:")
        output_lines = []
        menu_appeared = False
        
        for line in process.stdout:
            output_lines.append(line)
            print(line, end='')  # Echo output
            
            if "Enter choice:" in line:
                menu_appeared = True
                break
        
        if not menu_appeared:
            print("❌ ERROR: Menu did not appear")
            return False
        
        # Send Option 2
        print("\n>>> Sending Option 2 (Calculate ALL systems)...\n")
        process.stdin.write("2\n")
        process.stdin.flush()
        
        # Collect all output
        for line in process.stdout:
            output_lines.append(line)
            print(line, end='')
            
            # Stop after results shown
            if "=== CoAnQi MAIN MENU ===" in line:
                # Send Exit command
                process.stdin.write("18\n")
                process.stdin.flush()
                break
        
        # Wait for process to complete
        process.wait(timeout=60)
        
        print()
        print("=" * 80)
        print("✅ 121-SYSTEM VALIDATION COMPLETE")
        print("=" * 80)
        
        # Save output to log file
        log_path = Path(f"121_system_validation_{int(time.time())}.txt")
        with open(log_path, 'w', encoding='utf-8') as f:
            f.writelines(output_lines)
        
        print(f"📄 Output saved to: {log_path}")
        
        return True
        
    except subprocess.TimeoutExpired:
        print("❌ ERROR: Process timeout (>60s)")
        process.kill()
        return False
    except Exception as e:
        print(f"❌ ERROR: {e}")
        return False

if __name__ == "__main__":
    success = run_validation()
    sys.exit(0 if success else 1)
