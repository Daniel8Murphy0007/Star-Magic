"""
build_uqff_core.py - Build uqff_core Python module with Visual Studio 2022

Run from VS Developer Command Prompt or after running vcvars64.bat
"""

import os
import sys
import subprocess
import sysconfig

def main():
    # Get paths
    python_include = sysconfig.get_path('include')
    python_lib = sysconfig.get_config_var('LIBDIR') or os.path.dirname(sys.executable) + '/libs'
    python_version = f"python{sys.version_info.major}{sys.version_info.minor}"
    
    try:
        import pybind11
        pybind_include = pybind11.get_include()
    except ImportError:
        print("ERROR: pybind11 not installed. Run: pip install pybind11")
        return 1
    
    print(f"Python: {sys.version}")
    print(f"Python Include: {python_include}")
    print(f"Python Lib: {python_lib}")
    print(f"pybind11 Include: {pybind_include}")
    
    # Check for cl.exe
    try:
        subprocess.run(['cl'], capture_output=True)
    except FileNotFoundError:
        print("\nERROR: cl.exe not found!")
        print("Run this from VS Developer Command Prompt:")
        print('  "C:\\Program Files\\Microsoft Visual Studio\\2022\\Professional\\VC\\Auxiliary\\Build\\vcvars64.bat"')
        return 1
    
    # Compile
    source = "uqff_pybind.cpp"
    output = f"uqff_core.cp{sys.version_info.major}{sys.version_info.minor}-win_amd64.pyd"
    
    cmd = [
        'cl', '/O2', '/LD', '/MD', '/EHsc', '/std:c++17',
        f'/I{python_include}',
        f'/I{pybind_include}',
        source,
        '/link',
        f'/LIBPATH:{python_lib}',
        f'{python_version}.lib',
        f'/OUT:{output}'
    ]
    
    print(f"\nCompiling: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    
    if result.returncode == 0:
        print(f"\nSUCCESS: {output} created")
        
        # Copy to site-packages
        import site
        site_packages = site.getsitepackages()[0]
        dest = os.path.join(site_packages, "uqff_core.pyd")
        import shutil
        shutil.copy(output, dest)
        print(f"Installed to: {dest}")
        
        # Test import
        print("\nTesting import...")
        test_result = subprocess.run([sys.executable, '-c', 'import uqff_core; print(f"Version: {uqff_core.__version__}")'])
        return test_result.returncode
    else:
        print(f"\nFAILED: cl.exe returned {result.returncode}")
        return result.returncode

if __name__ == '__main__':
    sys.exit(main())
