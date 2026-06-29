"""
setup.py - Build script for uqff_core PyBind11 module

Build:
    pip install .
    # or
    pip install -e .  # Editable mode

Usage:
    import uqff_core
    sys = uqff_core.SystemParams("Test", M=1e30, r=1e10)
    result = uqff_core.compute_full(sys)
    print(result.F_U_Bi_i)

Author: Daniel T. Murphy
"""

import os
import sys
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """CMake-based extension for pybind11"""
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """Custom build command for CMake-based extensions"""
    
    def build_extension(self, ext):
        # Try pybind11 simple build first
        import subprocess
        
        # Get pybind11 include paths
        try:
            import pybind11
            pybind11_includes = pybind11.get_include()
        except ImportError:
            raise RuntimeError("pybind11 required: pip install pybind11")
        
        # Get Python includes
        from distutils import sysconfig
        python_includes = sysconfig.get_python_inc()
        
        # Compiler flags
        extra_compile_args = []
        if sys.platform == 'win32':
            extra_compile_args = ['/std:c++17', '/O2', '/EHsc']
        else:
            extra_compile_args = ['-std=c++17', '-O3', '-fPIC']
        
        # Use setuptools to build the extension directly
        from setuptools import Extension as SetuptoolsExtension
        from setuptools.command.build_ext import build_ext as SetuptoolsBuildExt
        
        # Fallback: Build with setuptools directly
        ext_module = SetuptoolsExtension(
            'uqff_core',
            sources=['uqff_pybind.cpp'],
            include_dirs=[pybind11_includes, python_includes],
            extra_compile_args=extra_compile_args,
            language='c++',
        )
        
        # Copy attributes from CMakeExtension
        ext_module.name = ext.name
        
        # Build using parent
        build_ext.build_extension(self, ext_module)


# Standard setuptools approach (simpler, works without CMake)
def build_ext_module():
    """Create extension module for pybind11"""
    try:
        import pybind11
        pybind11_includes = [pybind11.get_include(), pybind11.get_include(user=True)]
    except ImportError:
        pybind11_includes = []
    
    # Get Python include directories
    from distutils import sysconfig
    python_includes = sysconfig.get_python_inc()
    
    extra_compile_args = []
    extra_link_args = []
    
    if sys.platform == 'win32':
        extra_compile_args = ['/std:c++17', '/O2', '/EHsc', '/MD']
    elif sys.platform == 'darwin':
        extra_compile_args = ['-std=c++17', '-O3', '-fPIC', '-stdlib=libc++']
        extra_link_args = ['-stdlib=libc++']
    else:
        extra_compile_args = ['-std=c++17', '-O3', '-fPIC']
    
    return Extension(
        'uqff_core',
        sources=['uqff_pybind.cpp'],
        include_dirs=[python_includes] + pybind11_includes,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        language='c++',
    )



def _safe_ext_modules():
    """Return [Extension(...)] if pybind11 is available, else []. Lets
    `python -m build` produce the pure-Python wheel on CI runners without
    pybind11 installed. To force C++ build, install pybind11 first:
        pip install pybind11
        python setup.py build_ext --inplace
    Or set UQFF_SKIP_CPP=1 to skip even when pybind11 is present.
    """
    import os as _os
    if _os.environ.get('UQFF_SKIP_CPP') == '1':
        print("UQFF_SKIP_CPP=1 set, skipping C++ extension build")
        return []
    try:
        import pybind11  # noqa: F401
    except ImportError:
        print("pybind11 not installed, skipping C++ extension build "
              "(set up CI/build environment with `pip install pybind11` to enable)")
        return []
    if not _os.path.isfile('uqff_pybind.cpp'):
        print("uqff_pybind.cpp not found, skipping C++ extension build")
        return []
    return [build_ext_module()]


def read_readme():
    """Read README.md with fallback for encoding issues"""
    if not os.path.exists('README.md'):
        return ''
    try:
        with open('README.md', 'r', encoding='utf-8') as f:
            return f.read()
    except (UnicodeDecodeError, IOError):
        return 'UQFF Star-Magic C++ Core - Python bindings'

setup(
    name='uqff-core',
    version='3.0.0',
    author='Daniel T. Murphy',
    author_email='daniel.murphy00@gmail.com',
    description='UQFF Star-Magic C++ Core - Python bindings',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/Daniel8Murphy0007/Star-Magic',
    ext_modules=_safe_ext_modules(),
    python_requires='>=3.8',
    install_requires=[
        'pybind11>=2.10.0',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: C++',
    ],
    keywords='uqff physics gravity cosmology astronomy unified-field',
)
