"""
Star-Magic Plug/Play Module System
==================================

Provides dynamic loading, hot-reload, and sandboxed execution
for AI Clones, Gaming Modules, and Debug/Maintenance functions.

Usage:
    from modules import get_loader, get_debugger
    from modules.module_interface import ModuleInterface, AICloneModule, GamingModule
    
    # Load a module
    loader = get_loader()
    result = loader.load_module('ai_clones/encryption/aes_256_uqff.py')
    
    # Execute
    output = loader.execute_module(result.uuid, {'operation': 'encrypt', ...})
    
    # Debug
    debugger = get_debugger()
    health = debugger.health.get_health_report(module, path)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

from .module_interface import (
    ModuleInterface,
    ModuleType,
    ModuleFormat,
    ModuleMetadata,
    ModuleCapabilities,
    ModuleState,
    AICloneModule,
    GamingModule,
    DebugModule,
    CalculationResult,
    GameFrame,
    DiagnosticReport,
)

from .module_loader import (
    ModuleLoader,
    LoadResult,
    LoaderStatus,
    get_loader,
)

__version__ = "3.0.0"
__all__ = [
    # Interface classes
    'ModuleInterface',
    'AICloneModule',
    'GamingModule', 
    'DebugModule',
    # Enums
    'ModuleType',
    'ModuleFormat',
    'LoaderStatus',
    # Data classes
    'ModuleMetadata',
    'ModuleCapabilities',
    'ModuleState',
    'CalculationResult',
    'GameFrame',
    'DiagnosticReport',
    'LoadResult',
    # Functions
    'get_loader',
    'ModuleLoader',
]
