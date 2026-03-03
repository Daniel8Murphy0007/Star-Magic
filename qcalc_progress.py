#!/usr/bin/env python3
"""
qcalc_progress.py
=================
Phase 3: Progress tracking and status updates for long calculations

Features:
- Progress callback mechanism
- Status message streaming
- Incremental result updates
- Time estimation
- Cancellation support

Author: Phase 3 Polish
Date: March 3, 2026
"""

import sys
import json
import time
from typing import Callable, Optional, Dict, Any
from enum import Enum
from dataclasses import dataclass


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATION STAGES
# ═══════════════════════════════════════════════════════════════════════════════

class CalculationStage(Enum):
    """Stages of UQFF calculation for progress tracking"""
    INIT = "Initializing"
    IMPORT_MODULES = "Importing modules"
    PARSE_INPUT = "Parsing input parameters"
    VALIDATE_PARAMS = "Validating parameters"
    COMPUTE_FU = "Computing unified field (F_U)"
    COMPUTE_UG1 = "Computing Ug1 (magnetic dipole)"
    COMPUTE_UG2 = "Computing Ug2 (charge-reactivity)"
    COMPUTE_UG3 = "Computing Ug3 (string rotation)"
    COMPUTE_UG4 = "Computing Ug4 (vacuum concentration)"
    COMPUTE_UM = "Computing Um (magnetism)"
    COMPUTE_UBI = "Computing Ubi (buoyancy)"
    FORMAT_RESULTS = "Formatting results"
    COMPLETE = "Complete"


@dataclass
class ProgressUpdate:
    """Progress update message"""
    stage: CalculationStage
    progress_percent: float
    elapsed_seconds: float
    estimated_remaining_seconds: Optional[float]
    message: str
    current_result: Optional[Dict[str, Any]] = None
    
    def to_json(self) -> str:
        """Convert to JSON string for streaming"""
        return json.dumps({
            'type': 'PROGRESS_UPDATE',
            'stage': self.stage.value,
            'progress_percent': round(self.progress_percent, 1),
            'elapsed_seconds': round(self.elapsed_seconds, 2),
            'estimated_remaining_seconds': round(self.estimated_remaining_seconds, 2) if self.estimated_remaining_seconds else None,
            'message': self.message,
            'current_result': self.current_result
        })


# ═══════════════════════════════════════════════════════════════════════════════
# PROGRESS TRACKER
# ═══════════════════════════════════════════════════════════════════════════════

class ProgressTracker:
    """
    Track calculation progress and emit status updates.
    
    Usage:
        tracker = ProgressTracker(callback=print_progress)
        tracker.start()
        tracker.update(CalculationStage.IMPORT_MODULES, 10)
        tracker.update(CalculationStage.COMPUTE_FU, 50, result={'F_U': 1.23e-5})
        tracker.complete()
    """
    
    def __init__(self, callback: Optional[Callable[[ProgressUpdate], None]] = None,
                 enable_stderr_output: bool = True):
        """
        Initialize progress tracker.
        
        Args:
            callback: Function to call with progress updates
            enable_stderr_output: Whether to print to stderr
        """
        self.callback = callback
        self.enable_stderr_output = enable_stderr_output
        
        self.start_time = None
        self.last_update_time = None
        self.current_stage = None
        self.current_progress = 0.0
        
        # Stage weights for progress calculation
        self.stage_weights = {
            CalculationStage.INIT: 2,
            CalculationStage.IMPORT_MODULES: 15,
            CalculationStage.PARSE_INPUT: 3,
            CalculationStage.VALIDATE_PARAMS: 5,
            CalculationStage.COMPUTE_FU: 20,
            CalculationStage.COMPUTE_UG1: 10,
            CalculationStage.COMPUTE_UG2: 10,
            CalculationStage.COMPUTE_UG3: 10,
            CalculationStage.COMPUTE_UG4: 10,
            CalculationStage.COMPUTE_UM: 8,
            CalculationStage.COMPUTE_UBI: 7,
            CalculationStage.FORMAT_RESULTS: 5,
            CalculationStage.COMPLETE: 0
        }
        
        self.total_weight = sum(self.stage_weights.values())
        self.completed_weight = 0.0
    
    def start(self):
        """Start tracking"""
        self.start_time = time.time()
        self.last_update_time = self.start_time
        self.current_progress = 0.0
        self.completed_weight = 0.0
        
        self.update(CalculationStage.INIT, 0, "Starting calculation...")
    
    def update(self, stage: CalculationStage, progress_in_stage: float = 0.0,
               message: Optional[str] = None, current_result: Optional[Dict[str, Any]] = None):
        """
        Update progress.
        
        Args:
            stage: Current calculation stage
            progress_in_stage: Progress within stage (0-100)
            message: Optional custom message
            current_result: Optional partial results
        """
        if self.start_time is None:
            self.start()
        
        # Calculate overall progress
        if self.current_stage != stage:
            # Moved to new stage, add completed weight
            if self.current_stage:
                self.completed_weight += self.stage_weights[self.current_stage]
            self.current_stage = stage
        
        # Current stage contribution
        stage_weight = self.stage_weights[stage]
        stage_progress = (progress_in_stage / 100.0) * stage_weight
        
        # Overall progress
        overall_progress = ((self.completed_weight + stage_progress) / self.total_weight) * 100.0
        self.current_progress = min(overall_progress, 99.9)  # Never show 100% until complete
        
        # Time estimation
        elapsed = time.time() - self.start_time
        if self.current_progress > 0:
            estimated_total = elapsed / (self.current_progress / 100.0)
            estimated_remaining = estimated_total - elapsed
        else:
            estimated_remaining = None
        
        # Create update
        update = ProgressUpdate(
            stage=stage,
            progress_percent=self.current_progress,
            elapsed_seconds=elapsed,
            estimated_remaining_seconds=estimated_remaining,
            message=message or stage.value,
            current_result=current_result
        )
        
        # Emit update
        self._emit_update(update)
        
        self.last_update_time = time.time()
    
    def complete(self, final_result: Optional[Dict[str, Any]] = None):
        """Mark calculation as complete"""
        elapsed = time.time() - self.start_time if self.start_time else 0.0
        
        update = ProgressUpdate(
            stage=CalculationStage.COMPLETE,
            progress_percent=100.0,
            elapsed_seconds=elapsed,
            estimated_remaining_seconds=0.0,
            message="Calculation complete",
            current_result=final_result
        )
        
        self._emit_update(update)
    
    def _emit_update(self, update: ProgressUpdate):
        """Emit progress update"""
        # Call callback if provided
        if self.callback:
            self.callback(update)
        
        # Print to stderr if enabled
        if self.enable_stderr_output:
            print(f"[Progress {update.progress_percent:.1f}%] {update.message}", 
                  file=sys.stderr, flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# CONSOLE PROGRESS BAR
# ═══════════════════════════════════════════════════════════════════════════════

def print_progress_bar(update: ProgressUpdate, bar_length: int = 40):
    """
    Print ASCII progress bar to stderr.
    
    Args:
        update: Progress update
        bar_length: Length of progress bar in characters
    """
    progress = update.progress_percent / 100.0
    filled = int(bar_length * progress)
    bar = '=' * filled + '-' * (bar_length - filled)
    
    # Format time
    elapsed_str = f"{update.elapsed_seconds:.1f}s"
    if update.estimated_remaining_seconds is not None:
        remaining_str = f"{update.estimated_remaining_seconds:.1f}s"
        time_str = f"{elapsed_str} / ~{remaining_str}"
    else:
        time_str = elapsed_str
    
    # Print (carriage return to overwrite previous line)
    print(f"\r[{bar}] {update.progress_percent:.1f}% | {update.message} | {time_str}",
          end='', file=sys.stderr, flush=True)
    
    # Newline on completion
    if update.stage == CalculationStage.COMPLETE:
        print(file=sys.stderr)


# ═══════════════════════════════════════════════════════════════════════════════
# JSON PROGRESS STREAMER (for IPC)
# ═══════════════════════════════════════════════════════════════════════════════

def json_progress_callback(update: ProgressUpdate):
    """
    Stream progress as JSON to stdout (for IPC parsing).
    
    Prints JSON on separate line with newline for easy parsing.
    """
    print(update.to_json(), flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# DUMMY PROGRESS (for non-interactive mode)
# ═══════════════════════════════════════════════════════════════════════════════

class DummyProgressTracker:
    """No-op progress tracker for when progress is disabled"""
    
    def start(self):
        pass
    
    def update(self, stage, progress=0, message=None, current_result=None):
        pass
    
    def complete(self, final_result=None):
        pass


if __name__ == '__main__':
    print("Testing progress tracker...\n")
    
    # Test with ASCII progress bar
    tracker = ProgressTracker(callback=print_progress_bar)
    
    tracker.start()
    time.sleep(0.2)
    
    tracker.update(CalculationStage.IMPORT_MODULES, 50)
    time.sleep(0.3)
    
    tracker.update(CalculationStage.IMPORT_MODULES, 100)
    time.sleep(0.2)
    
    tracker.update(CalculationStage.COMPUTE_FU, 30, message="Computing F_U (SGR 1745+29)")
    time.sleep(0.3)
    
    tracker.update(CalculationStage.COMPUTE_FU, 100, 
                   current_result={'F_U': 1.234e-5})
    time.sleep(0.2)
    
    tracker.update(CalculationStage.COMPUTE_UG1, 100)
    time.sleep(0.2)
    
    tracker.update(CalculationStage.FORMAT_RESULTS, 100)
    time.sleep(0.1)
    
    tracker.complete(final_result={'success': True, 'F_U': 1.234e-5})
    
    print("\n\n✓ Progress module ready for Phase 3 integration!")
