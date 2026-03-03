#!/usr/bin/env python3
"""
qcalc_error_handler.py
======================
Phase 3: Robust error handling, retry logic, and circuit breaker

Features:
- Exponential backoff for transient failures
- Circuit breaker pattern (open/half-open/closed)
- Detailed error categorization
- Retry limits and timeout handling

Author: Phase 3 Polish
Date: March 3, 2026
"""

import time
import json
from typing import Callable, Any, Dict, Optional
from enum import Enum
from datetime import datetime, timedelta

# ═══════════════════════════════════════════════════════════════════════════════
# ERROR CATEGORIES
# ═══════════════════════════════════════════════════════════════════════════════

class ErrorCategory(Enum):
    """Error classification for appropriate handling"""
    TRANSIENT = "transient"           # Retry these (timeout, connection, temp failure)
    PERMANENT = "permanent"            # Don't retry (invalid input, missing module)
    RATE_LIMIT = "rate_limit"         # Backoff required (API limits, resource exhaustion)
    UNKNOWN = "unknown"                # Retry with caution


class CircuitState(Enum):
    """Circuit breaker states"""
    CLOSED = "closed"                  # Normal operation
    OPEN = "open"                      # Failing, reject requests
    HALF_OPEN = "half_open"           # Testing if recovered


# ═══════════════════════════════════════════════════════════════════════════════
# ERROR CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

def classify_error(error: Exception, stderr: str = "") -> ErrorCategory:
    """
    Classify error for retry decisions.
    
    Args:
        error: The exception raised
        stderr: Standard error output from subprocess
        
    Returns:
        ErrorCategory enum
    """
    error_str = str(error).lower()
    stderr_lower = stderr.lower()
    
    # Permanent errors (don't retry)
    permanent_indicators = [
        'modulenotfounderror',
        'importerror',
        'syntaxerror',
        'nameerror',
        'typeerror',
        'valueerror',
        'invalid json',
        'no such file',
        'file not found',
        'cannot open file'
    ]
    
    for indicator in permanent_indicators:
        if indicator in error_str or indicator in stderr_lower:
            return ErrorCategory.PERMANENT
    
    # Transient errors (retry)
    transient_indicators = [
        'timeout',
        'connection',
        'refused',
        'network',
        'temporary failure',
        'resource temporarily unavailable',
        'try again'
    ]
    
    for indicator in transient_indicators:
        if indicator in error_str or indicator in stderr_lower:
            return ErrorCategory.TRANSIENT
    
    # Rate limit errors
    if 'rate limit' in error_str or 'too many requests' in error_str:
        return ErrorCategory.RATE_LIMIT
    
    # Default: unknown (retry once)
    return ErrorCategory.UNKNOWN


# ═══════════════════════════════════════════════════════════════════════════════
# RETRY LOGIC WITH EXPONENTIAL BACKOFF
# ═══════════════════════════════════════════════════════════════════════════════

def retry_with_backoff(
    func: Callable,
    *args,
    max_retries: int = 3,
    initial_delay: float = 1.0,
    max_delay: float = 30.0,
    backoff_factor: float = 2.0,
    **kwargs
) -> Any:
    """
    Execute function with exponential backoff retry logic.
    
    Args:
        func: Function to execute
        *args: Positional arguments for func
        max_retries: Maximum retry attempts
        initial_delay: Initial retry delay (seconds)
        max_delay: Maximum retry delay (seconds)
        backoff_factor: Delay multiplier per retry
        **kwargs: Keyword arguments for func
        
    Returns:
        Function result or raises last exception
    """
    delay = initial_delay
    last_error = None
    
    for attempt in range(max_retries + 1):
        try:
            return func(*args, **kwargs)
            
        except Exception as e:
            last_error = e
            
            # Classify error
            category = classify_error(e)
            
            # Don't retry permanent errors
            if category == ErrorCategory.PERMANENT:
                raise
            
            # Last attempt - raise
            if attempt >= max_retries:
                raise
            
            # Log retry
            print(f"[Retry {attempt + 1}/{max_retries}] {type(e).__name__}: {str(e)}", 
                  file=sys.stderr)
            print(f"[Retry] Waiting {delay:.1f}s before retry...", file=sys.stderr)
            
            # Wait before retry
            time.sleep(delay)
            
            # Increase delay with backoff
            delay = min(delay * backoff_factor, max_delay)
    
    # Should never reach here, but just in case
    if last_error:
        raise last_error


# ═══════════════════════════════════════════════════════════════════════════════
# CIRCUIT BREAKER
# ═══════════════════════════════════════════════════════════════════════════════

class CircuitBreaker:
    """
    Circuit breaker pattern to prevent cascading failures.
    
    States:
    - CLOSED: Normal operation, requests pass through
    - OPEN: Failing, reject requests immediately
    - HALF_OPEN: Testing if service recovered
    """
    
    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: float = 60.0,
        half_open_max_calls: int = 3
    ):
        """
        Initialize circuit breaker.
        
        Args:
            failure_threshold: Failures before opening circuit
            recovery_timeout: Seconds before attempting recovery
            half_open_max_calls: Max requests in half-open state
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.half_open_max_calls = half_open_max_calls
        
        self.failure_count = 0
        self.success_count = 0
        self.state = CircuitState.CLOSED
        self.last_failure_time = None
        self.half_open_calls = 0
    
    def call(self, func: Callable, *args, **kwargs) -> Any:
        """
        Execute function through circuit breaker.
        
        Args:
            func: Function to execute
            *args, **kwargs: Function arguments
            
        Returns:
            Function result
            
        Raises:
            Exception: If circuit is open or function fails
        """
        # Check if we should attempt recovery
        if self.state == CircuitState.OPEN:
            if (self.last_failure_time and 
                time.time() - self.last_failure_time > self.recovery_timeout):
                # Attempt recovery
                self.state = CircuitState.HALF_OPEN
                self.half_open_calls = 0
                print("[Circuit Breaker] Entering HALF_OPEN state (testing recovery)", 
                      file=sys.stderr)
            else:
                raise Exception(
                    f"Circuit breaker is OPEN. Service unavailable "
                    f"(last failure: {time.time() - self.last_failure_time:.1f}s ago)"
                )
        
        # Limit calls in half-open state
        if self.state == CircuitState.HALF_OPEN:
            if self.half_open_calls >= self.half_open_max_calls:
                raise Exception("Circuit breaker HALF_OPEN limit reached")
            self.half_open_calls += 1
        
        # Execute function
        try:
            result = func(*args, **kwargs)
            self._on_success()
            return result
            
        except Exception as e:
            self._on_failure()
            raise
    
    def _on_success(self):
        """Handle successful call"""
        self.success_count += 1
        
        if self.state == CircuitState.HALF_OPEN:
            # Recovered! Close circuit
            if self.success_count >= 2:  # Need 2 successes to close
                self.state = CircuitState.CLOSED
                self.failure_count = 0
                self.success_count = 0
                print("[Circuit Breaker] Service recovered. Circuit CLOSED", 
                      file=sys.stderr)
        else:
            # Reset failure count on success
            self.failure_count = 0
    
    def _on_failure(self):
        """Handle failed call"""
        self.failure_count += 1
        self.success_count = 0
        self.last_failure_time = time.time()
        
        if self.state == CircuitState.HALF_OPEN:
            # Failed recovery attempt, reopen circuit
            self.state = CircuitState.OPEN
            print("[Circuit Breaker] Recovery failed. Circuit OPEN", file=sys.stderr)
        
        elif self.failure_count >= self.failure_threshold:
            # Too many failures, open circuit
            self.state = CircuitState.OPEN
            print(f"[Circuit Breaker] Failure threshold reached ({self.failure_count}). "
                  f"Circuit OPEN", file=sys.stderr)
    
    def reset(self):
        """Manually reset circuit breaker"""
        self.state = CircuitState.CLOSED
        self.failure_count = 0
        self.success_count = 0
        self.last_failure_time = None
        print("[Circuit Breaker] Manually reset to CLOSED", file=sys.stderr)
    
    def get_state(self) -> Dict[str, Any]:
        """Get circuit breaker state for monitoring"""
        return {
            'state': self.state.value,
            'failure_count': self.failure_count,
            'success_count': self.success_count,
            'last_failure_time': self.last_failure_time,
            'half_open_calls': self.half_open_calls if self.state == CircuitState.HALF_OPEN else None
        }


# ═══════════════════════════════════════════════════════════════════════════════
# USER-FRIENDLY ERROR MESSAGES
# ═══════════════════════════════════════════════════════════════════════════════

def format_error_response(error: Exception, category: ErrorCategory, 
                         retries_attempted: int = 0) -> Dict[str, Any]:
    """
    Format error into user-friendly JSON response.
    
    Args:
        error: The exception
        category: Error classification
        retries_attempted: Number of retries made
        
    Returns:
        JSON-serializable error response
    """
    # Base error info
    response = {
        'success': False,
        'error_type': type(error).__name__,
        'error_message': str(error),
        'error_category': category.value,
        'retries_attempted': retries_attempted,
        'timestamp': datetime.now().isoformat()
    }
    
    # User-friendly message
    if category == ErrorCategory.TRANSIENT:
        response['user_message'] = (
            "Temporary calculation failure. The system will automatically retry. "
            "If this persists, try again in a few moments."
        )
        response['suggestion'] = "Wait 30 seconds and retry your query"
        
    elif category == ErrorCategory.PERMANENT:
        response['user_message'] = (
            "Invalid request or missing dependencies. Please check your input parameters."
        )
        response['suggestion'] = "Verify mass (M), distance (r), and other parameters are valid"
        
    elif category == ErrorCategory.RATE_LIMIT:
        response['user_message'] = (
            "Rate limit exceeded. Too many requests in a short time."
        )
        response['suggestion'] = "Wait 60 seconds before submitting more queries"
        
    else:  # UNKNOWN
        response['user_message'] = (
            "Calculation failed due to unexpected error."
        )
        response['suggestion'] = "Check logs or contact support if this persists"
    
    return response


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL CIRCUIT BREAKER INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

# Global circuit breaker for QCalc/CP2 subprocess calls
CIRCUIT_BREAKER = CircuitBreaker(
    failure_threshold=5,      # Open after 5 failures
    recovery_timeout=60.0,    # Wait 60s before testing recovery
    half_open_max_calls=3     # Allow 3 test calls
)


if __name__ == '__main__':
    import sys
    
    # Test error classification
    print("Testing error classification...")
    
    test_cases = [
        (ImportError("No module named 'missing'"), ErrorCategory.PERMANENT),
        (TimeoutError("Connection timed out"), ErrorCategory.TRANSIENT),
        (Exception("Rate limit exceeded"), ErrorCategory.RATE_LIMIT),
    ]
    
    for error, expected in test_cases:
        category = classify_error(error)
        status = "✓" if category == expected else "✗"
        print(f"{status} {type(error).__name__}: {category.value} (expected: {expected.value})")
    
    print("\nError handler module ready for Phase 3 integration!")
