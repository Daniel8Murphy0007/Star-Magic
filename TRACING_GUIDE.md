# UQFF Tracing System Guide

## Overview

The UQFF Tracing System provides comprehensive operation tracking for the MAIN_1_CoAnQi C++ application. It records execution flows, performance metrics, and detailed operational data to help understand system behavior, identify bottlenecks, and debug complex physics calculations.

## Features

- **Automatic Duration Tracking**: Each traced operation automatically measures execution time
- **Hierarchical Spans**: Create nested operation traces for detailed execution flow
- **Thread-Safe**: Safe for use in multi-threaded environments (Windows CRITICAL_SECTION)
- **Performance Metrics**: Record custom metrics with units
- **Event Logging**: Log important events with severity levels
- **File Output**: All traces written to `uqff_trace.log` for analysis

## Trace Output File

By default, traces are written to: `uqff_trace.log`

This file contains:
- Session start/end markers
- Span start/end events with timestamps
- Duration measurements in microseconds (µs)
- Span attributes (custom key-value pairs)
- Performance metrics
- Event logs with severity levels

## Usage Examples

### 1. Initialize Tracing (in main())

```cpp
#include "uqff_tracing.h"

int main() {
    // Initialize tracing system
    TRACE_INIT("uqff_trace.log");
    
    // Log initialization event
    TRACE_EVENT("Application started", TraceLevel::INFO);
    
    // Your application code here...
    
    // Shutdown tracing before exit
    TRACE_SHUTDOWN();
    return 0;
}
```

### 2. Create a Traced Operation (Span)

```cpp
// Automatic duration tracking with RAII pattern
void calculateSystem(const SystemParams& params) {
    auto span = UQFFTracer::getInstance().createSpan(
        "calculate_system_" + params.name, 
        SpanType::SYSTEM_CALCULATION
    );
    
    if (span) {
        // Add custom attributes
        span->setAttribute("system_name", params.name);
        span->setAttribute("mass_kg", params.mass);
        span->setAttribute("radius_m", params.radius);
    }
    
    // Perform calculations...
    double result = complexCalculation(params);
    
    if (span) {
        span->setAttribute("result", result);
    }
    
    // Span automatically ends when it goes out of scope
}
```

### 3. Log Performance Metrics

```cpp
// Record a metric with units
TRACE_METRIC("F_U_Bi_i_sun", 1.234e56, "N");
TRACE_METRIC("calculation_time", 1250.5, "ms");
TRACE_METRIC("physics_terms_registered", 6785, "terms");
```

### 4. Log Events

```cpp
// Different severity levels
TRACE_EVENT("System initialization complete", TraceLevel::INFO);
TRACE_EVENT("Warning: Low memory available", TraceLevel::WARN);
TRACE_EVENT("Error: Physics term validation failed", TraceLevel::ERROR);
TRACE_EVENT("Debug: Entering loop iteration 42", TraceLevel::DEBUG);
```

### 5. Manual Span Control

```cpp
// Create span that persists across larger scope
auto mainSpan = UQFFTracer::getInstance().createSpan(
    "main_execution", 
    SpanType::SYSTEM_CALCULATION
);

if (mainSpan) {
    mainSpan->setAttribute("version", "2.0");
}

// ... lots of code ...

// Manually end the span when appropriate
if (mainSpan) {
    mainSpan->end();
}
```

## Span Types

The system supports different operation categories:

- `SpanType::SYSTEM_CALCULATION` - Full system physics calculations
- `SpanType::PHYSICS_TERM_EVAL` - Individual physics term evaluations
- `SpanType::MODULE_INIT` - Module registration and initialization
- `SpanType::WOLFRAM_CALL` - Wolfram WSTP integration calls
- `SpanType::OPTIMIZATION` - Self-optimization routines
- `SpanType::VALIDATION` - Observational data validation
- `SpanType::CROSS_MODULE_COMM` - Cross-module communication
- `SpanType::SIMULATION_STEP` - Individual simulation steps
- `SpanType::STATISTICAL_ANALYSIS` - Statistical analysis operations

## Trace Levels

Event logging supports standard severity levels:

- `TraceLevel::DEBUG` (0) - Detailed debugging information
- `TraceLevel::INFO` (1) - General informational messages
- `TraceLevel::WARN` (2) - Warning messages
- `TraceLevel::ERROR` (3) - Error messages
- `TraceLevel::FATAL` (4) - Fatal error messages

Set minimum level during initialization:
```cpp
UQFFTracer::getInstance().initialize("trace.log", TraceLevel::DEBUG);
```

## Trace File Format

Example trace file output:

```
========================================
UQFF TRACING SESSION STARTED
Time: 2025-12-01 15:30:45
========================================

[SPAN_START] 2025-12-01 15:30:45 | SYSTEM_CALC | main_execution
[INFO] 2025-12-01 15:30:45 | Application initialized
[SPAN_START] 2025-12-01 15:30:45 | MODULE_INIT | register_all_physics_terms
[METRIC] 2025-12-01 15:30:46 | physics_terms_registered = 6.785000e+03 terms
[SPAN_END] 2025-12-01 15:30:46 | MODULE_INIT | register_all_physics_terms | Duration: 850432 µs | Attributes: {expected_terms=6785, actual_terms=6785}
[SPAN_START] 2025-12-01 15:30:46 | SYSTEM_CALC | system_calculation_Sun
[METRIC] 2025-12-01 15:30:47 | F_U_Bi_i_Sun = 1.234567e+56 N
[SPAN_END] 2025-12-01 15:30:47 | SYSTEM_CALC | system_calculation_Sun | Duration: 523189 µs | Attributes: {system_name=Sun, mass_kg=1.989000e+30, F_U_Bi_i_result=1.234567e+56}
[INFO] 2025-12-01 15:35:20 | User shutdown initiated
[SPAN_END] 2025-12-01 15:35:20 | SYSTEM_CALC | main_execution | Duration: 275834921 µs
========================================
UQFF TRACING SESSION ENDED
Time: 2025-12-01 15:35:20
========================================
```

## Performance Analysis

### Using Trace Data

1. **Find Bottlenecks**: Look for spans with longest durations
2. **Verify Calculations**: Check attributes for unexpected values
3. **Track Metrics**: Monitor physics_terms_registered, system counts
4. **Debug Issues**: Review event logs with ERROR/WARN levels
5. **Analyze Flow**: Follow span hierarchy to understand execution flow

### Example Analysis Queries

Search trace file for specific patterns:

```powershell
# Find all errors
Select-String -Path uqff_trace.log -Pattern "\[ERROR\]"

# Find slowest operations (>1 second = 1,000,000 µs)
Select-String -Path uqff_trace.log -Pattern "Duration: [0-9]{7,} µs"

# Find specific system calculations
Select-String -Path uqff_trace.log -Pattern "system_calculation_NGC"

# Count total spans
(Select-String -Path uqff_trace.log -Pattern "\[SPAN_START\]").Count
```

## Integration with Existing Logging

The tracing system works alongside the existing `VerboseLogger` class:

- **VerboseLogger**: User-facing output, calculation details
- **UQFFTracer**: Performance tracking, debugging, analysis

Both can be used simultaneously for comprehensive logging.

## Thread Safety

The tracing system is thread-safe for Windows environments using CRITICAL_SECTION. All public methods acquire locks before accessing shared state.

## Best Practices

1. **Initialize Early**: Call `TRACE_INIT()` at the start of main()
2. **Shutdown Cleanly**: Call `TRACE_SHUTDOWN()` before exiting
3. **Use RAII**: Let spans auto-complete by scope exit
4. **Add Context**: Use attributes to record important parameters
5. **Log Metrics**: Record key performance indicators
6. **Categorize Spans**: Use appropriate SpanType for each operation
7. **Minimize Overhead**: Don't create spans in tight inner loops

## Future Enhancements

Potential additions:
- Export to OpenTelemetry format
- Real-time trace visualization
- Remote trace collection
- Statistical trace analysis tools
- Integration with profiling tools

## Troubleshooting

**Problem**: No trace file created
- **Solution**: Check file permissions, verify TRACE_INIT() was called

**Problem**: Spans missing duration
- **Solution**: Ensure span completes (goes out of scope or .end() called)

**Problem**: Trace file too large
- **Solution**: Reduce trace level, disable DEBUG traces, rotate log files

**Problem**: Performance degradation
- **Solution**: Reduce tracing in hot paths, disable tracing for production builds

## Author

Daniel T. Murphy  
December 1, 2025  
Star-Magic UQFF Framework
