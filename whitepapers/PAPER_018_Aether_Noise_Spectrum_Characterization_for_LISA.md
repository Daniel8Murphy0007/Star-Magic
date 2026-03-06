# PAPER_018: Aether Noise Spectrum Characterization for LISA

**Date:** 2026-03-06  
**Domain:** 1.2 — Gravitational Waves — LISA Future Detector  
**Primary Validation File:** `validate_lisa_extended.py`  
**Calibration constants:** κ = 0.0005/day, [SSq] = 0.57  
**Validation log:** `../paper18_validate_lisa_extended.out`

## Abstract
This paper captures the LISA aether-noise spectral characterization produced by the Star-Magic LISA extended validator. Key metrics (frequency range, power fraction, peak frequency, integrated SNR, detectability) are extracted directly from the validator log to ensure reproducibility.

## Validation
Run command (Windows-safe):
```bash
PYTHONUTF8=1 python3 ./validate_lisa_extended.py > paper18_validate_lisa_extended.out 2>&1
```

## Results (extracted from `paper18_validate_lisa_extended.out`)
| Metric | Value |
|---|---:|
| Frequency range | 0.10 - 100 mHz |
| Aether power fraction | 222.93% |
| Peak aether frequency | 0.99 mHz |
| Integrated SNR | 12695834.0 |
| Detectable (SNR > 5) | True |

### Test status
- PASSED: compute_aether_noise_spectrum
- ALL TESTS PASSED - LISA extended methods validated

## References
- `validate_lisa_extended.py`
- `paper18_validate_lisa_extended.out`
