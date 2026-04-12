# PAPER_960: VDS Polylogarithm Li_{26} Reference

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 215
**Source:** ramanujan_26d_summation.py (VDSPolylog26)
**Calculator:** VDSPolylog26Calc (CP4 #544)
**CVW:** v2.0.0 compliant

---

## Abstract

Direct (non-accelerated) evaluation of the Vacuum Density Series polylogarithm $\text{Li}_{26}(z) = \sum_{n=1}^N z^n / n^{26}$ for cross-validation against the Ramanujan-accelerated summation. Both must agree to working precision.

---

## 1. VDS Polylogarithm

$$\text{Li}_{26}(z) = \sum_{n=1}^{N} \frac{z^n}{n^{26}}$$

## 2. Cross-Validation

$$|S_{26}(z) - \text{Li}_{26}(z)| \xrightarrow{N \to \infty} 0$$

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
