# PAPER_373 — Morris-Thorne Wormhole Null Geodesics
## Star Magic UQFF Whitepaper Series
### Author: Daniel T. Murphy | Session 101 | Source: grok_share_11254865.txt (lines 2700–2800)
### Significance: FIRST wormhole physics integration in the Star Magic / UQFF CP pipeline.

---

## Abstract

This paper presents the implementation of Morris-Thorne traversable wormhole null geodesics
within the Star Magic UQFF framework. The wormhole metric, geodesic equations, traversal and
reflection conditions, and embedding functions are fully specified. This constitutes the
first wormhole physics to appear in the Condensed Physics calculator pipeline (CP1–CP4).
The wormhole coupling to MUGE resonance gravity via the a_worm term is derived in PAPER_375.

---

## 1. Wormhole Metric (Morris-Thorne 1988)

The static, spherically symmetric traversable wormhole metric:
$$
ds^2 = -dt^2 + dr^2 + (b^2 + r^2)(d\theta^2 + \sin^2\theta\, d\varphi^2)
$$

where $b = 1.0$ m is the throat radius. This is the simplest traversable wormhole shape function
$b(r) = b_0$ = const, corresponding to zero radial tidal forces at the throat.

---

## 2. Null Geodesic Equations

For null geodesics ($ds^2 = 0$) with conserved energy $E = dt/d\lambda$ and angular momentum
$L = (b^2 + r^2) d\varphi/d\lambda$:

$$
\frac{dr}{d\lambda} = \pm\sqrt{E^2 - \frac{L^2}{b^2 + r^2}}
$$

$$
\frac{d\varphi}{d\lambda} = \frac{L}{b^2 + r^2}
$$

$$
\frac{dt}{d\lambda} = E
$$

where $\lambda$ is the affine parameter.

---

## 3. Traversal and Reflection Cases

### Case 1: Traversal (L = 0.5, E = 1.0)
The effective potential $V_{\mathrm{eff}}(r) = L^2/(b^2+r^2) = 0.25/(1+r^2)$.
Since $V_{\mathrm{eff}} < E^2 = 1$ for all $r$, the null ray crosses the throat $r = 0$ and
continues to negative $r$ (the second universe).

### Case 2: Reflection (L = 1.5, E = 1.0)
$V_{\mathrm{eff}}(0) = L^2/b^2 = 2.25 > E^2 = 1$.  
The turning radius is:
$$
r_{\min} = \sqrt{\frac{L^2}{E^2} - b^2} = \sqrt{2.25 - 1} = \sqrt{1.25} \approx 1.118 \text{ m}
$$

| Case | L | E | Outcome | r_min (m) |
|------|---|---|---------|-----------|
| Traversal | 0.5 | 1.0 | Passes through throat | 0 |
| Reflection | 1.5 | 1.0 | Reflects at r_min | ≈1.118 |

---

## 4. Embedding Functions

The wormhole embedding in Euclidean 3-space (equatorial slice, $\theta = \pi/2$):
$$
z_{\mathrm{embed}}(r) = b \cdot \mathrm{arcsinh}(r/b)
$$
$$
\rho_{\mathrm{embed}}(r) = \sqrt{b^2 + r^2}
$$

These define the "funnel" shape visible in wormhole visualisations.

---

## 5. Numerical Integration Scheme

The C++ propagator uses a first-order Euler method for simplicity (upgradeable to RK4):
```
for step = 0..n_steps:
    dr = (±√(E²-L²/(b²+r²))) × dlambda
    dφ = (L/(b²+r²)) × dlambda
    dt = E × dlambda
    r += dr, φ += dφ, t += dt
```
Default: dlambda = 0.1 m, n_steps = 100.

---

## 6. Connection to UQFF

The wormhole throat geometry enters the MUGE resonance model via the **wormhole-MUGE coupling
term** (PAPER_375):
$$
a_{\mathrm{worm}} = \frac{f_{\mathrm{worm}} \cdot E_{\mathrm{vac,neb}}}{b^2 + r^2}
$$

where $f_{\mathrm{worm}} = 10^{-10}$ is the wormhole coupling constant and $b^2 + r^2$ is the
effective gravitational area of the wormhole funnel at radius r.

---

## 7. Implementation

**C++:** `STAR_MAGIC_09SEPT_UQFF_MODULE.cpp`, namespace `WormholeGeodesics`
- Struct: `GeodesicState {r, phi, t_coord}`
- Functions: `drdt()`, `dphidt()`, `z_embed()`, `rho_embed()`, `propagate()`, `selftest()`

**Python:** `CondensedPhysics4.py`, class `MorrisThorneWormholeNullGeodesicsCalculator` (CP4 #21)

**WOLFRAM_TERM:** `WOLFRAM_TERM_WORMHOLE_GEODESIC`

---

## References

- Morris, M.S. & Thorne, K.S. (1988). "Wormholes in spacetime and their use for interstellar
  travel: A tool for teaching general relativity." *Am. J. Phys.* 56(5): 395–412.
- Visser, M. (1995). *Lorentzian Wormholes*. Springer-Verlag.

---

*PAPER_373 | Session 101 | Star Magic UQFF Framework | ©2025 Daniel T. Murphy*
</content>
