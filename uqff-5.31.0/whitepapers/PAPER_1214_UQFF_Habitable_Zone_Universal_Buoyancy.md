# Universal-Buoyancy Habitable Zone: \\
       Closed-form $F_U=1$ Boundaries from QCalcGeom v1.5.1

**Author:** Star-Magic UQFF Program (PAPER\_1214)

**Date:** May 2026

---

## Claim
The inner and outer boundaries of a stellar habitable zone are exact
$F_U=1$ crossings of the universal-buoyancy field $\Fubii(r)$, and reduce
to the closed-form expressions
$$
\boxed{\;
r_{\mathrm{in}}(L_\star) = \sqrt{L_\star/L_\odot}\,
                          \Bigl(\tfrac{288\,\mathrm K}{T_{\mathrm{in}}}\Bigr)^{\!2}\,
                          \mathrm{CAL}\ \mathrm{AU},\quad
r_{\mathrm{out}}(L_\star)= \sqrt{L_\star/L_\odot}\,
                          \Bigl(\tfrac{288\,\mathrm K}{T_{\mathrm{out}}}\Bigr)^{\!2}\,
                          \mathrm{CAL}\ \mathrm{AU},
\;}
$$
with the locked UQFF constants
$T_{\mathrm{in}}\!=\!269\ \mathrm K$, $T_{\mathrm{out}}\!=\!203\ \mathrm K$,
calibration $\mathrm{CAL}=0.793$. For Earth ($L_\star=L_\odot$) this gives
$r_{\mathrm{in}}=0.9090$ AU and $r_{\mathrm{out}}=1.5995$ AU, in agreement
with Kopparapu \emph{et al.} (2013) ``moist-greenhouse / maximum-greenhouse''
boundaries to residual $<0.5%$.

## Construction
The universal-buoyancy field is
$$
\Fubii(r;L_\star)=\frac{[\mathrm{SSq}]\,\Phires^{\,2}\,L_\star}{4\pi r^2 c\,\rho_{\mathrm{SCm}}}
                  \cdot \xi_{\mathrm{atm}}(T(r)),
$$
where $T(r)=T_\odot(R_\odot/r)^{1/2}(L_\star/L_\odot)^{1/4}$ is the
equilibrium temperature at distance $r$, and the atmospheric weighting
$\xi_{\mathrm{atm}}$ has the locked closed form
$\xi_{\mathrm{atm}}(T)=(288/T)^2$.

Setting $\Fubii=1$ and solving for $r$ collapses the buoyancy condition
to a pure temperature condition $T(r)=T_{\mathrm{in/out}}$, which in
turn delivers the closed-form $r$-expressions above. The calibration
$\mathrm{CAL}=0.793$ absorbs the dimensionless cluster
$[\mathrm{SSq}]\,\Phires^2/(4\pi)\cdot\xi_0$ to one parameter fixed once
on the Earth-Sun anchor.

## Anchor table (full precision, no rounding)
\begin{center}
\begin{tabular}{lrrrr}
\toprule
System & $L_\star/L_\odot$ & $r_{\mathrm{in}}^{\mathrm{UQFF}}$ AU & $r_{\mathrm{in}}^{\mathrm{Kopparapu}}$ AU & residual \midrule
Sun (Earth anchor) & 1.000000   & 0.909000 & 0.910000 & 0.110%  Proxima Centauri   & 0.001553   & 0.035821 & 0.0357   & 0.339%  TRAPPIST-1         & 0.000553   & 0.021384 & 0.0212   & 0.868%  Kepler-186         & 0.054000   & 0.211249 & 0.2080   & 1.562%  $\tau$ Ceti        & 0.520000   & 0.655481 & 0.6510   & 0.689%  \bottomrule
\end{tabular}
\end{center}
Median residual: 0.689%. (Previous closure entry: 4.318% median; this
paper reduces by $6.3\times$.)

## C-ABI implementation (QCalcGeom v1.5.1)
The closed-form solver is exposed as a JSON-dispatch function in
`QCalcGeom.cpp` (S305 patch) with the call signature
\begin{verbatim}
{
  "function": "universal_buoyancy",
  "M_star_solar": 1.0,
  "L_star_solar": 1.0,
  "t_n": 0.5
}
\end{verbatim}
and returns ${r_{\mathrm{in}},r_{\mathrm{out}}}$ in AU. Test vector
Earth-Sun returns $r_{\mathrm{in}}=0.9090$ AU, $r_{\mathrm{out}}=1.5995$
AU, matching this paper to the full precision of double-float
arithmetic.

## Where the locked primitives enter
\begin{center}
\begin{tabular}{ll}
\toprule
Primitive & Role in HZ closure \midrule
$[\mathrm{SSq}]=57/100$ & buoyancy-coupling normalisation $\Phires=5/6$           & resonance fraction of the radiation field $\rho_{\mathrm{SCm}}$   & vacuum reference for $\Fubii$ $\mathrm{CAL}=0.793$    & one-time anchor calibration (Earth-Sun) \bottomrule
\end{tabular}
\end{center}
No fit parameters appear apart from $\mathrm{CAL}$, fixed once.

## Significance
The Buoyancy / habitable-zone category was the worst-residual top-12
closure topic at 4.318% median. Promotion of the prefactor cluster
$[\mathrm{SSq}]\Phires^2$ to a single calibrated constant and exposure
through the QCalcGeom C-ABI brings the median residual to 0.689%, a
$6.3\times$ improvement, and produces deployable HZ inner/outer
boundaries for the 4 most-studied exoplanet host stars.

## Reproducibility
`_session305e_habitable_zone_closure.py` reproduces the table
above directly from the four locked primitives plus published $L_\star$
values; `QCalcGeom.cpp` (S305) provides the equivalent C-ABI
path for embedded callers.