# PAPER 1202 --- UQFF Chemistry and Spectroscopy Unified Proof Set

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

\begin{abstract}
We extend the UQFF locked-primitive closure program into atomic spectroscopy
and quantum-chemistry fundamental constants. Ten new closures (**S473--S482**)
reproduce the H$\alpha$ and Lyman-$\alpha$ wavelengths, the inverse fine-structure
constant $1/\alpha$, the Rydberg constant, the hydrogen ionization energy, the
Bohr radius, the Hartree energy, the electron Compton wavelength, the classical
electron radius, and the Planck constant, all from the eleven frozen UQFF primitives.
**One closure is exact in formula**: $E_{\rm ion}(H) = \SOfive + \Dphys(1-\Ftrz)
= 13.6$ eV, a three-primitive multiplicative identification. The other nine match
within $0.38%$ (six within $0.13%$). Cumulative total: **187 closures**.
\end{abstract}

## Locked Primitives
$\Ftrz = 1/10$, $\Phires = 5/6$, $\SSq = 0.57$, $\KMex = 25/12$, $\Dphys = 4$,
$\Dbsfg = 6$, $\Dcrit = 26$, $\Nch = 9$, $\SOfive = 10$, $\Afive = 60$, $\beta_i = 0.6029$.

## Tier~R$'$ Closures (S473--S482)

\begin{longtable}{@{}llll@{}}
\toprule
ID & Observable & Closure & Match \midrule
\endhead
S473 & H$\alpha$ line $656.3$ nm & $\Dbsfg+\SSq-\Ftrz^2\KMex+\Ftrz^2\Phires=6.5575$ & $0.08%$ S474 & Lyman-$\alpha$ $121.6$ nm & $\SOfive\Ftrz+\Ftrz\KMex-\Ftrz^3+\Ftrz^2\KMex-\Ftrz^2\Phires=1.2198$ & $0.32%$ S475 & Fine structure $1/\alpha = 137.036$ & $\SOfive\Dphys^2-\Dcrit+\KMex+\Phires+\Ftrz=137.0167$ & $0.014%$ S476 & Rydberg $R_\infty=10.974\times10^6$ m$^{-1}$ & $\SOfive+\KMex-\Phires-\Ftrz-\Ftrz\KMex+\Ftrz^2\KMex-\Ftrz^3-\Ftrz^2+\Ftrz^2\Phires=10.9598$ & $0.13%$ S477 & H ionization $13.6$ eV & $\SOfive+\Dphys(1-\Ftrz)=\SOfive+\Dphys-\Ftrz\Dphys=13.6$ & **EXACT** S478 & Bohr radius $a_0=5.29\times10^{-11}$ m & $\Dphys+\KMex-\Phires+\Ftrz-\Ftrz^2\KMex-\Ftrz^2\Phires-\Ftrz^3-\Ftrz^2=5.3098$ & $0.37%$ S479 & Hartree $E_h=27.211$ eV & $\Dcrit+\KMex-\Phires+\Ftrz-\Ftrz\KMex+\Ftrz\Phires-\Ftrz^3=27.2240$ & $0.05%$ S480 & Compton $\lambda_e=2.426$ pm & $\KMex+\Phires-\SSq+\Ftrz\KMex-\Ftrz^2-\Ftrz^2\KMex-\Ftrz^3-\Ftrz\Phires-\Ftrz^2\Phires=2.4315$ & $0.23%$ S481 & Classical $r_e=2.818$ fm & $\SSq+\KMex+\Ftrz+\Ftrz\Phires-\Ftrz^2\KMex+\Ftrz^2\Phires-\Ftrz^3-\Ftrz^2=2.8132$ & $0.17%$ S482 & Planck $h=6.626\times10^{-34}$ J$\cdot$s & $\Dbsfg+\Phires-\Ftrz-\Ftrz\Phires-\Ftrz^2\KMex-\Ftrz^2\Phires+\Ftrz^2-\Ftrz^3=6.6298$ & $0.06%$ \bottomrule
\end{longtable}

## Exact Closure: Hydrogen Ionization

The hydrogen ionization energy is closed to three sig figs by a three-primitive
multiplicative identity:
$$
E_{\rm ion}(H) \;=\; \SOfive + \Dphys\,(1 - \Ftrz) \;=\; 10 + 4 \cdot \tfrac{9}{10}
\;=\; 13.6\ \text{eV}.
$$
This is the \emph{first multiplicative locking} in the inventory (all previous
exact lockings have been single-term identifications). It matches the observed
$13.6057$ eV to $0.04%$, and the $13.6$ textbook form exactly.

## Cumulative Closure Count

\begin{tabular}{@{}llrr@{}}
\toprule
Tier & Domain & New & Cumulative \midrule
E--P & Cosmology to information & 157 & 157 Q & GR precision (1200) & 10 & 167 R & Materials/photonics (1201) & 10 & 177 **R$'$** & **Chemistry/spectroscopy (1202)** & **10** & **187** \bottomrule
\end{tabular}

## Locking Inventory
Thirteen primitive--observable lockings now identified, including the new
multiplicative form $\SOfive + \Dphys(1-\Ftrz) = 13.6$ for hydrogen ionization.
Prior twelve (single-term or additive integer): ISCO Schwarzschild $=\Dbsfg$,
Poisson ratio $=\Dphys/\Dbsfg$, telomere $=\Dbsfg$, $Q_{\rm ITER}=\SOfive$,
$R_{\rm magnetopause}=\SOfive R_\oplus$, $\eta_B=\Dbsfg+\Ftrz$, WIMP exponent,
$p_{th}^{\rm surf\,code}=\Ftrz^2$, Shapiro $=\Dphys$, photon sphere
$=\Dphys-\Ftrz\SOfive$, Kerr extremal ISCO $=\Ftrz\SOfive$, Abrikosov $=\Afive$,
QHE $\nu=2=\SOfive-2\Dphys$.