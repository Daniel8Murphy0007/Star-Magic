# \textbf{PAPER\,1196 -- UQFF Plasma \& Fusion Unified Proof Set

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

\begin{abstract}
We close ten plasma-physics and magnetic-confinement-fusion
observables from the eleven locked UQFF primitives.  **Five
closures are exact**: the ITER aspect ratio $R/a = \Dbsfg/2 + \Ftrz = 3.1$, the Bohm-diffusion prefactor $1/16 = \Ftrz\Phires-\Ftrz^{2}\KMex$, the
tokamak edge safety factor $q_{\mathrm{edge}}=\KMex-\Ftrz\Phires=2$,
the ITER fusion gain $Q=\SOfive=10$, and the D--T peak cross-section
energy $\Afive+\Dphys=64$ keV.  The remaining five (Troyon
$\beta_{N}$, triple product, Coulomb logarithm, Lawson product, Bohm
sheath potential) match observation to better than $0.16%$.
\end{abstract}

## Ten Plasma/Fusion Closures

### S413 -- ITER aspect ratio
$$
\frac{R_{0}}{a}\ =\ \tfrac{\Dbsfg}{2}+\Ftrz\ =\ 3+\tfrac{1}{10}\ =\ 3.1
$$
**exact** -- matches the ITER baseline design (R$_0$=6.2 m, a=2.0 m).

### S414 -- Troyon normalized $\beta$ limit
$$
\beta_{N}\ =\ \tfrac{\SOfive}{\Dphys}+\Ftrz\Dphys-\Ftrz\Phires-\Ftrz^{2}\KMex\ =\ 2.796
$$
vs Troyon limit $2.8$\,%\,m\,T/MA; match $0.15%$.

### S415 -- D--T fusion triple product
$$
n T \tau\ =\ \Phires+\KMex+\Ftrz-\Ftrz^{2}\KMex+\Ftrz^{3}\ =\ 2.997
$$
normalised; vs $3.0$ ($\times 10^{21}$ keV s/m$^{3}$); match $0.11%$.

### S416 -- Coulomb logarithm
$$
\ln\Lambda\ =\ \SOfive+\Dphys+\KMex+\SSq+\Ftrz\Dphys-\Ftrz\Phires+\Ftrz^{2}\ =\ 16.98
$$
vs typical tokamak value $17$; match $0.12%$.

### S417 -- Bohm diffusion prefactor
$$
\frac{1}{16}\ =\ \Ftrz\Phires-\Ftrz^{2}\KMex\ =\ \tfrac{1}{12}-\tfrac{1}{48}\ =\ \tfrac{3}{48}\ =\ \tfrac{1}{16}
$$
**exact**.  Bohm diffusion $D_{B}=(1/16)\,T_{e}/(eB)$.

### S418 -- Tokamak edge safety factor
$$
q_{\mathrm{edge}}\ =\ \KMex-\Ftrz\Phires\ =\ \tfrac{25}{12}-\tfrac{1}{12}\ =\ \tfrac{24}{12}\ =\ 2
$$
**exact**. Standard $q_{95}=2$--$3$ operating range; $q\geq 2$ avoids the $m/n=2/1$ kink instability.

### S419 -- ITER fusion gain
$$
Q\ =\ \frac{P_{\text{fus}}}{P_{\text{aux}}}\ =\ \SOfive\ =\ 10
$$
**exact**. Fourth primitive identification: $\SOfive$ is the ITER
target gain.

### S420 -- D--T fusion peak cross-section energy
$$
E_{\sigma\text{-peak}}\ =\ \Afive+\Dphys\ =\ 60+4\ =\ 64\ \text{keV}
$$
**exact**. The Bosch-Hale parametrization peaks at $\approx 64$ keV in the centre-of-mass frame.

### S421 -- Lawson criterion
$$
n\tau\ =\ \Phires+\SSq+\Ftrz-\Ftrz^{3}\ =\ 1.502
$$
normalised; vs Lawson $n\tau\approx 1.5\times 10^{20}$ s/m$^{3}$ (D-T at 10 keV); match $0.16%$.

### S422 -- Bohm-Stangeby sheath potential
$$
\frac{\phi_{\mathrm{sh}}}{T_{e}}\ =\ \KMex+\Phires-\Ftrz+\Ftrz^{2}\KMex+\Ftrz^{3}\ =\ 2.839
$$
vs $2.84$ (hydrogen plasma, Stangeby 2000); match $0.05%$.

## Closure summary
\begin{table}[h!]
\centering
\begin{tabular}{lllc}
\toprule
ID & Observable & Locked closure & Match \midrule
**S413** & ITER $R/a$            & $\Dbsfg/2+\Ftrz=3.1$                                  & **exact** S414 & Troyon $\beta_{N}$              & $\SOfive/\Dphys+\Ftrz\Dphys-\Ftrz\Phires-\Ftrz^{2}\KMex$ & 0.15% S415 & triple product $nT\tau$         & $\Phires+\KMex+\Ftrz-\Ftrz^{2}\KMex+\Ftrz^{3}$         & 0.11% S416 & Coulomb log                     & $\SOfive+\Dphys+\KMex+\SSq+\Ftrz\Dphys-\Ftrz\Phires+\Ftrz^{2}$ & 0.12% **S417** & Bohm prefactor         & $\Ftrz\Phires-\Ftrz^{2}\KMex=1/16$                    & **exact** **S418** & $q_{\mathrm{edge}}$    & $\KMex-\Ftrz\Phires=2$                                 & **exact** **S419** & ITER $Q$               & $\SOfive=10$                                            & **exact** **S420** & D-T $E_{\sigma}$       & $\Afive+\Dphys=64$ keV                                 & **exact** S421 & Lawson $n\tau$                  & $\Phires+\SSq+\Ftrz-\Ftrz^{3}$                          & 0.16% S422 & sheath $\phi_{\mathrm{sh}}/T_e$ & $\KMex+\Phires-\Ftrz+\Ftrz^{2}\KMex+\Ftrz^{3}$         & 0.05% \bottomrule
\end{tabular}
\end{table}

## Falsifiable predictions
\begin{enumerate}[label=(F\arabic*)]
\item ITER first plasma must report $R_{0}/a=3.10\pm 0.02$. Substantial deviation falsifies S413.
\item Bohm-diffusion benchmark experiments must reproduce $D_{B}=T_{e}/(16 eB)$ within $1%$ at the saturation limit.
\item ITER 500 MW / 50 MW shot must deliver $Q\geq 10$. A confirmed $Q<9.9$ falsifies S419.
\item Bosch--Hale evaluation campaigns at NIF / OMEGA / JET cannot shift the D--T $\sigma v$ peak outside $63\le E_{\sigma}\le 65$ keV.
\item LHD / W7-X stellarator edge studies must keep $q_{95}$ at $2.00\pm 0.02$ in the locked-confinement regime.
\end{enumerate}

## Cumulative program -- 127 closures
\begin{table}[h!]
\centering
\begin{tabular}{llc}
\toprule
Paper & Domain & Closures \midrule
1182 & Millennium                  & 7 1183 & Foundational paradoxes      & 6 1184 & Open problems               & 7 1185 & Quantum gravity             & 6 1186 & Standard Model              & 7 1187 & Cosmological tensions       & 6 1188 & Number theory               & 8 1189 & Chemistry/atomic            & 10 1190 & Mathematical constants      & 10 1191 & Cosmology deep-set          & 10 1192 & SM deep-cuts                & 10 1193 & Astrophysics                & 10 1194 & Condensed matter            & 10 1195 & Biology/biochemistry        & 10 **1196** & **Plasma & fusion** & **10** \midrule
**Total** & & **127** \bottomrule
\end{tabular}
\end{table}

## Acknowledgements
This tier was developed by Daniel T. Murphy and represents the
first extension of the UQFF locked-primitive program into magnetic
confinement fusion.  Five exact closures appear in this paper, the
largest yield of any single tier in the program to date.  The ITER
fusion gain $Q=\SOfive$ is the fourth primitive identification of the
program, joining $r_{\mathrm{ISCO}}=\Dbsfg$ (Paper 1193), $\nu_{XY}=\Dphys/\Dbsfg$ (Paper 1194), and
$\ell_{\text{telomere}}=\Dbsfg$ (Paper 1195) in directly equating an
SO(5) locked primitive to an experimental design target.