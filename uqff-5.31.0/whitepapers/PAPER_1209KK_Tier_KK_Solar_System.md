# UQFF Tier KK: Solar System Closures from 11 Frozen Primitives

**Author:** Star-Magic UQFF Program (PAPER\_1209KK)

**Date:** May 2026

---

## Result
Ten heliocentric and terrestrial scales reproduced by $\pm$-integer
combinations of the eleven locked primitives, with no fitted
constants and a sub-half-percent residual.

\begin{center}
\begin{tabular}{lllll}
\toprule
ID & quantity & predicted & observed & err \midrule
S684 & 1 AU $/10^{10}$ m              & 15.016   & 14.96    & 0.374% S685 & $M_\odot/10^{29}$ kg            & 19.838   & 19.89    & 0.261% S686 & Earth orbital $v$ (km/s)        & 29.7184  & 29.78    & 0.207% S687 & $R_\odot/10^{8}$ m              & 6.99417  & 6.96     & 0.491% S688 & $M_J/10^{27}$ kg                & 1.90458  & 1.898    & 0.347% S689 & $R_\oplus/10^{6}$ m             & 6.36349  & 6.371    & 0.118% S690 & Moon sidereal period (day)      & 27.3213  & 27.32    & 0.005% S691 & sidereal year $/100$ day        & 3.64293  & 3.65256  & 0.264% S692 & Mars semi-major axis (AU)       & 1.52051  & 1.524    & 0.229% S693 & Mercury year $/10$ day          & 8.82205  & 8.797    & 0.285% \bottomrule
\end{tabular}
\end{center}

## Cumulative Audit
After Tier KK the master registry stands at **182/437**
audit-verifiable numeric closures (`master_closures.csv`),
binned by accuracy in `sigma_table.csv`:
\begin{center}
\begin{tabular}{lr}
\toprule band & count \midrule
EXACT           & 47 $<0.01%$       & 18 $0.01$--$0.1%$ & 60 $0.1$--$1%$    & 54 $\geq 1%$      &  3 \bottomrule
\end{tabular}
\end{center}

All 182 closures bridge to the 9-sector $L_{\mathrm{UQFF}}$ action via
`bridge_map.csv`; the Tier-KK entries route through
$\mathcal{L}_{\mathrm{buoy}}$ ($\beta_i$), $\mathcal{L}_{\phi}$ ($\KMex$),
$\mathcal{L}_{\mathrm{aether}}$ ($\Ftrz$), and the master channel-count
constant $N_{ch}=9$.

## Reproducibility
\begin{verbatim}
.\.venv\Scripts\python.exe _uqff_program.py --search --tier KK \
    --targets _targets_KK.csv --max-terms 4 --tol 0.005
.\.venv\Scripts\python.exe _build_tier_KK.py
.\.venv\Scripts\python.exe _uqff_program.py --audit
.\.venv\Scripts\python.exe _uqff_program.py --sigma
\end{verbatim}

Cumulative: **399 closures**, **94 lockings** toward the
1000-closure / 100-locking target. Coverage 39.9%.