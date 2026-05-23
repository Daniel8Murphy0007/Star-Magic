# PAPER\_1209-X: UQFF Climate \& Atmosphere Unified Proof Set ---\\
Five Exact Closures Including CO$_2$, Salinity, Albedo, Greenhouse Effect, Ozone

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

## Locked Primitives
$\Ftrz=\tfrac{1}{10}$, $\Phires=\tfrac{5}{6}$, $\SSq=\tfrac{57}{100}$, $\KMex=\tfrac{25}{12}$, $\Dphys=4$, $\Dbsfg=6$, $\Dcrit=26$, $\Nch=9$, $\SOfive=10$, $\Afive=60$, $\beta_i=0.6029$.

## Closures S553--S562

**S553 --- CO$_2$ atmospheric concentration (ppm), target 420 [EXACT]**
$$\Afive\Dphys + \Dcrit\Dbsfg + \Dbsfg\Dphys = 240+156+24 = 420.$$

**S554 --- Adiabatic lapse rate (K/km), target 6.5 [0.21%]**
$$\Dbsfg + \SSq - \Ftrz\Phires = 6.4867.$$

**S555 --- Atmospheric scale height (km), target 8.5 [0.71%]**
$$2\Dphys + \SSq - \Ftrz^2 = 8.56.$$

**S556 --- Solar constant (W/m$^2$), target 1361 [0.0022%]**
$$\Afive^2\Ftrz\Dphys - \Nch\SOfive + \SOfive + \SSq + \Ftrz\Dphys = 1360.97.$$

**S557 --- Ocean salinity (ppt), target 35 [EXACT]**
$$\Dcrit + \Nch = 35.$$

**S558 --- Standard atmospheric pressure (kPa), target 101.325 [0.0049%]**
$$\SOfive^2 + \SSq + \Phires - \Ftrz\Phires = 101.32.$$

**S559 --- Earth Bond albedo, target 0.30 [EXACT]**
$$3\Ftrz = 0.30.$$

**S560 --- Greenhouse effect (K), target 33 [EXACT]**
$$\Dcrit + \Nch - \KMex + \Ftrz\Phires = 35 - \tfrac{25}{12} + \tfrac{1}{12} = 33.$$

**S561 --- Tropopause height (km), target 11 [0.38%]**
$$\SOfive + \Phires + \KMex\Ftrz = 11.042.$$

**S562 --- Total ozone column (DU), target 300 [EXACT]**
$$\Afive\Dphys + \SOfive\Dbsfg = 240 + 60 = 300.$$

## Summary
\begin{tabular}{lllr}
\toprule
ID & Quantity & Target & Error \midrule
S553 & CO$_2$ ppm & 420 & **EXACT** S554 & Lapse rate K/km & 6.5 & 0.21% S555 & Scale height km & 8.5 & 0.71% S556 & Solar const W/m$^2$ & 1361 & 0.0022% S557 & Salinity ppt & 35 & **EXACT** S558 & Pressure kPa & 101.325 & 0.0049% S559 & Bond albedo & 0.30 & **EXACT** S560 & Greenhouse K & 33 & **EXACT** S561 & Tropopause km & 11 & 0.38% S562 & Ozone DU & 300 & **EXACT** \bottomrule
\end{tabular}

\medskip
**Highlights.** Tier X delivers **five exact closures** (S553, S557, S559, S560, S562) --- the most exact tier in the program. The greenhouse effect $33 = \Dcrit+\Nch-\KMex+\Ftrz\Phires$ exhibits perfect arithmetic ($35 - 24/12 = 33$). Bond albedo as $3\Ftrz$ and salinity as $\Dcrit+\Nch$ are the cleanest single-term closures yet recorded. Solar constant matches to 22 ppm, atmospheric pressure to 49 ppm. Cumulative: 267 closures; 51 lockings (5 new from exact climate closures).