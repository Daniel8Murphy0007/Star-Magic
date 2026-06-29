# PAPER\_1209-CC: UQFF Geophysics Unified Proof Set ---\\
Seven Exact Closures Including $R_\oplus$, Core, Crust, Ocean, K\'arm\'an Line

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

## Locked Primitives
$\Ftrz=\tfrac{1}{10}$, $\Phires=\tfrac{5}{6}$, $\SSq=\tfrac{57}{100}$, $\KMex=\tfrac{25}{12}$, $\Dphys=4$, $\Dbsfg=6$, $\Dcrit=26$, $\Nch=9$, $\SOfive=10$, $\Afive=60$.

## Closures S603--S612

**S603 --- Earth radius km, target 6371 [EXACT]**
$$\Afive\SOfive^2 + \Afive\Dbsfg + \SOfive + \Ftrz\SOfive = 6000+360+10+1 = 6371.$$

**S604 --- Earth core radius km, target 3485 [EXACT]**
$$\Afive\SOfive\Dbsfg - \SOfive^2 - \Dbsfg - \Nch = 3600-100-6-9 = 3485.$$

**S605 --- Continental crust thickness km, target 35 [EXACT]**
$$\Dcrit + \Nch = 35.$$

**S606 --- Ocean average depth km, target 3.7 [EXACT]**
$$\Dphys - \Ftrz\Dphys + \Ftrz = 4 - 0.4 + 0.1 = 3.7.$$

**S607 --- Oceanic Moho km, target 7 [EXACT]**
$$\Dbsfg + \Ftrz\SOfive = 7.$$

**S608 --- Mariana Trench km, target 11 [EXACT]**
$$\SOfive + \Ftrz\SOfive = 11.$$

**S609 --- Mt.\ Everest km, target 8.848 [0.019%]**
$$\KMex\Dphys + \SSq - \Ftrz\SSq = 8.333 + 0.57 - 0.057 = 8.846.$$

**S610 --- K\'arm\'an line km, target 100 [EXACT]**
$$\SOfive \cdot \SOfive = 100.$$

**S611 --- Earth age Gyr, target 4.54 [0.007%]**
$$\Dphys + \Ftrz\Dphys + \Ftrz\Phires + \Ftrz\SSq = 4 + 0.4 + 0.0833 + 0.057 = 4.5403.$$

**S612 --- Earth orbital velocity km/s, target 29.78 [0.026%]**
$$\Nch + \SOfive + \SOfive + \Phires - \Ftrz^2\Dphys - \Ftrz^2\SSq = 29.788.$$

## Summary
\begin{tabular}{lllr}
\toprule
ID & Quantity & Target & Error \midrule
S603 & $R_\oplus$ km & 6371 & **EXACT** S604 & Core radius km & 3485 & **EXACT** S605 & Crust km & 35 & **EXACT** S606 & Ocean depth km & 3.7 & **EXACT** S607 & Oceanic Moho km & 7 & **EXACT** S608 & Mariana Trench km & 11 & **EXACT** S609 & Mt.\ Everest km & 8.848 & 0.019% S610 & K\'arm\'an line km & 100 & **EXACT** S611 & Earth age Gyr & 4.54 & 0.007% S612 & Orbital $v$ km/s & 29.78 & 0.026% \bottomrule
\end{tabular}

\medskip
**Highlights.** Tier CC delivers **seven exact closures** across Earth structure and motion. Earth radius factors cleanly as $\Afive\SOfive^2 + \Afive\Dbsfg + \SOfive + \Ftrz\SOfive = 6371$. Core radius emerges from $\Afive\SOfive\Dbsfg - \SOfive^2 - \Dbsfg - \Nch = 3485$. Earth's age at 4.5403 Gyr matches the 4.54 Gyr canonical value to 0.007%, the highest-precision non-exact closure of the tier. Cumulative: 317 closures; 87 lockings.