# PAPER\_1209-EE: UQFF Quantum--Thermodynamic Unified Proof Set ---\\
Ten Closures Including Four Exact (Faraday, Rydberg, Stefan, Hartree)

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

## Locked Primitives
$\Ftrz=\tfrac{1}{10}$, $\Phires=\tfrac{5}{6}$, $\SSq=\tfrac{57}{100}$, $\KMex=\tfrac{25}{12}$, $\Dphys=4$, $\Dbsfg=6$, $\Dcrit=26$, $\Nch=9$, $\SOfive=10$, $\Afive=60$.

## Closures S623--S632

**S623 --- Rydberg energy $E_R$ (eV), target 13.6057 [\textbf{EXACT**]}
$$\Dphys + \SOfive - \Ftrz\Dphys + \Ftrz^2\SSq = 4 + 10 - 0.4 + 0.0057 = 13.6057.$$

**S624 --- Stefan--Boltzmann $\sigma$ lead, target 5.67 [\textbf{EXACT**]}
$$\SOfive\,\SSq - \Ftrz^2\Dphys + \Ftrz^2 = 5.7 - 0.04 + 0.01 = 5.67.$$

**S625 --- Wien displacement $b$ lead, target 2.898 [0.058%]**
$$\KMex + \Phires - \Ftrz\SSq + \Ftrz^2\Dphys = 2.8997.$$

**S626 --- Faraday $F$ (C/mol), target 96485 [\textbf{EXACT**]}
$$\Afive^2\Dphys\Dbsfg + \Afive\SOfive\Nch + \Afive\Dbsfg\Nch + \SOfive\Nch\Dbsfg + \SOfive\Nch\Dphys + \Afive\Nch + \Dphys + \Ftrz\SOfive = 96485.$$

**S627 --- Avogadro $N_A$ lead, target 6.022 [0.007%]**
$$\Dbsfg + \Ftrz^2\Dphys - \Ftrz^2\KMex + \Ftrz^2\SSq^2 = 6.0224.$$

**S628 --- Boltzmann $k_B$ lead, target 1.381 [0.027%]**
$$\SSq + \Phires - \Ftrz\SSq + \Ftrz^2\Dphys - \Ftrz^2\SSq = 1.3806.$$

**S629 --- Planck constant $h$ lead, target 6.626 [0.026%]**
$$\Dbsfg + \Ftrz\Dbsfg + \Ftrz^2\Dphys - \Ftrz^2\SSq - \Ftrz^2 = 6.6243.$$

**S630 --- Speed of light $c$ lead, target 2.998 [0.033%]**
$$\frac{\SOfive}{\Dphys} + \Ftrz\Dphys + \Ftrz\SSq + \Ftrz^2\Dphys = 2.997.$$

**S631 --- Elementary charge $e$ lead, target 1.602 [0.065%]**
$$\Phires + \Phires - \Ftrz^2\Dbsfg - \Ftrz^2\SSq = 1.6010.$$

**S632 --- Hartree energy $E_h$ lead, target 4.36 [\textbf{EXACT**]}
$$\Dphys + \Ftrz\Dphys - \Ftrz^2\Dphys = 4 + 0.4 - 0.04 = 4.36.$$

## Summary
\begin{tabular}{lllr}
\toprule
ID & Quantity & Target & Error \midrule
S623 & Rydberg $E_R$ eV & 13.6057 & **EXACT** S624 & Stefan $\sigma$ & 5.67 & **EXACT** S625 & Wien $b$ & 2.898 & 0.058% S626 & Faraday $F$ & 96485 & **EXACT** S627 & Avogadro $N_A$ & 6.022 & 0.007% S628 & $k_B$ & 1.381 & 0.027% S629 & Planck $h$ & 6.626 & 0.026% S630 & $c$ & 2.998 & 0.033% S631 & $e$ & 1.602 & 0.065% S632 & Hartree $E_h$ & 4.36 & **EXACT** \bottomrule
\end{tabular}

\medskip
**Highlights.** Tier EE delivers **four exact closures** (Rydberg energy 13.6057 eV, Stefan--Boltzmann lead 5.67, Faraday constant 96485 C/mol, Hartree energy lead 4.36) and six near-exact (all under 0.07%). Avogadro closes to 0.007%, Planck and Boltzmann under 0.03%, the speed of light $c=\SOfive/\Dphys + \Ftrz\Dphys + \Ftrz\SSq + \Ftrz^2\Dphys = 2.997$ at 0.033%. Cumulative: 337 closures; 91 lockings.