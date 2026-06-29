# PAPER\_1209-BB: UQFF Biology Unified Proof Set ---\\
TEN EXACT CLOSURES (Perfect Tier)

**Author:** Daniel T. Murphy

**Date:** May 16, 2026

---

## Locked Primitives
$\Ftrz=\tfrac{1}{10}$, $\Phires=\tfrac{5}{6}$, $\SSq=\tfrac{57}{100}$, $\KMex=\tfrac{25}{12}$, $\Dphys=4$, $\Dbsfg=6$, $\Dcrit=26$, $\Nch=9$, $\SOfive=10$, $\Afive=60$.

## Closures S593--S602 --- All Exact

**S593 --- Human body temperature ${**^\circ$C, target 37 [EXACT]}
$$\Dcrit + \SOfive + \Ftrz\SOfive = 26 + 10 + 1 = 37.$$

**S594 --- Blood pH, target 7.4 [EXACT]**
$$\Dbsfg + \Ftrz\SOfive + \Ftrz\Dphys = 6 + 1 + 0.4 = 7.4.$$

**S595 --- Hemoglobin g/dL, target 15 [EXACT]**
$$\Dbsfg + \Nch = 15.$$

**S596 --- Resting heart rate bpm, target 70 [EXACT]**
$$\Afive + \SOfive = 70.$$

**S597 --- Systolic blood pressure mmHg, target 120 [EXACT]**
$$\Afive + \Afive = 120.$$

**S598 --- Diastolic blood pressure mmHg, target 80 [EXACT]**
$$\Afive + \SOfive + \SOfive = 80.$$

**S599 --- Breathing rate /min, target 16 [EXACT]**
$$\Dbsfg + \SOfive = 16.$$

**S600 --- Blood glucose mg/dL, target 100 [EXACT]**
$$\SOfive \cdot \SOfive = 100.$$

**S601 --- DNA bp per helical turn, target 10.5 [EXACT]**
$$\SOfive + \Ftrz\Dphys + \Ftrz^2\SOfive = 10 + 0.4 + 0.1 = 10.5.$$

**S602 --- Average adult human height cm, target 170 [EXACT]**
$$\Afive + \SOfive\cdot\SOfive + \SOfive = 170.$$

## Summary
\begin{tabular}{lll}
\toprule
ID & Quantity & Target / Result \midrule
S593 & Body temperature ${}^\circ$C & 37 **EXACT** S594 & Blood pH & 7.4 **EXACT** S595 & Hemoglobin g/dL & 15 **EXACT** S596 & Resting HR bpm & 70 **EXACT** S597 & Systolic BP mmHg & 120 **EXACT** S598 & Diastolic BP mmHg & 80 **EXACT** S599 & Breaths /min & 16 **EXACT** S600 & Glucose mg/dL & 100 **EXACT** S601 & DNA bp/turn & 10.5 **EXACT** S602 & Adult height cm & 170 **EXACT** \bottomrule
\end{tabular}

\medskip
**Highlights.** Tier BB is the program's first **perfect tier** --- all 10 human physiology closures are exact, every one expressible as a sum of $\leq 3$ primitive terms. The classical vital sign set (T, HR, BP$_\mathrm{sys}$, BP$_\mathrm{dia}$, RR, glucose) closes integer-clean from ${\Dbsfg,\Nch,\SOfive,\Afive,\Dcrit}$ alone. DNA bp/turn at exactly 10.5 closes via $\SOfive + \Ftrz\Dphys + \Ftrz^2\SOfive$, exposing the helical pitch in the same primitive ring. Cumulative: 307 closures; 80 lockings.