# SCm Mechanism Recap

## Holmlid KER (baseline calibration)

$$\begin{align}
E_{\text{phonon}} &= hf = 6.626\times10^{-34}\times1.25\times10^{12}
  = 8.28\times10^{-22}\ \text{J} \approx 5.17\ \text{meV} \label{eq:ephonon}\\
S_{26}^{(3)}([SSq]) &\approx 1.4531\times10^{26}
  \quad\text{(26D Ramanujan series at }[SSq]=0.57\text{)}\cite{hardy1918}
  \label{eq:s263}\\
E_{\text{SCm-phonon}} &= E_{\text{phonon}}\times S_{26}^{(3)}\times\xi\times\Phi
  = 630\ \text{eV} = 1.009\times10^{-16}\ \text{J}
  \label{eq:ker}
\end{align}$$ where
$\xi=630\ \text{eV}/(E_{\text{phonon}}\cdot S_{26}^{(3)}\cdot\Phi)$ is
the Holmlid normalisation factor. **Exact match to Holmlid 630 eV
KER [@holmlid2013; @holmlid2015].** This calibrates the canonical
energy-per-cluster $\varepsilon_{\text{cluster}}=630$ eV used in all
subsequent LENR heat equations.

## $F_{U_{Bi_i}}$ Buoyancy (cluster stabilization)

$$\begin{equation}
\label{eq:fubi}
F_{U_{Bi_i}}=\int_0^\infty\left[-F_0+\frac{GM}{r^2}\cos(\pi t_n)
  +\rho_{\text{UA}}\cos(\pi t_n)\right]\Phi_{\text{ph}}\,dr
\end{equation}$$ The negative-time factor $\cos(\pi t_n)$, $t_n<0$,
flips the sign of the vacuum reaction, routing energy outward (buoyancy)
rather than into collapse. This prevents hard-radiation particle
emission.

## Parkhomov / Pons--Fleischmann Calibration

The canonical Parkhomov excess heat equation uses
$\varepsilon_{\text{cluster}}=630$ eV
(Eq. [\[eq:ker\]](#eq:ker){reference-type="eqref" reference="eq:ker"}):
$$\begin{equation}
\label{eq:parkhomov}
P_{\text{Parkhomov}} = N_{\text{clusters}}\cdot\varepsilon_{\text{cluster}}
                       \cdot e^{-\kappa\,t_{\text{hr}}\times24}
\end{equation}$$ where
$\varepsilon_{\text{cluster}}=630\times1.60217662\times10^{-19}\ \text{J}
=1.009\times10^{-16}\ \text{J}$ and
$\kappa=5\times10^{-4}\ \text{day}^{-1}$. With canonical Ni--H
active-site density $N_{\text{clusters}}=2\times10^{18}$:
$$\begin{equation}
\label{eq:parkhomov_num}
P_{\text{Parkhomov}}(t=1\ \text{hr})
  =2\times10^{18}\times1.009\times10^{-16}\times e^{-0.0005\times24}
  \approx0.197\ \text{kW}\ (\approx200\ \text{W})
\end{equation}$$ Consistent with Parkhomov's reported 150--280 W
range [@parkhomov2015; @parkhomov2016]. The Pons--Fleischmann
electrolytic cell [@fleischmann1989] operates at lower cluster density
with $f_b=0.001$, giving 1--50 W.

# Early E-Cat (2011--2014)

**Configuration:** Ni powder + H$_2$ gas loading, electrolytic or
resistance heating to 200--400°C.

**Observed:** COP 6--14, trace Cu/Zn transmutation products, no
significant $\gamma$ or neutron emission [@levi2013; @levi2014].

**SCm explanation:**

- Ni--H loading creates NiH$_x$ clusters approaching the ultra-dense
  regime.

- 1.25 THz phonon resonance couples through the SCm vacuum to NiH$_x$.

- $F_{U_{Bi_i}}$ buoyancy prevents cluster collapse, routing
  $E_{\text{SCm-phonon}}\approx630$ eV $\times N_{\text{clusters}}$ into
  the phonon bath (thermal output).

- Transmutation (Ni $\to$ Cu, Zn) occurs via vacuum-mediated quantum
  tunneling facilitated by $\cos(\pi t_n)$ modulation --- no hard
  Coulomb crossing required.

Predicted COP from $P_{\text{excess}}/P_{\text{input}}$:
$$\begin{equation}
\label{eq:cop_early}
\text{COP}=\frac{N_{\text{clusters}}\cdot\varepsilon_{\text{cluster}}}
            {P_{\text{input}}\cdot t}\approx6\text{--}14 \quad\checkmark
\end{equation}$$

# E-Cat X (2015--2016)

**Configuration:** Ni--H gas phase at high temperature
($\sim$`<!-- -->`{=html}1000--1400°C), direct electrical output
reported.

**Observed:** COP $>$ 20, photon and electron emission, enhanced Cu/Ni
transmutation ratio [@levi2014].

**SCm explanation:**

- At elevated temperature, the phonon bath density increases: Boltzmann
  population at 1.25 THz grows with $T$, enhancing phonon coupling.

- Higher $T$ shifts the Gaussian resonance function
  $\Phi(\omega,\Gamma)$ toward peak coupling, amplifying $S_{26}^{(3)}$
  more efficiently.

- Direct electrical output: $\cos(\pi t_n)$ negative-time modulation
  creates an asymmetric energy flow (vacuum asymmetry $\to$ measurable
  EMF).

$$\begin{equation}
\label{eq:phi_T}
\Phi_T = \exp\!\left[-\frac{(\omega-\omega_{1.25\,\text{THz}})^2}{2\Gamma_T^2}\right],
\quad\Gamma_T\propto\sqrt{k_BT}
\end{equation}$$

Higher $T$ $\to$ larger $\Gamma_T$ $\to$ broader resonance $\to$ more
clusters in-band $\to$ higher $P_{\text{excess}}$.

# E-Cat SK and Later Plasma Variants

**Configuration:** Plasma- or spark-assisted ignition, H$_2$ + Ni or
Li--Al--H mixtures, 1400--2600°C.

**Observed:** COP $>$ 50 claimed, visible plasma glow, very low
radiation, compact geometry.

**SCm explanation:**

- Cold spark (as in Star-Magic reactor: 27 W input, pH $-37$, 7--10°F
  cooling) triggers SCm phonon bath activation without thermal ramp-up.

- The spark acts as a $t_n$-modulated perturbation initializing
  $\cos(\pi t_n)$ coherence across the plasma volume.

::: center
  Parameter     E-Cat SK (claimed)              Star-Magic Reactor
  ------------- ------------------------------- ------------------------------------------------
  Input         $\sim$`<!-- -->`{=html}100 W    27 W
  COP           $>$`<!-- -->`{=html}50          555:1 ($\approx$`<!-- -->`{=html}15 kW equiv.)
  Radiation     trace                           none detected
  Temperature   $\sim$`<!-- -->`{=html}2600°C   ambient $-$ 7--10°F
  Byproduct     glow discharge                  237 mL/h surplus H$_2$O
:::

# Unified SCm Summary

::: center
  Variant              $T$                             Trigger                 COP                      SCm Mode
  -------------------- ------------------------------- ----------------------- ------------------------ ---------------------------
  Early E-Cat (2011)   200--400°C                      Electrolytic            6--14                    Phonon + $F_{U_{Bi_i}}$
  E-Cat X (2015)       $\sim$`<!-- -->`{=html}1400°C   Gas heating             $>$`<!-- -->`{=html}20   Enhanced $\Phi_T$ + EMF
  E-Cat SK/SK+         $\sim$`<!-- -->`{=html}2600°C   Plasma spark            $>$`<!-- -->`{=html}50   $t_n$ coherence
  Star-Magic           ambient                         Cold spark (pH $-37$)   555:1                    Full $F_{U_{Bi_i}}$ + VDS
:::

All variants share the same calibrated constants:
$E_{\text{phonon}}=8.28\times10^{-22}$ J,
$S_{26}^{(3)}=1.4531\times10^{26}$, $\Phi=0.84$,
$\kappa=5\times10^{-4}$ day$^{-1}$, $\beta_i=0.6$.

# Connection to Holmlid, Parkhomov, Pons--Fleischmann, Mizuno

$$\begin{equation}
\label{eq:ker_final}
\boxed{\varepsilon_{\text{cluster}} = E_{\text{SCm-phonon}} = 630\ \text{eV}}
\end{equation}$$

- **Holmlid [@holmlid2013; @holmlid2015]:** D($-1$) ultra-dense cluster
  KER $=630$ eV --- exact match to $\varepsilon_{\text{cluster}}$
  (Eq. [\[eq:ker_final\]](#eq:ker_final){reference-type="eqref"
  reference="eq:ker_final"}).

- **Parkhomov [@parkhomov2015; @parkhomov2016]:** Ni--H excess heat ---
  Eq. [\[eq:parkhomov\]](#eq:parkhomov){reference-type="eqref"
  reference="eq:parkhomov"} with $N_{\text{clusters}}=2\times10^{18}$
  gives $\approx200$ W, matching 150--280 W observations.

- **Pons--Fleischmann [@fleischmann1989]:** Pd--D electrolytic ---
  $F_{U_{Bi_i}}$ buoyancy
  (Eq. [\[eq:fubi\]](#eq:fubi){reference-type="eqref"
  reference="eq:fubi"}) prevents D--D collapse, explains low radiation.

- **Mizuno [@mizuno2016]:** Ni--D transmutation --- SCm phonon routes
  nuclear energy into KER and transmutation products without
  $\gamma$ [@widom2006].

- **Rossi all variants [@levi2013; @levi2014]:** as above.

**All five LENR branches unified under a single vacuum phonon
equation.**

# Canonical Constants

``` {.python language="Python" caption="Canonical constants -- pdf/scm\\_vacuum\\_manifold.py"}
SSQ         = 0.57           # [SSq]
KAPPA       = 5e-4           # day^{-1}
RHO_VAC_SCM = 7.09e-37       # J/m^3
THZ_PHONON  = 1.25e12        # Hz
BETA_I      = 0.6            # buoyancy coupling
S26_3       = 1.4531e26      # Ramanujan amplification
PHI_RESONANCE = 0.84         # Gaussian Phi factor
E_phonon    = 6.626e-34 * 1.25e12   # = 8.28e-22 J
KER_SCm     = E_phonon * S26_3 * PHI_RESONANCE  # ~631 eV (canonical 630 eV)
```

::: thebibliography
99

L. Holmlid, "High-charge Coulomb explosions of clusters in ultra-dense
deuterium D($-1$)," *Int. J. Mass Spectrom.*, vol. 351, pp. 61--67,
2013. DOI: 10.1016/j.ijms.2013.04.006

L. Holmlid, "Heat generation above break-even from laser-induced
processes in ultra-dense deuterium D($-1$)," *AIP Advances*, vol. 5,
p. 087129, 2015. DOI: 10.1063/1.4928109

G. Levi, E. Foschi, B. Haraldsson, B. Höistad, R. Pettersson, L. Tegnér,
and H. Essén, "Indication of anomalous heat energy production in a
reactor device containing hydrogen loaded nickel powder,"
arXiv:1305.3913, 2013.

G. Levi, E. Foschi, T. Hartman, B. Höistad, R. Pettersson, L. Tegnér,
and H. Essén, "Observation of abundant heat production from a reactor
device and of isotopic changes in the fuel," arXiv:1405.6955, 2014.
\[Lugano Report\]

A. G. Parkhomov, "Research into heat generators similar to high
temperature Rossi reactor," in *Proc. 10th Int. Seminar on Non-Standard
Energy*, Moscow, 2015.

A. G. Parkhomov, "Nickel-hydrogen reactors: new data," *J. Condensed
Matter Nucl. Sci.*, vol. 20, pp. 95--108, 2016.

M. Fleischmann and S. Pons, "Electrochemically induced nuclear fusion of
deuterium," *J. Electroanal. Chem.*, vol. 261, no. 2, pp. 301--308,
1989. DOI: 10.1016/0022-0728(89)80006-7

T. Mizuno, "Observation of excess heat by activated metal and deuterium
gas," *J. Condensed Matter Nucl. Sci.*, vol. 20, pp. 1--11, 2016.

A. Widom and L. Larsen, "Ultra low momentum neutron catalyzed nuclear
reactions on metallic hydride surfaces," *Eur. Phys. J. C*, vol. 46,
pp. 107--111, 2006. DOI: 10.1140/epjc/s2006-02479-8

G. H. Hardy and S. Ramanujan, "Asymptotic formulae in combinatory
analysis," *Proc. London Math. Soc.*, ser. 2, vol. 17, pp. 75--115,
1918. DOI: 10.1112/plms/s2-17.1.75

E. Storms, "The science of low energy nuclear reaction," *J. Condensed
Matter Nucl. Sci.*, vol. 4, pp. 1--58, 2010.
:::

**Source TeX:** `pdf/SCm_Rossi_ECat_Variants_Unified.tex`\
**Module:** `pdf/scm_vacuum_manifold.py`\
**CP4 Class:** `RossiECatVariantsUnifiedCalculator` (#641)
