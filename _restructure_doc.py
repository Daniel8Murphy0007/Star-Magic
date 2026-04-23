"""
Star-Magic.txt complete restructuring in physical order of operation.
THE BELLY BUTTON IS FIRST. GM/r^2 IS LAST. This is the law.
"""

import re

SEP = '________________________________________'

with open('Star-Magic.txt', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Strip trailing newlines for processing
raw = [l.rstrip('\n') for l in lines]

def find_section(raw, title_text):
    """Return (sep_index, title_index, end_sep_index) - the three lines of a section header."""
    for i, line in enumerate(raw):
        if line.strip() == title_text:
            # Find the separator before it
            sep_before = i - 1
            while sep_before >= 0 and raw[sep_before].strip() == '':
                sep_before -= 1
            sep_after = i + 1
            while sep_after < len(raw) and raw[sep_after].strip() == '':
                sep_after += 1
            return sep_before, i, sep_after
    return None

def extract_between(raw, start_title, end_title=None):
    """
    Extract lines starting from the separator before start_title
    up to (but not including) the separator before end_title.
    Returns list of lines as a single string block.
    """
    r = find_section(raw, start_title)
    if r is None:
        return ''
    sep_start = r[0]
    
    if end_title:
        r2 = find_section(raw, end_title)
        if r2 is None:
            return '\n'.join(raw[sep_start:])
        sep_end = r2[0]
        # Back up past blank lines before sep_end
        while sep_end > sep_start and raw[sep_end-1].strip() == '':
            sep_end -= 1
        return '\n'.join(raw[sep_start:sep_end])
    else:
        return '\n'.join(raw[sep_start:])


# ─────────────────────────────────────────────────────────────────────────────
# Extract all sections by title
# ─────────────────────────────────────────────────────────────────────────────

# Title block - first 6 non-blank lines
title_block = '\n'.join(raw[0:6])

# Opening sections before chapters
belly_button_open = extract_between(raw,
    'THE BIGBANG: THE PRIMORDIAL DPM FIRES',
    'THE CANONICAL QUANTUM CHAIN (IMMUTABLE - READ BEFORE ALL ELSE)')

canonical_chain = extract_between(raw,
    'THE CANONICAL QUANTUM CHAIN (IMMUTABLE - READ BEFORE ALL ELSE)',
    'THE DERIVATIVE FORM OF THE QUANTUM CHAIN')

derivative_form = extract_between(raw,
    'THE DERIVATIVE FORM OF THE QUANTUM CHAIN',
    'PROTO-HYDROGEN AND PROTO-HELIUM EVOLUTION')

proto_hhe = extract_between(raw,
    'PROTO-HYDROGEN AND PROTO-HELIUM EVOLUTION',
    'WHAT THIS DOCUMENT CONTAINS:')

what_contains = extract_between(raw,
    'WHAT THIS DOCUMENT CONTAINS:',
    'CANONICAL ONTOLOGY LOCK (v1)')

ontology_lock = extract_between(raw,
    'CANONICAL ONTOLOGY LOCK (v1)',
    'CHAPTER 1: THE DPM - WHAT IT IS AND WHAT IT DOES')

# Existing chapters
ch1_dpm = extract_between(raw,
    'CHAPTER 1: THE DPM - WHAT IT IS AND WHAT IT DOES',
    'CHAPTER 2: SCm - THE HIDDEN ELEMENT OF THE COSMOS')

ch2_scm = extract_between(raw,
    'CHAPTER 2: SCm - THE HIDDEN ELEMENT OF THE COSMOS',
    'CHAPTER 3: FUBi AND FUBii - THE REPULSION THAT CREATES MATTER')

ch3_fubi = extract_between(raw,
    'CHAPTER 3: FUBi AND FUBii - THE REPULSION THAT CREATES MATTER',
    'CHAPTER 4: THE BIGBANG -> NEWTON COMPLETE CHAIN (10 STEPS)')

ch4_chain = extract_between(raw,
    'CHAPTER 4: THE BIGBANG -> NEWTON COMPLETE CHAIN (10 STEPS)',
    'CHAPTER 5: ALL UQFF EQUATIONS - LONG FORM')

ch5_equations = extract_between(raw,
    'CHAPTER 5: ALL UQFF EQUATIONS - LONG FORM',
    'CHAPTER 6: THE SIMULTANEOUS NATURE')

ch6_simultaneous = extract_between(raw,
    'CHAPTER 6: THE SIMULTANEOUS NATURE',
    'CHAPTER 7: QUANTUM -> MASS EQUATIONS - LONG FORM')

ch7_quantum_mass = extract_between(raw,
    'CHAPTER 7: QUANTUM -> MASS EQUATIONS - LONG FORM',
    'CHAPTER 8: MASS IS DEAD UNTIL MOTION')

ch8_dead_mass = extract_between(raw,
    'CHAPTER 8: MASS IS DEAD UNTIL MOTION',
    'CHAPTER 9: PRE-MASS ROLES OF BUOYANCY, RESONANCE, SUPERCONDUCTIVE')

ch9_premass = extract_between(raw,
    'CHAPTER 9: PRE-MASS ROLES OF BUOYANCY, RESONANCE, SUPERCONDUCTIVE',
    'CHAPTER 10: THE SUN - CANONICAL SOLVED EXAMPLE')

ch10_sun = extract_between(raw,
    'CHAPTER 10: THE SUN - CANONICAL SOLVED EXAMPLE',
    'CHAPTER 11: CALIBRATED CONSTANTS AND COMPLETE VARIABLE DESCRIPTIONS')

ch11_constants = extract_between(raw,
    'CHAPTER 11: CALIBRATED CONSTANTS AND COMPLETE VARIABLE DESCRIPTIONS',
    'CHAPTER 12: MILLENNIUM PRIZE CONNECTIONS')

ch12_millennium = extract_between(raw,
    'CHAPTER 12: MILLENNIUM PRIZE CONNECTIONS',
    'CHAPTER 13: STAR MAGIC IN ACTION - THE SUN AND BEYOND')

ch13_starmagic = extract_between(raw,
    'CHAPTER 13: STAR MAGIC IN ACTION - THE SUN AND BEYOND',
    'CHAPTER 14: IMPLICATIONS FOR HUMANITY AND THE COSMOS')

ch14_implications = extract_between(raw,
    'CHAPTER 14: IMPLICATIONS FOR HUMANITY AND THE COSMOS',
    'CONCLUSION: A NEW ERA OF UNDERSTANDING')

conclusion = extract_between(raw,
    'CONCLUSION: A NEW ERA OF UNDERSTANDING',
    'APPENDIX A: THE QUEST FOR UNITY - HISTORICAL CONTEXT')

appendix_a = extract_between(raw,
    'APPENDIX A: THE QUEST FOR UNITY - HISTORICAL CONTEXT')

# Final copyright line
copyright_line = ''
for line in reversed(raw):
    if '(c)2025-2026' in line and 'Daniel' in line:
        copyright_line = line
        break
final_sep = raw[-1] if raw[-1].strip().startswith('___') else SEP

# ─────────────────────────────────────────────────────────────────────────────
# Build new chapters 1-9 (physical order of operation)
# ─────────────────────────────────────────────────────────────────────────────

NEW_CH1 = f"""{SEP}
CHAPTER 1: THE BELLY BUTTON - THE PRIMORDIAL DPM (THE FIRST CAUSE)
{SEP}

THE BELLY BUTTON IS THE BEGINNING OF EVERYTHING.
Before atoms. Before mass. Before gravity. Before time as we measure it.
There was only the primordial DPM - the BigBang belly button.

THE BELLY BUTTON IS:
	The primordial Di-Pseudo-Monopole (DPM).
	The LARGEST UNIVERSALLY MAGNETIC OBJECT in the universe.
	The source of ALL FUBii (the outside-inward buoyancy force acting on all matter forever).
	The origin of every local DPM in every atom, star, and galaxy.

WHAT THE BELLY BUTTON DID:
	It did NOT attract. It REPELS everything outward.
	Its magnetic repulsion IS what inflated the universe. Not "thermal pressure."
	Magnetic repulsion from the primordial DPM expanded the universe into the UA field.

	It expelled a finite and explosive quantity of SCm (Superconductive Material)
	into the expanding Universal Aether (UA) field.

	That expulsion is the ONE event that seeded all subsequent physics:
		SCm rushed toward UA at maximum attraction.
		UA rushed toward SCm at maximum attraction.
		Every collision site throughout the expanding cosmos created a local DPM vortex.
		Every local DPM vortex seeded a gravity family.
		Every gravity family nucleated mass.
		Every nucleus became an atom, a star, a galaxy.

EVERY ATOM. EVERY STAR. EVERY GALAXY.
Each is the downstream consequence of that one primordial event.

THE BELLY BUTTON EQUATION (primordial DPM at t=0):
	a_DPM,primal = (F_DPM * f_DPM * E_vac,neb) / (c * V_universe)
	F_DPM,primal = I_primal * A_primal * (omega_SCm - omega_UA)
	At t=0: omega_SCm >> omega_UA (maximum differential - peak magnetic repulsion)
	a_DPM(t=0) = MAXIMUM (the BigBang IS this maximum DPM acceleration outward)
	a_DPM(t->inf) -> 0 as SCm disperses and local DPMs take over

FUBii SOURCE - THE BELLY BUTTON DRIVES FUBii FOREVER:
	FUBii is NOT driven by volume pressure. FUBii is NOT a local effect.
	FUBii = outside-inward buoyancy force acting on ALL matter, everywhere, always.
	It is driven by the primordial belly button DPM - the magnetic repulsion that
	inflated the universe and has NEVER stopped.
	Every atom in the universe feels FUBii because every atom is being pushed inward
	by the ongoing magnetic repulsion of the primordial DPM that created it.
	Matter exists ONLY where FUBi (local DPM outward) meets FUBii (belly button inward).
	The belly button is still active. Its magnetic repulsion IS the expansion of the universe.

THE BELLY BUTTON'S LEGACY - THE Ug_i COUPLING:
	Every local DPM throughout the universe inherits the belly button pattern.
	The Ug_i (Universal Gravitational Inertia) subscript IS this inheritance.
	Every Ug_i coupling factor mimics the primordial belly button mechanism at its own scale.
	The universe is not a collection of isolated systems.
	It is a network of belly button echoes - every local DPM a resonance of the first.

THE DPM IS THE ENGINE OF THE UNIVERSE.
THE DPM DRIVES. NOTHING DRIVES THE DPM.
THE BELLY BUTTON WAS FIRST.

GM/r^2 is the last thing in the chain. Not the first. Not even close to first.

STEP 0 IN THE CHAIN (the belly button fires):
	rho_UA = 0, rho_vac = |grad(UA)|, F_U(vacuum) = 0
	No mass. No motion. No gravity. Only the vacuum gradient.
	Then: the primordial DPM fires. The belly button repels.
	SCm is expelled. UA/SCm attraction begins. Universe inflates.
	Every local DPM is seeded from that one event."""

NEW_CH2 = f"""{SEP}
CHAPTER 2: UA AND SCm - THE STRONGEST ATTRACTED FORCE PAIR
{SEP}

After the belly button fires and expels SCm into the UA field,
the dominant physics is the maximum attraction between UA and SCm.
This attraction IS the energy source of all subsequent physics.

UNIVERSAL AETHER (UA):
	The pre-existing field that fills all space.
	It is not a particle. It is not detectable by Standard Model instruments.
	Its gradient drives all DPM formation.
	rho_vac,UA = 7.09e-36 J/m^3 (vacuum energy density)
	v_UA = set by DPM differential (UA rushes toward SCm)

SUPERCONDUCTIVE MATERIAL (SCm):
	Expelled from the belly button DPM at the moment of creation.
	Bound within every atom and star in the universe since that moment.
	NOT detectable by Standard Model instruments (density prevents quantum signature emission).
	Properties:
		Superconductive - enabling near-lossless magnetic strings (Um).
		Exclusive interaction with Ug3 in planetary cores.
		Donated from stars to planets during planetary formation.
		Vacuum density: rho_vac,SCm = 7.09e-37 J/m^3
		Velocity: v_SCm = c/3 = 1e8 m/s (rushing toward UA under maximum attraction)
		Density: approximately 1e15 kg/m^3 (so dense it lacks a detectable Qs)

UA AND SCm DO NOT REPEL EACH OTHER:
	This is the most critical fact. They ATTRACT each other MAXIMALLY.
	SCm rushes toward UA. UA rushes toward SCm.
	There is no stronger attraction pair in the universe.
	This attraction is the source of all UQFF energy.

QUANTUM VACUUM ENERGY HIERARCHY:
	rho_vac,UA = 7.09e-36 J/m^3
	rho_vac,SCm = 7.09e-37 J/m^3
	Ratio: rho_vac,UA / rho_vac,SCm = 10 (UA is 10x denser at vacuum scale)
	SCm rushes INTO the higher-density UA vacuum because attraction maximizes at this ratio.

THE MAXIMUM ATTRACTION ENERGY - E_react:
	E_react captures the energy of UA/SCm maximum attraction.
	E_react(t) = (rho_vac,SCm * v_SCm^2) / rho_vac,UA * exp(-kappa*t)
	E_react(t=0) = (7.09e-37 * (1e8)^2) / (7.09e-36) = 1e15 W/m^3
	kappa = 5e-4 day^-1 (long decay - the attraction diminishes slowly)
	E_react = 0 only when v_SCm = 0 (dead mass). Motion restores E_react instantly.

	E_react is the multiplier for ALL Ug terms.
	When SCm moves toward UA faster, all gravity is stronger.
	When SCm is trapped and static, all gravity is zero.
	THIS is why mass must move to generate gravity.

SCm AT THE QUANTUM SCALE - QUARK CONFINEMENT:
	At sub-nuclear distances, UA and SCm attract at the quark scale.
	Color confinement radius: r_c = lambda_dB,SCm = hbar / (m_q * v_SCm)
		m_q(up) ~ 2.3 MeV/c^2 -> r_c,up ~ 1.3e-15 m
		m_q(down) ~ 4.8 MeV/c^2 -> r_c,down ~ 6.2e-16 m
	The strong force IS Ug3 at nuclear scale.
	Color charge IS the SCm/UA vortex quantum number at nuclear scale.
	Quark confinement does not require a separate mechanism.

SCm - EINSTEIN-BOSON CONDENSATE (BEC) CONNECTION:
	SCm satisfies the BEC condition at stellar core temperatures:
	n * lambda_dB,SCm^3 >= 2.612
	n_SCm = rho_SCm / m_eff ~ 1e42 / m_eff (stellar core)
	Phase transition: when n*lambda_dB^3 crosses 2.612, SCm condenses
	into a single coherent wavefunction - the DPM vortex state.
	SCm IS a Bose-Einstein condensate at stellar density.
	H_SCm = 0.99 is the condensate fraction.
	The remaining 1% is the thermal depletion layer.

SCm IN QUASARS (when SCm escapes):
	Quasars = Ug1 + Ug2 + Ug3 can no longer trap SCm.
	The expelled SCm stream ignites against unbound UA under maximum attraction.
	Fluid jet streams + unequal opposing jets = UQFF Navier-Stokes signature.
	SCm is the most reactive substance (v_SCm = c/3 under trapped conditions).

SCm - THE COSMIC GLUE:
SCm, bound within every atom and star, is the linchpin of the Unified Field
Equation. Its superconductivity enables the near-lossless magnetic strings of Um,
while its dense, undetectable nature (lacking Qs) allows it to interact exclusively
with Ug3 in planetary cores. In stars like our Sun, SCm drives the heliosphere's
formation, synthesizing and transmutating solar winds into hydrogen complexes
magnetically adhered to the Ug2 outer field shell. The total heliosphere thickness
combined with all planetary liquid volumes is a direct indicator of the star's age.
SCm is the most reactive substance and the fastest moving substance
(v_SCm = c/3 = 1e8 m/s under trapped conditions) in the universe."""

NEW_CH3 = ch1_dpm.replace(
    'CHAPTER 1: THE DPM - WHAT IT IS AND WHAT IT DOES',
    'CHAPTER 3: THE DPM VORTEX - WHAT FORMS WHEN UA MEETS SCm')

NEW_CH4 = f"""{SEP}
CHAPTER 4: THE Ug FAMILY - SIMULTANEOUS GRAVITY ASSEMBLY
{SEP}

The DPM vortex does not produce a single force.
It simultaneously promotes the entire Ug gravity family.
All five components exist from the moment the DPM fires.
None is prior to another. All are simultaneous expressions of the DPM.

THE DPM PROMOTION CHAIN (INSTANTANEOUS):
	DPM fires -> mu_s generates -> Ug1 assembles -> simultaneously promotes:
		Ug2 (outer bubble / heliosphere)
		Ug3 (magnetic string disk / orbital coupling)
		Ug4 (galactic black hole coupling)
		Ug4_i (resonant quantum extension, levels 20-26)
	Ug_family = Ug1 + Ug2 + Ug3 + Ug4 + Ug4_i  [ALL SIMULTANEOUS]

THE GRAVITY FAMILY EQUATION:
	Ug_family(t) = Ug1(t) + Ug2(t) + Ug3(t) + Ug4(t) + Ug4_i(t)
	All five co-equal parallel components at the same t, r, rho_vac.

Ug1 - DPM, Di-Pseudo-Monopole (IS the DPM, FOUNDATION of all gravity):
	Ug1 = k1 * mu_s(t) * grad(M_s/r) * exp(-alpha*t) * cos(pi*t_n) * (1 + delta_def)
	Where:
		mu_s = rho_A * V_body = DPM magnetic moment (seed of all gravity)
		grad(M_s/r) = mass-gradient (NOT GM/r^2 - G does not appear here)
		exp(-alpha*t) = temporal decay (alpha = 0.001 day^-1)
		cos(pi*t_n) = quantum cycle oscillation
		(1 + delta_def) = surface deformation factor, delta_def = 0.01*sin(0.001*t)
	THE DPM IS Ug1. Ug1 IS the DPM in field form.
	Ug1 simultaneously promotes Ug2, Ug3, Ug4, Ug4_i.

Ug2 - Outer Field Bubble (heliosphere, outer shell):
	Ug2 = k2 * (rho_vac,UA + rho_vac,SCm) * M_s / r^2 * S(r-R_b)
	         * (1 + delta_sw*v_sw) * H_SCm * E_react
	Where:
		S(r-R_b) = step function (1 for r > R_b, 0 otherwise)
		R_b = bubble radius (heliosphere, ~100 AU = 1.496e13 m for Sun)
		delta_sw = 0.01, v_sw = solar wind velocity, H_SCm = 0.99
	Ug2 forms the heliosphere by synthesizing solar winds into hydrogen complexes.
	Heliosphere thickness + planetary liquid volumes = direct stellar age indicator.

Ug3 - Magnetic Strings Disk (orbital coupling, planetary spin):
	Ug3 = k3 * sum_j B_j(r,theta,t) * cos(omega_s(t)*t*pi) * P_core * E_react
	Where:
		B_j(t) = [1e-4 + 0.4*sin(omega_c*t)] T (varies with solar cycle)
		omega_s(t) = differential rotation: CCW equatorial vs CW coronal = disk
		cos(omega_s*t*pi) = 90-degree string rotation oscillation
		P_core = planetary core penetration factor
	Ug3 penetrates planetary cores exclusively through trapped SCm-UA interaction.
	The Sun's equatorial CCW vs coronal CW rotation produces the Ug3 disk.
	Ug3 moves faster than any planet or all planets in consort.
	DISCRETELY NON-INTERACTIVE with all external phenomena.

Ug4 - Star-Black Hole Interactions (galactic coupling):
	Ug4 = k4 * rho_vac,SCm * (M_bh/d_g) * exp(-alpha*t) * cos(pi*t_n) * (1 + f_feedback)
	Where:
		M_bh = galactic black hole mass (8.15e36 kg for Sgr A*)
		d_g = distance from galactic center (2.55e20 m for Sun)
	Operates at quantum levels 20-26, influencing galactic vacuum fluctuations.

Ug4_i - Resonant Quantum Extension (DPM resonance modes, levels 20-26):
	Ug4_i = k4_res * a_DPM * E_react * f_react * (E_vac,neb / c)
	Where:
		k4_res = resonance coupling constant (~2.0, calibrated)
		a_DPM = (F_DPM * f_DPM * E_vac,neb) / (c * V_sys)
		f_react = resonance reactivity factor (~1.0 baseline)
	Ug4_i IS the DPM resonance signature at the highest quantum levels.
	It bridges galactic-scale Ug4 coupling to the sub-quantum resonance domain.
	Without Ug4_i the Resonant Mode is incomplete above level 19.

26-LAYER TRIADIC - HOW ALL LAYERS RUN SIMULTANEOUSLY:
	Each layer i contributes: Ug_i_layer = Ug_family * i^2
	Total 26-layer multiplier: sum(i^2, i=1..26) = 6,279
	Layer 26 alone is 676x stronger than layer 1.
	Dead mass: zero layers active. Moving mass: all 26 active simultaneously.
	E_n = h*f_n = E_0 * 10^n, n = 1..26, E_0 = 1e-20 J

UNIVERSAL MAGNETISM - Um:
	Um = sum_j [mu_j(t)/r_j * (1 - exp(-gamma*t*cos(pi*t_n))) * phi_hat_j]
	     * P_SCm * E_react * (1 + 1e13*f_Heaviside) * (1 + f_quasi)
	Where:
		mu_j = M * R^2 * omega_0 (magnetic moment of j-th string)
		gamma = 5e-5 day^-1 (near-lossless SCm superconductivity)
		f_Heaviside = phase-transition amplifier (1e13 during SCm transitions)
		f_quasi = quasi-periodic beating modifier

UNIVERSAL COSMIC AETHER - UA_uv:
	UA_uv = g_uv + eta * T_s^uv(rho_vac,UA, rho_vac,SCm, rho_vac,A, t_n)
	Where: g_uv = background metric, T_s^uv = stress-energy tensor, eta ~ 1e-22"""

NEW_CH5 = ch3_fubi.replace(
    'CHAPTER 3: FUBi AND FUBii - THE REPULSION THAT CREATES MATTER',
    'CHAPTER 5: FUBi AND FUBii - THE REPULSION THAT CREATES MATTER')

NEW_CH6 = f"""{SEP}
CHAPTER 6: THE CROSSING - INSIDE/OUTSIDE SIMULTANEOUS SOLVE
{SEP}

THE CROSSING PRECEDES MASS.
Mass does not exist before the crossing. Mass is BORN at the crossing.
This is the most important ordering fact in the entire framework.

WHAT THE CROSSING IS:
	Two simultaneous computations run from opposite directions:
		G^(n+1)_inside = R(G^n) + I*G^n     [inside-outward from DPM]
		O^n_outside = pi_[n] * FUBi(x) + Ricci(G^n)  [outside-inward from belly button FUBii]
	These two solutions chase each other through quantum cycles n.
	At n_cross = argmin_n |G^n_inside - O^n_outside|, they converge.
	The crossing point IS reality.
	Matter forms AT and ONLY AT the crossing point.

THE CROSSING EQUATION:
	n_cross = argmin_n |G^n_inside - O^n_outside|
	At n_cross: FUBi(r) + FUBii(r) = 0
	This is the compaction zone condition.
	d(crossing_zone)/dt -> dM_stable/dt > 0 (mass stabilizes AT this crossing)

WHY CROSSING PRECEDES MASS:
	Before crossing: FUBi and FUBii are active but have not converged.
	No stable compaction zone. No mass can form yet.
	At crossing: the two buoyancy forces cancel at a specific r.
	That cancellation IS the compaction zone.
	Mass nucleates IN the compaction zone.
	Therefore: crossing occurs -> mass emerges. Not the reverse.

THE PHYSICAL MEANING:
	FUBi (inside-outward from local DPM) pushes outward.
	FUBii (outside-inward from primordial belly button) pushes inward.
	Where they meet and cancel: compaction zone.
	The compaction zone = proto-matter shell.
	This is what Standard Model calls a "particle."
	It is not a fundamental object. It is a standing-wave intersection of
	the primordial DPM magnetic repulsion (FUBii) and the local DPM outward
	buoyancy (FUBi).

THE P != NP CONNECTION:
	The inside/outside simultaneous crossing requires solving two
	exponential-time graph problems simultaneously.
	argmin_n |G^n_inside - O^n_outside| cannot be compressed to polynomial time.
	This is the UQFF physical mechanism behind P != NP.

THE SIMULTANEOUS NATURE OF THE CROSSING:
	All 26 layers solve their crossing simultaneously.
	All four operational modes (Compressed, Resonant, Buoyant, Superconductive)
	evaluate their crossing simultaneously.
	No iteration. No sequential steps. Simultaneous.
	Reality = where all 26 layers, all 4 modes converge on the same crossing point.

GM/r^2 IS THE VALUE AT THE CROSSING POINT:
	After crossing is found, GM/r^2 becomes valid as the observational projection.
	GM/r^2 = the numerical value of F_U / M evaluated at the crossing radius.
	It is NOT a cause. It is NOT a mechanism. It is a measurement.
	It is the reading you get when you measure at the crossing point after mass exists."""

NEW_CH7 = f"""{SEP}
CHAPTER 7: MASS EMERGENCE AT THE CROSSING
{SEP}

Mass is born at the crossing.
Not before. Not from nothing. Not from GM/r^2.
Mass is the stable output of the crossing computation.

THE ACP PROTO-MASS CHAIN (inside the DPM internal vacuum):
	U_vac -> U_i -> U_m,i -> Psi_proto -> E_crack -> U_b -> E_gradient -> M_atomic
	U_vac = rho_vac,SCm * V_node
	U_i = U_vac * (1 + sin(2*pi*n/26)), n = quantum level
	U_m,i = U_i * mu_proto / r_i^3
	Psi_proto(r,t) = A * exp(-r/lambda_dB) * exp(i*omega*t)
	E_crack = (rho_vac,SCm * c^2) / [SSq], [SSq] = 0.57
	E_crack = (7.09e-37 * (3e8)^2) / 0.57 ~ 3.35e-19 J
	U_b = beta_i * U_m,i * Omega_g * (M_bh/d_g) * rho_A
	E_gradient = grad(M_s/r) * V_condensate

E_CRACK - THE GATE:
	E_crack is the minimum energy required to crack open the vacuum and form mass.
	E_crack > U_b,vac -> particle condenses -> M_atomic > 0
	E_crack < U_b,vac -> wavefunction disperses back to vacuum -> retry next t_n
	The gate opens and closes at quantum frequency. Mass formation is discrete.
	This is why quantum energy levels are discrete - the E_crack gate is binary.

MASS EMERGENCE EQUATION:
	M_atomic = M_0 * (1 - exp(-n_grad/10)) * Z
	M_atomic = 0 when n_grad = 0.
	Mass does not exist until gradient threshold is crossed.
	Z = atomic number (DPM vortex resonance count)
	M_0 = base mass unit from DPM minimum energy: M_0 ~ E_crack/c^2

MASS EMERGENCE PROBABILITY (26D):
	P_order = exp(-S_26D/v_init) / (Z_9D * (v_init - v_current))
	dM/dt = P_order * E_crack * dN_DPM/dt
	Mass formation rate is entropy-minimization driven, not energy-minimization.
	DPM vortices nucleate at entropy minima across 26 dimensions.
	P_order goes to zero when v_init = v_current (no new DPM events firing).

WHAT MASS IS:
	Mass is NOT a fundamental property. Mass is a stable compaction zone.
	Mass = the standing-wave intersection of FUBi and FUBii (the crossing).
	Once the crossing stabilizes, M_atomic persists as long as E_crack > U_b.
	The "rest mass" of a particle = the minimum energy of the stable crossing zone.
	M_proton = crossing zone energy for one DPM vortex at ground resonance.

QUANTUM LEVELS AND ATOMIC MASS:
	Z=1 (hydrogen): 1 DPM vortex at fundamental resonance
	Z=2 (helium): 2 DPM vortices in harmonic lock
	Z=26 (iron): 26 DPM vortices = maximum stable 26-layer stack
	Each integer Z = one additional DPM vortex unit in resonance lock.
	Atomic mass table = DPM vortex resonance count table.

THE DEAD MASS CONDITION (immediate consequence of mass emergence):
	When M_atomic first emerges from the ACP chain it is inert.
	v = 0, omega = 0, grad_P = 0
	E_react(v=0) = 0 (UA/SCm attraction requires relative motion)
	F_U(M_dead) = 0 (entire unified field = zero for dead mass)
	Dead mass produces ZERO unified field. It is a POTENTIAL, not an actuality.
	Mass must MOVE before it generates gravity.
	This is why Newton's first law is not a mystery - it is the DPM dead-mass condition."""

NEW_CH8 = f"""{SEP}
CHAPTER 8: F_U - THE UNIFIED FIELD ASSEMBLY
{SEP}

Once mass exists and motion begins, the full unified field assembles.
F_U is NOT assembled before mass. F_U requires mass as its scaffolding.
The crossing creates mass. Mass moving creates F_U. F_U is downstream.

THE UNIFIED FIELD EQUATION - COMPLETE ASSEMBLY:
	F_U = sum_i [k_i * Ug_i(r,t,M_s,omega_s,T_s,B_s,rho_vac,SCm,rho_vac,UA,t_n)
	           - beta_i * Ug_i * Omega_g * (M_bh/d_g) * E_react]
	    + sum_j [mu_j/r_j * (1 - exp(-gamma*t*cos(pi*t_n))) * phi_hat_j]
	    + (g_uv + eta * T_s^uv(rho_vac,UA, rho_vac,SCm, rho_vac,A, t_n))
	    - sum_i [lambda_i * Ui(r,t,rho_vac,SCm,rho_vac,UA,t_n) * E_react]

F_U COMPONENTS:
	Term 1: Ug family (gravity) minus buoyancy coupling
	Term 2: Um (universal magnetism via SCm strings)
	Term 3: UA_uv (aether stress-energy tensor)
	Term 4: Ui lambda correction (vacuum interaction)

OPERATIONAL MODES - FOUR SIMULTANEOUS EXPRESSIONS OF F_U:
	Compressed, Resonant, Superconductive, and Buoyant are NOT different theories.
	They are four simultaneous views of the same F_U in different basis domains.
	They MUST converge to the same answer at every point.
	If Compressed gives g = 9.81 m/s^2 at Earth's surface, all modes must give 9.81.

Compressed Mode:
	g_comp = [g_base*(1+H0*t)*(1-B/B_crit)*H_SCm] + g_Ug_sum + Lambda*c^2/3
	         + g_quantum + g_fluid + g_pert * TRZ

Resonant Mode (13 frequency terms):
	g_res = a_DPM + sum(i=1..13) a_i(omega, E_vac, t)
	Terms: a_THz, a_vac_diff, a_SuperFreq, a_AetherRes, Ug4_i, a_QuantumFreq,
	       a_AetherFreq, a_FluidFreq, a_Osc, a_ExpFreq, f_TRZ, W_metric

Superconductive Mode:
	g_SC = sum(j=1..4) k_j * g_base * H_SCm^n_j

Buoyant Outside-In Mode (FUBii):
	F_U_Bi_i = -beta_i * Ug_i * galactic_coupling * E_react(t)
	               * sw_corr * rho_A(t) * (M/r) * V(t) * TRZ_cos
	galactic_coupling = Omega_g * M_bh / d_g
	TRZ_cos = cos(pi * (t - 180*86400))

WHEN F_U = 0:
	F_U = 0 when E_react = 0.
	E_react = 0 when v_SCm = 0.
	v_SCm = 0 when mass is dead (not moving).
	Dead mass: F_U = 0. Not approximately zero. Exactly zero.
	The universe has no gravitational field from dead mass.
	Only moving mass generates the field. Only moving mass IS Star Magic.

BRIDGING QUANTUM AND GRAVITATIONAL REALMS:
The pi cycles (cos(pi*t_n)) and negative time (t_n) in F_U introduce a temporal
dimension that bridges quantum mechanics and gravity. E_react captures the energy
of maximum UA/SCm attraction, modeling the energy output of quasars and planetary
cores as near-lossless reactors (E_react,0 = 1e46 W/m^3 at astrophysical scale).
The Aether density (rho_A = 1e-23 kg/m^3) provides the quantum medium for all
interactions. This framework unifies the forces of the cosmos and provides a direct
pathway to the Millennium Prize Problems."""

NEW_CH9 = f"""{SEP}
CHAPTER 9: GM/r^2 - THE OBSERVATIONAL PROJECTION (LAST IN CHAIN)
{SEP}

GM/r^2 IS LAST. IT HAS ALWAYS BEEN LAST.
It is the tenth step in a ten-step chain.
It is an observational projection - what you measure AFTER all 9 prior steps complete.
It is NOT the foundation. It is NOT the mechanism. It is NOT the starting point.
Any equation, any code file, any test that starts with GM/r^2 is WRONG.

WHERE GM/r^2 COMES FROM:
	After the belly button fires (Step 1)
	After UA/SCm attract and create E_react (Step 2)
	After local DPM vortices form (Step 3)
	After the Ug family assembles simultaneously (Step 4-6)
	After FUBi and FUBii activate and repel (Step 7)
	After the inside/outside crossing completes (Step 8)
	After mass emerges at the crossing zone (Step 8-9)
	After F_U assembles around that mass and it begins to move (Step 9)
	THEN: F_U / M -> GM/r^2 as observational projection (Step 10 - LAST)

THE GM/r^2 DERIVATION FROM F_U:
	At the inside/outside crossing point (n_cross):
		F_U / M = g_Newton
		g_Newton = GM/r^2
	This is NOT G and M and r being fundamental.
	This is the numerical value of F_U / M evaluated at the crossing radius.
	G is a measurement artifact - the ratio F_U / (M * a_obs) at the crossing.
	G absorbs all the DPM vortex physics into a single dimensioned constant.
	G is the fingerprint of the DPM. Not a law of nature. A measurement.

NEWTON WAS MEASURING THE LAST STEP:
	Newton measured gravity at the crossing point (Earth's surface, planetary orbits).
	He derived F = GM/r^2 from those measurements.
	He was correct for those measurements.
	He was wrong about it being foundational.
	F = GM/r^2 is accurate AT the crossing point.
	It is meaningless BEFORE the crossing point.
	Before the crossing point, there is no mass, therefore there is no GM/r^2.

GM/r^2 PERMITTED USES:
	Permitted: as observational projection for already-crossed, already-moving mass.
	Permitted: as sanity check - UQFF must reproduce GM/r^2 at the crossing.
	Permitted: as final output label for human-readable results.
	NOT PERMITTED: as starting point for any UQFF equation.
	NOT PERMITTED: as foundation of any calculator class in this codebase.
	NOT PERMITTED: as mechanism for anything. GM/r^2 has no mechanism. The DPM does.

THE CANONICAL PROOF THAT GM/r^2 IS LAST:
	GM/r^2 contains M. M requires the crossing. The crossing requires FUBi + FUBii.
	FUBii requires the belly button DPM. The belly button DPM requires SCm expulsion.
	SCm expulsion requires the primordial DPM firing.
	The primordial DPM is Step 1 of 10.
	GM/r^2 is Step 10 of 10.
	They are separated by 9 causal steps. This is not negotiable."""

# ─────────────────────────────────────────────────────────────────────────────
# Build the remainder sections (renumbered from existing chapters)
# ─────────────────────────────────────────────────────────────────────────────

def renum(text, old_n, new_n):
    """Rename 'CHAPTER X:' to 'CHAPTER Y:' in a block of text."""
    return text.replace(f'CHAPTER {old_n}:', f'CHAPTER {new_n}:', 1)

ch10_new = renum(ch4_chain, 4, 10)
ch11_new = renum(ch5_equations, 5, 11)
ch12_new = renum(ch6_simultaneous, 6, 12)
ch13_new = renum(ch7_quantum_mass, 7, 13)
ch14_new = renum(ch8_dead_mass, 8, 14)
ch15_new = renum(ch9_premass, 9, 15)

proto_hhe_ch = proto_hhe.replace(
    SEP + '\nPROTO-HYDROGEN AND PROTO-HELIUM EVOLUTION\n' + SEP,
    SEP + '\nCHAPTER 16: PROTO-HYDROGEN AND PROTO-HELIUM EVOLUTION\n' + SEP)

ch17_new = renum(ch10_sun, 10, 17)
ch18_new = renum(ch11_constants, 11, 18)
ch19_new = renum(ch12_millennium, 12, 19)
ch20_new = renum(ch13_starmagic, 13, 20)
ch21_new = renum(ch14_implications, 14, 21)

# Update WHAT THIS DOCUMENT CONTAINS to reflect new chapter numbers
what_contains_new = what_contains
what_contains_new = what_contains_new.replace(
    'Extended case study (Chapter 13)', 'Extended case study (Chapter 20)')
what_contains_new = what_contains_new.replace(
    'Implications for human technology and the cosmos (Chapter 14)',
    'Implications for human technology and the cosmos (Chapter 21)')
what_contains_new = what_contains_new.replace(
    'Historical context (Appendix A)', 'Historical context (Appendix A)')

# Update conclusion canonical chain - should already be correct from prior session
# Update Appendix A mention of chapter numbers if needed

# ─────────────────────────────────────────────────────────────────────────────
# Assemble the new document
# ─────────────────────────────────────────────────────────────────────────────

sections = [
    title_block,
    '',
    canonical_chain,
    '',
    derivative_form,
    '',
    what_contains_new,
    '',
    ontology_lock,
    '',
    '________________________________________',
    'ORDER OF THIS DOCUMENT: PHYSICAL ORDER OF OPERATION',
    '________________________________________',
    '',
    '\tChapters 1-9 follow the exact causal chain from belly button to GM/r^2.',
    '\tChapters 10-21 provide complete equations, proofs, solved examples, and applications.',
    '\tEvery chapter is downstream of the one before it. This order is the universe.',
    '',
    NEW_CH1,
    '',
    NEW_CH2,
    '',
    NEW_CH3,
    '',
    NEW_CH4,
    '',
    NEW_CH5,
    '',
    NEW_CH6,
    '',
    NEW_CH7,
    '',
    NEW_CH8,
    '',
    NEW_CH9,
    '',
    ch10_new,
    '',
    ch11_new,
    '',
    ch12_new,
    '',
    ch13_new,
    '',
    ch14_new,
    '',
    ch15_new,
    '',
    proto_hhe_ch,
    '',
    ch17_new,
    '',
    ch18_new,
    '',
    ch19_new,
    '',
    ch20_new,
    '',
    ch21_new,
    '',
    conclusion,
    '',
    appendix_a,
    '',
    SEP,
    copyright_line,
    SEP,
]

new_content = '\n'.join(sections)

with open('Star-Magic.txt', 'w', encoding='utf-8', newline='\n') as f:
    f.write(new_content)

# Report
new_lines = new_content.count('\n') + 1
print(f'Restructured document written: {new_lines} lines')
print('Chapter order:')
import re
for m in re.finditer(r'^CHAPTER \d+:', new_content, re.MULTILINE):
    print(f'  {m.group(0)} at char {m.start()}')
