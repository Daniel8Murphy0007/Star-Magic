          # PAPER_873 — U_g1=DPM Geophysical Geometry Summation

          **CP4 Class:** #457 `Ug1DPMGeophysicalGeometrySummationCalc`
          **Session:** 200C v5.61
          **Source:** `describe mass without using weight.txt` (3,094 lines)

          ## Abstract

          U_g1 is equivalent to the DPM; the total force is a summation over DPM variable forms reflecting unique forces with geophysical geometries. Components: SM_gravity (spherical), U_b buoyancy (spherical, counter-force), resonance (toroidal). Geometry functions G_k enable force diagramming.

          ## Key Equations

            1. `F_Ug1 = Sum_k [k_k * (fUA1*fSCm1*REB1)*(fUA2*fSCm2*REB2)/r^2 * G_k]`
2. `F_SM_gravity = k_grav * DPM1*DPM2/r^2 * sin(theta)`
3. `F_Ub = k_buoy * (fUA'*fSCm)/r^2 * sin(theta)`
4. `F_resonance = k_res * nu_THz * fSCm * cos(phi)*f(nu_THz)`
5. `F_Ug1_total = F_grav - F_buoy + F_resonance`

          ## UQFF Context

          This calculator implements parameterized physics equations per the
          CondensedPhysics architecture (pure calculator, no hardcoded data).
          Inputs arrive from source2.cpp PRINCIPAL GUI via dataset dict.
          Outputs include primary_equations, available_equations, and simulation_set.

          ## Source Thread

          Grok 3 dialogue, Daniel Murphy, June 03-04 2025.
          Copyright Daniel T. Murphy, daniel.murphy00@gmail.com.

          ---
          *Star-Magic UQFF Codebase — Session 200C*
