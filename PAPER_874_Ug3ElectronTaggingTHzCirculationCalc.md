          # PAPER_874 — U_g3=U_i+U_m Electron Tagging and THz Circulation

          **CP4 Class:** #458 `Ug3ElectronTaggingTHzCirculationCalc`
          **Session:** 200C v5.61
          **Source:** `describe mass without using weight.txt` (3,094 lines)

          ## Abstract

          U_g3 = U_i + U_m in motion. U_i individually tags each electron via the THz hole pipeline from nucleus (Point A) to electron (Point B), counteracting the strong-force to position electrons at projected quantum energy shell balance. Electron DPM uses limited U_i for coherent circulation in U_g2 imaginary orbit shell. U_g2 monitors and adjusts shell position.

          ## Key Equations

            1. `F_Ug3 = (F_Ui + F_Um + F_DPM_e) * G_geo / r_shell^2`
2. `F_Ui = k_i * fUA'_nuc * nu_THz * R_EB   (repulsive tagging)`
3. `F_Um = k_m * fSCm_nuc * nu_res   (SCm-driven motion)`
4. `F_DPM_e = k_e * (fUA'_e * fSCm_e) * nu_THz   (electron circulation)`
5. `G_geo = sin(theta)*cos(phi)*f(nu_THz)`
6. `E_shell = c * nu_res * hbar * fSCm * sin(theta)   (U_g2 monitor)`

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
