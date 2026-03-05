[33mcommit 53ae55ba52644519507da48ad4e9c406f20bbc5a[m
Author: Daniel T Murphy <daniel.murphy00@gmail.com>
Date:   Mon Feb 16 12:03:54 2026 -0500

    Add 26D quantum state framework parameters and 19-system constants (Feb 16, 2026)
    
    - Add 26D framework params: T_0_lock, tau_lock, r_THz_base, f_Um_scale, B_crit_UQFF, hbar_c
    - Add vacuum density scales: rho_UA_scale, rho_SCm_scale, rho_UA_SCm_ratio
    - Add angular params: theta_base, theta_step, v_char, SCm_B_scale, f_TRZ_scale
    - Add 19 astrophysical systems from source172.cpp (SOURCE115):
      NGC 2264, UGC 10214, NGC 4676, Red Spider, NGC 3372, AG Carinae, M42,
      Tarantula, NGC 2841, Mystic Mountain, NGC 6217, Stephan's Quintet,
      NGC 7049, NGC 3324, M74, NGC 1672, NGC 5866, M82, Spirograph
    - Add document-variant parameters: NGC2264_M_doc, NGC2264_r_doc, wind velocities
    - Total CONSTANTS: 799 keys
    - All 34 tests passing

 README.md | 41 [32m+++++++++++++++++++++++++++++++++[m[31m--------[m
 1 file changed, 33 insertions(+), 8 deletions(-)
