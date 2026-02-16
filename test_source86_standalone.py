from source86_extended_extract import Source86_Extended, SystemType

# Test compressed MUGE (Magnetar SGR 1745-2900)
print("Testing SOURCE86: Extended Fields MUGE")
print("=" * 80)
print("\n1. Compressed MUGE (Magnetar SGR 1745-2900):")
print("-" * 80)

result_comp = Source86_Extended.calculate_muge_compressed(t=3.799e10)

print(f"g_total         = {result_comp['g_total']:.6e} m/s²")
print(f"g_base          = {result_comp['g_base']:.6e} m/s² (dominant)")
print(f"expansion       = {result_comp['expansion_factor']:.6f}")
print(f"sc_correction   = {result_comp['sc_correction']:.6f}")
print(f"ug_sum          = {result_comp['ug_sum']:.6e} m/s²")
print(f"lambda_term     = {result_comp['lambda_term']:.6e} m/s²")
print(f"quantum_term    = {result_comp['quantum_term']:.6e} m/s²")
print(f"fluid_term      = {result_comp['fluid_term']:.6e} m/s²")
print(f"sys_term        = {result_comp['sys_term']:.6e} m/s²")

print("\n2. Resonance MUGE (Nebula mode):")
print("-" * 80)

result_res = Source86_Extended.calculate_muge_resonance(t=3.799e10)

print(f"g_total         = {result_res['g_total']:.6e} m/s²")
print(f"adpm (base)     = {result_res['adpm']:.6e} m/s² (foundation)")
print(f"athz            = {result_res['athz']:.6e} m/s²")
print(f"avac_diff       = {result_res['avac_diff']:.6e} m/s²")
print(f"asuper_freq     = {result_res['asuper_freq']:.6e} m/s²")
print(f"aaether_res     = {result_res['aaether_res']:.6e} m/s²")
print(f"afluid_freq     = {result_res['afluid_freq']:.6e} m/s²")

print("\n3. Sagittarius A* system-specific term:")
print("-" * 80)

sys_term_sgra = Source86_Extended.calculate_system_specific_term(
    SystemType.SAGITTARIUS_A,
    t=3.799e10,
    M=4.3e6 * 1.989e30,  # 4.3 million M_sun
    r=2.5e10,            # ~25 AU Schwarzschild radius
    dOmega_dt=1e-10
)
print(f"SgrA* frame-dragging = {sys_term_sgra:.6e} m/s²")

print("\n✅ SOURCE86 Extended Fields standalone test PASSED")
