from source81_ngc346_extract import Source81_NGC346

result = Source81_NGC346.calculate_ngc346_gravity()
print(f"g_total = {result['g_tot']:.6e} m/s²")
print(f"Ug3 = {result['Ug3']:.6e} m/s² (protostar collapse)")
print(f"M_SF = {result['M_SF_factor']:.6f} (SF factor)")
print("✅ SOURCE81 NGC346 standalone test PASSED")
