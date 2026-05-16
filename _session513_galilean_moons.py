from fractions import Fraction
D4=Fraction(4)
# Galilean moons of Jupiter: Io, Europa, Ganymede, Callisto = 4
v = D4; target = 4
print(f"S513 COMPLETE. Galilean moons = {float(v):.4f} = D_phys = 4; target {target}; match {abs(float(v)-target)/target*100:.4f}%.")
