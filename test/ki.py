from milcapy.frame import local_stiffnes_matrix
import pandas as pd

# Parametros:
L = 4

# Materiales y secciones
hr = 0.5
br = 0.3
E: float = 2.1e6
v: float = 0.2
G: float = E / (2 * (1 + v))
A: float = hr * br
Jx: float = hr*br**3*(16/3 - 3.36*br/hr*(1-(br/hr)**4/12))/16
Iy: float = hr * br**3 / 12
Iz: float = br * hr**3 / 12
Asy: float = 5/6 * A
Asz: float = 5/6 * A
phiY = 12.0 * E * Iz / (L**2 * G * Asy)
phiZ = 12.0 * E * Iy / (L**2 * G * Asz)


kl = local_stiffnes_matrix(E, G, A, L, Jx, Iy, Iz, phiY, phiZ)


df = pd.DataFrame(kl.round(3))
print(df.to_string(index=False, header=False))

print(Jx)