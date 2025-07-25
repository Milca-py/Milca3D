from milcapy.model import Model
from math import tanh, pi
portico = Model(3)

# Parametros:
L = 4

# Materiales y secciones
h = 0.5
b = 0.3
E = 2.1e6
v = 0.2

G = E / (2 * (1 + v))
A = h * b
Jxx = (1/3)*b**3*h*(1-192*b/(h*pi**5)*(tanh(h*pi/(2*b))))
Iy = h * b**3 / 12
Iz = b * h**3 / 12
Asy = 5/6 * A
Asz = 5/6 * A

coords = {
    1: (0, 0, 0),
    2: (L, 0, 0),
    3: (2*L, 0, 0),
}

for i, coord in coords.items():
    portico.add_node(i, *coord)

portico.def_geom_transf(2, *[0, 0, 1])

portico.add_member(1, 1, 2, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(2, 2, 3, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)


fixed = (False, False, True, True, True, True)
for i in [1, 2, 3]:
    portico.add_restraint(i, fixed)


portico.add_point_load(2, fz=100)

portico.solve()

portico.show_deformed(5, 100)

portico.print.node_results()
