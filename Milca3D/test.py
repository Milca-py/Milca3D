from milcapy.core.model import Model
from math import tanh, pi

portico = Model()

# Parametros:
H = 4
B = 4
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
    3: (L, B, 0),
    4: (0, B, 0),
    5: (0, 0, H),
    6: (L, 0, H),
    7: (L, B, H),
    8: (0, B, H),
}

for i, coord in coords.items():
    portico.add_node(i, *coord)


portico.def_geom_transf(1, *[1, 0, 0])
portico.def_geom_transf(2, *[0, 0, 1])

portico.add_member(1, 1, 5, E, G, A, Jxx, Iy, Iz, Asy, Asz, 1)
portico.add_member(2, 2, 6, E, G, A, Jxx, Iy, Iz, Asy, Asz, 1)
portico.add_member(3, 3, 7, E, G, A, Jxx, Iy, Iz, Asy, Asz, 1)
portico.add_member(4, 4, 8, E, G, A, Jxx, Iy, Iz, Asy, Asz, 1)
portico.add_member(5, 5, 6, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(6, 6, 7, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(7, 8, 7, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(8, 5, 8, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)

portico.add_member(9, 2, 5, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(10, 1, 8, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)

portico.add_member(11, 8, 3, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)
portico.add_member(12, 7, 2, E, G, A, Jxx, Iy, Iz, Asy, Asz, 2)


fixed = (True, True, True, True, True, True)
for i in [1, 2, 3, 4]:
    portico.add_restraint(i, fixed)

portico.add_point_load(5, fx=100, fy=100, fz=100)
portico.add_distributed_load(5, direction=2, wi=-100, wj=-100)
portico.add_distributed_load(6, direction=2, wi=-100, wj=-100)
portico.add_distributed_load(7, direction=2, wi=-100, wj=-100)
portico.add_distributed_load(8, direction=2, wi=-100, wj=-100)

portico.solve()

portico.show_model()
portico.show_deformed(10, 100)
portico.print_node_results()
