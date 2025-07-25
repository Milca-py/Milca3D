from milcapy.model import Model

portico = Model(3)

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

vec_col = [1, 0, 0]
vec_viga = [0, 0, 1]

# Tranformacion de coordenadas
portico.def_geom_transf(1, *vec_col)
portico.def_geom_transf(2, *vec_viga)

# defimnicion de los nodos
portico.add_node(1, 0, 0, 0)
portico.add_node(2, 0, 0, L)
portico.add_node(3, L, 0, L)
portico.add_node(4, 0, -L, L)

# defimnicion de los elementos
el1 = portico.add_member(1, 1, 2, E, G, A, Jx, Iy, Iz, Asy, Asz, 1)
el2 = portico.add_member(2, 2, 3, E, G, A, Jx, Iy, Iz, Asy, Asz, 2)
el3 = portico.add_member(3, 2, 4, E, G, A, Jx, Iy, Iz, Asy, Asz, 2)

fixed = (True, True, True, True, True, True)
portico.add_restraint(1, fixed)
portico.add_restraint(3, fixed)
portico.add_restraint(4, fixed)

portico.add_point_load(2, fx=5, fy=10, fz=20)

u, r =portico.solve()

portico.show_model()
portico.show_deformed(1000)
# portico.show_deformed(1000)
import pandas as pd
disp = portico.nodes[2].displacements
print(pd.DataFrame(disp).round(6))