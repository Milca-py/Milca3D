from milcapy.model import Model
portico = Model()
# Materiales y secciones
h = 0.5
b = 0.3
E = 2.1e6
v = 0.2
G = E / (2 * (1 + v))
A = h * b
Iy = h * b**3 / 12
Iz = b * h**3 / 12
Jx = 0.002822 # h*b**3*(16/3 - 3.36*b/h*(1-(b/h)**4/12))/16
Asy = 5/6 * A
Asz = 5/6 * A

import numpy as np
np.set_printoptions(formatter={'float': '{: 0.6f}'.format})
portico.add_node(1, 0, 0, 0)
portico.add_node(2, 0, 4, 0)
portico.def_geom_transf(1, *[0, 0, 1])
portico.add_member(1, 1, 2, E, G, A, Jx, Iy, Iz, Asy, Asz, 1)
portico.add_restraint(1, (True, True, True, True, True, True))
# portico.add_restraint(2, (True, True, True, True, True, True))
portico.add_distributed_load(1, direction=3, wi=-10000, wj=-10000)
portico.add_distributed_load(1, direction=2, wi=-10000, wj=-10000)
portico.solve()

# portico.show_model()
portico.show_deformed(0.01)
portico.print.node_results()
