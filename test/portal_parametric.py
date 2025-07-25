from milcapy.model import Model, PortalParametric

# Parametros:
m = 1   # numero de bahias en la direccion X
n = 1   # numero de bahias en la direccion Y
p = 4   # numero de pisos de la estructura (Z)
l = 7   # ancho de las bahias en la direccion X
b = 7   # ancho de las bahias en la direccion Y
h = 4   # altura de los pisos en la direccion Z

# Materiales y secciones
hr = 0.5
br = 0.3
E = 2.1e6
v = 0.2
G = E / (2 * (1 + v))
A = hr * br
Jx = 0.002822 # hr*br**3*(16/3 - 3.36*br/hr*(1-(br/hr)**4/12))/16
Iy = hr * br**3 / 12
Iz = br * hr**3 / 12
Asy = 5/6 * A
Asz = 5/6 * A


portico = Model(3)

data = PortalParametric(m, n, p, l, b, h)
nodes = data.nodes()
members = data.members()
restrained_nodes = data.nodes_from_story(0)
forces = data.nodes_from_story(p)


for i, coord in nodes.items():
    portico.add_node(i, *coord)


portico.def_geom_transf(1, *[1, 0, 0]) # columnas
portico.def_geom_transf(2, *[0, 0, 1]) # vigas
for i, (node_i, node_j) in members["columna"].items():
    portico.add_member(i, node_i, node_j, E, G, A, Jx, Iy, Iz, Asy, Asz, 1)
for i, (node_i, node_j) in members["viga"].items():
    portico.add_member(i, node_i, node_j, E, G, A, Jx, Iy, Iz, Asy, Asz, 2)


fixed = (True, True, True, True, True, True)
for i in restrained_nodes:
    portico.add_restraint(i, fixed)


for mbr in members["viga"].keys():
    portico.add_distributed_load(mbr, direction=1, wi=-10, wj=-10)
for mbr in members["viga"].keys():
    portico.add_distributed_load(mbr, direction=2, wi=-100, wj=-100)
for mbr in members["viga"].keys():
    portico.add_distributed_load(mbr, direction=3, wi=-10, wj=-10)


# portico.mesh(d_max=1)
portico.solve()

print("Analisis completado")

portico.show_model()
portico.show_deformed(1)
