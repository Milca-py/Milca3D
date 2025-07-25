from milcapy.model import Model

cercha = Model(2)
cercha.set_dofs(False, False, False, True, True, True)

cercha.add_node(1, 0, 0, 0)
cercha.add_node(2, 4, 0, 0)
cercha.add_node(3, 4, 4, 0)
cercha.add_node(4, 0, 4, 0)
cercha.add_node(5, 2, 2, 4)


fixed = (True, True, True, True, True, True)
for i in [1, 2, 3, 4]:
    cercha.add_restraint(i, fixed)


cercha.add_truss_member(1, 1, 5, 2.1e6, 0.5*0.3)
cercha.add_truss_member(2, 2, 5, 2.1e6, 0.5*0.3)
cercha.add_truss_member(3, 3, 5, 2.1e6, 0.5*0.3)
cercha.add_truss_member(4, 4, 5, 2.1e6, 0.5*0.3)

cercha.add_point_load(5, fz=-10)

cercha.solve()

cercha.show_model()
cercha.show_deformed(0.01, 100)

cercha.print.node_results()



