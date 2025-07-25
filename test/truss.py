from milcapy.model import Model

cercha = Model(3)
cercha.set_dofs(False, False, False, True, True, True)

cercha.add_node(1, 0, 0, 0)
cercha.add_node(2, 0, 0, 4)
cercha.add_node(3, 4, 0, 4)
cercha.add_node(4, 4, 0, 0)
cercha.add_node(5, 0, 4, 0)
cercha.add_node(6, 0, 4, 4)
cercha.add_node(7, 4, 4, 4)
cercha.add_node(8, 4, 4, 0)
cercha.add_node(9, 0, 8, 0)
cercha.add_node(10, 0, 8, 4)
cercha.add_node(11, 4, 8, 4)
cercha.add_node(12, 4, 8, 0)


fixed = (True, True, True, True, True, True)
for i in [1, 2, 3, 4]:
    cercha.add_restraint(i, fixed)


cercha.add_truss_member(1, 1, 2, 2.1e6, 0.5*0.3)
cercha.add_truss_member(2, 2, 3, 2.1e6, 0.5*0.3)
cercha.add_truss_member(3, 3, 4, 2.1e6, 0.5*0.3)
cercha.add_truss_member(4, 4, 1, 2.1e6, 0.5*0.3)

cercha.add_truss_member(5, 5, 6, 2.1e6, 0.5*0.3)
cercha.add_truss_member(6, 6, 7, 2.1e6, 0.5*0.3)
cercha.add_truss_member(7, 7, 8, 2.1e6, 0.5*0.3)
cercha.add_truss_member(8, 8, 5, 2.1e6, 0.5*0.3)

cercha.add_truss_member(9, 9, 10, 2.1e6, 0.5*0.3)
cercha.add_truss_member(10, 10, 11, 2.1e6, 0.5*0.3)
cercha.add_truss_member(11, 11, 12, 2.1e6, 0.5*0.3)
cercha.add_truss_member(12, 12, 9, 2.1e6, 0.5*0.3)

cercha.add_truss_member(13, 1, 5, 2.1e6, 0.5*0.3)
cercha.add_truss_member(14, 2, 6, 2.1e6, 0.5*0.3)
cercha.add_truss_member(15, 3, 7, 2.1e6, 0.5*0.3)
cercha.add_truss_member(16, 4, 8, 2.1e6, 0.5*0.3)

cercha.add_truss_member(17, 5, 9, 2.1e6, 0.5*0.3)
cercha.add_truss_member(18, 6, 10, 2.1e6, 0.5*0.3)
cercha.add_truss_member(19, 7, 11, 2.1e6, 0.5*0.3)
cercha.add_truss_member(20, 8, 12, 2.1e6, 0.5*0.3)

cercha.add_truss_member(21, 1, 6, 2.1e6, 0.5*0.3)
cercha.add_truss_member(22, 4, 7, 2.1e6, 0.5*0.3)
cercha.add_truss_member(23, 6, 9, 2.1e6, 0.5*0.3)
cercha.add_truss_member(24, 7, 12, 2.1e6, 0.5*0.3)

cercha.add_point_load(10, fz=-10)
cercha.add_point_load(11, fz=-10)

cercha.solve()

cercha.show_model()
cercha.show_deformed(0.01, 100)

cercha.print.node_results()



