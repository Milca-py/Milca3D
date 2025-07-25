from milcapy.geom_transf import LinearCrdTransf3D
from milcapy.node import Node
import numpy as np


node1 = Node(1, np.array([0, 0, 0]))
node2 = Node(2, np.array([0, 0, 4]))
node3 = Node(3, np.array([4, 0, 4]))
node4 = Node(4, np.array([0, -4, 4]))
node5 = Node(5, np.array([1, 1, 0]))

t_col = LinearCrdTransf3D(tag=1, v_xz=np.array([1, 0, 0]))
t_col.connect(node1, node2)


t_vigaY = LinearCrdTransf3D(tag=2, v_xz=np.array([0, 0, 1]))
t_vigaY.connect(node2, node3)


t_vigaX = LinearCrdTransf3D(tag=3, v_xz=np.array([0, 0, 1]))
t_vigaX.connect(node2, node4)

vigaXY = LinearCrdTransf3D(tag=4, v_xz=np.array([0, 0, 1]))
vigaXY.connect(node1, node3)

print(t_col.get_rotation_matrix())
print(t_vigaY.get_rotation_matrix())
print(t_vigaX.get_rotation_matrix())
print(vigaXY.get_rotation_matrix())