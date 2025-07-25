from numpy import array, sqrt, zeros, concatenate, vstack, float64
from numpy.linalg import inv, det
from numpy.typing import NDArray
from shells.integration import NumericalIntegration2D

class Node:
    _ID = 1
    def __init__(self, tag: int, x: float, y: float) -> None:
        self.id = Node._ID
        Node._ID += 1
        self.tag = tag
        self.x = x
        self.y = y
        self.coords = array([x, y], dtype=float64)
        self.forces = array([0, 0, 0], dtype=float64)
        self.restraints = (False, False, False)
        self.dofs = array([
            self.id*3-2,
            self.id*3-1,
            self.id*3,
        ])

        self.reactions = array([0, 0, 0], dtype=float64)
        self.displacements = array([0, 0, 0], dtype=float64)

    def add_load(self, fz: float, mx: float, my: float) -> None:
        self.forces += array([fz, mx, my], dtype=float64)

    def set_restraints(self, uz: bool, rx: bool, ry: bool) -> None:
        self.restraints = (uz, rx, ry)

    def get_load_vector(self) -> NDArray[float64]:
        return self.forces


class PlateQ6:
    def __init__(self, tag: int, node1: Node, node2: Node, node3: Node, node4: Node, E: float, v: float, t: float) -> None:
        self.tag = tag
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.E = E          # Modulo de elasticidad
        self.v = v          # Cueficiente de poisson
        self.t = t          # Espesor de la placa
        self.dofs = concatenate((node1.dofs, node2.dofs, node3.dofs, node4.dofs))

    def D(self) -> NDArray[float64]:
        v = self.v
        t = self.t
        k = self.E * t**3 / (12 * (1 - v**2))
        D = k * array([
                        [1, v, 0,         0,                0],
                        [v, 1, 0,         0,                0],
                        [0, 0, (1 - v)/2, 0,                0],
                        [0, 0, 0,         5/(t**2*(1 + v)), 0],
                        [0, 0, 0,         0,                5/(t**2*(1 + v))]
                    ])
        return D

    def J(self, xi, eta):
        (x1, x2, x3, x4), (y1, y2, y3, y4) = self.__coords()
        J = array([
            [
                -x1*(1 - eta)/4 + x2*(1 - eta)/4 + x3*(eta + 1)/4 - x4*(eta + 1)/4, # diff(x, xi)
                -y1*(1 - eta)/4 + y2*(1 - eta)/4 + y3*(eta + 1)/4 - y4*(eta + 1)/4  # diff(y, xi)
            ],
            [
                -x1*(1 - xi)/4 - x2*(xi + 1)/4 + x3*(xi + 1)/4 + x4*(1 - xi)/4,     # diff(x, eta)
                -y1*(1 - xi)/4 - y2*(xi + 1)/4 + y3*(xi + 1)/4 + y4*(1 - xi)/4      # diff(y, eta)
            ]
        ])
        return J

    def __Mpx(self, xi, eta):
        _, _, sina = self.__parameters()
        s1, s2, s3, s4 = sina
        Mpx = array([[-s1*xi*(1 - eta), s2*(1 - eta**2)/2, -s3*xi*(eta + 1), -s4*(1 - eta**2)/2],
                     [-s1*(1 - xi**2)/2, -eta*s2*(xi + 1), s3*(1 - xi**2)/2, -eta*s4*(1 - xi)]])
        return Mpx

    def __Mpy(self, xi, eta):
        _, cosa, _ = self.__parameters()
        c1, c2, c3, c4 = cosa
        Mpy = array([[c1*xi*(1 - eta), -c2*(1 - eta**2)/2, c3*xi*(eta + 1), c4*(1 - eta**2)/2],
                     [c1*(1 - xi**2)/2, c2*eta*(xi + 1), -c3*(1 - xi**2)/2, c4*eta*(1 - xi)]])
        return Mpy

    def B(self, xi, eta):
        J = self.J(xi, eta)
        (x1, x2, x3, x4), (y1, y2, y3, y4) = self.__coords()
        L, cosa, sina = self.__parameters()
        L1, L2, L3, L4 = L
        cosa1, cosa2, cosa3, cosa4 = cosa
        sina1, sina2, sina3, sina4 = sina

        Mpx = self.__Mpx(xi, eta)
        Mpy = self.__Mpy(xi, eta)
        dMx = inv(J) @ Mpx
        dMy = inv(J) @ Mpy
        dNN = inv(J) @ self.__dN(xi, eta)

        B = zeros((3, 16))
        for i in range(4):
            B[0, 4*(i+1) - 3] = dNN[0, i]
            B[0, 4*(i+1)-1]   = dMy[0, i]
            B[1, 4*(i+1) - 4] = -dNN[1, i]
            B[1, 4*(i+1)-1]   = -dMx[1, i]
            B[2, 4*(i+1) - 3] = dNN[1, i]
            B[2, 4*(i+1) - 4] = -dNN[0, i]
            B[2, 4*(i+1)-1]   = dMy[1, i] - dMx[0, i]

        V1 = array([[-sina1/2, cosa1/2, -1/L1, -2/3],
                    [-sina2/2, cosa2/2, -1/L2, -2/3],
                    [-sina3/2, cosa3/2, -1/L3, -2/3],
                    [-sina4/2, cosa4/2, -1/L4, -2/3]])
        V2 = array([[-sina1/2, cosa1/2, 1/L1, 0],
                    [-sina2/2, cosa2/2, 1/L2, 0],
                    [-sina3/2, cosa3/2, 1/L3, 0],
                    [-sina4/2, cosa4/2, 1/L4, 0]])

        M = [0, 0, 0, 0]
        for i in range(4):
            N = zeros((2, 16))
            if i == 0:
                for j in range(4):
                    N[0, j] = V1[i, j]
                for j in range(4):
                    N[0, j + 4] = V2[i, j]
                for j in range(4):
                    N[1, j] = V2[3, j]
                for j in range(4):
                    N[1, j + 12] = V1[3, j]
                CS = array([
                    [cosa1, sina1],
                    [cosa4, sina4]
                ])
                N = inv(CS) @ N
            if i == 1:
                for j in range(4):
                    N[0, j + 4] = V1[i, j]
                for j in range(4):
                    N[0, j + 8] = V2[i, j]
                for j in range(4):
                    N[1, j] = V1[i-1, j]
                for j in range(4):
                    N[1, j + 4] = V2[i-1, j]
                CS = array([
                    [cosa2, sina2],
                    [cosa1, sina1]
                ])
                N = inv(CS) @ N
            if i == 2:
                for j in range(4):
                    N[0, j + 8] = V1[i, j]
                for j in range(4):
                    N[0, j + 12] = V2[i, j]
                for j in range(4):
                    N[1, j + 4] = V1[i-1, j]
                for j in range(4):
                    N[1, j + 8] = V2[i-1, j]
                CS = array([
                    [cosa3, sina3],
                    [cosa2, sina2]
                ])
                N = inv(CS) @ N
            if i == 3:
                for j in range(4):
                    N[0, j + 12] = V1[i, j]
                for j in range(4):
                    N[0, j] = V2[i, j]
                for j in range(4):
                    N[1, j + 8] = V1[i-1, j]
                for j in range(4):
                    N[1, j + 12] = V2[i-1, j]
                CS = array([
                    [cosa4, sina4],
                    [cosa3, sina3]
                ])
                N = inv(CS) @ N
            M[i] = N

        NShape = self.__Nshape(xi, eta)
        Mxyz = zeros((2, 16))
        for i in range(4):
            Mxyz += M[i]*NShape[i]

        Bt = vstack((B, Mxyz))

        return Bt

    def Kgauss(self):
        r3 = 3**0.5
        xii = [-1/r3, -1/r3, 1/r3, 1/r3]
        yii = [-1/r3, 1/r3, -1/r3, 1/r3]
        w = [1, 1, 1, 1]
        K = zeros((16,16), dtype=float64)
        for i in range(4):
            B = self.B(xii[i], yii[i])
            J = self.J(xii[i], yii[i])
            D = self.D()
            W = w[i]
            detJ = det(J)

            phi = B.T @ D @ B * detJ * W
            K += phi

        return K

    def Kwilson(self):
        integr = NumericalIntegration2D(rule='5-point', params={'w0': 224/81})
        D = self.D()

        phi = lambda xi, eta: self.B(xi, eta).T @ D @ self.B(xi, eta) * det(self.J(xi, eta))

        return integr.integrate(phi)

    def static_condensation(self):
        Ke = self.Kgauss() # Cambiar el tipo de integracion
        # indices = [0, 4, 8, 12]
        indices = [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]
        nonindi = [3, 7, 11, 15]
        k11 = Ke[indices, :][:, indices]
        k12 = Ke[indices, :][:, nonindi]
        k21 = Ke[nonindi, :][:, indices]
        k22 = Ke[nonindi, :][:, nonindi]
        Ke_red = k11 - k12 @ inv(k22) @ k21
        return Ke_red

    def K(self):
        return self.static_condensation()

    def __Nshape(self, xi, eta):
        return array([
            [(1 - eta)*(1 - xi)/4],
            [(1 - eta)*(xi + 1)/4],
            [(eta + 1)*(xi + 1)/4],
            [(1 - xi)*(eta + 1)/4],
            [(1 - eta)*(1 - xi**2)/2],
            [(1 - eta**2)*(xi + 1)/2],
            [(1 - xi**2)*(eta + 1)/2],
            [(1 - eta**2)*(1 - xi)/2]
        ])

    def __parameters(self):
        (x1, x2, x3, x4), (y1, y2, y3, y4) = self.__coords()
        L1 = sqrt((x2-x1)**2+(y2-y1)**2)
        L2 = sqrt((x3-x2)**2+(y3-y2)**2)
        L3 = sqrt((x4-x3)**2+(y4-y3)**2)
        L4 = sqrt((x1-x4)**2+(y1-y4)**2)
        L = array([L1, L2, L3, L4])

        cosa1 = (x2-x1)/L1
        cosa2 = (x3-x2)/L2
        cosa3 = (x4-x3)/L3
        cosa4 = (x1-x4)/L4
        cosa = array([cosa1, cosa2, cosa3, cosa4])

        sina1 = (y2-y1)/L1
        sina2 = (y3-y2)/L2
        sina3 = (y4-y3)/L3
        sina4 = (y1-y4)/L4
        sina = array([sina1, sina2, sina3, sina4])

        return L, cosa, sina

    def __coords(self):
        x1, y1 = self.node1.coords
        x2, y2 = self.node2.coords
        x3, y3 = self.node3.coords
        x4, y4 = self.node4.coords
        return (x1, x2, x3, x4), (y1, y2, y3, y4)

    def __dN(self, xi, eta):
        return array([
            [
                eta/4 - 1/4, 1/4 - eta/4,
                eta/4 + 1/4, -eta/4 - 1/4
            ],
            [
                xi/4 - 1/4, -xi/4 - 1/4,
                xi/4 + 1/4, 1/4 - xi/4
            ]
        ])


if __name__ == "__main__":
    import pandas as pd
    import time
    import numpy as np
    from numpy import linalg
    t_ini = time.time()
    pd.set_option('display.float_format', '{:.5f}'.format)
    # Units: tonf, m
    E = 2.5356  # tonf/m2
    v = 0.2
    t = 250  # m
    nd1 = Node(1, 0, 0)
    nd2 = Node(2, 3000, 0)
    nd3 = Node(3, 3000, 3000)
    nd4 = Node(4, 0, 3000)
    el = PlateQ6(1, nd1, nd2, nd3, nd4, E, v, t)
    K = el.K()
    k_red = K[6:12, 6:12]
    Q = array([0, 0, 3, 0, 0, 3])
    u = linalg.solve(k_red, Q)
    t_end = time.time()
    print("El analisis duro: ", t_end-t_ini)
    # print(pd.DataFrame(K))
    print(pd.DataFrame(u))