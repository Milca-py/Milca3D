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



class Quad6:
    def __init__(self, tag: int, node1: Node, node2: Node, node3: Node, node4: Node, E: float, v: float, t: float) -> None:
        self.tag = tag
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.E = E
        self.v = v
        self.t = t
        self.dofs = concatenate(
            (node1.dofs, node2.dofs, node3.dofs, node4.dofs))

        r3 = sqrt(3)
        self.xi = [-1/r3, -1/r3, 1/r3, 1/r3]
        self.eta = [-1/r3, 1/r3, -1/r3, 1/r3]
        self.w = [1, 1, 1, 1]

    def set_gauss_points(self, xi: list[float], eta: list[float], w: list[float]):
        """Setea los puntos de gauss"""
        self.xi = xi
        self.eta = eta
        self.w = w

    def shape_function(self, xi: float, eta: float):
        """Funciones de forma"""
        N1 = (1-xi)*(1-eta)/4
        N2 = (1+xi)*(1-eta)/4
        N3 = (1+xi)*(1+eta)/4
        N4 = (1-xi)*(1+eta)/4
        N5 = (1-xi**2)*(1-eta)/2
        N6 = (1-eta**2)*(1+xi)/2
        N7 = (1-xi**2)*(1+eta)/2
        N8 = (1-eta**2)*(1-xi)/2
        return N1, N2, N3, N4, N5, N6, N7, N8

    def coordinates(self, xi: float, eta: float):
        """Interpolacion de coordenadas del elemento"""
        N1, N2, N3, N4, *_ = self.shape_function(xi, eta)
        x1, y1 = self.node1.coords
        x2, y2 = self.node2.coords
        x3, y3 = self.node3.coords
        x4, y4 = self.node4.coords
        x = N1*x1 + N2*x2 + N3*x3 + N4*x4
        y = N1*y1 + N2*y2 + N3*y3 + N4*y4
        return x, y

    def Ia(self) -> np.ndarray:
        """Matriz de unitaria ampliada"""
        Ia = np.array([
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 1, 1, 0]
        ])
        return Ia

    def J(self, xi: float, eta: float) -> np.ndarray:
        """Matriz de jacobiano"""
        x1, y1 = self.node1.coords
        x2, y2 = self.node2.coords
        x3, y3 = self.node3.coords
        x4, y4 = self.node4.coords
        J11 = -x1*(1-eta) + x2*(1-eta) + x3*(1+eta) - x4*(1+eta)
        J12 = -y1*(1-eta) + y2*(1-eta) + y3*(1+eta) - y4*(1+eta)
        J21 = -x1*(1-xi) - x2*(1+xi) + x3*(1+xi) + x4*(1-xi)
        J22 = -y1*(1-xi) - y2*(1+xi) + y3*(1+xi) + y4*(1-xi)
        return np.array([[J11, J12], [J21, J22]]) * (1/4)

    def Gamma(self, xi: float, eta: float) -> np.ndarray:
        """Matriz de jacobiano invertida"""
        J = self.J(xi, eta)
        gamma = np.linalg.inv(J)
        zero = np.zeros((2, 2))
        Gamma = np.concatenate((
            np.concatenate((gamma, zero), axis=1),
            np.concatenate((zero, gamma), axis=1)
        ),
            axis=0)
        return Gamma

    def dNdisp(self, xi: float, eta: float) -> np.ndarray:
        """Derivada de las funciones de forma (dof: desplazamientos)"""
        dN = (1/4)*np.array([
            [eta-1, 0, 1-eta, 0, 1+eta, 0, -1-eta, 0],
            [xi-1,  0, -1-xi, 0, 1+xi,   0,  1-xi, 0],
            [0,     eta-1, 0, 1-eta, 0, 1+eta, 0, -1-eta],
            [0,     xi-1,  0, -1-xi, 0, 1+xi,   0,  1-xi]
        ])
        return dN

    def Bdisp(self, xi: float, eta: float) -> np.ndarray:
        """Matriz de deformacion de desplazamientos"""
        dN = self.dNdisp(xi, eta)
        Gamma = self.Gamma(xi, eta)
        Ia = self.Ia()
        B = Ia @ Gamma @ dN
        return B

    def dNrot(self, xi: float, eta: float) -> np.ndarray:
        """Derivada de las funciones de forma modificadas (dof: rotaciones)"""
        x1, y1 = self.node1.coords
        x2, y2 = self.node2.coords
        x3, y3 = self.node3.coords
        x4, y4 = self.node4.coords

        dN = (1/16)*np.array([
            [
                2 * xi * (eta - 1) * (y1 - y2)
                + (eta**2 - 1) * (y1 - y4),

                -2 * xi * (eta - 1) * (y1 - y2)
                - (eta**2 - 1) * (y2 - y3),

                -2 * xi * (eta + 1) * (y3 - y4)
                + (eta**2 - 1) * (y2 - y3),

                2 * xi * (eta + 1) * (y3 - y4)
                - (eta**2 - 1) * (y1 - y4)
            ],
            [
                2 * eta * (xi - 1) * (y1 - y4)
                + (xi**2 - 1) * (y1 - y2),

                -2 * eta * (xi + 1) * (y2 - y3)
                - (xi**2 - 1) * (y1 - y2),

                2 * eta * (xi + 1) * (y2 - y3)
                - (xi**2 - 1) * (y3 - y4),

                -2 * eta * (xi - 1) * (y1 - y4)
                + (xi**2 - 1) * (y3 - y4)
            ],
            [
                -2 * xi * (eta - 1) * (x1 - x2)
                - (eta**2 - 1) * (x1 - x4),

                2 * xi * (eta - 1) * (x1 - x2)
                + (eta**2 - 1) * (x2 - x3),

                2 * xi * (eta + 1) * (x3 - x4)
                - (eta**2 - 1) * (x2 - x3),

                -2 * xi * (eta + 1) * (x3 - x4)
                + (eta**2 - 1) * (x1 - x4)
            ],
            [
                -2 * eta * (x1 - x4) * (xi - 1)
                - (x1 - x2) * (xi**2 - 1),

                2 * eta * (x2 - x3) * (xi + 1)
                + (x1 - x2) * (xi**2 - 1),

                -2 * eta * (x2 - x3) * (xi + 1)
                + (x3 - x4) * (xi**2 - 1),

                2 * eta * (x1 - x4) * (xi - 1)
                - (x3 - x4) * (xi**2 - 1)
            ]
        ])
        return dN

    def Brot(self, xi: float, eta: float) -> np.ndarray:
        """Matriz de deformacion de rotaciones"""
        dN = self.dNrot(xi, eta)
        Gamma = self.Gamma(xi, eta)
        Ia = self.Ia()
        B = Ia @ Gamma @ dN
        return B

    def B(self, xi: float, eta: float) -> np.ndarray:
        """Matriz de deformacion"""
        B = np.hstack((self.Bdisp(xi, eta), self.Brot(xi, eta)))

        # # order = [0, 1, 8, 2, 3, 9, 4, 5, 10, 6, 7, 11]
        # order = [1, 0, 8, 3, 2, 9, 5, 4, 10, 7, 6, 11]

        # B_new = B[:, order]

        return B

    def D(self):
        """Matriz constitutiva de esfuerzo plano"""
        return (self.E/(1-self.v**2))*np.array([
            [1, self.v, 0],
            [self.v, 1, 0],
            [0, 0, (1-self.v)/2],
        ])

    def phiKi(self, xi: float, eta: float):
        """Matriz de rigidez local"""
        detJ = np.linalg.det(self.J(xi, eta))
        B = self.B(xi, eta)
        D = self.D()
        t = self.t
        return B.T @ D @ B * detJ * t

    def Ki(self):
        """Calcula la matriz rigidez atravez de la intergracion numerica cuadratura de gauss"""
        Ki = np.zeros((12, 12))
        for i in range(len(self.xi)):
            Ki += self.phiKi(self.xi[i], self.eta[i]) * self.w[i]


        # # Identidad de 12x12
        # P = np.eye(12)

        # # Nuevo orden
        # order = [1, 0, 8, 3, 2, 9, 5, 4, 10, 7, 6, 11]

        # # Matriz de permutación
        # P_new = P[:, order]

        # # Reordenar
        # Ki = P_new @ Ki @ P_new.T


        return Ki


if __name__ == "__main__":
    import pandas as pd
    pd.set_option('display.float_format', '{:.2f}'.format)
    E = 25  # MPa
    v = 0.20
    t = 400  # mm
    l = 2000  # mm
    h = 600  # mm
    nd1 = Node(1, 0, 0)
    nd2 = Node(2, h, 0)
    nd3 = Node(3, h, l)
    nd4 = Node(4, 0, l)
    el = Quad6(1, nd1, nd2, nd3, nd4, E, v, t)
    print(pd.DataFrame(el.Ki()))
