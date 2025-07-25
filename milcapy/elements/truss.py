from milcapy.types import ElementModelType
from milcapy.core.node import Node

from numpy.typing import NDArray
import numpy as np

from milcapy.elements.frame import Element

class ElasticTruss3D(Element):
    """
    Clase que representa un elemento de trus en 3D.
    """
    def __init__(
        self,
        tag: int,
        node_i: Node,
        node_j: Node,
        E: float,
        A: float,
    ):
        super().__init__(tag, ElementModelType.ElasticTruss3d)

        self.node_i: Node = node_i
        self.node_j: Node = node_j

        # Propiedades del material y sección
        self.E: float = E               # Módulo de elasticidad (Young)
        self.A: float = A               # Área de la sección transversal

        # Propiedades geometricas
        self.length: float | None = None   # Longitud del elemento

        # Vectores y matrices
        # Resultados ==========================
        self.ul:  NDArray | None = None   # Vector de desplazamientos en sistema local
        self.ql:  NDArray | None = None   # Vector de fuerzas en sistema local

        # Precondicionamiento =================
        self.ql0: NDArray | None = None   # Vector de fuerzas fijas en sistema local
        self.kl:  NDArray | None = None   # Matriz de rigidez en sistema local

        self.Tlg: NDArray | None = None   # Matriz de transformación local a global
        self.Ki:  NDArray | None = None   # Matriz de rigidez inicial en sistema global

        # Ensamblaje =========================
        self.dofs: NDArray | None = None   # [ux, uy, uz, rx, ry, rz] * 2

    def set_up(self):
        """Configura las matrices y parámetros iniciales del elemento."""

        self.dofs = np.concatenate((self.node_i.dofs[0:3], self.node_j.dofs[0:3]))

        self.length = np.linalg.norm(self.node_j.coords - self.node_i.coords)

        self.Tlg = self.transform_matrix()

        self.kl = self.E*self.A/self.length*np.array([[1, -1], [-1, 1]])

        self.Ki = self.Tlg.T @ self.kl @ self.Tlg
        # # rellenar con ceros los grados de libertd de rotacion [rx, ry, rz]
        # self.Ki = np.zeros((12, 12))
        # self.Ki[0:3, 0:3] = Ki[0:3, 0:3]
        # self.Ki[6:9, 0:3] = Ki[3:6, 0:3]
        # self.Ki[0:3, 6:9] = Ki[0:3, 3:6]
        # self.Ki[6:9, 6:9] = Ki[3:6, 3:6]
        # print(self.Ki)

    def transform_matrix(self):
        L = self.length
        cx = (self.node_j.coords[0] - self.node_i.coords[0]) / L
        cy = (self.node_j.coords[1] - self.node_i.coords[1]) / L
        cz = (self.node_j.coords[2] - self.node_i.coords[2]) / L
        T = np.array([
            [cx, cy, cz, 0, 0, 0],
            [0, 0, 0, cx, cy, cz],
        ])
        return T

    def get_global_stiffness_matrix(self) -> NDArray:
        return self.Ki

    def get_global_load_vector(self) -> None:
        return None

    def set_load(self, load: 'DistributedLoad'):
        Warning("No se puede aplicar cargas a elemetos tipo Trus")
        print("La carga aplicada de despreciara para este analisis")

    def add_load(self, load: 'DistributedLoad'):
        """Agrega cargas aplicadas al elemento."""
        Warning("No se puede aplicar cargas a elemetos tipo Trus")
        print("La carga aplicada de despreciara para este analisis")