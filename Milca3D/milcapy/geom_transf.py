from typing import Optional, Tuple
import numpy as np
from numpy.typing import NDArray

from milcapy.core.node import Node


class LinearCrdTransf3D:
    """
    Transformador lineal de coordenadas para elementos en 3D.
    Convierte entre sistemas de coordenadas global y local.
    """

    def __init__(self, tag: int, v_xz: NDArray[np.float64]) -> None:
        """
        Inicializa el transformador de coordenadas.

        :param tag: Identificador del sistema de transformación.
        :param v_xz: Vector auxiliar que define el plano local XZ (no paralelo al eje del elemento).
        """
        self.tag: int = tag
        self.v_xz: NDArray[np.float64] = np.asarray(v_xz, dtype=np.float64)

        self.node_i: Optional[Node] = None
        self.node_j: Optional[Node] = None

    def connect(self, node_i: Node, node_j: Node) -> None:
        """
        Conecta dos nodos al transformador.

        :param node_i: Nodo inicial del elemento.
        :param node_j: Nodo final del elemento.
        """
        self.node_i = node_i
        self.node_j = node_j

    @property
    def length(self) -> float:
        """
        Longitud del elemento.
        """
        self._validate_connection()
        return np.linalg.norm(self.node_j.coords - self.node_i.coords)

    @property
    def local_axes(self) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """
        Ejes locales X', Y', Z' del sistema de coordenadas del elemento.

        :return: Tupla con los vectores unitarios (x', y', z').
        """
        self._validate_connection()

        x_axis = (self.node_j.coords - self.node_i.coords) / self.length
        z_axis = np.cross(x_axis, self.v_xz)
        z_axis /= np.linalg.norm(z_axis)
        y_axis = np.cross(z_axis, x_axis)

        return x_axis, y_axis, z_axis

    @property
    def rotation_matrix(self) -> NDArray[np.float64]:
        """
        Matriz de rotación de coordenadas locales a globales (3x3).

        :return: Matriz de rotación R.
        """
        x_axis, y_axis, z_axis = self.local_axes
        return np.vstack((x_axis, y_axis, z_axis))

    @property
    def transformation_matrix(self) -> NDArray[np.float64]:
        """
        Matriz de transformación (rotación ampliada a 12x12) de coordenadas locales a globales.

        :return: Matriz de transformación T.
        """
        R = self.rotation_matrix
        T = np.zeros((12, 12), dtype=np.float64)
        for i in range(4):
            T[3 * i:3 * i + 3, 3 * i:3 * i + 3] = R
        return T

    def _validate_connection(self) -> None:
        """
        Lanza un error si los nodos aún no han sido conectados.
        """
        if self.node_i is None or self.node_j is None:
            raise ValueError("Los nodos aún no han sido conectados al transformador.")
