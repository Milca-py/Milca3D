from typing import ClassVar
import numpy as np
from numpy.typing import NDArray
from milcapy.types import RestraintType, TagType
from milcapy.core.load import PointLoad


class Node:
    """
    Representa un nodo 3D con 6 grados de libertad (DOF) para análisis FEM.
    """

    _ID: ClassVar[int] = 1
    _DEFAULT_RESTRAINT: ClassVar[RestraintType] = (False, False, False, False, False, False)

    def __init__(self, tag: TagType, coords: NDArray[np.float64]) -> None:
        """
        Inicializa un nodo con coordenadas y etiqueta.

        :param tag: Identificador del nodo (str o int).
        :param coords: Coordenadas 3D del nodo como array [x, y, z].
        """
        self.id: int = Node._ID
        Node._ID += 1

        self.tag:    TagType = tag
        self.coords: NDArray[np.float64] = np.asarray(coords, dtype=np.float64)

        # Grados de libertad (6 por nodo)
        self.dofs: NDArray[np.int32] = np.array([
            self.id * 6 - 5, self.id * 6 - 4, self.id * 6 - 3,
            self.id * 6 - 2, self.id * 6 - 1, self.id * 6
        ], dtype=np.int32)

        # Restricciones y cargas
        self.restraint: RestraintType = Node._DEFAULT_RESTRAINT
        self.the_load:  PointLoad = PointLoad()

        # Resultados del análisis (se asignan luego)
        self.displacements: NDArray[np.float64] | None = None
        self.reactions:     NDArray[np.float64] | None = None

    def set_restraint(self, restraint: RestraintType) -> None:
        """
        Establece o acumula las restricciones del nodo.
        Las nuevas restricciones se combinan con las existentes (OR lógico).

        :param restraint: Tupla de 6 booleanos indicando restricciones en UX, UY, UZ, RX, RY, RZ.
        """
        self.restraint = tuple(
            r_old or r_new for r_old, r_new in zip(self.restraint, restraint)
        )

    def set_load(self, load: PointLoad) -> None:
        """
        Reemplaza la carga actual aplicada al nodo.

        :param load: Objeto PointLoad con vector de fuerzas y momentos.
        """
        self.the_load = load

    def add_load(self, load: PointLoad) -> None:
        """
        Suma una carga adicional al nodo.

        :param load: Carga adicional como objeto PointLoad.
        """
        self.the_load += load

    def get_load_vector(self) -> NDArray[np.float64]:
        """
        Retorna el vector de fuerzas y momentos aplicado al nodo.

        :return: Vector numpy de 6 componentes.
        """
        return self.the_load.vector

    def __repr__(self) -> str:
        return f"Node(tag={self.tag}, coords={self.coords.tolist()})"
