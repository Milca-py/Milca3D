from enum import Enum, auto
from typing import Tuple, Union

# Tipo para restricciones de 6 grados de libertad: (UX, UY, UZ, RX, RY, RZ)
RestraintType = Tuple[bool, bool, bool, bool, bool, bool]

# Etiqueta de un nodo o elemento: puede ser número o string (flexible para el usuario)
TagType = Union[int, float, str]


class ElementModelType(Enum):
    """
    Tipos de modelos de elementos estructurales admitidos por la librería FEM.
    """
    ElasticTimoshenkoBeam3d = auto()
    ElasticTruss3d = auto()


class DirectionType(Enum):
    """
    Direcciones espaciales principales, para identificar ejes, cargas o desplazamientos.
    """
    X = auto()
    Y = auto()
    Z = auto()
