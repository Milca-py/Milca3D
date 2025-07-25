from typing import Union, Dict
import numpy as np


class PointLoad:
    __slots__ = ('fx', 'fy', 'fz', 'mx', 'my', 'mz')

    def __init__(self, fx: float = 0.0, fy: float = 0.0, fz: float = 0.0,
                mx: float = 0.0, my: float = 0.0, mz: float = 0.0) -> None:
        """
        Inicializa una nueva carga puntual.

        Args:
            fx: Fuerza en la direccion X GLOBAL.
            fy: Fuerza en la direccion Y GLOBAL.
            fz: Fuerza en la direccion Z GLOBAL.
            mx: Momento en la direccion X GLOBAL.
            my: Momento en la direccion Y GLOBAL.
            mz: Momento en la direccion Z GLOBAL.
        """
        self.fx: float = float(fx)
        self.fy: float = float(fy)
        self.fz: float = float(fz)
        self.mx: float = float(mx)
        self.my: float = float(my)
        self.mz: float = float(mz)

    @property
    def vector(self) -> np.ndarray:
        """
        Obtiene los componentes de la carga como un array NumPy.

        Returns:
            np.ndarray: Array con [fx, fy, fz, mx, my, mz].
        """
        return np.array([self.fx, self.fy, self.fz, self.mx, self.my, self.mz], dtype=np.float64)

    def to_dict(self) -> Dict[str, Union[int, float, None]]:
        """
        Convierte la carga puntual a un diccionario serializable.

        Returns:
            Dict[str, Union[int, float, None]]: Diccionario con los componentes de la carga.
        """
        return {"fx": self.fx, "fy": self.fy, "fz": self.fz, "mx": self.mx, "my": self.my, "mz": self.mz}

    def __add__(self, other: "PointLoad") -> "PointLoad":
        """
        Suma dos cargas puntuales.

        Args:
            other: Otra carga puntual a sumar.

        Returns:
            PointLoad: Nueva carga puntual resultante de la suma.

        Raises:
            TypeError: Si other no es una instancia de PointLoad.
        """
        if not isinstance(other, PointLoad):
            return NotImplemented
        return PointLoad(
            self.fx + other.fx,
            self.fy + other.fy,
            self.fz + other.fz,
            self.mx + other.mx,
            self.my + other.my,
            self.mz + other.mz,
        )

    def __sub__(self, other: "PointLoad") -> "PointLoad":
        """
        Resta dos cargas puntuales.

        Args:
            other: Carga puntual a restar.

        Returns:
            PointLoad: Nueva carga puntual resultante de la resta.

        Raises:
            TypeError: Si other no es una instancia de PointLoad.
        """
        if not isinstance(other, PointLoad):
            return NotImplemented
        return PointLoad(
            self.fx - other.fx,
            self.fy - other.fy,
            self.fz - other.fz,
            self.mx - other.mx,
            self.my - other.my,
            self.mz - other.mz,
        )

    def __mul__(self, scalar: Union[float, int]) -> "PointLoad":
        """
        Multiplica la carga por un escalar.

        Args:
            scalar: Factor de escala.

        Returns:
            PointLoad: Nueva carga puntual escalada.

        Raises:
            TypeError: Si scalar no es un número.
        """
        if not isinstance(scalar, (float, int)):
            return NotImplemented
        return PointLoad(
            self.fx * scalar,
            self.fy * scalar,
            self.fz * scalar,
            self.mx * scalar,
            self.my * scalar,
            self.mz * scalar,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: Union[float, int]) -> "PointLoad":
        """
        Divide la carga por un escalar.

        Args:
            scalar: Divisor.

        Returns:
            PointLoad: Nueva carga puntual dividida.

        Raises:
            TypeError: Si scalar no es un número.
            ZeroDivisionError: Si scalar es cero.
        """
        if not isinstance(scalar, (float, int)):
            return NotImplemented
        if scalar == 0:
            raise ZeroDivisionError("No se puede dividir por cero.")
        return PointLoad(
            self.fx / scalar,
            self.fy / scalar,
            self.fz / scalar,
            self.mx / scalar,
            self.my / scalar,
            self.mz / scalar,
        )


class FullLinearDistributedLoad:
    """
    Representa una carga distribuida lineal en un elemento estructural 2D.

    Esta clase modela cargas distribuidas con valores iniciales (i) y finales (j):
    - p: Carga distribuida axial a lo largo del eje de la viga en la direccion 1 LOCAL
    - q: Carga distribuida perpendicular al eje de la viga en la direccion 2 LOCAL
    - w: Carga distribuida perpendicular al eje de la viga en la direccion 3 LOCAL

    Attributes:
        p_i (float): Carga axial inicial, a lo largo del eje x de la viga.
        p_j (float): Carga axial final, a lo largo del eje x de la viga.
        q_i (float): Carga distribuida inicial, perpendicular al eje de la viga, eje y.
        q_j (float): Carga distribuida final, perpendicular al eje de la viga, eje y.
        w_i (float): Carga distribuida inicial, perpendicular al eje de la viga, eje z.
        w_j (float): Carga distribuida final, perpendicular al eje de la viga, eje z.
    """

    __slots__ = ('p_i', 'p_j', 'q_i', 'q_j', 'w_i', 'w_j')

    def __init__(self, p_i: float = 0.0, p_j: float = 0.0,
                q_i: float = 0.0, q_j: float = 0.0,
                w_i: float = 0.0, w_j: float = 0.0) -> None:

        self.p_i: float = p_i
        self.p_j: float = p_j
        self.q_i: float = q_i
        self.q_j: float = q_j
        self.w_i: float = w_i
        self.w_j: float = w_j

    @property
    def vector(self) -> np.ndarray:
        """
        Obtiene los componentes de la carga como un array NumPy.

        Returns:
            np.ndarray: Array con [p_i, p_j, q_i, q_j, w_i, w_j].
        """
        return np.array([self.p_i, self.p_j, self.q_i, self.q_j, self.w_i, self.w_j],
                       dtype=np.float64)

    def to_dict(self) -> Dict[str, Union[int, float, None]]:
        """
        Convierte la carga distribuida a un diccionario serializable.

        Returns:
            Dict[str, Union[int, float, None]]: Diccionario con los componentes de la carga.
        """
        return {
            "p_i": self.p_i, "p_j": self.p_j,
            "q_i": self.q_i, "q_j": self.q_j,
            "w_i": self.w_i, "w_j": self.w_j,
        }

    def __add__(self, other: "DistributedLoad") -> "DistributedLoad":
        """
        Suma dos cargas distribuidas.

        Args:
            other: Otra carga distribuida a sumar.

        Returns:
            DistributedLoad: Nueva carga distribuida resultante de la suma.

        Raises:
            TypeError: Si other no es una instancia de DistributedLoad.
        """
        if not isinstance(other, DistributedLoad):
            return NotImplemented
        return DistributedLoad(
            self.p_i + other.p_i,
            self.p_j + other.p_j,
            self.q_i + other.q_i,
            self.q_j + other.q_j,
            self.w_i + other.w_i,
            self.w_j + other.w_j,
        )

    def __sub__(self, other: "DistributedLoad") -> "DistributedLoad":
        """
        Resta dos cargas distribuidas.

        Args:
            other: Carga distribuida a restar.

        Returns:
            DistributedLoad: Nueva carga distribuida resultante de la resta.

        Raises:
            TypeError: Si other no es una instancia de DistributedLoad.
        """
        if not isinstance(other, DistributedLoad):
            return NotImplemented
        return DistributedLoad(
            self.p_i - other.p_i,
            self.p_j - other.p_j,
            self.q_i - other.q_i,
            self.q_j - other.q_j,
            self.w_i - other.w_i,
            self.w_j - other.w_j,
        )

    def __mul__(self, scalar: Union[float, int]) -> "DistributedLoad":
        """
        Multiplica la carga por un escalar.

        Args:
            scalar: Factor de escala.

        Returns:
            DistributedLoad: Nueva carga distribuida escalada.

        Raises:
            TypeError: Si scalar no es un número.
        """
        if not isinstance(scalar, (float, int)):
            return NotImplemented
        return DistributedLoad(
            self.p_i * scalar,
            self.p_j * scalar,
            self.q_i * scalar,
            self.q_j * scalar,
            self.w_i * scalar,
            self.w_j * scalar,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: Union[float, int]) -> "DistributedLoad":
        """
        Divide la carga por un escalar.

        Args:
            scalar: Divisor.

        Returns:
            DistributedLoad: Nueva carga distribuida dividida.

        Raises:
            TypeError: Si scalar no es un número.
            ZeroDivisionError: Si scalar es cero.
        """
        if not isinstance(scalar, (float, int)):
            return NotImplemented
        if scalar == 0:
            raise ZeroDivisionError("No se puede dividir por cero.")
        return DistributedLoad(
            self.p_i / scalar,
            self.p_j / scalar,
            self.q_i / scalar,
            self.q_j / scalar,
            self.w_i / scalar,
            self.w_j / scalar,
        )


