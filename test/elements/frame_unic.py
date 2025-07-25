from milcapy.types import ElementModelType
from typing import TYPE_CHECKING
from milcapy.geom_transf import LinearCrdTransf3D
from milcapy.types import ElementModelType, DirectionType
from milcapy.core.node import Node

from milcapy.core.load import DistributedLoad

from numpy.typing import NDArray
from abc import ABC
import numpy as np

# from __future__ import annotations
if TYPE_CHECKING:
    from milcapy.core.model import Model


class Element(ABC):
    """Clase base abstracta para elementos. Define la interfaz para elementos estructurales."""

    def __init__(self, tag: int, model_type: 'ElementModelType'):
        """Inicializa el elemento con su etiqueta y tipo de modelo usado.

        Args:
            tag (int): Identificador único del elemento.
            model_type (ElementModelType): Tipo de modelo asociado al elemento.
        """
        self.tag: int = tag
        self.model_type: ElementModelType = model_type


class ElasticTimoshenkoBeam3D(Element):
    """
    Clase que representa un elemento de viga Timoshenko elástica en 3D.
    Modela deformaciones por cortante y efectos de inercia rotacional.
    """

    def __init__(
        self,
        tag: int,
        node_i: Node,
        node_j: Node,
        E: float,
        G: float,
        A: float,
        Jx: float,
        Iy: float,
        Iz: float,
        Asy: float,
        Asz: float,
        the_transf: LinearCrdTransf3D,
    ):
        super().__init__(tag, ElementModelType.ElasticTimoshenkoBeam3d)

        self.node_i: Node = node_i
        self.node_j: Node = node_j

        # Propiedades del material y sección
        self.E: float = E               # Módulo de elasticidad (Young)
        self.G: float = G               # Módulo de cortante (Shear)
        self.A: float = A               # Área de la sección transversal
        self.Jx: float = Jx             # Constante torsional de la sección
        self.Iy: float = Iy             # Momento de inercia respecto al eje y
        self.Iz: float = Iz             # Momento de inercia respecto al eje z
        self.Asy: float = Asy           # Área efectiva para cortante en dirección y
        self.Asz: float = Asz           # Área efectiva para cortante en dirección z

        # Propiedades geometricas
        self.the_transf: LinearCrdTransf3D = the_transf
        self.phiY: float | None = None   # Relación de rigidez flexión-cortante en y
        self.phiZ: float | None = None   # Relación de rigidez flexión-cortante en z
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

        # Cargas =============================
        self.load: DistributedLoad | None = None   # Cargas aplicadas

        # Ensamblaje =========================
        self.dofs: NDArray | None = None   # [ux, uy, uz, rx, ry, rz] * 2

    def set_up(self):
        """Configura las matrices y parámetros iniciales del elemento."""
        # self.the_transf.connect(self.node_i, self.node_j)

        self.v_xz = self.the_transf.v_xz

        self.dofs = np.concatenate((self.node_i.dofs, self.node_j.dofs))

        self.length = np.linalg.norm(self.node_j.coords - self.node_i.coords)

        self.Tlg = self.transform_matrix()

        self.phiY = 12.0 * self.E * self.Iz / \
            (self.length**2 * self.G * self.Asy)
        self.phiZ = 12.0 * self.E * self.Iy / \
            (self.length**2 * self.G * self.Asz)

        self.kl = self.local_stiffnes_matrix()

        if self.load is not None:
            self.ql0 = self.q_phi()

    def get_global_stiffness_matrix(self) -> NDArray:
        return self.Tlg.T @ self.kl @ self.Tlg

    def get_global_load_vector(self) -> NDArray:
        if self.ql0 is None:
            return np.zeros(12)
        # return self.the_transf.get_global_load_vector(self.ql0, self.Tlg)
        return self.Tlg.T @ self.ql0

    def set_load(self, load: DistributedLoad):
        """Establece las cargas aplicadas al elemento."""
        self.load = load

    def add_load(self, load: DistributedLoad):
        """Agrega cargas aplicadas al elemento."""
        if self.load is None:
            self.load = load
        else:
            self.load += load

    def __del__(self):
        """Elimina el elemento."""
        self.node_i = None
        self.node_j = None
        self.the_transf = None
        self.kl = None
        self.Tlg = None
        self.Ki = None
        self.load = None
        self.dofs = None

    def deformed_shape(self, n: int = 100, escale: float = 1) -> None:
        """
        Plotea la forma deformada del miembro.
        """
        deformed_shape(self, n, escale)

    def local_stiffnes_matrix(self) -> NDArray:
        """Calcula la matriz de rigidez local del elemento siguiendo la teoría de Timoshenko."""
        # Constantes
        E = self.E
        G = self.G
        A = self.A
        L = self.length
        Jx = self.Jx
        Iy = self.Iy
        Iz = self.Iz
        phiY = self.phiY
        phiZ = self.phiZ

        k11 = E * A / L
        k22 = 12 * E * Iz / (L**3 * (1 + phiY))
        k33 = 12 * E * Iy / (L**3 * (1 + phiZ))
        k44 = G * Jx / L
        k55 = (4 + phiZ) * E * Iy / (L * (1 + phiZ))
        k66 = (4 + phiY) * E * Iz / (L * (1 + phiY))
        k62 = 6 * E * Iz / (L**2 * (1 + phiY))
        k53 = 6 * E * Iy / (L**2 * (1 + phiZ))
        k115 = (2 - phiZ) * E * Iy / (L * (1 + phiZ))
        k116 = (2 - phiY) * E * Iz / (L * (1 + phiY))

        kii = np.array([
            [k11,   0,      0,      0,      0,      0],
            [0,     k22,    0,      0,      0,      k62],
            [0,     0,      k33,    0,     -k53,    0],
            [0,     0,      0,      k44,    0,      0],
            [0,     0,     -k53,    0,      k55,    0],
            [0,     k62,    0,      0,      0,      k66],
        ])

        kij = np.array([
            [-k11,  0,      0,      0,      0,      0],
            [0,    -k22,    0,      0,      0,     -k62],
            [0,     0,     -k33,    0,      k53,    0],
            [0,     0,      0,     -k44,    0,      0],
            [0,     0,     -k53,    0,      k115,   0],
            [0,     k62,    0,      0,      0,      k116],
        ])

        kjj = np.array([
            [k11,   0,      0,      0,      0,      0],
            [0,     k22,    0,      0,      0,     -k62],
            [0,     0,      k33,    0,      k53,    0],
            [0,     0,      0,      k44,    0,      0],
            [0,     0,      k53,    0,      k55,    0],
            [0,    -k62,    0,      0,      0,      k66],
        ])

        kl = np.zeros((12, 12))
        kl[0:6, 0:6] = kii
        kl[6:12, 0:6] = kij
        kl[0:6, 6:12] = kij.T
        kl[6:12, 6:12] = kjj

        return kl

    def transform_matrix(self) -> NDArray:
        """Transforma la matriz de rigidez local a global."""
        node_i = self.node_i
        node_j = self.node_j

        # Vector de dirección del elemento (x')
        vec_elem = (node_j.coords - node_i.coords)
        x_axis = vec_elem / np.linalg.norm(vec_elem)

        # z' = x' × v_aux
        z_axis = np.cross(x_axis, self.v_xz)
        z_axis /= np.linalg.norm(z_axis)

        # y' = z' × x'
        y_axis = np.cross(z_axis, x_axis)

        R = np.vstack((x_axis, y_axis, z_axis))

        T = np.zeros((12, 12))
        for i in range(4):
            T[3 * i:3 * i + 3, 3 * i:3 * i + 3] = R

        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis

        return T

    def q_phi(self) -> NDArray:

        # constantes:
        L = self.length
        phiY = self.phiY
        phiZ = self.phiZ
        pi = self.load.p_i
        pj = self.load.p_j
        qi = self.load.q_i
        qj = self.load.q_j
        wi = self.load.w_i
        wj = self.load.w_j

        # Acciones axiales
        FXi = (2 * pi + pj) * L / 6
        FXj = (pi + 2 * pj) * L / 6

        # Acciones torsionales
        MXi = 0  # Torsión
        MXj = 0  # Torsión

        # Acciones en dirección LOCAL Y
        ay = (qj - qi) / L
        by = qi

        FYi = by * L / 2 + ay * L**2 / 60 * ((10 * phiZ + 9) / (phiZ + 1))
        FYj = by * L / 2 + ay * L**2 / 60 * ((20 * phiZ + 21) / (phiZ + 1))

        MZi = by * L**2 / 12 + ay * L**3 / 120 * ((5 * phiZ + 4) / (phiZ + 1))
        MZj = -(by * L**2 / 12 + ay * L**3 / 120 *
                ((5 * phiZ + 6) / (phiZ + 1)))

        # Acciones en dirección LOCAL Z
        az = (wj - wi) / L
        bz = wi

        FZi = bz * L / 2 + az * L**2 / 60 * ((10 * phiY + 9) / (phiY + 1))
        FZj = bz * L / 2 + az * L**2 / 60 * ((20 * phiY + 21) / (phiY + 1))

        MYi = -bz * L**2 / 12 + az * L**3 / 120 * ((5 * phiY + 4) / (phiY + 1))
        MYj = bz * L**2 / 12 + az * L**3 / 120 * ((5 * phiY + 6) / (phiY + 1))

        return np.array([FXi, FYi, FZi, MXi, MYi, MZi, FXj, FYj, FZj, MXj, MYj, MZj])


    def coefficients(self):
        L = self.Length()
        q1 = self.q1
        q2 = self.q2
        E = self.E
        I = self.I
        phi = self.phi
        v1 = self.v1
        v2 = self.v2
        theta1 = self.theta1
        theta2 = self.theta2

        A = -(q2 - q1) / L
        B = -q1

        M = np.array([
            [0,                     0,          0, 1],
            [0,                     0,          1, 0],
            [L**3 * (2 - phi) / 12, L**2 / 2,   L, 1],
            [L**2 / 2,              L,          1, 0],
        ])

        N = np.array([
            E * I * v1,
            E * I * theta1,
            E * I * v2 + A * L**5 * (0.6 - phi) / 72 +
            B * L**4 * (1 - phi) / 24,
            E * I * theta2 + A * L**4 / 24 + B * L**3 / 6,
        ])

        C = np.linalg.solve(M, N)

        self.C1 = C[0]
        self.C2 = C[1]
        self.C3 = C[2]
        self.C4 = C[3]

        return tuple(map(float, C))

    def Length(self):
        return self.x2 - self.x1

    def shear(self, x):
        q1 = self.q1
        q2 = self.q2
        L = self.Length()
        A = (q2 - q1)/L
        B = q1

        V = self.C1 - B*x - (A/2)*x**2

        return V

    def moment(self, x):
        q1 = self.q1
        q2 = self.q2
        L = self.Length()
        A = (q2 - q1)/L
        B = q1

        M = self.C2 + self.C1*x - (B/2)*x**2 - (A/6)*x**3

        return M

    def slope(self, x):
        q1 = self.q1
        q2 = self.q2
        L = self.Length()
        EI = self.E * self.I
        As = self.As
        G = self.G
        A = (q2 - q1)/L
        B = q1

        C1 = self.C1
        C2 = self.C2
        C3 = self.C3

        theta_x = 1/EI * (C3 + C2*x + (C1/2)*x**2 - (B/6)*x**3 -
                          (A/24)*x**4 + ((EI/G*As)*((A/2)*x**2 + B*x)))

        return theta_x

    def deflection(self, x):

        q1 = self.q1
        q2 = self.q2
        L = self.Length()
        EI = self.E * self.I
        A = (q2 - q1)/L
        B = q1

        C1 = self.C1
        C2 = self.C2
        C3 = self.C3
        C4 = self.C4

        phi = self.phi

        term1 = A * L**2 * x**3 * (0.6 * (x / L)**2 - phi) / 72
        term2 = B * L**2 * x**2 * ((x / L)**2 - phi) / 24
        term3 = C1 * x * L**2 * (2 * (x / L)**2 - phi) / 12
        term4 = C2 * x**2 / 2
        term5 = C3 * x
        term6 = C4

        u = (term1 + term2 + term3 + term4 + term5 + term6) / (EI)

        return u

    def process_builder(self, direction: DirectionType) -> None:
        if self.load is None:
            self.load = DistributedLoad()
        self.x1 = 0
        self.x2 = self.length
        self.pi = self.load.p_i
        self.pj =   self.load.p_j
        self.P1 = self.ql[0]
        self.P2 = self.ql[6]
        self.T1 = self.ql[3]
        self.T2 = self.ql[9]
        self.u1 = self.ul[0]
        self.u2 = self.ul[6]
        self.E = self.E
        self.A = self.A

        if direction == DirectionType.Y:
            self.x1 = 0
            self.x2 = self.length
            self.q1 = self.load.q_i
            self.q2 = self.load.q_j
            self.V1 = self.ql[1]
            self.V2 = self.ql[7]
            self.M1 = self.ql[5]
            self.M2 = self.ql[11]
            self.v1 = self.ul[1]
            self.v2 = self.ul[7]
            self.theta1 = self.ul[5]
            self.theta2 = self.ul[11]
            self.E = self.E
            self.G = self.G
            self.I = self.Iz
            self.A = self.A
            self.As = self.Asy
            self.phi = self.phiY
        elif direction == DirectionType.Z:
            self.x1 = 0
            self.x2 = self.length
            self.q1 = self.load.w_i
            self.q2 = self.load.w_j
            self.V1 = self.ql[2]
            self.V2 = self.ql[8]
            self.M1 = -self.ql[4]
            self.M2 = -self.ql[10]
            self.v1 = self.ul[2]
            self.v2 = self.ul[8]
            self.theta1 = -self.ul[4]
            self.theta2 = -self.ul[10]
            self.E = self.E
            self.G = self.G
            self.I = self.Iy
            self.A = self.A
            self.As = self.Asz
            self.phi = self.phiZ
        self.coefficients()

    def axial(self, x):
        a = (self.pj - self.pi)/self.length
        b = self.pi

        EA = self.A * self.E
        L = self.Length()

        C1 = EA/L * (self.u2 - self.u1) + (a/6)*L**2 + (b/2)*L

        P = - (a/2)*x**2 - b*x + C1

        return P

    def torsion(self, x):

        t1 = 0
        t2 = 0
        L = self.Length()
        A = (t2 - t1)/L
        B = t1
        C1 = self.T1

        T = - (A/2)*x**2 - B*x + C1

        return T

    def displacement(self, x):

        a = (self.pj - self.pi)/self.Length()
        b = self.pi

        EA = self.A * self.E
        L = self.Length()

        C1 = EA/L * (self.u2 - self.u1) + (a/6)*L**2 + (b/2)*L
        C2 = EA * self.u1

        u = 1/EA * (-(a/6)*x**3 - (b/2)*x**2 + C1*x + C2)

        return u




def deformed_shape(model: 'Model', n: int, escale: float, ax: 'Axes3D') -> None:
    """
    Plotea la forma deformada del miembro.
    """
    # plotear:
    for member in model.members.values():
        if member.model_type == ElementModelType.ElasticTimoshenkoBeam3d:
            x = np.linspace(0, member.length, n)
            # Obtener las coordenadas en sistema local
            ux = np.array([member.displacement(xi)
                          for xi in x]) * escale
            uy = np.array([member.deflection(xi)
                           for xi in x]) * escale
            uz = np.array([member.deflection(xi)
                           for xi in x]) * escale

            # PLOTEAR EN COORDENADAS GLOBALES:
            Tlg = member.Tlg[:3, :3]    # matriz de rotacion

            nd1 = member.node_i.coords
            disp_shap = np.array([x + ux, uy, uz])
            disp_shap = np.dot(Tlg.T, disp_shap)
            disp_shap = disp_shap + nd1[:, np.newaxis]

            ax.plot(disp_shap[0], disp_shap[1], disp_shap[2], color="#0000ff")

        elif member.model_type == ElementModelType.ElasticTruss3d:

            nd1 = member.node_i.coords
            nd2 = member.node_j.coords
            # forma rigida de la deformada con escala
            coord_rig_i = nd1 + [member.node_i.displacements[0]*escale,
                                 member.node_i.displacements[1]*escale, member.node_i.displacements[2]*escale]
            coord_rig_j = nd2 + [member.node_j.displacements[0]*escale,
                                 member.node_j.displacements[1]*escale, member.node_j.displacements[2]*escale]

            ax.plot([coord_rig_i[0], coord_rig_j[0]], [coord_rig_i[1], coord_rig_j[1]], [
                    coord_rig_i[2], coord_rig_j[2]], color="#0000ff")

    ax.axis('equal')
