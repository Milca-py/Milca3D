from __future__ import annotations
from typing import TYPE_CHECKING, Tuple
import numpy as np

if TYPE_CHECKING:
    from milcapy.core.model import Model
    from numpy.typing import NDArray

class LinearStaticAnalysis:
    """
    Solucionador basado en el Método de Rigidez Directa.
    Resuelve sistemas K·u = F usando solución directa con numpy.
    """

    def __init__(self, model: Model) -> None:
        self.model = model
        self.nn = len(self.model.nodes)  # Número de nodos

    def _dof_map(self) -> Tuple[NDArray[np.int32], NDArray[np.int32]]:
        """Mapea grados de libertad libres y restringidos."""
        restraints = np.zeros(self.nn * 6, dtype=np.bool_)

        for node in self.model.nodes.values():
            restraints[node.dofs - 1] = node.restraint
            # restraints[node.dofs[0] - 1] = node.restraint[0]
            # restraints[node.dofs[1] - 1] = node.restraint[1]
            # restraints[node.dofs[2] - 1] = node.restraint[2]
            # restraints[node.dofs[3] - 1] = node.restraint[3]
            # restraints[node.dofs[4] - 1] = node.restraint[4]
            # restraints[node.dofs[5] - 1] = node.restraint[5]

        free_dofs = np.where(~restraints)[0].astype(np.int32)
        restrained_dofs = np.where(restraints)[0].astype(np.int32)

        return free_dofs, restrained_dofs

    def assemble_global_load_vector(self) -> NDArray[float64]:
        """Calcula el vector de carga global del sistema.

        Returns:
            NDArray: Vector de carga global.
        """
        F = np.zeros(self.nn * 6, dtype=np.float64)

        # Asignar fuerzas nodales almacenadas en los nodos
        for node in self.model.nodes.values():
            F[node.dofs - 1] += node.get_load_vector()

        # Agregar el vector de fuerzas globales almacenadas en los miembros
        for member in self.model.members.values():
            f = member.get_global_load_vector()
            if f is not None:
                F[member.dofs - 1] += f

        return F

    def assemble_global_stiffness_matrix(self) -> NDArray[np.float64]:
        """Ensamblaje de la matriz global de rigidez."""
        K = np.zeros((self.nn * 6, self.nn * 6), dtype=np.float64)

        for member in self.model.members.values():
            k = member.get_global_stiffness_matrix()
            dofs = member.dofs
            rows = np.repeat(dofs - 1, 12)
            cols = np.tile(dofs - 1, 12)
            K[rows, cols] += k.flatten()

        return K

    def apply_boundary_conditions(
        self,
        K: NDArray[np.float64],
        F: NDArray[np.float64],
        free_dofs: NDArray[np.int32],
        restrained_dofs: NDArray[np.int32]
    ) -> Tuple[
        NDArray[np.float64], NDArray[np.float64], NDArray[np.float64],
        NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
    ]:
        """Aplica condiciones de frontera y retorna particiones reducidas."""
        # Reducir la matriz de rigidez (d: desconocidos, c: conocidos)
        #     | Kd   Kdc |
        # K = |          |
        #     | Kcd   Kc |
        K_d = K[np.ix_(free_dofs, free_dofs)]
        K_dc = K[np.ix_(free_dofs, restrained_dofs)]
        K_cd = K[np.ix_(restrained_dofs, free_dofs)]
        K_c = K[np.ix_(restrained_dofs, restrained_dofs)]

        # Reducir el vector de fuerzas
        #     | Fd |
        # F = |    |
        #     | Fc |
        F_d = F[free_dofs]
        F_c = F[restrained_dofs]

        return K_d, K_dc, K_cd, K_c, F_d, F_c

    def solve(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Resuelve el sistema estructural K·u = F."""
        K = self.assemble_global_stiffness_matrix()
        F = self.assemble_global_load_vector()
        free_dofs, restrained_dofs = self._dof_map()
        K_d, K_dc, K_cd, K_c, F_d, F_c = self.apply_boundary_conditions(K, F, free_dofs, restrained_dofs)

        # Resolver el sistema de ecuaciones
        #     |Ud|
        # U = |  | = displacements
        #     |Uc|

        # Estabilidad estructural
        if np.linalg.det(K_d) == 0:
            print("Advertencia: estructura posiblemente inestable.")
            K_d += np.eye(K_d.shape[0]) * 1e-12

        # Resolver desplazamientos
        U_d = np.linalg.solve(K_d, F_d)

        # Vector completo de desplazamientos
        displacements = np.zeros(self.nn * 6, dtype=np.float64)
        displacements[free_dofs] = U_d

        # Reacciones en apoyos
        reactions = np.zeros(self.nn * 6, dtype=np.float64)
        reactions[restrained_dofs] = K_cd @ U_d - F_c

        # forma compacta de calcular las reacciones
        # R  =  K_global * U_global - F_global
        # reactions = K @ displacements - F

        return displacements, reactions