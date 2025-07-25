import pandas as pd
from milcapy.elements.frame import ElasticTimoshenkoBeam3D
from milcapy.geom_transf import LinearCrdTransf3D
from milcapy.core.node import Node
from milcapy.core.load import PointLoad, FullLinearDistributedLoad
from milcapy.analysis import LinearStaticAnalysis
from typing import Dict, Tuple
import numpy as np
from milcapy.types import ElementModelType
from numpy.typing import NDArray


class Model:
    """Clase que representa un modelo estructural."""

    def __init__(self) -> None:
        """Inicializa un nuevo modelo estructural."""
        self.transf: Dict[int, LinearCrdTransf3D] = {}
        self.nodes: Dict[int, Node] = {}
        self.members: Dict[int, ElasticTimoshenkoBeam3D] = {}

    def set_dofs(self, ux: bool, uy: bool, uz: bool, rx: bool, ry: bool, rz: bool):
        if self.nodes != {}:
            raise ValueError(
                "La eleccion de DOF's debe realizarse antes de crear nodos")
        Node._RESTRAINTS = (ux, uy, uz, rx, ry, rz)

    def def_geom_transf(self, tag: int, vx: float, vy: float, vz: float) -> LinearCrdTransf3D:
        """Define la transformación geométrica para el modelo."""
        if tag in self.transf:
            raise ValueError(
                f"Ya existe una transformación geométrica con el tag {tag}.")

        vec_xz = np.array([vx, vy, vz])
        transf = LinearCrdTransf3D(tag, vec_xz)
        self.transf[tag] = transf

        return transf

    def add_node(self, tag: int, x: float, y: float, z: float) -> Node:
        """Agrega un nodo al modelo."""
        if tag in self.nodes:
            raise ValueError(f"Ya existe un nodo con el tag {tag}.")
        node = Node(tag, np.array([x, y, z]))
        self.nodes[tag] = node

        return node

    def add_member(self, tag: int, node_i_tag: int, node_j_tag: int, E: float,
                   G: float, A: float, Jx: float, Iy: float, Iz: float,
                   Asy: float, Asz: float, transf_tag: int = 1,
                   ) -> ElasticTimoshenkoBeam3D:
        """ Agrega un miembro estructural al modelo."""
        if tag in self.members:
            raise ValueError(f"Ya existe un miembro con el tag {tag}.")
        if node_i_tag not in self.nodes:
            raise ValueError(f"No existe el nodo {node_i_tag}.")
        if node_j_tag not in self.nodes:
            raise ValueError(f"No existe el nodo {node_j_tag}.")

        # Obtener las coordenadas de los nodos
        node_i = self.nodes[node_i_tag]
        node_j = self.nodes[node_j_tag]

        coords_i = node_i.coords
        coords_j = node_j.coords

        # Comparar con prioridad: Z > Y > X
        if (coords_j[2], coords_j[1], coords_j[0]) < (coords_i[2], coords_i[1], coords_i[0]):
            node_i_tag, node_j_tag = node_j_tag, node_i_tag

        member = ElasticTimoshenkoBeam3D(tag,
                                         self.nodes[node_i_tag], self.nodes[node_j_tag], E, G, A,
                                         Jx, Iy, Iz, Asy, Asz, self.transf[transf_tag])
        self.members[tag] = member

        return member


    def add_restraint(self, node_tag: int, restraint: 'Restraint') -> None:
        """Agrega un restricción al nodo."""
        if node_tag not in self.nodes:
            raise ValueError(f"No existe el nodo {node_tag}.")

        self.nodes[node_tag].set_restraint(restraint)

    def add_point_load(self, node_tag: int, fx: float = 0.0, fy: float = 0.0, fz: float = 0.0,
                       mx: float = 0.0, my: float = 0.0, mz: float = 0.0, replace: bool = False) -> None:
        """Agrega una carga puntual al modelo."""
        if node_tag not in self.nodes:
            raise ValueError(f"No existe el nodo {node_tag}.")

        load = PointLoad(fx, fy, fz, mx, my, mz)

        if replace:
            self.nodes[node_tag].set_load(load)
        else:
            self.nodes[node_tag].add_load(load)

    def add_distributed_load(self, member_tag: int, direction: int, wi: float, wj: float,
                             replace: bool = False) -> None:
        """Agrega una carga distribuida al modelo."""
        if member_tag not in self.members:
            raise ValueError(f"No existe el miembro {member_tag}.")

        if direction not in [1, 2, 3]:
            raise ValueError(f"La dirección debe ser '1', '2' o '3'.")

        if direction == 1:
            p_i = wi
            p_j = wj
            q_i = 0.0
            q_j = 0.0
            w_i = 0.0
            w_j = 0.0
        elif direction == 2:
            p_i = 0.0
            p_j = 0.0
            q_i = wi
            q_j = wj
            w_i = 0.0
            w_j = 0.0
        elif direction == 3:
            p_i = 0.0
            p_j = 0.0
            q_i = 0.0
            q_j = 0.0
            w_i = wi
            w_j = wj

        load = FullLinearDistributedLoad(p_i, p_j, q_i, q_j, w_i, w_j)

        if replace:
            self.members[member_tag].set_load(load)
        else:
            self.members[member_tag].add_load(load)

    def solve(self) -> Tuple[NDArray, NDArray]:
        """Resuelve el modelo."""

        Node._ID = 1
        for member in self.members.values():
            member.set_up()

        analysis = LinearStaticAnalysis(self)
        displacements, reactions = analysis.solve()

        for node in self.nodes.values():
            node.displacements = displacements[node.dofs - 1]

        for member in self.members.values():
            member.ul = member.Tlg @ displacements[member.dofs - 1]
            if member.ql0 is None:
                member.ql0 = np.zeros(12)
            if len(member.dofs) == 12:
                member.ql = member.kl @ member.ul - member.ql0
            else:
                continue

        return displacements, reactions

    def print_node_results(self):
        data = {}
        for i in self.nodes.values():
            data[f"Node {i.tag}"] = i.displacements
        df = pd.DataFrame(
            data, index=['UX', 'UY', 'UZ', 'RX', 'RY', 'RZ']).round(6)
        print(df)

    def show_model(self):
        """Plotea el modelo."""
        import matplotlib.pyplot as plt

        self.graph_model = plt.figure(figsize=(8, 8))
        ax = self.graph_model.add_subplot(111, projection='3d')

        for node in self.nodes.values():
            ax.plot(node.coords[0], node.coords[1],
                    node.coords[2], 'o', color='k', markersize=5)
            ax.text(node.coords[0], node.coords[1],
                    node.coords[2], str(node.tag))
            if node.restraint == (True, True, True, True, True, True):
                ax.scatter(node.coords[0], node.coords[1],
                           node.coords[2], marker='s', s=50, c='#a1179d')
            else:
                continue

        for member in self.members.values():

            member.set_up()

            ax.plot([member.node_i.coords[0], member.node_j.coords[0]],
                    [member.node_i.coords[1], member.node_j.coords[1]],
                    [member.node_i.coords[2], member.node_j.coords[2]], color='b')

            if member.model_type == ElementModelType.ElasticTimoshenkoBeam3d:
                # Graficamos los ejes locales:
                p_med = (member.node_i.coords + member.node_j.coords) / 2
                # vectores unitarios de los ejes locales:
                x_axis = member.x_axis
                y_axis = member.y_axis
                z_axis = member.z_axis
                # graficamos los ejes locales:
                ax.quiver(p_med[0], p_med[1], p_med[2], x_axis[0],
                          x_axis[1], x_axis[2], color='#ff0000', length=0.5)
                ax.quiver(p_med[0], p_med[1], p_med[2], y_axis[0],
                          y_axis[1], y_axis[2], color='#40ff40', length=0.5)
                ax.quiver(p_med[0], p_med[1], p_med[2], z_axis[0],
                          z_axis[1], z_axis[2], color='#0000ff', length=0.5)

            # Graficamos los ejes globales:
            ax.quiver(0, 0, 0, 1, 0, 0, color='#ff0000', length=0.5)
            ax.quiver(0, 0, 0, 0, 1, 0, color='#40ff40', length=0.5)
            ax.quiver(0, 0, 0, 0, 0, 1, color='#0000ff', length=0.5)
        plt.tight_layout()
        plt.axis('equal')
        plt.show()

    def show_deformed(self, escale: float = 1.0, ndp: int = 100):
        """Plotea el modelo deformado."""
        import matplotlib.pyplot as plt

        self.graph_deformed = plt.figure(figsize=(8, 8))
        ax = self.graph_deformed.add_subplot(111, projection='3d')

        for member in self.members.values():

            # undeformed:
            ax.plot([member.node_i.coords[0], member.node_j.coords[0]],
                    [member.node_i.coords[1], member.node_j.coords[1]],
                    [member.node_i.coords[2], member.node_j.coords[2]], color='#742e9b', linewidth=1.5, linestyle=':')

            # deformed:
            xi = member.node_i.coords[0] + \
                member.node_i.displacements[0] * escale
            yi = member.node_i.coords[1] + \
                member.node_i.displacements[1] * escale
            zi = member.node_i.coords[2] + \
                member.node_i.displacements[2] * escale

            xj = member.node_j.coords[0] + \
                member.node_j.displacements[0] * escale
            yj = member.node_j.coords[1] + \
                member.node_j.displacements[1] * escale
            zj = member.node_j.coords[2] + \
                member.node_j.displacements[2] * escale
            ax.plot([xi, xj], [yi, yj], [zi, zj],
                    color='b', linewidth=1.5, linestyle=':')

        # deformed interpolate:
        from milcapy.elements.frame import deformed_shape
        deformed_shape(self, n=ndp, escale=escale, ax=ax)

        # restraints:
        for node in self.nodes.values():
            if node.restraint == (True, True, True, True, True, True):
                ax.scatter(node.coords[0], node.coords[1],
                           node.coords[2], marker='s', s=50, c='#a1179d')

        plt.tight_layout()
        plt.axis('equal')
        plt.show()

    def mesh(self, d_max: float = 1.0):
        """Discretiza los miembros con una longitud máxima d_max.

        - Mantiene los nodos originales y crea nodos intermedios interpolados.
        - Divide los elementos según la longitud máxima, creando nuevos elementos con las mismas propiedades.
        - Elimina el elemento original y almacena los nuevos nodos y elementos en el modelo.
        - Maneja cargas distribuidas ajustándolas proporcionalmente para cada sub-elemento.
        """

        # Funciones auxiliares para generar tags únicos
        def _get_new_node_tag():
            return max(self.nodes.keys(), default=0) + 1

        def _get_new_member_tag():
            return max(self.members.keys(), default=0) + 1

        # Procesar una copia de los miembros para evitar modificar el diccionario durante la iteración
        members_to_process = list(self.members.items())

        for member_tag, member in members_to_process:
            if member.length is None:
                member.length = np.linalg.norm(
                    member.node_j.coords - member.node_i.coords)

            n = max(1, int(np.ceil(member.length / d_max)))
            if n == 1:
                continue  # No se necesita discretizar

            # Calcular posiciones de los nodos intermedios
            node_i = member.node_i
            node_j = member.node_j
            coords_i = node_i.coords
            coords_j = node_j.coords

            # Lista de tags: [tag_i, tag_inter1, ..., tag_j]
            new_node_tags = [node_i.tag]

            # Crear nodos intermedios
            for i in range(1, n):
                # Interpolar coordenadas (fracción a lo largo del elemento)
                t = i / n
                new_coords = (1 - t) * coords_i + t * coords_j

                # Crear nuevo nodo
                new_tag = _get_new_node_tag()
                self.add_node(new_tag, *new_coords)
                new_node_tags.append(new_tag)

            new_node_tags.append(node_j.tag)

            # Crear nuevos elementos conectando los nodos secuencialmente
            for i in range(len(new_node_tags) - 1):
                new_member_tag = _get_new_member_tag()

                # Crear elemento con las mismas propiedades que el original
                new_member = self.add_member(
                    tag=new_member_tag,
                    node_i_tag=new_node_tags[i],
                    node_j_tag=new_node_tags[i + 1],
                    E=member.E,
                    G=member.G,
                    A=member.A,
                    Jx=member.Jx,
                    Iy=member.Iy,
                    Iz=member.Iz,
                    Asy=member.Asy,
                    Asz=member.Asz,
                    transf_tag=member.the_transf.tag
                )

                # Configurar el elemento (necesario para cálculos posteriores)
                new_member.set_up()

                # Copiar y ajustar cargas distribuidas
                if member.load is not None:
                    # Calcular fracción de carga para este sub-elemento
                    load = member.load
                    t_start = i / n
                    t_end = (i + 1) / n

                    # Ajustar cargas lineales (interpolación)
                    new_p_i = load.p_i + (load.p_j - load.p_i) * t_start
                    new_p_j = load.p_i + (load.p_j - load.p_i) * t_end
                    new_q_i = load.q_i + (load.q_j - load.q_i) * t_start
                    new_q_j = load.q_i + (load.q_j - load.q_i) * t_end
                    new_w_i = load.w_i + (load.w_j - load.w_i) * t_start
                    new_w_j = load.w_i + (load.w_j - load.w_i) * t_end

                    new_load = FullLinearDistributedLoad(
                        new_p_i, new_p_j,
                        new_q_i, new_q_j,
                        new_w_i, new_w_j
                    )
                    new_member.add_load(new_load)

            # Eliminar el elemento original
            del self.members[member_tag]

