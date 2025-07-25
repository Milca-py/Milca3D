

class PortalParametric:
    def __init__(self, m: int, n: int, p: int, l: float, b: float, h: float):
        self.m = m
        self.n = n
        self.p = p
        self.l = l
        self.b = b
        self.h = h

    def nodes_from_story(self, story: int) -> list[int]:
        # devuele los nodos de una piso
        nodes = [x for x in range(
            story * (self.n + 1) * (self.m + 1) + 1, (story + 1) * (self.n + 1) * (self.m + 1) + 1)]
        return nodes

    def nodes(self) -> dict[int, tuple[float, float, float]]:
        # Nodos
        nodes = {}  # {id, (x, y, z)}
        for i in range(self.p + 1):
            for j in range(self.n + 1):
                for k in range(self.m + 1):
                    nodes[i * (self.n + 1) * (self.m + 1) + j * (self.m + 1) +
                          k + 1] = (k * self.l, j * self.b, i * self.h)

        return nodes

    def members(self) -> dict[str, dict[int, tuple[int, int]]]:
        # Miembros
        members = {"columna": {}, "viga": {}}  # {id, (node_id_i, node_id_j)}
        id = 1
        for i in range(self.p):
            # columnas
            for j in range(self.n + 1):
                for k in range(self.m + 1):
                    # columnas en el plano XY en todos los pisos
                    members["columna"][id] = (i * (self.n + 1) * (self.m + 1) + j * (
                        self.m + 1) + k + 1, (i + 1) * (self.n + 1) * (self.m + 1) + j * (self.m + 1) + k + 1)
                    id += 1

        for i in range(self.p + 1):
            # vigas
            if i == 0:
                pass
            else:
                for j in range(self.n + 1):
                    for k in range(self.m + 1):
                        # vigas en el plano XY en todos los pisos
                        if k == self.m:
                            pass
                        else:
                            members["viga"][id] = (i * (self.n + 1) * (self.m + 1) + j * (
                                self.m + 1) + k + 1, (i * (self.n + 1) * (self.m + 1) + j * (self.m + 1) + k + 2))
                            id += 1

                for j in range(self.m + 1):
                    for k in range(self.n + 1):
                        # vigas en el plano XY en todos los pisos
                        if k == self.n:
                            pass
                        else:
                            members["viga"][id] = (i * (self.n + 1) * (self.m + 1) + k * (
                                self.m + 1) + j + 1, (i * (self.n + 1) * (self.m + 1) + (k + 1) * (self.m + 1) + j + 1))
                            id += 1

        return members
