import numpy as np



nodes = {
    1: np.array([0, 0, 0]),
    2: np.array([0, 0, 4]),
    3: np.array([4, 0, 4]),
    4: np.array([4, 0, 0]),
    5: np.array([0, 4, 0]),
    6: np.array([0, 4, 4]),
    7: np.array([4, 4, 4]),
    8: np.array([4, 4, 0]),
    9: np.array([0, 8, 0]),
    10: np.array([0, 8, 4]),
    11: np.array([4, 8, 4]),
    12: np.array([4, 8, 0]),
}

members = {
    1: (1, 2, 2.1e6, 0.5*0.3),
    2: (2, 3, 2.1e6, 0.5*0.3),
    3: (3, 4, 2.1e6, 0.5*0.3),
    4: (4, 1, 2.1e6, 0.5*0.3),
    5: (5, 6, 2.1e6, 0.5*0.3),
    6: (6, 7, 2.1e6, 0.5*0.3),
    7: (7, 8, 2.1e6, 0.5*0.3),
    8: (8, 5, 2.1e6, 0.5*0.3),
    9: (9, 10, 2.1e6, 0.5*0.3),
    10: (10, 11, 2.1e6, 0.5*0.3),
    11: (11, 12, 2.1e6, 0.5*0.3),
    12: (12, 9, 2.1e6, 0.5*0.3),
    13: (1, 5, 2.1e6, 0.5*0.3),
    14: (2, 6, 2.1e6, 0.5*0.3),
    15: (3, 7, 2.1e6, 0.5*0.3),
    16: (4, 8, 2.1e6, 0.5*0.3),
    17: (5, 9, 2.1e6, 0.5*0.3),
    18: (6, 10, 2.1e6, 0.5*0.3),
    19: (7, 11, 2.1e6, 0.5*0.3),
    20: (8, 12, 2.1e6, 0.5*0.3),
    21: (1, 6, 2.1e6, 0.5*0.3),
    22: (4, 7, 2.1e6, 0.5*0.3),
    23: (6, 9, 2.1e6, 0.5*0.3),
    24: (7, 12, 2.1e6, 0.5*0.3),
    25: (3, 6, 2.1e6, 0.5*0.3),
    26: (6, 11, 2.1e6, 0.5*0.3),
    27: (4, 5, 2.1e6, 0.5*0.3),
    28: (5, 12, 2.1e6, 0.5*0.3),
}

loads = {#id: fx, fy, fz
    10: (0, 0, -10),
    11: (0, 0, -10),
}

# =====================================================
def ke(E, A, L):
    return E*A/L*np.array([[1, -1], [-1, 1]])

def transform_matrix(coords_i, coords_j):
    L = np.linalg.norm(coords_j - coords_i)
    cx = (coords_j[0] - coords_i[0]) / L
    cy = (coords_j[1] - coords_i[1]) / L
    cz = (coords_j[2] - coords_i[2]) / L
    T = np.array([
        [cx, cy, cz, 0, 0, 0],
        [0, 0, 0, cx, cy, cz],
    ])
    return T

def global_stiffness_matrix(ke, transform_matrix):
    return transform_matrix.T @ ke @ transform_matrix

def global_load_vector(ql, transform_matrix):
    return transform_matrix.T @ ql

# ASSEMBLY
ndf = 3
K = np.zeros((ndf*len(nodes), ndf*len(nodes)))
for i, (node_i, node_j, E, A) in members.items():
    L = np.linalg.norm(nodes[node_j] - nodes[node_i])
    k_e = ke(E, A, L)
    T = transform_matrix(nodes[node_i], nodes[node_j])
    K_e = global_stiffness_matrix(k_e, T)
    dof_i = np.array([(node_i - 1) * ndf, (node_i - 1) * ndf + 1, (node_i - 1) * ndf + 2])
    dof_j = np.array([(node_j - 1) * ndf, (node_j - 1) * ndf + 1, (node_j - 1) * ndf + 2])
    for i in range(ndf):
        for j in range(ndf):
            K[dof_i[i], dof_j[j]] += K_e[i, j]

print(K)

# FORCES
F = np.zeros((ndf*len(nodes), 1))
for i, (fx, fy, fz) in loads.items():
    # completar el vector de fuerzas
    F[i*ndf] = fx
    F[i*ndf+1] = fy
    F[i*ndf+2] = fz

print(F)

# RESTRAINTS
soportes = {
    1: (True, True, True),
    2: (True, True, True),
    3: (True, True, True),
    4: (True, True, True),
}

# Aplicar condiciones de restringido a la matriz K, F
restraints = np.zeros(len(nodes) * ndf, dtype=bool)

for node in nodes.keys():
    restraints[(node - 1) * ndf] = soportes.get(node, (False, False, False))[0]
    restraints[(node - 1) * ndf + 1] = soportes.get(node, (False, False, False))[1]
    restraints[(node - 1) * ndf + 2] = soportes.get(node, (False, False, False))[2]

free_dofs = np.where(~restraints)[0]
restrained_dofs = np.where(restraints)[0]


# SOLVE
K_d = K[np.ix_(free_dofs, free_dofs)]
F_d = F[free_dofs]

u_d = np.linalg.solve(K_d, F_d)
u = np.zeros(len(nodes) * ndf)
u[free_dofs] = u_d

print(u)

# REACTION FORCES




