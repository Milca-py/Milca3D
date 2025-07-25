import numpy as np
import math
from scipy.linalg import solve, lu_factor, lu_solve

# Reutilizamos la clase Vec3 de coordtrans.py
class Vec3:
    """Estructura para representar coordenadas 3D."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

# Importamos funciones de coordtrans.py (asumimos que están en el mismo módulo)
from coordtrans import coord_trans, atma

def assemble_K(K, DoF, nE, xyz, r, L, Le, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, p, shear, geom, Q, debug):
    """Ensambla la matriz de rigidez global a partir de los elementos individuales.
    
    Args:
        K: Matriz de rigidez global (DoF x DoF, modificada en lugar).
        DoF: Número de grados de libertad.
        nE: Número de elementos de viga.
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        r: Radios rígidos de los nodos.
        L, Le: Longitudes de los elementos (total y efectiva).
        N1, N2: Conectividad de nodos (índices).
        Ax, Asy, Asz: Áreas de sección de los elementos.
        Jx, Iy, Iz: Inercias de sección de los elementos.
        E, G: Módulos elástico y de corte.
        p: Ángulos de rotación (radianes).
        shear: 1 para incluir deformación por corte, 0 para no incluir.
        geom: 1 para incluir rigidez geométrica, 0 para no incluir.
        Q: Fuerzas en los extremos de los elementos.
        debug: 1 para guardar matrices de rigidez de elementos, 0 para no guardar.
    """
    k = np.zeros((12, 12))  # Matriz de rigidez del elemento
    ind = np.zeros((12, nE), dtype=int)  # Tabla de índices de grados de libertad

    # Inicializar matriz de rigidez global
    K.fill(0.0)

    # Calcular índices de grados de libertad para cada elemento
    for i in range(nE):
        ind[0, i] = 6 * N1[i] - 6
        ind[6, i] = 6 * N2[i] - 6
        for j in range(1, 6):
            ind[j, i] = ind[0, i] + j
            ind[j + 6, i] = ind[6, i] + j

    # Ensamblar matriz de rigidez global
    for i in range(nE):
        elastic_K(k, xyz, r, L[i], Le[i], N1[i], N2[i], Ax[i], Asy[i], Asz[i], Jx[i], Iy[i], Iz[i], E[i], G[i], p[i], shear)
        if geom:
            geometric_K(k, xyz, r, L[i], Le[i], N1[i], N2[i], Ax[i], Asy[i], Asz[i], Jx[i], Iy[i], Iz[i], E[i], G[i], p[i], -Q[i, 0], shear)
        
        # Ensamblar en matriz global
        for l in range(12):
            ii = ind[l, i]
            for ll in range(12):
                jj = ind[ll, i]
                K[ii, jj] += k[l, ll]

def elastic_K(k, xyz, r, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, shear):
    """Calcula la matriz de rigidez elástica del elemento en coordenadas globales.
    
    Args:
        k: Matriz de rigidez del elemento (12x12, modificada en lugar).
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        r: Radios rígidos de los nodos.
        L, Le: Longitudes del elemento (total y efectiva).
        n1, n2: Índices de nodos.
        Ax, Asy, Asz: Áreas de sección.
        J, Iy, Iz: Inercias de sección.
        E, G: Módulos elástico y de corte.
        p: Ángulo de rotación (radianes).
        shear: 1 para incluir deformación por corte, 0 para no incluir.
    """
    t1 = [0.0]; t2 = [0.0]; t3 = [0.0]; t4 = [0.0]; t5 = [0.0]
    t6 = [0.0]; t7 = [0.0]; t8 = [0.0]; t9 = [0.0]
    
    coord_trans(xyz, L, n1, n2, t1, t2, t3, t4, t5, t6, t7, t8, t9, p)
    
    k.fill(0.0)
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0
    
    k[0, 0] = k[6, 6] = E * Ax / Le
    k[1, 1] = k[7, 7] = 12.0 * E * Iz / (Le**3 * (1.0 + Ksy))
    k[2, 2] = k[8, 8] = 12.0 * E * Iy / (Le**3 * (1.0 + Ksz))
    k[3, 3] = k[9, 9] = G * J / Le
    k[4, 4] = k[10, 10] = (4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz))
    k[5, 5] = k[11, 11] = (4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy))
    
    k[4, 2] = k[2, 4] = -6.0 * E * Iy / (Le**2 * (1.0 + Ksz))
    k[5, 1] = k[1, 5] = 6.0 * E * Iz / (Le**2 * (1.0 + Ksy))
    k[6, 0] = k[0, 6] = -k[0, 0]
    
    k[11, 7] = k[7, 11] = k[7, 5] = k[5, 7] = -k[5, 1]
    k[10, 8] = k[8, 10] = k[8, 4] = k[4, 8] = -k[4, 2]
    k[9, 3] = k[3, 9] = -k[3, 3]
    k[10, 2] = k[2, 10] = k[4, 2]
    k[11, 1] = k[1, 11] = k[5, 1]
    
    k[7, 1] = k[1, 7] = -k[1, 1]
    k[8, 2] = k[2, 8] = -k[2, 2]
    k[10, 4] = k[4, 10] = (2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz))
    k[11, 5] = k[5, 11] = (2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy))
    
    atma(t1[0], t2[0], t3[0], t4[0], t5[0], t6[0], t7[0], t8[0], t9[0], k, r[n1], r[n2])
    
    # Verificar y forzar simetría
    for i in range(12):
        for j in range(i + 1, 12):
            if abs(k[i, j] - k[j, i]) > 1e-6 and (abs(k[i, j] / k[i, i]) > 1e-6 or abs(k[j, i] / k[i, i]) > 1e-6):
                print(f"elastic_K: Matriz de rigidez no simétrica en k[{i}][{j}] = {k[i,j]:.6e}, k[{j}][{i}] = {k[j,i]:.6e}")
                k[i, j] = k[j, i] = 0.5 * (k[i, j] + k[j, i])

def geometric_K(k, xyz, r, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, T, shear):
    """Calcula la matriz de rigidez geométrica del elemento en coordenadas globales.
    
    Args:
        k: Matriz de rigidez del elemento (12x12, modificada en lugar).
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        r: Radios rígidos de los nodos.
        L, Le: Longitudes del elemento (total y efectiva).
        n1, n2: Índices de nodos.
        Ax, Asy, Asz: Áreas de sección.
        J, Iy, Iz: Inercias de sección.
        E, G: Módulos elástico y de corte.
        p: Ángulo de rotación (radianes).
        T: Fuerza axial para rigidez geométrica.
        shear: 1 para incluir deformación por corte, 0 para no incluir.
    """
    t1 = [0.0]; t2 = [0.0]; t3 = [0.0]; t4 = [0.0]; t5 = [0.0]
    t6 = [0.0]; t7 = [0.0]; t8 = [0.0]; t9 = [0.0]
    
    coord_trans(xyz, L, n1, n2, t1, t2, t3, t4, t5, t6, t7, t8, t9, p)
    
    kg = np.zeros((12, 12))
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0
    Dsy = (1 + Ksy)**2 if shear else 1.0
    Dsz = (1 + Ksz)**2 if shear else 1.0
    
    kg[0, 0] = kg[6, 6] = 0.0  # T/L
    kg[1, 1] = kg[7, 7] = T / L * (1.2 + 2.0 * Ksy + Ksy**2) / Dsy
    kg[2, 2] = kg[8, 8] = T / L * (1.2 + 2.0 * Ksz + Ksz**2) / Dsz
    kg[3, 3] = kg[9, 9] = T / L * J / Ax
    kg[4, 4] = kg[10, 10] = T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz
    kg[5, 5] = kg[11, 11] = T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy
    
    kg[0, 6] = kg[6, 0] = 0.0  # -T/L
    kg[4, 2] = kg[2, 4] = kg[10, 2] = kg[2, 10] = -T / 10.0 / Dsz
    kg[8, 4] = kg[4, 8] = kg[10, 8] = kg[8, 10] = T / 10.0 / Dsz
    kg[5, 1] = kg[1, 5] = kg[11, 1] = kg[1, 11] = T / 10.0 / Dsy
    kg[7, 5] = kg[5, 7] = kg[11, 7] = kg[7, 11] = -T / 10.0 / Dsy
    kg[9, 3] = kg[3, 9] = -kg[3, 3]
    
    kg[7, 1] = kg[1, 7] = -kg[1, 1]
    kg[8, 2] = kg[2, 8] = -kg[2, 2]
    kg[10, 4] = kg[4, 10] = -T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz
    kg[11, 5] = kg[5, 11] = -T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy
    
    atma(t1[0], t2[0], t3[0], t4[0], t5[0], t6[0], t7[0], t8[0], t9[0], kg, r[n1], r[n2])
    
    # Verificar y forzar simetría
    for i in range(12):
        for j in range(i + 1, 12):
            if abs(kg[i, j] - kg[j, i]) > 1e-6 and (abs(kg[i, j] / kg[i, i]) > 1e-6 or abs(kg[j, i] / kg[i, i]) > 1e-6):
                print(f"geometric_K: Matriz de rigidez no simétrica en kg[{i}][{j}] = {kg[i,j]:.6e}, kg[{j}][{i}] = {kg[j,i]:.6e}")
                kg[i, j] = kg[j, i] = 0.5 * (kg[i, j] + kg[j, i])
    
    # Sumar matriz geométrica a la matriz elástica
    k += kg

def frame_element_force(s, xyz, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, f_t, f_m, D, shear, geom, axial_strain):
    """Calcula las fuerzas en los extremos del elemento en coordenadas locales.
    
    Args:
        s: Vector de fuerzas en los extremos (12 elementos, modificado en lugar).
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        L, Le: Longitudes del elemento (total y efectiva).
        n1, n2: Índices de nodos.
        Ax, Asy, Asz: Áreas de sección.
        J, Iy, Iz: Inercias de sección.
        E, G: Módulos elástico y de corte.
        p: Ángulo de rotación (radianes).
        f_t, f_m: Fuerzas equivalentes por temperatura y cargas mecánicas.
        D: Vector de desplazamientos globales.
        shear: 1 para incluir deformación por corte, 0 para no incluir.
        geom: 1 para incluir rigidez geométrica, 0 para no incluir.
        axial_strain: Deformación axial del elemento (modificada en lugar).
    """
    t1 = [0.0]; t2 = [0.0]; t3 = [0.0]; t4 = [0.0]; t5 = [0.0]
    t6 = [0.0]; t7 = [0.0]; t8 = [0.0]; t9 = [0.0]
    
    coord_trans(xyz, L, n1, n2, t1, t2, t3, t4, t5, t6, t7, t8, t9, p)
    
    n1 = 6 * (n1 - 1)
    n2 = 6 * (n2 - 1)
    
    d = D[n1:n1+6].tolist() + D[n2:n2+6].tolist()
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0
    Dsy = (1 + Ksy)**2 if shear else 1.0
    Dsz = (1 + Ksz)**2 if shear else 1.0
    
    delta = (d[6] - d[0]) * t1[0] + (d[7] - d[1]) * t2[0] + (d[8] - d[2]) * t3[0]
    axial_strain[0] = delta / Le
    
    s[0] = -(Ax * E / Le) * delta
    T = -s[0] if geom else 0.0
    
    s[1] = -(12.0 * E * Iz / (Le**3 * (1.0 + Ksy)) + T / L * (1.2 + 2.0 * Ksy + Ksy**2) / Dsy) * \
           ((d[6] - d[0]) * t4[0] + (d[7] - d[1]) * t5[0] + (d[8] - d[2]) * t6[0]) + \
           (6.0 * E * Iz / (Le**2 * (1.0 + Ksy)) + T / 10.0 / Dsy) * \
           ((d[3] + d[9]) * t7[0] + (d[4] + d[10]) * t8[0] + (d[5] + d[11]) * t9[0])
    
    s[2] = -(12.0 * E * Iy / (Le**3 * (1.0 + Ksz)) + T / L * (1.2 + 2.0 * Ksz + Ksz**2) / Dsz) * \
           ((d[6] - d[0]) * t7[0] + (d[7] - d[1]) * t8[0] + (d[8] - d[2]) * t9[0]) - \
           (6.0 * E * Iy / (Le**2 * (1.0 + Ksz)) + T / 10.0 / Dsz) * \
           ((d[3] + d[9]) * t4[0] + (d[4] + d[10]) * t5[0] + (d[5] + d[11]) * t6[0])
    
    s[3] = -(G * J / Le) * ((d[9] - d[3]) * t1[0] + (d[10] - d[4]) * t2[0] + (d[11] - d[5]) * t3[0])
    
    s[4] = (6.0 * E * Iy / (Le**2 * (1.0 + Ksz)) + T / 10.0 / Dsz) * \
           ((d[6] - d[0]) * t7[0] + (d[7] - d[1]) * t8[0] + (d[8] - d[2]) * t9[0]) + \
           ((4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz)) + T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz) * \
           (d[3] * t4[0] + d[4] * t5[0] + d[5] * t6[0]) + \
           ((2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz)) - T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz) * \
           (d[9] * t4[0] + d[10] * t5[0] + d[11] * t6[0])
    
    s[5] = -(6.0 * E * Iz / (Le**2 * (1.0 + Ksy)) + T / 10.0 / Dsy) * \
           ((d[6] - d[0]) * t4[0] + (d[7] - d[1]) * t5[0] + (d[8] - d[2]) * t6[0]) + \
           ((4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy)) + T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy) * \
           (d[3] * t7[0] + d[4] * t8[0] + d[5] * t9[0]) + \
           ((2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy)) - T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy) * \
           (d[9] * t7[0] + d[10] * t8[0] + d[11] * t9[0])
    
    s[6] = -s[0]
    s[7] = -s[1]
    s[8] = -s[2]
    s[9] = -s[3]
    
    s[10] = (6.0 * E * Iy / (Le**2 * (1.0 + Ksz)) + T / 10.0 / Dsz) * \
            ((d[6] - d[0]) * t7[0] + (d[7] - d[1]) * t8[0] + (d[8] - d[2]) * t9[0]) + \
            ((4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz)) + T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz) * \
            (d[9] * t4[0] + d[10] * t5[0] + d[11] * t6[0]) + \
            ((2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz)) - T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz**2 / 12.0) / Dsz) * \
            (d[3] * t4[0] + d[4] * t5[0] + d[5] * t6[0])
    
    s[11] = -(6.0 * E * Iz / (Le**2 * (1.0 + Ksy)) + T / 10.0 / Dsy) * \
            ((d[6] - d[0]) * t4[0] + (d[7] - d[1]) * t5[0] + (d[8] - d[2]) * t6[0]) + \
            ((4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy)) + T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy) * \
            (d[9] * t7[0] + d[10] * t8[0] + d[11] * t9[0]) + \
            ((2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy)) - T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy**2 / 12.0) / Dsy) * \
            (d[3] * t7[0] + d[4] * t8[0] + d[5] * t9[0])
    
    # Agregar fuerzas fijas en los extremos
    f = f_t + f_m
    s[0] -= (f[0] * t1[0] + f[1] * t2[0] + f[2] * t3[0])
    s[1] -= (f[0] * t4[0] + f[1] * t5[0] + f[2] * t6[0])
    s[2] -= (f[0] * t7[0] + f[1] * t8[0] + f[2] * t9[0])
    s[3] -= (f[3] * t1[0] + f[4] * t2[0] + f[5] * t3[0])
    s[4] -= (f[3] * t4[0] + f[4] * t5[0] + f[5] * t6[0])
    s[5] -= (f[3] * t7[0] + f[4] * t8[0] + f[5] * t9[0])
    s[6] -= (f[6] * t1[0] + f[7] * t2[0] + f[8] * t3[0])
    s[7] -= (f[6] * t4[0] + f[7] * t5[0] + f[8] * t6[0])
    s[8] -= (f[6] * t7[0] + f[7] * t8[0] + f[8] * t9[0])
    s[9] -= (f[9] * t1[0] + f[10] * t2[0] + f[11] * t3[0])
    s[10] -= (f[9] * t4[0] + f[10] * t5[0] + f[11] * t6[0])
    s[11] -= (f[9] * t7[0] + f[10] * t8[0] + f[11] * t9[0])

def element_end_forces(Q, nE, xyz, L, Le, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, p, eqF_temp, eqF_mech, D, shear, geom, axial_strain_warning):
    """Evalúa las fuerzas en los extremos de todos los elementos.
    
    Args:
        Q: Fuerzas en los extremos de los elementos (nE x 12, modificada en lugar).
        nE: Número de elementos.
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        L, Le: Longitudes de los elementos (total y efectiva).
        N1, N2: Conectividad de nodos (índices).
        Ax, Asy, Asz: Áreas de sección.
        Jx, Iy, Iz: Inercias de sección.
        E, G: Módulos elástico y de corte.
        p: Ángulos de rotación (radianes).
        eqF_temp, eqF_mech: Fuerzas equivalentes por temperatura y cargas mecánicas.
        D: Vector de desplazamientos globales.
        shear: 1 para incluir deformación por corte, 0 para no incluir.
        geom: 1 para incluir rigidez geométrica, 0 para no incluir.
        axial_strain_warning: Contador de advertencias de deformación axial (modificado en lugar).
    """
    s = np.zeros(12)
    axial_strain = [0.0]
    axial_strain_warning[0] = 0
    
    for m in range(nE):
        frame_element_force(s, xyz, L[m], Le[m], N1[m], N2[m], Ax[m], Asy[m], Asz[m], Jx[m], Iy[m], Iz[m], E[m], G[m], p[m], eqF_temp[m], eqF_mech[m], D, shear, geom, axial_strain)
        Q[m, :] = s
        if abs(axial_strain[0]) > 0.001:
            axial_strain_warning[0] += 1

def solve_system(K, D, F, R, DoF, q, r, ok, verbose, rms_resid):
    """Resuelve el sistema {F} = [K]{D} mediante descomposición LDL.
    
    Args:
        K: Matriz de rigidez (DoF x DoF).
        D: Vector de desplazamientos (modificado en lugar).
        F: Vector de cargas externas.
        R: Vector de reacciones (modificado en lugar).
        DoF: Número de grados de libertad.
        q: 1 para coordenadas libres, 0 para reacciones.
        r: 0 para coordenadas libres, 1 para reacciones.
        ok: Indicador de matriz definida positiva (modificado en lugar).
        verbose: 1 para salida detallada, 0 para ninguna.
        rms_resid: Error RMS del residuo (modificado en lugar).
    """
    # Nota: ldl_dcmp_pm y ldl_mprove_pm no están definidas, usamos scipy.linalg.solve
    try:
        free_dofs = [i for i in range(DoF) if q[i]]
        K_qq = K[np.ix_(free_dofs, free_dofs)]
        F_q = F[free_dofs]
        D_q = solve(K_qq, F_q, assume_a='pos')
        
        for i, idx in enumerate(free_dofs):
            D[idx] = D_q[i]
        
        # Calcular reacciones
        for i in range(DoF):
            if r[i]:
                R[i] = -F[i]
                for j in range(DoF):
                    R[i] += K[i, j] * D[j]
        
        ok[0] = 1
        rms_resid[0] = equilibrium_error(np.zeros(DoF), F, K, D, DoF, q, r)
    except np.linalg.LinAlgError:
        ok[0] = -1
        print("Error: Matriz de rigidez no definida positiva. Asegúrese de que todas las traslaciones rígidas estén restringidas.")

def equilibrium_error(dF, F, K, D, DoF, q, r):
    """Calcula el error de equilibrio {dF} = {F} - [K_qq]{D_q} - [K_qr]{D_r}.
    
    Args:
        dF: Vector de error de equilibrio (modificado en lugar).
        F: Vector de cargas externas.
        K: Matriz de rigidez.
        D: Vector de desplazamientos.
        DoF: Número de grados de libertad.
        q: 1 para coordenadas libres, 0 para reacciones.
        r: 0 para coordenadas libres, 1 para reacciones.
    
    Returns:
        float: Norma del error dividido por la norma de F.
    """
    ss_dF = 0.0
    ss_F = 0.0
    
    for i in range(DoF):
        if q[i]:
            errF = F[i]
            for j in range(DoF):
                if q[j]:
                    errF -= K[i, j] * D[j] if i <= j else K[j, i] * D[j]
                if r[j]:
                    errF -= K[i, j] * D[j]
            dF[i] = errF
            ss_dF += dF[i] ** 2
            ss_F += F[i] ** 2
    
    return math.sqrt(ss_dF) / math.sqrt(ss_F) if ss_F > 0 else 0.0

def compute_reaction_forces(R, F, K, D, DoF, r):
    """Calcula las fuerzas de reacción: R(r) = [K(r,q)]*{D(q)} + [K(r,r)]*{D(r)} - F(r).
    
    Args:
        R: Vector de reacciones (modificado en lugar).
        F: Vector de cargas externas.
        K: Matriz de rigidez.
        D: Vector de desplazamientos.
        DoF: Número de grados de libertad.
        r: 0 para coordenadas libres, 1 para reacciones.
    """
    for i in range(DoF):
        if r[i]:
            R[i] = -F[i]
            for j in range(DoF):
                R[i] += K[i, j] * D[j]

def assemble_M(M, DoF, nN, nE, xyz, r, L, N1, N2, Ax, Jx, Iy, Iz, p, d, EMs, NMs, NMx, NMy, NMz, lump, debug):
    """Ensambla la matriz de masa global a partir de masas e inercias de elementos.
    
    Args:
        M: Matriz de masa global (DoF x DoF, modificada en lugar).
        DoF: Número de grados de libertad.
        nN, nE: Número de nodos y elementos.
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        r: Radios rígidos de los nodos.
        L: Longitudes de los elementos.
        N1, N2: Conectividad de nodos (índices).
        Ax, Jx, Iy, Iz: Áreas e inercias de sección.
        p: Ángulos de rotación (radianes).
        d: Densidad de los elementos.
        EMs: Masa adicional de los elementos.
        NMs, NMx, NMy, NMz: Masas e inercias nodales.
        lump: 1 para matriz de masa concentrada, 0 para consistente.
        debug: 1 para guardar matrices de masa de elementos, 0 para no guardar.
    """
    m = np.zeros((12, 12))
    ind = np.zeros((12, nE), dtype=int)
    
    M.fill(0.0)
    
    for i in range(nE):
        ind[0, i] = 6 * N1[i] - 6
        ind[6, i] = 6 * N2[i] - 6
        for j in range(1, 6):
            ind[j, i] = ind[0, i] + j
            ind[j + 6, i] = ind[6, i] + j
    
    for i in range(nE):
        if lump:
            lumped_M(m, xyz, L[i], N1[i], N2[i], Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i])
        else:
            consistent_M(m, xyz, r, L[i], N1[i], N2[i], Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i])
        
        for l in range(12):
            ii = ind[l, i]
            for ll in range(12):
                jj = ind[ll, i]
                M[ii, jj] += m[l, ll]
    
    for j in range(nN):
        i = 6 * j
        M[i, i] += NMs[j]
        M[i + 1, i + 1] += NMs[j]
        M[i + 2, i + 2] += NMs[j]
        M[i + 3, i + 3] += NMx[j]
        M[i + 4, i + 4] += NMy[j]
        M[i + 5, i + 5] += NMz[j]
    
    for i in range(DoF):
        if M[i, i] <= 0.0:
            print(f"Error: Matriz de masa no definida positiva en M[{i}][{i}] = {M[i,i]}")

def lumped_M(m, xyz, L, n1, n2, Ax, J, Iy, Iz, p, d, EMs):
    """Calcula la matriz de masa concentrada del elemento en coordenadas globales.
    
    Args:
        m: Matriz de masa del elemento (12x12, modificada en lugar).
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        L: Longitud del elemento.
        n1, n2: Índices de nodos.
        Ax, J, Iy, Iz: Áreas e inercias de sección.
        p: Ángulo de rotación (radianes).
        d: Densidad del elemento.
        EMs: Masa adicional del elemento.
    """
    t1 = [0.0]; t2 = [0.0]; t3 = [0.0]; t4 = [0.0]; t5 = [0.0]
    t6 = [0.0]; t7 = [0.0]; t8 = [0.0]; t9 = [0.0]
    
    coord_trans(xyz, L, n1, n2, t1, t2, t3, t4, t5, t6, t7, t8, t9, p)
    
    t = (d * Ax * L + EMs) / 2.0
    ry = d * Iy * L / 2.0
    rz = d * Iz * L / 2.0
    po = d * L * J / 2.0
    
    m.fill(0.0)
    
    m[0, 0] = m[1, 1] = m[2, 2] = m[6, 6] = m[7, 7] = m[8, 8] = t
    m[3, 3] = m[9, 9] = po * t1[0]**2 + ry * t4[0]**2 + rz * t7[0]**2
    m[4, 4] = m[10, 10] = po * t2[0]**2 + ry * t5[0]**2 + rz * t8[0]**2
    m[5, 5] = m[11, 11] = po * t3[0]**2 + ry * t6[0]**2 + rz * t9[0]**2
    
    m[3, 4] = m[4, 3] = m[9, 10] = m[10, 9] = po * t1[0] * t2[0] + ry * t4[0] * t5[0] + rz * t7[0] * t8[0]
    m[3, 5] = m[5, 3] = m[9, 11] = m[11, 9] = po * t1[0] * t3[0] + ry * t4[0] * t6[0] + rz * t7[0] * t9[0]
    m[4, 5] = m[5, 4] = m[10, 11] = m[11, 10] = po * t2[0] * t3[0] + ry * t5[0] * t6[0] + rz * t8[0] * t9[0]

def consistent_M(m, xyz, r, L, n1, n2, Ax, J, Iy, Iz, p, d, EMs):
    """Calcula la matriz de masa consistente del elemento en coordenadas globales.
    
    Args:
        m: Matriz de masa del elemento (12x12, modificada en lugar).
        xyz: Lista de objetos Vec3 con coordenadas de nodos.
        r: Radios rígidos de los nodos.
        L: Longitud del elemento.
        n1, n2: Índices de nodos.
        Ax, J, Iy, Iz: Áreas e inercias de sección.
        p: Ángulo de rotación (radianes).
        d: Densidad del elemento.
        EMs: Masa adicional del elemento.
    """
    t1 = [0.0]; t2 = [0.0]; t3 = [0.0]; t4 = [0.0]; t5 = [0.0]
    t6 = [0.0]; t7 = [0.0]; t8 = [0.0]; t9 = [0.0]
    
    coord_trans(xyz, L, n1, n2, t1, t2, t3, t4, t5, t6, t7, t8, t9, p)
    
    t = d * Ax * L
    ry = d * Iy
    rz = d * Iz
    po = d * J * L
    
    m.fill(0.0)
    
    m[0, 0] = m[6, 6] = t / 3.0
    m[1, 1] = m[7, 7] = 13.0 * t / 35.0 + 6.0 * rz / (5.0 * L)
    m[2, 2] = m[8, 8] = 13.0 * t / 35.0 + 6.0 * ry / (5.0 * L)
    m[3, 3] = m[9, 9] = po / 3.0
    m[4, 4] = m[10, 10] = t * L**2 / 105.0 + 2.0 * L * ry / 15.0
    m[5, 5] = m[11, 11] = t * L**2 / 105.0 + 2.0 * L * rz / 15.0
    
    m[4, 2] = m[2, 4] = -11.0 * t * L / 210.0 - ry / 10.0
    m[5, 1] = m[1, 5] = 11.0 * t * L / 210.0 + rz / 10.0
    m[6, 0] = m[0, 6] = t / 6.0
    
    m[7, 5] = m[5, 7] = 13.0 * t * L / 420.0 - rz / 10.0
    m[8, 4] = m[4, 8] = -13.0 * t * L / 420.0 + ry / 10.0
    m[9, 3] = m[3, 9] = po / 6.0
    m[10, 2] = m[2, 10] = 13.0 * t * L / 420.0 - ry / 10.0
    m[11, 1] = m[1, 11] = -13.0 * t * L / 420.0 + rz / 10.0
    
    m[10, 8] = m[8, 10] = 11.0 * t * L / 210.0 + ry / 10.0
    m[11, 7] = m[7, 11] = -11.0 * t * L / 210.0 - rz / 10.0
    
    m[7, 1] = m[1, 7] = 9.0 * t / 70.0 - 6.0 * rz / (5.0 * L)
    m[8, 2] = m[2, 8] = 9.0 * t / 70.0 - 6.0 * ry / (5.0 * L)
    m[10, 4] = m[4, 10] = -L**2 * t / 140.0 - ry * L / 30.0
    m[11, 5] = m[5, 11] = -L**2 * t / 140.0 - rz * L / 30.0
    
    for i in range(3):
        m[i, i] += 0.5 * EMs
    for i in range(6, 9):
        m[i, i] += 0.5 * EMs
    
    atma(t1[0], t2[0], t3[0], t4[0], t5[0], t6[0], t7[0], t8[0], t9[0], m, r[n1], r[n2])
    
    # Verificar y forzar simetría
    for i in range(12):
        for j in range(i + 1, 12):
            if abs(m[i, j] - m[j, i]) > 1e-6 and (abs(m[i, j] / m[i, i]) > 1e-6 or abs(m[j, i] / m[i, i]) > 1e-6):
                print(f"consistent_M: Matriz de masa no simétrica en m[{i}][{j}] = {m[i,j]:.6e}, m[{j}][{i}] = {m[j,i]:.6e}")
                m[i, j] = m[j, i] = 0.5 * (m[i, j] + m[j, i])

def static_condensation(A, N, c, n, Ac, verbose):
    """Condensación estática de la matriz de rigidez de NxN a nxn.
    
    Args:
        A: Matriz cuadrada (NxN).
        N: Dimensión de la matriz original.
        c: Lista de índices a retener.
        n: Dimensión de la matriz condensada.
        Ac: Matriz condensada (nxn, modificada en lugar).
        verbose: 1 para salida detallada, 0 para ninguna.
    """
    r = [i for i in range(1, N + 1) if i not in c[:n]]
    N_r = N - n
    
    Arr = np.zeros((N_r, N_r))
    Arc = np.zeros((N_r, n))
    
    for i in range(N_r):
        for j in range(i, N_r):
            ri, rj = r[i] - 1, r[j] - 1
            Arr[j, i] = Arr[i, j] = A[ri, rj] if ri <= rj else A[rj, ri]
    
    for i in range(N_r):
        for j in range(n):
            ri, cj = r[i] - 1, c[j] - 1
            Arc[i, j] = A[ri, cj] if ri < cj else A[cj, ri]
    
    # Calcular Ac = Arc^T * inv(Arr) * Arc
    try:
        inv_Arr_Arc = solve(Arr, Arc, assume_a='sym')
        Ac_temp = np.dot(Arc.T, inv_Arr_Arc)
        
        for i in range(n):
            for j in range(i, n):
                ci, cj = c[i] - 1, c[j] - 1
                Ac[j, i] = Ac[i, j] = A[ci, cj] - Ac_temp[i, j] if ci <= cj else A[cj, ci] - Ac_temp[i, j]
    except np.linalg.LinAlgError:
        if verbose:
            print("Error en static_condensation: Matriz Arr no invertible.")

def paz_condensation(M, K, N, c, n, Mc, Kc, w2, verbose):
    """Condensación dinámica de Paz para matrices de masa y rigidez.
    
    Args:
        M, K: Matrices de masa y rigidez (NxN).
        N: Dimensión de las matrices.
        c: Lista de grados de libertad a retener.
        n: Dimensión de las matrices condensadas.
        Mc, Kc: Matrices condensadas (nxn, modificadas en lugar).
        w2: Cuadrado de la frecuencia objetivo.
        verbose: 1 para salida detallada, 0 para ninguna.
    """
    w2 = 4.0 * math.pi * math.pi * w2 * w2
    r = [i for i in range(1, N + 1) if i not in c[:n]]
    N_r = N - n
    
    Drr = np.zeros((N_r, N_r))
    Drc = np.zeros((N_r, n))
    T = np.zeros((N, n))
    
    for i in range(N_r):
        for j in range(N_r):
            ri, rj = r[i] - 1, r[j] - 1
            Drr[j, i] = Drr[i, j] = (K[ri, rj] - w2 * M[ri, rj]) if ri <= rj else (K[rj, ri] - w2 * M[rj, ri])
    
    for i in range(N_r):
        for j in range(n):
            ri, cj = r[i] - 1, c[j] - 1
            Drc[i, j] = (K[ri, cj] - w2 * M[ri, cj]) if ri < cj else (K[cj, ri] - w2 * M[cj, ri])
    
    try:
        invDrrDrc = solve(Drr, Drc, assume_a='sym')
        for i in range(n):
            T[c[i] - 1, i] = 1.0
        for i in range(N_r):
            for j in range(n):
                T[r[i] - 1, j] = -invDrrDrc[i, j]
        
        Kc[:] = np.dot(np.dot(T.T, K), T)
        Mc[:] = np.dot(np.dot(T.T, M), T)
    except np.linalg.LinAlgError:
        if verbose:
            print("Error en paz_condensation: Matriz Drr no invertible.")

def modal_condensation(M, K, N, R, p, n, Mc, Kc, V, f, m, verbose):
    """Condensación dinámica para igualar respuesta en frecuencias y modos específicos.
    
    Args:
        M, K: Matrices de masa y rigidez (NxN).
        N: Dimensión de las matrices.
        R: 1 para grados de libertad fijos, 0 para libres.
        p: Lista de grados de libertad primarios.
        n: Dimensión de las matrices condensadas.
        Mc, Kc: Matrices condensadas (nxn, modificadas en lugar).
        V: Formas modales.
        f: Frecuencias naturales.
        m: Lista de modos a igualar.
        verbose: 1 para salida detallada, 0 para ninguna.
    """
    P = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            P[i, j] = V[p[i] - 1, m[j] - 1]
    
    try:
        invP = np.linalg.pinv(P, rcond=1e-9)
        traceM = sum(M[i, i] for i in range(N) if not R[i])
        
        Mc[:] = np.dot(np.dot(invP.T, np.eye(n)), invP)
        traceMc = np.trace(Mc)
        Mc *= traceM / traceMc if traceMc != 0 else 1.0
        
        Kc[:] = np.dot(np.dot(invP.T, np.diag([4.0 * math.pi**2 * f[m[k] - 1]**2 for k in range(n)])), invP)
        Kc *= traceM / traceMc if traceMc != 0 else 1.0
    except np.linalg.LinAlgError:
        if verbose:
            print("Error en modal_condensation: Matriz P no invertible.")

def deallocate(nN, nE, nL, nF, nU, nW, nP, nT, DoF, modes, xyz, rj, L, Le, N1, N2, q, r, Ax, Asy, Asz, J, Iy, Iz, E, G, p, U, W, P, T, Dp, F_mech, F_temp, eqF_mech, eqF_temp, F, dF, K, Q, D, dD, R, dR, d, EMs, NMs, NMx, NMy, NMz, M, f, V, c, m, pkNx, pkVy, pkVz, pkTx, pkMy, pkMz, pkDx, pkDy, pkDz, pkRx, pkSy, pkSz):
    """Libera la memoria asignada para las estructuras de datos.
    
    Args:
        nN, nE, nL: Número de nodos, elementos y casos de carga.
        nF, nU, nW, nP, nT: Contadores de varias estructuras.
        DoF: Número de grados de libertad.
        modes: Número de modos.
        xyz, rj, L, Le, N1, N2, q, r, Ax, Asy, Asz, J, Iy, Iz, E, G, p, U, W, P, T, Dp, F_mech, F_temp, eqF_mech, eqF_temp, F, dF, K, Q, D, dD, R, dR, d, EMs, NMs, NMx, NMy, NMz, M, f, V, c, m, pkNx, pkVy, pkVz, pkTx, pkMy, pkMz, pkDx, pkDy, pkDz, pkRx, pkSy, pkSz: Estructuras de datos a liberar.
    """
    # En Python, la recolección de basura maneja la memoria automáticamente
    pass