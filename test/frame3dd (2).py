```python
import numpy as np
from math import cos, sin, sqrt, fabs

# Constante para el número máximo de casos de carga
_NL_ = 32

class Vec3:
    """Estructura para vectores cartesianos en 3D."""
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

def coord_trans(xyz, L, n1, n2, p):
    """
    Calcula los 9 elementos de la matriz de transformación de coordenadas 3x3.
    
    Args:
        xyz: Lista de objetos Vec3 con las coordenadas de los nodos.
        L: Longitud del elemento.
        n1, n2: Índices de los nodos inicial y final.
        p: Ángulo de rotación (roll) en radianes.
    
    Returns:
        tuple: Coeficientes t1, t2, t3, t4, t5, t6, t7, t8, t9 de la matriz de transformación.
    """
    # Calcular cosenos directores
    Cx = (xyz[n2].x - xyz[n1].x) / L
    Cy = (xyz[n2].y - xyz[n1].y) / L
    Cz = (xyz[n2].z - xyz[n1].z) / L

    # Inicializar coeficientes
    t1 = t2 = t3 = t4 = t5 = t6 = t7 = t8 = t9 = 0.0

    Cp = cos(p)
    Sp = sin(p)

    # Configuración para eje Z vertical (Zvert == 1)
    Zvert = True  # Asumimos Z vertical; ajustar según configuración

    if Zvert:
        if fabs(Cz) == 1.0:
            t3 = Cz
            t4 = -Cz * Sp
            t5 = Cp
            t7 = -Cz * Cp
            t8 = -Sp
        else:
            den = sqrt(1.0 - Cz * Cz)
            t1 = Cx
            t2 = Cy
            t3 = Cz
            t4 = (-Cx * Cz * Sp - Cy * Cp) / den
            t5 = (-Cy * Cz * Sp + Cx * Cp) / den
            t6 = Sp * den
            t7 = (-Cx * Cz * Cp + Cy * Sp) / den
            t8 = (-Cy * Cz * Cp - Cx * Sp) / den
            t9 = Cp * den
    else:
        # Configuración para eje Y vertical
        if fabs(Cy) == 1.0:
            t2 = Cy
            t4 = -Cy * Cp
            t6 = Sp
            t7 = Cy * Sp
            t9 = Cp
        else:
            den = sqrt(1.0 - Cy * Cy)
            t1 = Cx
            t2 = Cy
            t3 = Cz
            t4 = (-Cx * Cy * Cp - Cz * Sp) / den
            t5 = den * Cp
            t6 = (-Cy * Cz * Cp + Cx * Sp) / den
            t7 = (Cx * Cy * Sp - Cz * Cp) / den
            t8 = -den * Sp
            t9 = (Cy * Cz * Sp + Cx * Cp) / den

    return t1, t2, t3, t4, t5, t6, t7, t8, t9

def atma(t1, t2, t3, t4, t5, t6, t7, t8, t9, m, r1, r2):
    """
    Realiza la transformación de coordenadas de una matriz de locales a globales.
    
    Args:
        t1, t2, t3, t4, t5, t6, t7, t8, t9: Coeficientes de la matriz de transformación.
        m: Matriz 12x12 a transformar (rigidez o masa).
        r1, r2: Radios rígidos de los nodos inicial y final.
    """
    # Construir la matriz de transformación 12x12
    a = np.zeros((13, 13))  # Índices basados en 1
    for i in range(4):
        a[3*i+1][3*i+1] = t1
        a[3*i+1][3*i+2] = t2
        a[3*i+1][3*i+3] = t3
        a[3*i+2][3*i+1] = t4
        a[3*i+2][3*i+2] = t5
        a[3*i+2][3*i+3] = t6
        a[3*i+3][3*i+1] = t7
        a[3*i+3][3*i+2] = t8
        a[3*i+3][3*i+3] = t9

    # Efecto de los radios rígidos (desactivado en el código original)
    # Descomentar si se desea incluir
    """
    a[5][1] = r1 * t7
    a[5][2] = r1 * t8
    a[5][3] = r1 * t9
    a[6][1] = -r1 * t4
    a[6][2] = -r1 * t5
    a[6][3] = -r1 * t6
    a[11][7] = -r2 * t7
    a[11][8] = -r2 * t8
    a[11][9] = -r2 * t9
    a[12][7] = r2 * t4
    a[12][8] = r2 * t5
    a[12][9] = r2 * t6
    """

    # Calcular ma = m * a
    ma = np.zeros((13, 13))
    for j in range(1, 13):
        for i in range(1, 13):
            for k in range(1, 13):
                ma[i][j] += m[i][k] * a[k][j]

    # Calcular m = a^T * ma
    for i in range(1, 13):
        for j in range(i, 13):
            m[i][j] = m[j][i] = 0.0
            for k in range(1, 13):
                m[i][j] += a[k][i] * ma[k][j]

# Resto del código de frame3dd.py (incluye las funciones previas)

def assemble_K(K, DoF, nE, xyz, r, L, Le, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, p, shear, geom, Q, debug):
    """Forma la matriz de rigidez global a partir de elementos individuales."""
    k = np.zeros((13, 13))  # Matriz de rigidez del elemento
    ind = np.zeros((13, nE + 1), dtype=int)  # Tabla de índices de DoF

    K[1:DoF+1, 1:DoF+1] = 0.0

    for i in range(1, nE + 1):
        ind[1][i] = 6 * N1[i] - 5
        ind[7][i] = 6 * N2[i] - 5
        ind[2][i] = ind[1][i] + 1
        ind[8][i] = ind[7][i] + 1
        ind[3][i] = ind[1][i] + 2
        ind[9][i] = ind[7][i] + 2
        ind[4][i] = ind[1][i] + 3
        ind[10][i] = ind[7][i] + 3
        ind[5][i] = ind[1][i] + 4
        ind[11][i] = ind[7][i] + 4
        ind[6][i] = ind[1][i] + 5
        ind[12][i] = ind[7][i] + 5

    for i in range(1, nE + 1):
        elastic_K(k, xyz, r, L[i], Le[i], N1[i], N2[i], Ax[i], Asy[i], Asz[i], Jx[i], Iy[i], Iz[i], E[i], G[i], p[i], shear)
        
        if geom:
            geometric_K(k, xyz, r, L[i], Le[i], N1[i], N2[i], Ax[i], Asy[i], Asz[i], Jx[i], Iy[i], Iz[i], E[i], G[i], p[i], -Q[i][1], shear)
        
        if debug:
            stiffness_fn = f"k_{i:03d}"
            save_dmatrix(stiffness_fn, k, 1, 12, 1, 12, 0, "w")
        
        for l in range(1, 13):
            ii = ind[l][i]
            for ll in range(1, 13):
                jj = ind[ll][i]
                K[ii][jj] += k[l][ll]

def elastic_K(k, xyz, r, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, shear):
    """Calcula la matriz de rigidez elástica en coordenadas globales."""
    t1, t2, t3, t4, t5, t6, t7, t8, t9 = coord_trans(xyz, L, n1, n2, p)
    
    k[1:13, 1:13] = 0.0
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0

    k[1][1] = k[7][7] = E * Ax / Le
    k[2][2] = k[8][8] = 12.0 * E * Iz / (Le * Le * Le * (1.0 + Ksy))
    k[3][3] = k[9][9] = 12.0 * E * Iy / (Le * Le * Le * (1.0 + Ksz))
    k[4][4] = k[10][10] = G * J / Le
    k[5][5] = k[11][11] = (4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz))
    k[6][6] = k[12][12] = (4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy))

    k[5][3] = k[3][5] = -6.0 * E * Iy / (Le * Le * (1.0 + Ksz))
    k[6][2] = k[2][6] = 6.0 * E * Iz / (Le * Le * (1.0 + Ksy))
    k[7][1] = k[1][7] = -k[1][1]

    k[12][8] = k[8][12] = k[8][6] = k[6][8] = -k[6][2]
    k[11][9] = k[9][11] = k[9][5] = k[5][9] = -k[5][3]
    k[10][4] = k[4][10] = -k[4][4]
    k[11][3] = k[3][11] = k[5][3]
    k[12][2] = k[2][12] = k[6][2]

    k[8][2] = k[2][8] = -k[2][2]
    k[9][3] = k[3][9] = -k[3][3]
    k[11][5] = k[5][11] = (2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz))
    k[12][6] = k[6][12] = (2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy))

    atma(t1, t2, t3, t4, t5, t6, t7, t8, t9, k, r[n1], r[n2])

    for i in range(1, 13):
        for j in range(i + 1, 13):
            if abs(k[i][j] - k[j][i]) > 1e-6 and (abs(k[i][j] / k[i][i]) > 1e-6 or abs(k[j][i] / k[i][i]) > 1e-6):
                print(f"elastic_K: matriz de rigidez no simétrica ... k[{i}][{j}] = {k[i][j]:.6e}, k[{j}][{i}] = {k[j][i]:.6e}")
                k[i][j] = k[j][i] = 0.5 * (k[i][j] + k[j][i])

def geometric_K(k, xyz, r, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, T, shear):
    """Calcula la matriz de rigidez geométrica en coordenadas globales."""
    t1, t2, t3, t4, t5, t6, t7, t8, t9 = coord_trans(xyz, L, n1, n2, p)
    
    kg = np.zeros((13, 13))
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0
    Dsy = (1 + Ksy) ** 2 if shear else 1.0
    Dsz = (1 + Ksz) ** 2 if shear else 1.0

    kg[1][1] = kg[7][7] = 0.0
    kg[2][2] = kg[8][8] = T / L * (1.2 + 2.0 * Ksy + Ksy * Ksy) / Dsy
    kg[3][3] = kg[9][9] = T / L * (1.2 + 2.0 * Ksz + Ksz * Ksz) / Dsz
    kg[4][4] = kg[10][10] = T / L * J / Ax
    kg[5][5] = kg[11][11] = T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz
    kg[6][6] = kg[12][12] = T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy

    kg[1][7] = kg[7][1] = 0.0
    kg[5][3] = kg[3][5] = kg[11][3] = kg[3][11] = -T / 10.0 / Dsz
    kg[9][5] = kg[5][9] = kg[11][9] = kg[9][11] = T / 10.0 / Dsz
    kg[6][2] = kg[2][6] = kg[12][2] = kg[2][12] = T / 10.0 / Dsy
    kg[8][6] = kg[6][8] = kg[12][8] = kg[8][12] = -T / 10.0 / Dsy
    kg[4][10] = kg[10][4] = -kg[4][4]
    kg[8][2] = kg[2][8] = -T / L * (1.2 + 2.0 * Ksy + Ksy * Ksy) / Dsy
    kg[9][3] = kg[3][9] = -T / L * (1.2 + 2.0 * Ksz + Ksz * Ksz) / Dsz
    kg[11][5] = kg[5][11] = -T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz
    kg[12][6] = kg[6][12] = -T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy

    atma(t1, t2, t3, t4, t5, t6, t7, t8, t9, kg, r[n1], r[n2])

    for i in range(1, 13):
        for j in range(i + 1, 13):
            if abs(kg[i][j] - kg[j][i]) > 1e-6 and (abs(kg[i][j] / kg[i][i]) > 1e-6 or abs(kg[j][i] / kg[i][i]) > 1e-6):
                print(f"geometric_K: matriz de rigidez no simétrica ... kg[{i}][{j}] = {kg[i][j]:.6e}, kg[{j}][{i}] = {kg[j][i]:.6e}")
                kg[i][j] = kg[j][i] = 0.5 * (kg[i][j] + kg[j][i])

    k[1:13, 1:13] += kg[1:13, 1:13]

def solve_system(K, D, F, R, DoF, q, r, ok, verbose, rms_resid):
    """Resuelve el sistema {F} = [K]{D} mediante descomposición LDL."""
    diag = np.zeros(DoF + 1)
    verbose = 0

    ldl_dcmp_pm(K, DoF, diag, F, D, R, q, r, 1, 0, ok)
    if ok[0] >= 0:
        ldl_dcmp_pm(K, DoF, diag, F, D, R, q, r, 0, 1, ok)
        rms_resid[0] = 1
        ok[0] = 1
        while ok[0]:
            ldl_mprove_pm(K, DoF, diag, F, D, R, q, r, rms_resid, ok)

def equilibrium_error(dF, F, K, D, DoF, q, r):
    """Calcula el error de equilibrio {dF} = {F} - [K]{D}."""
    ss_dF = 0.0
    ss_F = 0.0

    for i in range(1, DoF + 1):
        errF = 0.0
        if q[i]:
            errF = F[i]
            for j in range(1, DoF + 1):
                if q[j]:
                    errF -= K[i][j] if i <= j else K[j][i]
                if r[j]:
                    errF -= K[i][j] * D[j]
        dF[i] = errF

    for i in range(1, DoF + 1):
        if q[i]:
            ss_dF += dF[i] * dF[i]
            ss_F += F[i] * F[i]

    return np.sqrt(ss_dF) / np.sqrt(ss_F) if ss_F > 0 else 0.0

def element_end_forces(Q, nE, xyz, L, Le, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, p, eqF_temp, eqF_mech, D, shear, geom, axial_strain_warning):
    """Evalúa las fuerzas en los extremos de todos los elementos."""
    s = np.zeros(13)
    axial_strain_warning[0] = 0

    for m in range(1, nE + 1):
        frame_element_force(s, xyz, L[m], Le[m], N1[m], N2[m], Ax[m], Asy[m], Asz[m], Jx[m], Iy[m], Iz[m], E[m], G[m], p[m], eqF_temp[m], eqF_mech[m], D, shear, geom, axial_strain)
        Q[m, 1:13] = s[1:13]
        if abs(axial_strain) > 0.001:
            axial_strain_warning[0] += 1

def frame_element_force(s, xyz, L, Le, n1, n2, Ax, Asy, Asz, J, Iy, Iz, E, G, p, f_t, f_m, D, shear, geom, axial_strain):
    """Evalúa las fuerzas en los extremos en coordenadas locales."""
    t1, t2, t3, t4, t5, t6, t7, t8, t9 = coord_trans(xyz, L, n1, n2, p)
    
    n1 = 6 * (n1 - 1)
    n2 = 6 * (n2 - 1)
    
    d1, d2, d3, d4, d5, d6 = D[n1 + 1:n1 + 7]
    d7, d8, d9, d10, d11, d12 = D[n2 + 1:n2 + 7]
    
    Ksy = 12.0 * E * Iz / (G * Asy * Le * Le) if shear else 0.0
    Ksz = 12.0 * E * Iy / (G * Asz * Le * Le) if shear else 0.0
    Dsy = (1 + Ksy) ** 2 if shear else 1.0
    Dsz = (1 + Ksz) ** 2 if shear else 1.0

    delta = (d7 - d1) * t1 + (d8 - d2) * t2 + (d9 - d3) * t3
    axial_strain[0] = delta / Le

    T = 0.0
    s[1] = -(Ax * E / Le) * delta
    if geom:
        T = -s[1]

    s[2] = -(12.0 * E * Iz / (Le * Le * Le * (1.0 + Ksy)) + T / L * (1.2 + 2.0 * Ksy + Ksy * Ksy) / Dsy) * ((d7 - d1) * t4 + (d8 - d2) * t5 + (d9 - d3) * t6) + (6.0 * E * Iz / (Le * Le * (1.0 + Ksy)) + T / 10.0 / Dsy) * ((d4 + d10) * t7 + (d5 + d11) * t8 + (d6 + d12) * t9)
    s[3] = -(12.0 * E * Iy / (Le * Le * Le * (1.0 + Ksz)) + T / L * (1.2 + 2.0 * Ksz + Ksz * Ksz) / Dsz) * ((d7 - d1) * t7 + (d8 - d2) * t8 + (d9 - d3) * t9) - (6.0 * E * Iy / (Le * Le * (1.0 + Ksz)) + T / 10.0 / Dsz) * ((d4 + d10) * t4 + (d5 + d11) * t5 + (d6 + d12) * t6)
    s[4] = -(G * J / Le) * ((d10 - d4) * t1 + (d11 - d5) * t2 + (d12 - d6) * t3)
    s[5] = (6.0 * E * Iy / (Le * Le * (1.0 + Ksz)) + T / 10.0 / Dsz) * ((d7 - d1) * t7 + (d8 - d2) * t8 + (d9 - d3) * t9) + ((4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz)) + T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz) * (d4 * t4 + d5 * t5 + d6 * t6) + ((2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz)) - T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz) * (d10 * t4 + d11 * t5 + d12 * t6)
    s[6] = -(6.0 * E * Iz / (Le * Le * (1.0 + Ksy)) + T / 10.0 / Dsy) * ((d7 - d1) * t4 + (d8 - d2) * t5 + (d9 - d3) * t6) + ((4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy)) + T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy) * (d4 * t7 + d5 * t8 + d6 * t9) + ((2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy)) - T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy) * (d10 * t7 + d11 * t8 + d12 * t9)
    s[7] = -s[1]
    s[8] = -s[2]
    s[9] = -s[3]
    s[10] = -s[4]
    s[11] = (6.0 * E * Iy / (Le * Le * (1.0 + Ksz)) + T / 10.0 / Dsz) * ((d7 - d1) * t7 + (d8 - d2) * t8 + (d9 - d3) * t9) + ((4.0 + Ksz) * E * Iy / (Le * (1.0 + Ksz)) + T * L * (2.0 / 15.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz) * (d10 * t4 + d11 * t5 + d12 * t6) + ((2.0 - Ksz) * E * Iy / (Le * (1.0 + Ksz)) - T * L * (1.0 / 30.0 + Ksz / 6.0 + Ksz * Ksz / 12.0) / Dsz) * (d4 * t4 + d5 * t5 + d6 * t6)
    s[12] = -(6.0 * E * Iz / (Le * Le * (1.0 + Ksy)) + T / 10.0 / Dsy) * ((d7 - d1) * t4 + (d8 - d2) * t5 + (d9 - d3) * t6) + ((4.0 + Ksy) * E * Iz / (Le * (1.0 + Ksy)) + T * L * (2.0 / 15.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy) * (d10 * t7 + d11 * t8 + d12 * t9) + ((2.0 - Ksy) * E * Iz / (Le * (1.0 + Ksy)) - T * L * (1.0 / 30.0 + Ksy / 6.0 + Ksy * Ksy / 12.0) / Dsy) * (d4 * t7 + d5 * t8 + d6 * t9)

    f = f_t + f_m
    s[1] -= (f[1] * t1 + f[2] * t2 + f[3] * t3)
    s[2] -= (f[1] * t4 + f[2] * t5 + f[3] * t6)
    s[3] -= (f[1] * t7 + f[2] * t8 + f[3] * t9)
    s[4] -= (f[4] * t1 + f[5] * t2 + f[6] * t3)
    s[5] -= (f[4] * t4 + f[5] * t5 + f[6] * t6)
    s[6] -= (f[4] * t7 + f[5] * t8 + f[6] * t9)
    s[7] -= (f[7] * t1 + f[8] * t2 + f[9] * t3)
    s[8] -= (f[7] * t4 + f[8] * t5 + f[9] * t6)
    s[9] -= (f[7] * t7 + f[8] * t8 + f[9] * t9)
    s[10] -= (f[10] * t1 + f[11] * t2 + f[12] * t3)
    s[11] -= (f[10] * t4 + f[11] * t5 + f[12] * t6)
    s[12] -= (f[10] * t7 + f[11] * t8 + f[12] * t9)

def compute_reaction_forces(R, F, K, D, DoF, r):
    """Calcula las fuerzas de reacción que satisfacen el equilibrio."""
    for i in range(1, DoF + 1):
        R[i] = 0
        if r[i]:
            R[i] = -F[i]
            for j in range(1, DoF + 1):
                R[i] += K[i][j] * D[j]

def assemble_M(M, DoF, nN, nE, xyz, r, L, N1, N2, Ax, Jx, Iy, Iz, p, d, EMs, NMs, NMx, NMy, NMz, lump, debug):
    """Ensambla la matriz de masa global a partir de masas e inercias de elementos."""
    m = np.zeros((13, 13))
    ind = np.zeros((13, nE + 1), dtype=int)

    M[1:DoF+1, 1:DoF+1] = 0.0

    for i in range(1, nE + 1):
        ind[1][i] = 6 * N1[i] - 5
        ind[7][i] = 6 * N2[i] - 5
        ind[2][i] = ind[1][i] + 1
        ind[8][i] = ind[7][i] + 1
        ind[3][i] = ind[1][i] + 2
        ind[9][i] = ind[7][i] + 2
        ind[4][i] = ind[1][i] + 3
        ind[10][i] = ind[7][i] + 3
        ind[5][i] = ind[1][i] + 4
        ind[11][i] = ind[7][i] + 4
        ind[6][i] = ind[1][i] + 5
        ind[12][i] = ind[7][i] + 5

    for i in range(1, nE + 1):
        if lump:
            lumped_M(m, xyz, L[i], N1[i], N2[i], Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i])
        else:
            consistent_M(m, xyz, r, L[i], N1[i], N2[i], Ax[i], Jx[i], Iy[i], Iz[i], p[i], d[i], EMs[i])
        
        if debug:
            mass_fn = f"m_{i:3d}"
            save_dmatrix(mass_fn, m, 1, 12, 1, 12, 0, "w")
        
        for l in range(1, 13):
            ii = ind[l][i]
            for ll in range(1, 13):
                jj = ind[ll][i]
                M[ii][jj] += m[l][ll]

    for j in range(1, nN + 1):
        i = 6 * (j - 1)
        M[i + 1][i + 1] += NMs[j]
        M[i + 2][i + 2] += NMs[j]
        M[i + 3][i + 3] += NMs[j]
        M[i + 4][i + 4] += NMx[j]
        M[i + 5][i + 5] += NMy[j]
        M[i + 6][i + 6] += NMz[j]

    for i in range(1, DoF + 1):
        if M[i][i] <= 0.0:
            print(f"error: Matriz de masa no definida positiva. M[{i}][{i}] = {M[i][i]}")

def lumped_M(m, xyz, L, n1, n2, Ax, J, Iy, Iz, p, d, EMs):
    """Calcula la matriz de masa concentrada en coordenadas globales."""
    t1, t2, t3, t4, t5, t6, t7, t8, t9 = coord_trans(xyz, L, n1, n2, p)
    
    t = (d * Ax * L + EMs) / 2.0
    ry = d * Iy * L / 2.0
    rz = d * Iz * L / 2.0
    po = d * L * J / 2.0

    m[1:13, 1:13] = 0.0
    m[1][1] = m[2][2] = m[3][3] = m[7][7] = m[8][8] = m[9][9] = t
    m[4][4] = m[10][10] = po * t1 * t1 + ry * t4 * t4 + rz * t7 * t7
    m[5][5] = m[11][11] = po * t2 * t2 + ry * t5 * t5 + rz * t8 * t8
    m[6][6] = m[12][12] = po * t3 * t3 + ry * t6 * t6 + rz * t9 * t9
    m[4][5] = m[5][4] = m[10][11] = m[11][10] = po * t1 * t2 + ry * t4 * t5 + rz * t7 * t8
    m[4][6] = m[6][4] = m[10][12] = m[12][10] = po * t1 * t3 + ry * t4 * t6 + rz * t7 * t9
    m[5][6] = m[6][5] = m[11][12] = m[12][11] = po * t2 * t3 + ry * t5 * t6 + rz * t8 * t9

def consistent_M(m, xyz, r, L, n1, n2, Ax, J, Iy, Iz, p, d, EMs):
    """Calcula la matriz de masa consistente en coordenadas globales."""
    t1, t2, t3, t4, t5, t6, t7, t8, t9 = coord_trans(xyz, L, n1, n2, p)
    
    t = d * Ax * L
    ry = d * Iy
    rz = d * Iz
    po = d * J * L

    m[1:13, 1:13] = 0.0
    m[1][1] = m[7][7] = t / 3.0
    m[2][2] = m[8][8] = 13.0 * t / 35.0 + 6.0 * rz / (5.0 * L)
    m[3][3] = m[9][9] = 13.0 * t / 35.0 + 6.0 * ry / (5.0 * L)
    m[4][4] = m[10][10] = po / 3.0
    m[5][5] = m[11][11] = t * L * L / 105.0 + 2.0 * L * ry / 15.0
    m[6][6] = m[12][12] = t * L * L / 105.0 + 2.0 * L * rz / 15.0
    m[5][3] = m[3][5] = -11.0 * t * L / 210.0 - ry / 10.0
    m[6][2] = m[2][6] = 11.0 * t * L / 210.0 + rz / 10.0
    m[7][1] = m[1][7] = t / 6.0
    m[8][6] = m[6][8] = 13.0 * t * L / 420.0 - rz / 10.0
    m[9][5] = m[5][9] = -13.0 * t * L / 420.0 + ry / 10.0
    m[10][4] = m[4][10] = po / 6.0
    m[11][3] = m[3][11] = 13.0 * t * L / 420.0 - ry / 10.0
    m[12][2] = m[2][12] = -13.0 * t * L / 420.0 + rz / 10.0
    m[11][9] = m[9][11] = 11.0 * t * L / 210.0 + ry / 10.0
    m[12][8] = m[8][12] = -11.0 * t * L / 210.0 - rz / 10.0
    m[8][2] = m[2][8] = 9.0 * t / 70.0 - 6.0 * rz / (5.0 * L)
    m[9][3] = m[3][9] = 9.0 * t / 70.0 - 6.0 * ry / (5.0 * L)
    m[11][5] = m[5][11] = -L * L * t / 140.0 - ry * L / 30.0
    m[12][6] = m[6][12] = -L * L * t / 140.0 - rz * L / 30.0

    for i in range(1, 4):
        m[i][i] += 0.5 * EMs
    for i in range(7, 10):
        m[i][i] += 0.5 * EMs

    atma(t1, t2, t3, t4, t5, t6, t7, t8, t9, m, r[n1], r[n2])

    for i in range(1, 13):
        for j in range(i + 1, 13):
            if abs(m[i][j] - m[j][i]) > 1e-6 and (abs(m[i][j] / m[i][i]) > 1e-6 or abs(m[j][i] / m[i][i]) > 1e-6):
                print(f"consistent_M: matriz de masa no simétrica ... m[{i}][{j}] = {m[i][j]:.6e}, m[{j}][{i}] = {m[j][i]:.6e}")
                m[i][j] = m[j][i] = 0.5 * (m[i][j] + m[j][i])

def static_condensation(A, N, c, n, Ac, verbose):
    """Condensa estáticamente la matriz de rigidez de NxN a nxn."""
    r = np.zeros(N - n + 1, dtype=int)
    Arr = np.zeros((N - n + 1, N - n + 1))
    Arc = np.zeros((N - n + 1, n + 1))

    k = 1
    for i in range(1, N + 1):
        ok = True
        for j in range(1, n + 1):
            if c[j] == i:
                ok = False
                break
        if ok:
            r[k] = i
            k += 1

    for i in range(1, N - n + 1):
        for j in range(i, N - n + 1):
            ri = r[i]
            rj = r[j]
            if ri <= rj:
                Arr[j][i] = Arr[i][j] = A[ri][rj]

    for i in range(1, N - n + 1):
        for j in range(1, n + 1):
            ri = r[i]
            cj = c[j]
            Arc[i][j] = A[ri][cj] if ri < cj else A[cj][ri]

    xtinvAy(Arc, Arr, Arc, N - n, n, Ac, verbose)

    for i in range(1, n + 1):
        for j in range(i, n + 1):
            ci = c[i]
            cj = c[j]
            if ci <= cj:
                Ac[j][i] = Ac[i][j] = A[ci][cj] - Ac[i][j]

def paz_condensation(M, K, N, c, n, Mc, Kc, w2, verbose):
    """Realiza condensación dinámica de Paz para matrices de masa y rigidez."""
    r = np.zeros(N - n + 1, dtype=int)
    Drr = np.zeros((N - n + 1, N - n + 1))
    Drc = np.zeros((N - n + 1, n + 1))
    invDrrDrc = np.zeros((N - n + 1, n + 1))
    T = np.zeros((N + 1, n + 1))

    w2 = 4.0 * pi * pi * w2 * w2

    k = 1
    for i in range(1, N + 1):
        ok = True
        for j in range(1, n + 1):
            if c[j] == i:
                ok = False
                break
        if ok:
            r[k] = i
            k += 1

    for i in range(1, N - n + 1):
        for j in range(1, N - n + 1):
            ri = r[i]
            rj = r[j]
            if ri <= rj:
                Drr[j][i] = Drr[i][j] = K[ri][rj] - w2 * M[ri][rj]
            else:
                Drr[j][i] = Drr[i][j] = K[rj][ri] - w2 * M[rj][ri]

    for i in range(1, N - n + 1):
        for j in range(1, n + 1):
            ri = r[i]
            cj = c[j]
            Drc[i][j] = K[ri][cj] - w2 * M[ri][cj] if ri < cj else K[cj][ri] - w2 * M[cj][ri]

    ok = [0]
    invAB(Drr, Drc, N - n, n, invDrrDrc, ok, verbose)

    for i in range(1, n + 1):
        T[c[i]][i] = 1.0
    for i in range(1, N - n + 1):
        for j in range(1, n + 1):
            T[r[i]][j] = -invDrrDrc[i][j]

    xtAx(K, T, Kc, N, n)
    xtAx(M, T, Mc, N, n)

def modal_condensation(M, K, N, R, p, n, Mc, Kc, V, f, m, verbose):
    """Realiza condensación dinámica para frecuencias y modos específicos."""
    P = np.zeros((n + 1, n + 1))
    invP = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            P[i][j] = V[p[i]][m[j]]

    pseudo_inv(P, invP, n, n, 1e-9, verbose)

    traceM = sum(M[i][i] for i in range(1, N + 1) if not R[i])
    traceMc = 0.0

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            Aij = sum(invP[k][i] * invP[k][j] for k in range(1, n + 1))
            Mc[i][j] = Aij
        traceMc += Mc[i][i]

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            Aij = sum(invP[k][i] * 4.0 * pi * pi * f[m[k]] * f[m[k]] * invP[k][j] for k in range(1, n + 1))
            Kc[i][j] = Aij

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            Mc[i][j] *= (traceM / traceMc)
            Kc[i][j] *= (traceM / traceMc)

def deallocate(nN, nE, nL, nF, nU, nW, nP, nT, DoF, nM, xyz, rj, L, Le, N1, N2, q, r, Ax, Asy, Asz, J, Iy, Iz, E, G, p, U, W, P, T, Dp, F_mech, F_temp, eqF_mech, eqF_temp, F, dF, K, Q, D, dD, R, dR, d, EMs, NMs, NMx, NMy, NMz, c, m, pkNx, pkVy, pkVz, pkTx, pkMy, pkMz, pkDx, pkDy, pkDz, pkRx, pkSy, pkSz):
    """Libera la memoria asignada para todas las estructuras de datos."""
    xyz = None
    rj = None
    L = None
    Le = None
    N1 = None
    N2 = None
    q = None
    r = None
    Ax = None
    Asy = None
    Asz = None
    J = None
    Iy = None
    Iz = None
    E = None
    G = None
    p = None
    U = None
    W = None
    P = None
    T = None
    Dp = None
    F_mech = None
    F_temp = None
    eqF_mech = None
    eqF_temp = None
    F = None
    dF = None
    K = None
    Q = None
    D = None
    dD = None
    R = None
    dR = None
    d = None
    EMs = None
    NMs = None
    NMx = None
    NMy = None
    NMz = None
    c = None
    m = None
    pkNx = None
    pkVy = None
    pkVz = None
    pkTx = None
    pkMy = None
    pkMz = None
    pkDx = None
    pkDy = None
    pkDz = None
    pkRx = None
    pkSy = None
    pkSz = None
```