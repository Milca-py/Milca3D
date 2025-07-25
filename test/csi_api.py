import comtypes.client


def ultima_version():
    helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    helper_interface = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    # Crea un nuevo modelo con la ultima version SAP2000
    mySapObject = helper_interface.CreateObjectProgID(
        "CSI.SAP2000.API.SapObject")
    mySapObject.ApplicationStart()  # Inicia la app de SAP2000
    SapModel = mySapObject.SapModel
    SapModel.File.NewBlank()  # Prepara SAP2000 para crear un nuevo modelo
    return SapModel


def especifica_version(version):
    program_path = rf'C:\Program Files\Computers and Structures\SAP2000 {version}\SAP2000.exe'
    helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    helper_interface = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    # Especifica la ruta al ejecutable de SAP2000
    mySapObject = helper_interface.CreateObject(program_path)
    mySapObject.ApplicationStart()  # Inicia la app de SAP2000
    SapModel = mySapObject.SapModel
    SapModel.InitializeNewModel()  # Prepara SAP2000 para crear un nuevo modelo
    SapModel.File.NewBlank()
    return SapModel


def crear_nuevo_model():
    """
    Conecta con una instancia abierta de SAP2000 y devuelve el objeto SapModel.

    Retorna:
    SapModel : objeto
        El objeto SapModel vinculado a la instancia abierta de SAP2000.

    Excepciones:
    Lanza una excepción si no se puede establecer la conexión.
    """
    try:
        # Crear objeto helper para SAP2000
        helper = comtypes.client.CreateObject('SAP2000v1.Helper')
        # Obtener la interfaz de cHelper
        helper_interface = helper.QueryInterface(
            comtypes.gen.SAP2000v1.cHelper)
        # Obtener el objeto SapObject del modelo abierto
        mySapObject = helper_interface.GetObject("CSI.SAP2000.API.SapObject")
        # Obtener y devolver el objeto SapModel
        SapModel = mySapObject.SapModel
        print("Conexión exitosa con SAP2000.")
        return SapModel
    except Exception as e:
        print(f"Error al conectar con SAP2000: {e}")
        raise


def vincular_model():
    helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    helper_interface = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    mySapObject = helper_interface.GetObject(
        "CSI.SAP2000.API.SapObject")  # Vincula a un modelo abierto
    SapModel = mySapObject.SapModel
    return SapModel

    # helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    # helper_interface = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    # mySapObject = helper_interface.GetObject(
    #     "CSI.SAP2000.API.SapObject")
    # SapModel = mySapObject.SapModel
    # return SapModel


def vincular_and_new_model():
    # TOMAR OBJETO SAP2000 Y CREAR NUEVO MODELO EN BLANCO
    helper = comtypes.client.CreateObject('SAP2000v1.Helper')
    helper_interface = helper.QueryInterface(comtypes.gen.SAP2000v1.cHelper)
    mySapObject = helper_interface.GetObject("CSI.SAP2000.API.SapObject")
    SapModel = mySapObject.SapModel
    SapModel.File.NewBlank()  # Prepara SAP2000 para crear un nuevo modelo
    return SapModel


def StarSAP2000(option):
    """
    Inicia una instancia de SAP2000 según la opción seleccionada.

    Parámetros:
    option (str): Opción para iniciar SAP2000. Puede ser:
        - 'crear_nuevo': Vincular o crear un nuevo modelo.
        - 'ultimo': Vincular a la última versión abierta.
        - 'especifico': Vincular a una versión específica.
        - 'vincular': Vincular a una instancia existente.

    Retorna:
    SapModel (obj): Objeto vinculado de SAP2000.
    """
    if option == 'crear_nuevo':
        SapModel = vincular_and_new_model()
    elif option == 'ultimo':
        SapModel = ultima_version()
    elif option == 'especifico':
        version = 22  # Aquí puedes pedir o definir la versión deseada
        SapModel = especifica_version(version)
    elif option == 'vincular':
        SapModel = vincular_model()
    else:
        raise ValueError(
            "Opción no válida. Elija entre 'crear_nuevo', 'ultimo', 'especifico' o 'vincular'.")

    return SapModel


def cambiar_unidades(SapModel, unidad_code):
    """
    Cambia las unidades de SAP2000 según el código proporcionado.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    unidad_code : int
        El código que representa la unidad deseada:
        1 - lb_in_F
        2 - lb_ft_F
        3 - kip_in_F
        4 - kip_ft_F
        5 - kN_mm_C
        6 - kN_m_C
        7 - kgf_mm_C
        8 - kgf_m_C
        9 - N_mm_C
        10 - N_m_C
        11 - Ton_mm_C
        12 - Ton_m_C
        13 - kN_cm_C
        14 - kgf_cm_C
        15 - N_cm_C
        16 - Ton_cm_C
    """
    # Asignar las unidades de SAP2000 según el código proporcionado
    SapModel.SetPresentUnits(unidad_code)

    # Diccionario opcional de unidades para mostrar al usuario (en caso de ser necesario)
    unidades = {
        1: "lb/in (Force)",
        2: "lb/ft (Force)",
        3: "kip/in (Force)",
        4: "kip/ft (Force)",
        5: "kN/mm (C)",
        6: "kN/m (C)",
        7: "kgf/mm (C)",
        8: "kgf/m (C)",
        9: "N/mm (C)",
        10: "N/m (C)",
        11: "Ton/mm (C)",
        12: "Ton/m (C)",
        13: "kN/cm (C)",
        14: "kgf/cm (C)",
        15: "N/cm (C)",
        16: "Ton/cm (C)"
    }

    print(f"Unidades cambiadas a: {unidades[unidad_code]}.")


def crear_material_concreto(SapModel,
                            Elas,
                            nombre_material,
                            f_c=280,
                            gamma_c=2.4,
                            poisson=0.15,
                            deltaT=1e-5,
                            factor_resistencia=0.85,
                            curva_tipo=1,
                            hist_resistencia=0,
                            deformacion_resistencia=0.0021,
                            deformacion_ultima=0.003):
    """
    Crea un material de concreto en SAP2000 con valores predeterminados para concreto 280 (si no se especifica el módulo de elasticidad).

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    nombre_material : str
        Nombre del material.
    f_c : float, opcional
        Resistencia a compresión del concreto (kg/cm²), por defecto 280 kg/cm².
    gamma_c : float, opcional
        Peso específico del concreto (kg/cm³), por defecto 2.4 kg/cm³.
    poisson : float, opcional
        Coeficiente de Poisson del concreto, por defecto 0.2.
    deltaT : float, opcional
        Coeficiente térmico (1/°C), por defecto 1e-5.
    Elas : float, opcional
        Módulo de elasticidad (kg/cm²), si no se especifica se calcula como 15000 * (f_c ** 0.5).
    factor_resistencia : float, opcional
        Factor de resistencia según las normas del material (ACI), por defecto 0.85.
    curva_tipo : int, opcional
        Tipo de curva tensión-deformación, por defecto 1 (Mander).
    hist_resistencia : float, opcional
        Histéresis del material, por defecto 0 (Takeda).
    deformacion_resistencia : float, opcional
        Deformación a resistencia (adimensional), por defecto 0.0021.
    deformacion_ultima : float, opcional
        Deformación última (adimensional), por defecto 0.003.
    """

    # Si no se proporciona un módulo de elasticidad, se calcula basado en f_c
    if Elas is None:
        # Módulo de elasticidad (kg/cm²) basado en f_c
        Elas = 15100 * (f_c ** 0.5)

    # Crear el material en SAP2000
    SapModel.PropMaterial.SetMaterial(nombre_material, 2)  # 2 = Concreto

    # Asignar propiedades mecánicas isotrópicas
    SapModel.PropMaterial.SetMPIsotropic(
        nombre_material,    # Nombre del material
        Elas,               # Módulo de elasticidad (kg/cm²)
        poisson,            # Coeficiente de Poisson
        deltaT              # Coeficiente térmico (1/°C)
    )

    # Asignar peso específico
    SapModel.PropMaterial.SetWeightAndMass(
        nombre_material,    # Nombre del material
        1,                  # 1 = Peso específico
        gamma_c             # Peso específico (kg/cm³)
    )

    # Definir propiedades adicionales del concreto según ACI
    SapModel.PropMaterial.SetOConcrete(
        nombre_material,       # Nombre del material
        f_c,                   # Resistencia a compresión (kg/cm²)
        False,                 # Concreto normal (False = No es liviano)
        factor_resistencia,    # Factor de resistencia (ACI)
        curva_tipo,            # Tipo de curva tensión-deformación
        hist_resistencia,      # Histéresis
        deformacion_resistencia,  # Deformación a resistencia
        deformacion_ultima     # Deformación última
    )

    print("Material de concreto creado exitosamente.")


def crear_seccion_rectangular(SapModel, nombre_seccion, material, peralte, ancho):
    """
    Crea una sección rectangular en SAP2000 para un elemento tipo 'Frame'.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    nombre_seccion : str
        Nombre de la sección.
    material : str
        Nombre del material asignado.
    peralte : float
        Peralte (altura) de la sección en cm.
    ancho : float
        Ancho de la sección en cm.
    """

    # Asignar la sección rectangular en SAP2000
    SapModel.PropFrame.SetRectangle(
        nombre_seccion,  # Nombre de la sección
        material,        # Material asignado
        peralte,         # Peralte (altura) en cm
        ancho            # Ancho en cm
    )
    print("Sección rectangular creada exitosamente.")


def crear_seccion_circular(SapModel, nombre_seccion, material, diametro):
    """
    Crea una sección circular en SAP2000 para un elemento tipo 'Frame'.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    nombre_seccion : str
        Nombre de la sección.
    material : str
        Nombre del material asignado.
    diametro : float
        Diámetro de la sección en cm.
    """

    # Asignar la sección circular en SAP2000
    SapModel.PropFrame.SetCircle(
        nombre_seccion,  # Nombre de la sección
        material,        # Material asignado
        diametro         # Diámetro en cm
    )
    print("Sección circular creada exitosamente.")


def crear_seccion_I(SapModel, nombre_seccion, material, t3, t2, tf, tw, t2b, tfb):
    """
    Crea una sección I en SAP2000 para un elemento tipo 'Frame'.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    nombre_seccion : str
        Nombre de la sección.
    material : str
        Nombre del material asignado.
    t3 : float
        Profundidad de la sección. [L]
    t2 : float
        Ancho de la brida superior. [L]
    tf : float
        Espesor de la brida superior. [L]
    tw : float
        Espesor del alma. [L]
    t2b : float
        Ancho de la brida inferior. [L]
    tfb : float
        Espesor de la brida inferior. [L]
    """

    # Asignar la sección I en SAP2000
    SapModel.PropFrame.SetISection(
        nombre_seccion,  # Nombre de la sección
        material,        # Material asignado
        t3,              # Profundidad de la sección. [L]
        t2,              # Ancho de la brida superior. [L]
        tf,              # Espesor de la brida superior. [L]
        tw,              # Espesor del alma. [L]
        t2b,             # Ancho de la brida inferior. [L]
        tfb              # Espesor de la brida inferior. [L]
    )
    print("Sección I creada exitosamente.")

def crear_nodos(SapModel, coord_list, node_names):
    """
    Crea nodos en SAP2000 con las coordenadas especificadas.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    coord_list : list
        Lista de coordenadas [X, Y, Z] de los nodos a crear.
    node_names : list
        Lista de nombres de los nodos a crear.
    """
    for coord, name in zip(coord_list, node_names):
        # Crear nodo en SAP2000
        SapModel.PointObj.AddCartesian(
            coord[0],   # Coordenada X
            coord[1],   # Coordenada Y
            coord[2],   # Coordenada Z
            name        # Nombre del nodo
        )
        print(f"Nodo {name} creado en las coordenadas: {coord}")





def crear_frame_por_coord(SapModel, Xi, Yi, Zi, Xj, Yj, Zj, PropName, FrameName=""):
    """
    Crea un frame en SAP2000 usando las coordenadas iniciales y finales.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    Xi, Yi, Zi : float
        Coordenadas iniciales del frame (en cm).
    Xj, Yj, Zj : float
        Coordenadas finales del frame (en cm).
    PropName : str
        Nombre de la propiedad de la sección (material y dimensiones).
    FrameName : str, opcional
        Nombre del frame (si se deja vacío, SAP2000 asignará un nombre automáticamente).
    """

    # Agregar el frame en SAP2000
    FrameName1, ret = SapModel.FrameObj.AddByCoord(
        Xi, Yi, Zi,  # Coordenadas iniciales
        Xj, Yj, Zj,  # Coordenadas finales
        FrameName,    # Nombre del frame (vacío para que SAP lo asigne)
        PropName      # Propiedad de la sección
    )

    if ret == 0:
        print(f"Frame creado exitosamente: {FrameName1}")
    else:
        print("Error al crear el frame.")


def crear_frame_por_nodos(SapModel, NodeI, NodeJ, PropName, FrameName=""):
    """
    Crea un frame en SAP2000 usando las etiquetas de los nodos inicial y final.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    NodeI : str
        Etiqueta del nodo inicial.
    NodeJ : str
        Etiqueta del nodo final.
    PropName : str
        Nombre de la propiedad de la sección (material y dimensiones).
    FrameName : str, opcional
        Nombre del frame (si se deja vacío, SAP2000 asignará un nombre automáticamente).
    """
    FrameName1, ret = SapModel.FrameObj.AddByPoint(
        NodeI,    # Etiqueta del nodo inicial
        NodeJ,    # Etiqueta del nodo final
        FrameName,  # Nombre del frame (vacío para que SAP lo asigne)
        PropName    # Propiedad de la sección
    )

    if ret == 0:
        print(f"Frame creado exitosamente: {FrameName1}")
    else:
        print("Error al crear el frame.")



def asignar_restricciones(SapModel, nodes, restraints_list):
    """
    Asigna restricciones a los nodos especificados en SAP2000.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    nodes : list
        Lista de etiquetas de nodos a los que se les asignarán las restricciones.
    restraints_list : list
        Lista de booleans con las restricciones para cada grado de libertad:
        [X, Y, Z, RotX, RotY, RotZ]
        True significa restringido, False significa libre.
    """
    for idnode, typeRest in zip(nodes, restraints_list):
        # Asignar las restricciones correspondientes a cada nodo
        SapModel.PointObj.SetRestraint(str(idnode), typeRest)
        print(f"Restricciones asignadas al nodo {idnode}: {typeRest}")

    # # Ejemplo de uso:
    # nodes = ['1', '2', '3']  # Lista de nodos
    # restraints_list = [True, True, True, True, True, True]  # Fijar todos los grados de libertad

    # # Llamada a la función para asignar restricciones
    # asignar_restricciones(SapModel, nodes, restraints_list)


def agregar_patrones_carga(SapModel, names, types, factors):
    """
    Agrega patrones de carga al modelo de SAP2000.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    names : list
        Lista con los nombres de los patrones de carga. Ejemplos: ["DEAD", "SUPERDEAD", "LIVE", ...]
    types : list
        Lista con los tipos de carga. Ejemplos: [1, 2, 3, ...] (1 = Dead Load, 2 = Super Dead Load, 3 = Live Load, etc.)
    factors : list
        Lista con los factores de peso propio. Ejemplos: [1.0, 0.0, 1.0, ...] (1.0 para Dead Load, 0.0 para otros tipos de carga)

    Ejemplo de uso:
    names = ["DEAD", "SUPERDEAD", "LIVE", "REDUCELIVE", "QUAKE", "WIND", "SNOW", "OTHER", "MOVE", "TEMPERATURE", "ROOFLIVE", "NOTIONAL", "PATTERNLIVE", "WAVE", "BRAKING", "CENTRIFUGAL", "FRICTION", "ICE", "WINDONLIVELOAD", "HORIZONTALEARTHPRESSURE", "VERTICALEARTHPRESSURE", "EARTHSURCHARGE", "DOWNDRAG", "VEHICLECOLLISION", "VESSELCOLLISION", "TEMPERATUREGRADIENT", "SETTLEMENT", "SHRINKAGE", "CREEP", "WATERLOADPRESSURE", "LIVELOADSURCHARGE", "LOCKEDINFORCES", "PEDESTRIANLL", "PRESTRESS", "HYPERSTATIC", "BOUYANCY", "STREAMFLOW", "IMPACT", "CONSTRUCTION"]
    types = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
    factors = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Llamada a la función para agregar los patrones de carga
    agregar_patrones_carga(SapModel, names, types, factors)
    """
    # Iterar y agregar los patrones de carga en SAP2000
    for i in range(len(names)):
        SapModel.LoadPatterns.Add(
            names[i],     # Nombre del patrón de carga
            types[i],     # Tipo de carga (número según SAP2000)
            factors[i],   # Factor de peso propio
            True          # True = Crear Load Case
        )
        print(f"Patrón de carga '{names[i]}' agregado con éxito.")
    # # Ejemplo de uso:
    # names = ["Dead_Load", "Live_Load"]
    # types = [1, 3]
    # factors = [1.0, 0.0]

    # # Llamada a la función para agregar los patrones de carga
    # agregar_patrones_carga(SapModel, names, types, factors)


def asignar_fuerzas_a_nodos(SapModel, tags, forces, load_pattern):
    """
    Asigna fuerzas a los nodos en SAP2000.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    tags : list
        Lista con las etiquetas de los nodos a los que se les asignará la carga.
    load_pattern : str
        Nombre del patrón de carga que se aplicará en los nodos.
    forces : list
        Lista de listas con las fuerzas a aplicar en los nodos [F1, F2, F3, M1, M2, M3].
    """
    for i in range(len(tags)):
        SapModel.PointObj.SetLoadForce(
            str(tags[i]),          # Etiqueta del nodo
            load_pattern,          # Nombre del patrón de carga
            forces[i],             # Fuerzas [F1, F2, F3, M1, M2, M3]
            Replace=True           # True = Reemplazar fuerzas si ya existen
        )
        print(f"Fuerzas aplicadas al nodo {tags[i]} con éxito.")


def asignar_carga_distribuida_a_frames(SapModel, tagsl, load_pattern, forces, direction=4, dist1=0, dist2=1):
    """
    Asigna una carga distribuida a varios frames en SAP2000.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    tagsl : list
        Lista con las etiquetas de los frames a los que se les asignará la carga distribuida.
    load_pattern : str
        Nombre del patrón de carga que se aplicará en los frames.
    forces : list
        Lista de fuerzas a aplicar en los frames [F1, F2].
    direction : int, opcional
        Dirección de la carga (predeterminado es 10).
    dist1 : float, opcional
        Distancia de inicio de la carga distribuida (predeterminado es 0).
    dist2 : float, opcional
        Distancia final de la carga distribuida (predeterminado es 1).
    """
    # Iterar sobre las etiquetas de los frames y asignar las cargas distribuidas
    for i in range(len(tagsl)):
        SapModel.FrameObj.SetLoadDistributed(
            tagsl[i],           # Etiqueta del frame
            load_pattern,       # Nombre del patrón de carga
            1,                  # 1: Fuerza distribuida / Momento distribuido
            direction,          # Dirección de la carga, [1, 2, 3, X, Y, Z, X_PROY, Y_PROY, Z_PROY, GRAV, GRAV_PROY]
            dist1,              # Distancia inicial
            dist2,              # Distancia final
            forces[i][0],          # Fuerza distribuida en el primer segmento
            forces[i][1],          # Fuerza distribuida en el segundo segmento
            CSys="Global",      # Sistema de coordenadas (Global o Local)
            RelDist=True,       # Distancias relativas (dl/L) o absoluta (dl)
            Replace=True        # Reemplazar la carga en el frame, false para agregar
        )
        print(f"Carga distribuida aplicada al frame {tagsl[i]} con éxito.")

    # # Ejemplo de uso:
    # tagsl = ['1']  # Etiquetas de los frames
    # load_pattern = "CV"  # Tipo de carga (CV: Carga distribuida uniforme)
    # forces = [100] * len(tagsl)  # Fuerzas distribuidas de 100 kgf para todos los frames

    # # Llamada a la función para asignar la carga distribuida
    # asignar_carga_distribuida_a_frames(SapModel, tagsl, load_pattern, forces)


def guardar_y_analizar_modelo(SapModel, model_path):
    """
    Guarda el modelo de SAP2000 en la ruta especificada y luego ejecuta el análisis del modelo.

    Parámetros:
    SapModel : objeto de SAP2000
        El objeto del modelo de SAP2000.
    model_path : str
        Ruta donde se guardará el modelo de SAP2000 (por ejemplo, "C:\\Users\\User\\Desktop\\Modelo.sdb").

    Retorna:
    ret : int
        Código de retorno de la operación (0 si exitosa, otro valor en caso de error).
    """
    # Guardar el modelo
    ret = SapModel.File.Save(model_path)
    if ret != 0:
        print("Error al guardar el modelo.")
        return ret

    # Ejecutar el análisis del modelo
    ret = SapModel.Analyze.RunAnalysis()
    if ret != 0:
        print("Error al ejecutar el análisis del modelo.")
        return ret

    print("Modelo guardado y análisis ejecutado con éxito.")
    return ret

    # # Ejemplo de uso
    # model_path = r"C:\Users\AMILCAR\Desktop\MODELOS\modelo.sdb"  # Ruta para guardar el modelo
    # guardar_y_analizar_modelo(SapModel, model_path)


# # parametros de modelo
# # ===========================================================================================
# option = 'ultimo'  # 'ultimo', 'especifico', 'vincular'

# unidades = 14  # 14: kg, cm, C

# nombre_material = 'Hormigon'
# Elast = 250000
# poisson = 0.15

# nombre_seccion = 'Rectangular'
# material = nombre_material
# peralte = 30
# ancho = 20

# xi = 0
# yi = 0
# zi = 0
# xj = 0
# yj = 0
# zj = 3
# nombre_seccion_frame = nombre_seccion

# nodes = ['1']
# restraints_list = [True, True, True, True, True, True]

# names_Patrones = ['Cargas en nudos']
# types = [3]
# factors = [0]

# tags = ['2']
# load_pattern = names_Patrones[0]
# forces = [[1000, 0, 0, 0, 0, 0]]  # [F1, F2, F3, M1, M2, M3]


# tagsl = ['1']
# load_pattern = names_Patrones[0]
# forces = [100] * len(tagsl)

# # Ruta para guardar el modelo
# model_path = r"C:\Users\AMILCAR\Desktop\MODELOS\modelo.sdb"


# # modelado de la estructura en SAP2000 22
# # ===========================================================================================
# if option == 'ultimo':
#     SapModel = api.ultima_version()

# elif option == 'especifico':
#     program_path = r"C:\Program Files\Computers and Structures\SAP2000 22\SAP2000.exe"
#     SapModel = api.especifica_version(program_path)

# elif option == 'vincular':
#     SapModel = api.vincular_model()
# else:
#     raise ValueError("Opción no válida.")


# api.cambiar_unidades(SapModel, unidades)  # 14: kg, cm, C

# api.crear_material_concreto(SapModel,
#                             nombre_material,
#                             f_c=280,
#                             gamma_c=2.4,
#                             poisson=poisson,
#                             Elas=Elast)


# api.crear_seccion_rectangular(SapModel,
#                               nombre_seccion,
#                               material,
#                               peralte, ancho)


# api.crear_frame_por_coord(SapModel, xi, yi, zi, xj,
#                           yj, zj, nombre_seccion_frame)


# api.asignar_restricciones(SapModel, nodes, restraints_list)


# api.agregar_patrones_carga(SapModel, names_Patrones, types, factors)


# api.asignar_fuerzas_a_nodos(SapModel, tags, load_pattern, forces)


# api.asignar_carga_distribuida_a_frames(SapModel,
#                                        tagsl,
#                                        load_pattern,
#                                        forces)


# api.guardar_y_analizar_modelo(SapModel, model_path)
