#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Filtrar Recursivamente.

Created on Fri Aug 16 10:36:32 2024

@author: rgrimson
"""

from pathlib import Path
from pprint import pprint


from packaging import version
#import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
#import numpy as np
#import seaborn as sns
#from tqdm import tqdm

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
#from shapely import wkt
#from shapely import intersects
#from shapely import union

# Importar condicionalmente unary_union desde shapely (v >= 2.0.0) o desde
#  shapely.ops (v < 2.0.0).
import shapely
if version.parse(shapely.__version__) >= version.parse("2.0.0"):
    from shapely import unary_union  # pylint: disable=no-name-in-module
else:
    from shapely.ops import unary_union  # pylint: disable=no-name-in-module
#import itertools


# guardar_intermedios = True


#%%
# Listar Polígonos
def lpolys(PMP):
    """Generar una lista de polígonos singlepart."""
    # Devuelve una lista con los polígonos simples que componen el P o MP.
    if type(PMP)==Polygon:
        return [PMP]
    elif type(PMP)==MultiPolygon:
        return list(PMP.geoms)
    elif type(PMP) == GeometryCollection:
        return [p for p in list(PMP.geoms) if type(p)==Polygon]
    else:
        a = 0
        b = a/a
        return b


#%%
#  Convertir GeometryCollection a MultiPolygon
def GC2MP(P):
    """Generar un MultiPolyogn a partir de una GeometryCollection."""
    if type(P) == GeometryCollection:
        return MultiPolygon([p for p in list(P.geoms) if type(p)==Polygon])
    else:
        return P


#%%
# Topologia de Polygon, MultiPolygon o GeometryCollection
def topo(P):
    """Generar la topología de un polígono o multipolígono."""
    # Devuelve una descripcion de la topología del polígono o del multip.
    # Es recursiva para multipolígonos
    if P:
        if type(P)==Polygon:
            return [len(P.interiors)]
        elif type(P)==MultiPolygon:
            return [topoP(Q) for Q in list(P.geoms)]
        elif type(P) == GeometryCollection:
            return [p for p in list(P.geoms) if type(p)==Polygon]
        else:
            print("Tipo no reconocido.")
    else:
        return []


#%%
# Topologia de Polygon
def topoP(P):
    """Generar la topologia de un polígono."""
    if P and type(P)==Polygon:
        return len(P.interiors)
    else:
        return []

#%%
# Limpiar interiores de Polygon
def clean_interiors_P(P, eps=0.0001):
    """Limpiar interiores de un polígono."""
    list_interiors = []

    for interior in P.interiors:
        p = Polygon(interior)

        if p.area > eps:
            list_interiors.append(interior)

    return Polygon(P.exterior.coords, holes=list_interiors)


# Limpiar interiores de MultiPolygon
def clean_interiors_MP(P, eps=0.0001):
    """Limpiar interiores de un multipolígono."""
    list_parts = []

    for polygon in P.geoms:
        list_parts.append(clean_interiors_P(polygon,eps))

    return MultiPolygon(list_parts)


# Limpiar interiores
def clean_interiors(PMP, eps=0.0001):
    """Limpiar interiores de un polígono o multipolígono."""
    if type(PMP)==Polygon:
        return clean_interiors_P(PMP,eps)
    elif type(PMP)==MultiPolygon:
        return clean_interiors_MP(PMP,eps)


#%%
# Imprimir Polígono
def plot_polygon(polygon):
    """Imprimir un polígono o multipolígono."""
    # Función auxiliar para imprimir polígonos singlepart.
    def _plot_single_polygon(polygon):
        x, y = polygon.exterior.xy
        plt.plot(x, y)
        for interior in polygon.interiors:
            x, y = interior.xy
            plt.plot(x, y)

    if type(polygon) == Polygon:
        _plot_single_polygon(polygon)
    elif type(polygon) == MultiPolygon:
        for p in polygon.geoms:
            _plot_single_polygon(p)

    plt.axis('scaled')
    plt.show()
    return None


#%%
# Calcular Filtración Recursiva
def calcular_filtracion_recursiva(P, r, cod, r_step=1, verb=0):
    """Calcular la filtración recursiva de P.

    Returns:
        Un Polygon() si P es Polygon() o si se absorbe completamente en
        r + rstep. O una lista `rta` con al menos dos elementos, donde cada
        elemento puede ser un polígono núcleo indivisible de la descomposición
        de P, una tupla de dos elementos donde el primero es un polígono
        núcleo divisible de la descomposición de P, y el segundo elemento es la
        lista `rta` de su división, en forma recursiva.
    """
    # P tiene que ser un polígono singlepart, si no eleva un error.
    if type(P)!=Polygon:
        a = 0
        b = a/a
        return b

    ta = topo(P)

    # Si la topología es vacía (si P es Polygon()) sale de la recursión
    #  devolviendo Polygon().
    if ta==[]:
        return Polygon()

    if verb > 0:
        print(cod,r,ta)


    # Distancia de filtración (no analiza r en el primer paso).
    d = r + r_step

    # Q es la intersección entre P y su buffereado
    #  (prácticamente el mismo buffereado, que debería estar contenido en P
    #  excepto porque tiene vértices que P no tiene).
    Q = ((P.intersection(P.buffer(-d).buffer(d))).buffer(0)).normalize()

    t = topo(Q)

    # Mientras que P y Q tengan misma cantidad de partes:
    while len(t)==len(ta):
        d += r_step
        Q = ((P.intersection(P.buffer(-d).buffer(d))).buffer(0)).normalize()
        t = topo(Q)

    # En este momento, Q tiene más partes que P, o es un polígono vacío (si P
    #  se absorbió sin haberse descompuesto).
    # LQ es una lista de Polygons que componen a Q o [Polygon()].
    # Estos polígonos son los núcleos componentes de P.
    LQ = lpolys(Q)

    # Si Q tiene más de una componente Polygon, descomponer P:
    if len(LQ)>1:
        # `rta`` es una lista con al menos dos elementos.
        # Cada elemento puede ser Q si es indivisible, o una tupla (Q, frQ),
        #  donde frQ es el resultado de filtrar recursivamente a Q.
        rta = []
        # Por cada parte de Q:
        for i, Q in enumerate(LQ):
            # Calcular esta misma filtración para cada núcleo, llevando un
            #  código que concatena índices de cada parte en forma recursiva, y
            #  una distancia inicial que es la distancia en la que se generó
            #  esta parte.
            # La filtración recursiva no analiza la descomposición en d sino en
            #  d + r_step (Q es una componente del buffereado de P en d, por lo
            #  que bufferear a Q en d no debería alterarlo).
            frQ = calcular_filtracion_recursiva(Q, d, cod+str(i), r_step, verb)

            # La filtración recursiva puede devolver un polígono vacío o
            #  una lista `rta`.
            # Si devolvió una lista (Q es divisible a partir de d + r_step):
            if frQ:
                # Agregar un elemento tupla a rta, formada por esta parte
                #  núcleo componente y el resultado de filtrarla recursivamente.
                rta.append((Q, frQ))

            # Si frQ devolvió Polygon() (falsy value)
            else:
                # Agregar esta parte Q (componente núcleo indivisible de P)
                #  a `rta`.
                if verb > 1:
                    print(f'{frQ = }')
                rta.append(Q)

        return rta

    # Si no, LQ tiene una sola parte.
    # Si esa parte no es un polígono vacío (quizás Q tenía partes no poligonales
    #  que se perdieron al convertir a lista?)
    elif LQ[0]:
        if verb > 1:
            print(f'{LQ[0] = }')
        # Filtrar esa parte en forma recursiva
        return calcular_filtracion_recursiva(LQ[0], d, cod, r_step, verb)

    # Si no, LQ == [Polygon()] y sale de la recursión devolviendo Polygon().
    else:
        return Polygon()


#%%
# Calcular Filtración de Polygon
def calcular_filtracion(P, r=0, r_step=1, verb=0):
    """Calcular la filtración de un polígono singlepart.

    Returns:
        Tupla del poligono original y el resultado de llamar a
        `calcular_filtracion_recursiva` (Polygon() si `P` se absorbe
        completamente sin descomponerse, o una lista conteniendo la
        descomposición de P).
    """
    # Limpiar la topología de P y calcular su filtración
    P_limpio = (P.buffer(0)).normalize()
    F = calcular_filtracion_recursiva(P_limpio, r, '', r_step, verb)

    # Si calcular_filtracion_recursiva devuelve Polygon(), entonces P es
    #  indivisible para los parámetros de r, r_step.
    # Devolver P:
    if not F:
        return P

    # Si no, devolver el P limpio y su filtración.
    else:
        return (P_limpio, F)


#%%
# Imprimir resultados de filtración recursiva
def antirecursion(F, cod):
    """Descomponer la recursión de `F` con código inicial `cod`."""
    # Lista de polígonos que se van a plotear al final de las impresiones.
    poligonos = []

    def _antirecursion(F, cod):
        """Función privada que imprime recursivamente la descomposición de F.

        F puede ser un polígono núcleo indivisible de descomposición o una lista
        de dos elementos o más, cada uno correspondiendo a una parte componente.
        Estos elementos pueden ser polígonos si la parte es núcleo indivisible
        de descomposición, o una tupla donde el primer elemento es la parte
        núcleo divisible y el segundo es una lista de dos elementos o más,
        cada uno correspondiendo a una parte componente.
        """
        # Crear una lista base de valores a imprimir.
        base = [f'{cod = }', f'{type(F) = }']
        # Si F es polígono agregar también su topología.
        if type(F) == Polygon:
            base.append(f'{topo(F)}')

        print("; ".join(base) + ".")

        if type(F) == Polygon:
            poligonos.append(F)
            # La única salida de recursión es agregando el polígono que trae F
            #  a la lista de polígonos a imprimir.
            return None

        else:
            for i, f in enumerate(F):
                _antirecursion(f, cod + str(i))

    _antirecursion(F, cod)
    plot_polygon(MultiPolygon(poligonos))
    return None


#%%
# Diccionario de Filtración
def dicc_filtracion(FP, cod=''):
    """Generar un diccionario con los polígonos de una filtración.

    Función recursiva que recorre la filtración y en cada elemento que no
     sea solo un polígono, agrega a un diccionario el polígono núcleo divisible
     e itera en su descomposición.
    """
    # Si el elemento es sólo un polígono (es indivisible), crear un diccionario
    #  con él.
    if type(FP)==Polygon:
        d = {cod:FP}

    # Si no, el elemento es una tupla de dos elementos, un polígono divisible y
    #  su descomposición.
    else:
        # Iniciar el diccionario con el polígono núcleo divisible.
        d = {cod:FP[0]}
        # El segundo elemento es una lista con su descomposición.
        # Iterarla y recurrir por cada elemento.
        for i, p in enumerate(FP[1]):
            d.update(dicc_filtracion(p,cod+str(i)))

    # Salir de la recursión devolviendo el diccionario creado.
    return d


#%%
# Diccionario de Patriarcas
def dicc_patriarcas(FP,cod=''):
    """Generar un diccionario con los polígonos indivisibles de una filtración.

    Función recursiva que recorre una filtración y agrega a un diccionario
     únicamente los polígonos núcleo indivisibles, iterando sobre la lista de la
     descomposición de los que no.
    """
    # Si el elemento es sólo un polígono (es indivisible), crear un diccionario
    #  con él y salir de esa recursión devolviéndolo.
    if type(FP)==Polygon:
        d = {cod:FP}

    # Si no, el elemento es una tupla de dos elementos, un polígono divisible y
    #  su descomposición.
    else:
        #d = {cod:FP[0]}
        # Iniciar un diccionario vacío
        d = {}
        # Si la lista de descomposición tiene más de un elemento (siempre?):
        if len(FP[1])>1:
            # Iterarla y recurrir por cada elemento
            for i, p in enumerate(FP[1]):
                d.update(dicc_patriarcas(p,cod+str(i)))

    # Salir de la recursión devolviendo el diccionario creado.
    return d


#%%
#Agrupar Filtración
def agrupar_filtracion(F, verb=0):
    """Agrupar una filtración.

    Función recursiva que del árbol de filtración va tomando los últimos
     polígonos núcleos indivisibles (hojas) y uniéndolos a las partes de
     diferencias (respecto de su núcleo divisible de origen) que sólo son
     adyacentes con cada uno de ellos y bajando un nivel del árbol hasta
     analizar los núcleos de descomposición del polígono original.

    Returns:
        Una lista de polígonos con las partes significativas del polígono
         original.
    """
    # Si el elemento es sólo un polígono (es indivisible), crear una lista
    #  con él y salir de esa recursión devolviéndola.
    if type(F)==Polygon:
        return [F]

    # Si no, el elemento es una tupla de dos elementos, un polígono divisible y
    #  su descomposición.
    else:
        # PG es tipo Polygon.
        PG = F[0].buffer(0)  # Polígono grande

        # DescPG es tipo lista, iniciada con la lista de la descomposición de PG
        DescPG = F[1]  # Descomposición del PG
        # Reasignar a DescPG una lista con cada elemento una recursión de la
        #  descomposición de PG
        DescPG = [agrupar_filtracion(P) for P in DescPG]

        # A partir de acá empiezan a salir los elementos de la descomposición
        #  empezando por el último polígono divisible asignado en PG y sus
        #  componentes (indivisibles [F] en la primer salida y divisibles
        #  unidos con sus adyacentes después) listados en DescPG.
        if verb > 1:
            print(f'{PG = }')
            print('DescPG:')
            pprint(DescPG)

        # Lista de polígonos chicos.
        LPC = []
        for LP in DescPG:
            # Extender la lista LPC con todos los núcleos que componen a PG.
            LPC.extend(LP)

        ## Diferencia entre el polígono divisible y todas sus componentes.
        # Crear union de todas las componentes y limpiar la topología.
        unidos = unary_union(LPC).buffer(0)  # pylint: disable=unused-variable
        # Asegurar que PG y unidos tengan los vértices de su intersección.
        PG_n_uni = (PG & unidos).buffer(0)
        PG_dif_uni = (PG - unidos).buffer(0)
        uni_dif_PG = (unidos - PG).buffer(0)

        PG_v = (PG_n_uni | PG_dif_uni).buffer(0)
        unidos_v = (PG_n_uni | uni_dif_PG).buffer(0)

        # Calcular la diferencia
        D = (PG_v - unidos_v).buffer(0)
        # Listar todos los polígonos simples que componen la diferencia
        LD = lpolys(D)

        ## Recalcular LPC de forma que todas las componentes compartan los
        ##  vértices de las diferencias.
        D_unidos = unary_union(D).buffer(0)
        LPC_v = []
        for c in LPC:
            c_limpio = c.buffer(0)
            c_n_duni = (c_limpio & D_unidos).buffer(0)
            c_dif_duni = (c_limpio - D_unidos).buffer(0)
            c_v = (c_n_duni | c_dif_duni).buffer(0)
            LPC_v.append(c_v)

        # Crear una lista para almacenar cuellos (polígonos de la diferencia
        #  que se intersecan con más de una parte componente.
        cuellos = []

        # Iterar sobre los polígonos de la diferencia.
        for i, p in enumerate(LD):  # Las miro una a una

            # Iniciar una lista para componentes adyacentes.
            J=[]
            if verb > 0:
                print(i,end=': ')

            # Agregar a la lista los índices de las componentes adyacentes.
            for j, q in enumerate(LPC_v):  # Me fijo que PChicos tocas
                # p es un polígono de la diferencia, q es una componente.
                if p.intersects(q):
                    if verb > 1:
                        print(i, j)
                    J.append(j)

            # Si la lista de componentes adyacentes tiene más de uno (la
            #  diferencia es un cuello, que es una parte significativa y no
            #  se une a nada).
            if len(J)>1:
                if verb > 0:
                    print(J)
                cuellos.append(p)


            # Si la lista tiene un sólo núcleo (PC) adyacente, unirse
            elif len(J)==1:
                j = J[0]
                q = LPC_v[j]
                LPC_v[j] = (q.union(p)).buffer(0)

            # Si la lista de núcleos adyacentes está vacía, imprimir advertencia.
            else:
                print("PROBLEMA")

        # Extender la lista LPC con los cuellos.
        LPC_v.extend(cuellos)

        # La recursión devuelve la lista de componentes con sus
        #  diferencias adyacentes unidas y los cuellos.
        return LPC_v


#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%

home_dir = Path.home()
wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
fn = wdir / 'laguito' # 'saladito_muy_corto'

gdf = gpd.read_file(str(fn) + '.shp')

# Miro solo el primer polígono.
R = gdf.iloc[0].geometry # union(gdf.iloc[0]['geometry']...)

coords = ((0., 0.), (0., 3.), (3., 3.), (3., 2.), (10.,2.), (10.,4), (15.,4.),(15.,-1.),(10.,-1.), (10.,1.),  (3., 1.), (3., 0.), (0., 0.))
S = Polygon(coords)

#%%

#gdf = gpd.geodataframe.GeoDataFrame(pd.DataFrame(D.items(),columns=['cod','geometry']))

P = R
FP = calcular_filtracion(P)

#%%

# Guardar filtración.
D = dicc_filtracion(FP, cod='')
gdfo = gpd.geodataframe.GeoDataFrame(D.items(), columns=['cod', 'geometry'])
gdfo = gdfo.set_crs(gdf.crs)
gdfo.to_file(str(fn) + '_filt.shp')

#%%

# Guardar patriarcas.
DP = dicc_patriarcas(FP, cod='')
gdfo = gpd.geodataframe.GeoDataFrame(DP.items(), columns=['cod', 'geometry'])
gdfo = gdfo.set_crs(gdf.crs)
gdfo.to_file(str(fn) + '_leaves.shp')

#%%

FA = agrupar_filtracion(FP)
gdfo = gpd.geodataframe.GeoDataFrame(FA, columns=['geometry'])
gdfo = gdfo.set_crs(gdf.crs)
gdfo.to_file(str(wdir / 'agrupada.shp'))

#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
