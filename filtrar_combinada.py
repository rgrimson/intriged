#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Filtrar para detectar radios y descomponer a partir de cuellos."""

# %% Librerías
from pathlib import Path
from pprint import pprint

import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
from shapely import unary_union

from aux import helpers
from aux import exceptions


# %% Calcular Filtración Recursiva
def calcular_filtracion_recursiva(P, r, cod, r_step=1, verb=0, eps=0.001):
    """Calcular la filtración recursiva de P.

    Returns:
        Un Polygon() si P es Polygon() o si se absorbe completamente en
        r + rstep. O una lista `rta` con al menos dos elementos, donde cada
        elemento puede ser un polígono núcleo indivisible de la descomposición
        de P, o una tupla de tres elementos donde el primero es un polígono
        núcleo divisible de la descomposición de P, el segundo es una distancia
        `d` de descomposición y el tercer elemento es la lista `rta` de su
        descomposición, en forma recursiva.
    """
    # P tiene que ser un polígono singlepart, si no eleva un error.
    if type(P)!=Polygon:
        a = 0
        b = a/a
        return b

    ta = helpers.topo(P)

    # Si la topología es vacía (si P es Polygon()) sale de la recursión
    #  devolviendo Polygon().
    if ta==[]:
        return Polygon()

    if verb > 0:
        print(f'\n{cod = }\n{r = }\n{ta = }')

    # Distancia de filtración (no analiza r en el primer paso).
    d = r + r_step

    # Q es la intersección entre P y su buffereado
    #  (prácticamente el mismo buffereado, que debería estar contenido en P
    #  excepto porque tiene vértices que P no tiene).
    buffered = helpers.reinflado(P, d, eps)        
    Q = (P.intersection(buffered)).normalize()
    t = helpers.topo(Q)

    # Mientras que P y Q tengan misma cantidad de partes:
    while len(t)==len(ta):
        d += r_step
        buffered = helpers.reinflado(P,d, eps)        
        Q = (P.intersection(buffered)).normalize()
        t = helpers.topo(Q)

    # En este momento, Q tiene más partes que P, o es un polígono vacío (si P
    #  se absorbió sin haberse descompuesto).
    # LQ es una lista de Polygons que componen a Q o [Polygon()].
    # Estos polígonos son los núcleos componentes de P.
    LQ = helpers.lpolys(Q)

    # Si Q tiene más de una componente Polygon, descomponer P:
    if len(LQ)>1:
        # `rta` es una lista con al menos dos elementos.
        # Cada elemento puede ser Q si es indivisible, o una tupla (Q, d, frQ),
        #  donde frQ es el resultado de filtrar recursivamente a Q en d.
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
            frQ = calcular_filtracion_recursiva(Q, d, cod+str(i), r_step, verb, eps)

            # La filtración recursiva puede devolver un polígono vacío o
            #  una lista `rta`.
            # Si devolvió una lista (Q es divisible a partir de d + r_step):
            if frQ:
                # Agregar un elemento tupla a rta, formada por esta parte
                #  núcleo componente, d y el resultado de filtrarla recursivamente.
                rta.append(((Q, d), frQ))

            # Si frQ devolvió Polygon() (falsy value)
            else:
                # Agregar esta parte Q (componente núcleo indivisible de P)
                #  a `rta`.
                if verb > 1:
                    print('Hoja')
                rta.append((Q, d))

        return rta

    # Si no, LQ tiene una sola parte.
    # Si esa parte no es un polígono vacío (quizás Q tenía partes no poligonales
    #  que se perdieron al convertir a lista?)
    elif LQ[0]:
        if verb > 1:
            print('Se creó una lista de una sola parte?')
            print(f'{LQ[0] = }')
        # Filtrar esa parte en forma recursiva
        return calcular_filtracion_recursiva(LQ[0], d, cod, r_step, verb, eps)

    # Si no, LQ == [Polygon()] y sale de la recursión devolviendo Polygon().
    else:
        return Polygon()


# %% Calcular filtración a partir de un polígono origen
def calcular_filtracion(P, r=0, r_step=1, verb=0, eps = 0.001):
    """Calcular la filtración de un polígono singlepart.

    Returns:
        Tupla de (poligono original, r) y el resultado de llamar a
        `calcular_filtracion_recursiva` (Polygon() si `P` se absorbe
        completamente sin descomponerse, o una lista conteniendo la
        descomposición de P).
    """
    # Limpiar P y calcular su filtración
    P_limpio = P.buffer(0)
    F = calcular_filtracion_recursiva(P_limpio, r, '', r_step, verb, eps)

    # Si calcular_filtracion_recursiva devuelve Polygon(), entonces P es
    #  indivisible para los parámetros de r, r_step.
    # Devolver P:
    if not F:
        return (P, 0)

    # Si no, devolver el P limpio y su filtración.
    else:
        return ((P_limpio, 0), F)


# %% Imprimir resultados de filtración recursiva
def antirecursion(F, cod, verb=0):
    """Descomponer la recursión de `F` con código inicial `cod`."""
    # Lista de polígonos que se van a plotear al final de las impresiones.
    poligonos = []

    def _antirecursion(F, cod, verb=0):
        """Función privada que imprime recursivamente la descomposición de F.

        F puede ser una tupla (P, d) de polígono núcleo indivisible de
        descomposición en d, o una lista de dos elementos o más, cada uno
        correspondiendo a una parte componente.
        """

        f_0 = F[0]
        d = F[1]

        # Crear una lista base de textos a imprimir.
        textos = [f'{cod = }', f'{d = }', f'{type(f_0) = }']

        # Si f_0 es polígono (si no, F es una lista de tuplas):
        if type(f_0) == Polygon:
            # Agregar la topología e imprimir los textos.
            textos.append(f'{helpers.topo(f_0)}')
            if verb > 0:
                print("; ".join(textos), end='.\n')

            # Agregar el polígono a la lista de polígonos a plotear y
            #  salir de la recursión.
            poligonos.append(f_0)

            return None

        # Si no, F es una lista de tuplas, recurrir:
        for i, f in enumerate(F):
            _antirecursion(f, cod + str(i), verb=verb)

    _antirecursion(F, cod, verb=verb)
    helpers.plot_polygon(MultiPolygon(poligonos))

    return None


# %% Extraer distancias
def extraer_distancias(F):
    """Extraer las distancias de una filtración."""
    distancias = []

    def _antirecursion(F):
        """Función privada que extrae las distancias en forma recursiva."""
        f_0 = F[0]
        d = F[1]
        if type(f_0) == Polygon:
            distancias.append(d)
        else:
            for f in F:
                _antirecursion(f)

    _antirecursion(F)

    return distancias

#%%

#%%
#Agrupar Filtración
def agrupar_filtracion(F, verb=0, eps = 0.001):
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
    P = F[0]
    d = F[1]


    if type(P)==Polygon:
        print('o',end='')
        return [P]
    

    # Si no, el elemento es una tupla de dos elementos, un polígono divisible y
    #  su descomposición.
    else:
        P = F[0][0]
        PED = F[1][0]
        if type(PED[0]) == Polygon:
            d = PED[1]
            print('.',end='')
        else:
            d = PED[0][1]
            print('-',end='')
        # DescPG es tipo lista, iniciada con la lista de la descomposición de PG
        DescPG = F[1]  # Descomposición del PG

        # PG es tipo Polygon.
        PG = P.buffer(0)  # Polígono grande

        
        # Reasignar a DescPG una lista con cada elemento una recursión de la
        #  descomposición de PG
        DescPG = [agrupar_filtracion(F) for F in DescPG]

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
            
        #R = helpers.reinflado(PG, d, eps) 

        ## Diferencia entre el polígono divisible y todas sus componentes.
        # Crear union de todas las componentes y limpiar la topología.
        unidos = unary_union(LPC).buffer(0)  # pylint: disable=unused-variable
        
        
        PG_v, unidos_v = helpers.agregar_vertices(PG, unidos)#, R)

        # Calcular la diferencia
        D = (PG_v - unidos_v).buffer(0)
        # Listar todos los polígonos simples que componen la diferencia
        LD = helpers.lpolys(D)

        ## Recalcular LPC de forma que todas las componentes compartan los
        ##  vértices de las diferencias.
        D_unidos = unary_union(D).buffer(0)
        LPC_v = []
        for c in LPC:
            c_limpio = c.buffer(0)
            c_v,  frula = helpers.agregar_vertices(c_limpio, D_unidos)
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




# %% Main
def main():
    """Leer un shapefile, filtrarlo y verificar los radios."""
    home_dir = Path.home()
    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
    fn = wdir / 'tramo' #'laguito' # 'saladito_muy_corto'
    nombre = str(fn) + '.shp'
    nombre_salida = str(fn) + '_desc.shp'

    #R = helpers.gen_poly(tipo='sintetico', nombre='pol_single_hole')
    R = helpers.gen_poly(tipo='fn', nombre=nombre)

    # helpers.plot_polygon(R)

    F = calcular_filtracion(R, verb=3)

    #distancias = extraer_distancias(F)
    #print(f'{distancias = }')

    A = agrupar_filtracion(F, verb=0, eps = 0.001)
    helpers.save_plist(A,nombre_salida)

    return None


if __name__ == '__main__':
    main()
