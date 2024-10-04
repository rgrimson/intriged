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
        print(cod,r,ta)

    # Distancia de filtración (no analiza r en el primer paso).
    d = r + r_step

    # Q es la intersección entre P y su buffereado
    #  (prácticamente el mismo buffereado, que debería estar contenido en P
    #  excepto porque tiene vértices que P no tiene).
    Q = (P.intersection(P.buffer(-d, quad_segs=16).buffer(d+eps , quad_segs=16))).normalize()
    t = helpers.topo(Q)

    # Mientras que P y Q tengan misma cantidad de partes:
    while len(t)==len(ta):
        d += r_step
        buffered = P.buffer(-d, quad_segs=16).buffer(d+eps, quad_segs=16)
        if not buffered.is_valid:
            textos = [
                'Buffered es invalido.',
                f'{P = }',
                f'{d = }',
                f'{buffered = }',
            ]
            msg = '\n'.join(textos)
            raise exceptions.InvalidGeometryError(msg)

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
                    print(f'{frQ = }')
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
def antirecursion(F, cod):
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
            _antirecursion(f, cod + str(i))

    _antirecursion(F, cod)
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


# %% Main
def main():
    """Leer un shapefile, filtrarlo y verificar los radios."""
    home_dir = Path.home()
    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
    fn = wdir / 'laguito' # 'saladito_muy_corto'
    nombre = str(fn) + '.shp'

    # R = helpers.gen_poly(tipo='sintetico', nombre='pol_single_hole')
    R = helpers.gen_poly(tipo='fn', nombre=nombre)

    # helpers.plot_polygon(R)

    F = calcular_filtracion(R, verb=0)

    distancias = extraer_distancias(F)
    print(f'{distancias = }')

    return None


if __name__ == '__main__':
    main()
