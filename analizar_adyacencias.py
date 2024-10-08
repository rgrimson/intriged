#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analizar adyacencias en una filtración."""

# %% Librerías
from pathlib import Path
from pprint import pprint
from rtree import index

import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
from shapely import unary_union

from aux import helpers
from aux import exceptions


# %% Calcular Filtración Recursiva
def calcular_filtracion_recursiva(P, cod='0', r=0, r_step=1, verb=0, eps=0.001):
    """Calcular la filtración recursiva de P.

    Returns:
        Un diccionario {'cod': cod, 'd': d, 'P': P, 'F': F}, donde F puede ser
        None si P es indivisible y se absorbe en d, o una lista de dos o más
        elementos, donde cada elemento es un diccionario correspondiente a cada
        parte de la descomposición de P en d.
    """
    # Empezar limpiando los vértices de P.
    P = P.buffer(0)

    # P tiene que ser un polígono singlepart, si no eleva un error.
    if not isinstance(P, Polygon):
        textos = ['P no es un poligono singlepart.',
                  f'{cod = }',
                  f'{r = }',
                  f'{P = }']
        msg = '\n'.join(textos)
        raise exceptions.NotAPolygonError(msg)

    # P tiene que ser un polígono válido, si no eleva un error.
    if not P.is_valid:
        textos = ['P es invalido.',
                  f'{cod = }',
                  f'{r = }',
                  f'{P = }']
        msg = '\n'.join(textos)
        raise exceptions.InvalidGeometryError(msg)

    # P no puede ser un polígono vacío, si no eleva un error.
    if P.is_empty:
        textos = ['P es vacío.',
                  f'{cod = }',
                  f'{r = }',
                  f'{P = }']
        msg = '\n'.join(textos)
        raise exceptions.EmptyGeometryError(msg)

    # Topología de P.
    topo_P = helpers.topo(P)

    if verb > 0:
        print(f'\n{cod = }\n{r = }\n{topo_P = }')

    # Distancia de filtración.
    d = r + r_step

    # Q es la intersección entre P y su buffereado
    #  (prácticamente el mismo buffereado, que debería estar contenido en P
    #  excepto porque tiene vértices que P no tiene).
    buffered = helpers.reinflado(P, d, eps)
    Q = (P.intersection(buffered)).normalize()
    topo_Q = helpers.topo(Q)

    # Mientras que P y Q tengan misma cantidad de partes:
    while len(topo_P) == len(topo_Q):
        d += r_step
        buffered = helpers.reinflado(P, d, eps)
        Q = (P.intersection(buffered)).normalize()
        topo_Q = helpers.topo(Q)

    # En este momento, Q tiene más partes que P, o es un polígono vacío (si P
    #  se absorbió sin haberse descompuesto); y d es la distancia para la que
    #  cambió la topología de P.

    # Iniciar el diccionario de respuesta.
    rta = {
        'cod': cod,
        'P': P,
        'd': d,
        'F': None
    }

    # LQ es una lista de Polygons partes de Q o [Polygon()].
    # Estos polígonos son los núcleos componentes de P.
    LQ = helpers.lpolys(Q)

    # Si Q tiene más de una componente Polygon, descomponer P:
    if len(LQ) > 1:
        # rta['F'] pasa a ser una lista con al menos dos elementos.
        # Cada elemento es el diccionario de respuesta de la filtración de Q.
        rta['F'] = []
        # Por cada parte de Q:
        for i, Q in enumerate(LQ):
            # Calcular esta misma filtración para cada núcleo, llevando un
            #  código que concatena índices de cada parte en forma recursiva, y
            #  una distancia inicial que es la distancia en la que se generó
            #  esta parte.
            frQ = calcular_filtracion_recursiva(Q, cod+str(i), d, r_step, verb,
                                                 eps)
            rta['F'].append(frQ)

        return rta

    # Si no, LQ tiene una sola parte.
    # Si esa parte no es un polígono vacío (falsy value) quizás Q tenía partes
    #  no poligonales que se perdieron al convertir a lista? Analizarlo.
    if LQ[0]:
        textos = ['Filtración de un solo elemento no vacío.',
                  f'{cod = }',
                  f'{d = }',
                  f'{LQ = }']
        msg = '\n'.join(textos)
        raise exceptions.FiltrationError(msg)

    # Si no, LQ == [Polygon()] y sale de la recursión devolviendo el
    #  diccionario de respuesta construído (rta['F'] == None).
    return rta


# %% Imprimir resultados de filtración recursiva
def antirecursion(F, verb=0):
    """Descomponer la recursión de una filtración."""
    # Lista de polígonos que se van a plotear al final de las impresiones.
    poligonos = []

    def _antirecursion(F, verb=0):
        """Función privada que imprime recursivamente la descomposición de F.

        F es un diccionario {'cod': cod, 'd': d, 'P': P, 'F': F}.
        """
        # Agregar F['P'] a la lista de polígonos.
        poligonos.append(F['P'])

        if verb > 0:
            textos = [f'{F["cod"] = }',
                    f'{helpers.topo(F["P"]) = }',
                    f'{F["d"] = }']
            print('\n'.join(textos), end='.\n')

        # Si F['F'] no es None, entonces es una nueva filtración.
        if F['F']:
            for f in F['F']:
                _antirecursion(f, verb)

    _antirecursion(F, verb)
    print()
    helpers.plot_polygon(MultiPolygon(poligonos))


# %% Extraer distancias
def extraer_distancias(F):
    """Extraer las distancias de una filtración."""
    distancias = []

    def _antirecursion(F):
        """Función privada que extrae las distancias en forma recursiva."""
        distancias.append(F['d'])
        if F['F']:
            for f in F['F']:
                _antirecursion(f)

    _antirecursion(F)

    return distancias


# %% Crear lista de diccionarios
def crear_lista_de_diccionarios(F):
    """Crear una lista de diccionarios con cada elemento de la filtración F."""
    diccionarios = []

    def _antirecursion(F):
        """Función privada que extrae las geometrías en forma recursiva."""
        dicc = {'cod': F['cod'], 'd': F['d'], 'geometry': F['P']}
        diccionarios.append(dicc)
        if F['F']:
            for f in F['F']:
                _antirecursion(f)

    _antirecursion(F)

    return diccionarios


# %% Crear lista de diccionarios de hojas
def crear_lista_de_hojas(F):
    """Crear una lista de diccionarios con cada hoja de la filtración F."""
    diccionarios = []

    def _antirecursion(F):
        """Función privada que extrae únicamente las hojas."""
        if F['F']:
            for f in F['F']:
                _antirecursion(f)
        else:
            dicc = {'cod': F['cod'], 'd': F['d'], 'geometry': F['P']}
            diccionarios.append(dicc)

    _antirecursion(F)

    return diccionarios


# %% Obtener cuellos
def obtener_cuellos(P, hojas, eps=0.001, verb=0):
    """Obtener una lista de cuellos de un polígono a partir de sus hojas."""
    # inter son las hojas, diff son cuellos y medialunas.
    inter, diff = helpers.get_inter_diff(P, hojas, eps)
    if verb > 0:
        print(f'{len(inter) = }')
        print(f'{len(diff) = }')

    # Inicializar una lista para almacenar cuellos.
    cuellos = []

    # Crear un índice R-tree
    idx = index.Index()

    # Completar el indice usando el bbox de cada hoja
    if verb > 0:
        print("Armando indice...")
    for i, hoja in enumerate(inter):
        idx.insert(i, hoja.bounds)

    # Analizar cada diferencia para detectar si es cuello o medialuna
    if verb > 0:
        print("Calculando intersecciones...")

    # D es una geometría de la diferencia, que puede ser un cuello o una
    #  medialuna, dependiendo de la cantidad de hojas con las que interseca.
    for i, D in enumerate(diff):
        # TODO: Decidir por ínidice de Miller y por ratio interseccion/area.
        # n lleva la cuenta de intersecciones (si n > 1, D es cuello).
        n = 0
        # j es el índice de las hojas cuyo bbox interseca el de D.
        for j in idx.intersection(D.bounds):  # pylint: disable=not-an-iterable
            if D.intersects(inter[j]):
                n += 1

        # Si n == 0, hubo un problema con la filtración. Analizarlo.
        if n == 0:
            textos = ['Parte de la diferencia no intersecta ninguna hoja.',
                      f'{i = }',
                      f'{D = }']
            msg = '\n'.join(textos)
            raise exceptions.FiltrationError(msg)

        # Si no, si n == 1, D es medialuna, si no es cuello.
        if n > 1:
            # TODO: Analizar el d de las hojas que intersecan
            cuellos.append(D)

    return cuellos


# %% Main
def main():
    """Leer un shapefile, filtrarlo y verificar los radios."""
    home_dir = Path.home()
    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
    fn = wdir / 'laguito' #'saladito_muy_corto'  # 'laguito' # 'saladito_muy_corto'
    nombre = str(fn) + '.shp'

    # R = helpers.gen_poly(tipo='sintetico', nombre='pol_single_hole')
    R = helpers.gen_poly(tipo='fn', nombre=nombre)

    # helpers.plot_polygon(R)

    print('Calculando filtración recursiva...')
    F = calcular_filtracion_recursiva(R, r_step=1)

    print('Guardando shapefile de filtración...')
    D = crear_lista_de_diccionarios(F)
    nombre_filtr = str(fn) + '_filtr.shp'
    helpers.shapefile_from_data(D, crs='EPSG:32721', fn=nombre_filtr)

    # antirecursion(F, verb=0)

    distancias = extraer_distancias(F)
    print(f'{distancias = }')

    print('Creando lista de hojas...')
    L = crear_lista_de_hojas(F)
    # pprint(L)

    print('Guardando shapefile de hojas...')
    nombre_hojas = str(fn) + '_hojas.shp'
    helpers.shapefile_from_data(L, crs='EPSG:32721', fn=nombre_hojas)
    # A = agrupar_filtracion(F, verb=0, eps = 0.001)
    # helpers.save_plist(A,nombre_salida)

    hojas = [hoja['geometry'] for hoja in L]
    # pprint(hojas)
    # helpers.plot_polygon(MultiPolygon(hojas))
    print('Obteniendo cuellos...')
    cuellos = obtener_cuellos(R, hojas)

    print('Guardando shapefile de cuellos...')
    nombre_cuellos = str(fn) + '_cuellos.shp'
    helpers.shapefile_from_geom(cuellos, crs='EPSG:32721', fn=nombre_cuellos)
    # helpers.plot_polygon(MultiPolygon(cuellos))

    # Una vez que se obtienen los cuellos, las partes significativas son
    #  la intersección y la diferencia entre R y cuellos.
    print('Obteniendo partes significativas...')
    inter, diff = helpers.get_inter_diff(R, cuellos, 0.001)
    partes = inter + diff

    print('Guardando shapefile de descomposición...')
    nombre_desc = str(fn) + '_desc.shp'
    helpers.shapefile_from_geom(partes, crs='EPSG:32721', fn=nombre_desc)
    # helpers.plot_polygon(MultiPolygon(partes))


if __name__ == '__main__':
    main()
