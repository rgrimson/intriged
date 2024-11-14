#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analizar adyacencias en una filtración."""

# %% Librerías
from pathlib import Path
from pprint import pprint
from rtree import index

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np


from sklearn.mixture import GaussianMixture

from shapely.geometry.linestring import LineString
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
from shapely.ops import split
from shapely import (
    unary_union,
    line_merge,
)

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


# %% Etiquetar distancias
def etiquetar_distancias(distancias, k=3):
    """Etiquetar las distancias en k categorías."""
    A = np.array(distancias)
    GMM = GaussianMixture(n_components=k, covariance_type='full', max_iter=20, random_state=0)
    GMM.fit(A.reshape(-1,1))
    E=GMM.predict(A.reshape(-1,1))
    P = np.vstack([A,E]).T
    print(P)
    return 0





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


# %% Obtener diferencias
def obtener_diferencias(P, hojas, eps=0.001, verb=0):
    """Obtener una lista de diferencias de un polígono a partir de sus hojas.

    hojas es una lista de diccionarios de hojas, extraído con
    crear_lista_de_hojas.
    """
    # inter son las hojas, diff son cuellos y medialunas.
    inter, diff = helpers.get_inter_diff(P, hojas, eps)
    if verb > 0:
        print(f'{len(inter) = }')
        print(f'{len(diff) = }')

    # Crear un índice R-tree
    idx = index.Index()

    # Completar el indice usando el bbox de cada hoja
    if verb > 0:
        print("Armando indice...")
    for i, hoja in enumerate(inter):
        idx.insert(i, hoja['geometry'].bounds)

    # Analizar cada diferencia para detectar si es cuello o medialuna
    if verb > 0:
        print("Calculando intersecciones...")

    # D es un diccionario {'geometry': diff} con la geometría de la diferencia,
    #  que puede ser un cuello o una medialuna.
    for i, D in enumerate(diff):
        # len_inter lleva el largo de las intersecciones.
        len_inter = 0

        # j es el índice de las hojas cuyo bbox interseca el de D.
        for j in idx.intersection(D['geometry'].bounds):  # pylint: disable=not-an-iterable
            if D['geometry'].intersects(inter[j]['geometry']):

                # Sumar una adyacencia.
                D['n'] += 1
                # Agregar el código de la hoja a la lista de adyacencias de D.
                D['cods'].append(inter[j]['cod'])
                # Agregar la distancia de filtración de la hoja a la lista de
                #  distancias.
                D['dists'].append(inter[j]['d'])

                # Calcular la intersección y agregarla a la lista de líneas.
                line = inter[j]['geometry'].intersection(D['geometry'])
                D['lines'].append(line)
                # Calcular el largo de la intersección y sumarlo al total.
                len_inter += helpers.get_length(line)

        # Si n == 0, hubo un problema con la filtración. Analizarlo.
        if D['n'] == 0:
            textos = ['Parte de la diferencia no intersecta ninguna hoja.',
                      f'{i = }',
                      f'{D = }']
            msg = '\n'.join(textos)
            raise exceptions.FiltrationError(msg)

        # Calcular el area y el ratio len_inter^2 / area.
        D['ratio'] = (len_inter * len_inter) / D['geometry'].area

        # Agregar el índice de Miller.
        D['miller'] = helpers.get_miller(D['geometry'])

    return diff


# %% Etiquetar cuellos
def etiquetar_cuellos(diff, umbrales, max_ratio=1.5, max_miller=0.5):
    """Etiquetar cuellos en la lista de diferencias."""

    for d in diff:
        # Analizar por separado según cantidad de adyacencias.
        if d['n'] == 1:
            # Solo analizar si tiene ratio y miller chicos.
            if d['ratio'] <= max_ratio and d['miller'] <= max_miller:
                d['es_cuello'] = True

        # Si no, analizar que no todas las distancias sean menores ni todas
        #  mayores que cada umbral.
        else:
            umbrales = np.array(sorted(set(umbrales)))
            distancias = np.array(d['dists'])

            # Comprobar si todas las distancias se encuentran dentro de un rango
            for umbral in umbrales:
                if not (all(distancias < umbral) or
                        all(distancias >= umbral)):
                    d['es_cuello'] = True
                    break

    return diff


# %% Extraer lineas
def extraer_lineas(cuellos, hojas, umbrales, eps=0.001, min_length=1.0):
    """Extraer las lineas de intersección entre cuellos y hojas."""
    # Quitar de cada cuello la línea, código de hoja y distancia,
    #  de la mayor distancia (el cuello se unirá a esa hoja)
    cuellos_filtrados = []
    for cuello in cuellos:
        # Si hay una sola adyacencia, mantener esa adyacencia.
        if cuello['n'] == 1:
            cuellos_filtrados.append(cuello)
        else:
            # Encontrar el índice del elemento de mayor distancia
            max_index = cuello['dists'].index(max(cuello['dists']))
            # Encontrar el máximo umbral por debajo de la mayor distancia
            umbral_max = max(umb
                             for umb in umbrales
                             if umb < cuello['dists'][max_index])

            # Mantener los elementos de distancia menor al máximo umbral
            indices_mantener = [i
                                for i, dist in enumerate(cuello['dists'])
                                if dist < umbral_max]
            cuello['cods'] = [cuello['cods'][i] for i in indices_mantener]
            cuello['dists'] = [cuello['dists'][i] for i in indices_mantener]
            cuello['lines'] = [cuello['lines'][i] for i in indices_mantener]

            cuellos_filtrados.append(cuello)

    # Extraer primero cada línea de cada cuello, con los atributos de
    #  las adyacencias.
    lineas = [{'geometry': line_merge(line),
               'cuello_id': cuello['id'],
               'cod': cuello['cods'][i],
               'dist': cuello['dists'][i]}
               for cuello in cuellos_filtrados
               for i, line in enumerate(cuello['lines'])]

    # Extender las líneas un epsilon en cada extremo y agregarles un id.
    lineas = [{'geometry': helpers.extender_linea(line['geometry'], eps),
               'id': i,
               'cuello_id': line['cuello_id'],
               'cod': line['cod'],
               'dist': line['dist']}
               for i, line in enumerate(lineas)
               if line['geometry'].length >= min_length]

    return lineas


# %% Rectificar líneas
def rectificar_lineas(P, lineas, eps=0.001):
    """Rectificar las líneas en sus intersecciones con P."""
    # Obtener los puntos de intersección de cada línea con los anillos de P.
    anillos = unary_union(helpers.extraer_anillos(P))
    intersecciones = []
    for linea in lineas:
        puntos = linea['geometry'].intersection(anillos)
        inicio, fin = puntos.geoms[0], puntos.geoms[-1]
        rectificada = LineString([inicio, fin])
        extendida = helpers.extender_linea(rectificada, eps)
        intersecciones.append({'id': linea['id'],
                               'cuello_id': linea['cuello_id'],
                               'cod': linea['cod'],
                               'dist': linea['dist'],
                               'geometry': extendida})

    return intersecciones


# %% Dividir polígono
def dividir_poligono(P, lineas):
    """Dividir un polígono por una lista de líneas."""
    # Iniciar lista de subpolígonos con el polígono original.
    subpoligonos = [{'geometry': P}]

    for linea in lineas:
        # Iniciar lista de supolígonos generados por esta línea.
        nuevos_subpol = []
        for subpol in subpoligonos:
            if subpol['geometry'].intersects(linea['geometry']):
                division = split(subpol['geometry'], linea['geometry'])
                nuevos_subpol.extend(list(division.geoms))
            else:
                nuevos_subpol.append(subpol['geometry'])
        # Una vez recorridos todos los polígonos, nuevos_subpol tiene la
        #  división de los subpoligonos anteriores por la linea actual.
        # Reasignar subpoligonos.
        subpoligonos = [{'geometry': poligono} for poligono in nuevos_subpol]

    return subpoligonos


# %% Obtener partes significativas
def obtener_partes_significativas(P, cuellos, eps=0.001, quad_segs=16):
    """Obtener las partes significativas de P.

    `cuellos` es una lista de geometrías de cuellos.
    """
    unidos = unary_union(cuellos)
    buffered = unidos.buffer(eps, quad_segs=quad_segs)

    if not buffered.is_valid:
        raise exceptions.InvalidGeometryError

    intersecciones = helpers.lpolys(P.intersection(buffered))
    diferencias = helpers.lpolys(P.difference(buffered))

    partes = intersecciones + diferencias

    return partes


# %% Main
def main():
    """Leer un shapefile, filtrarlo y verificar los radios."""
    home_dir = Path.home()
    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
    fn = wdir / 'laguito' #'saladito_muy_corto'  # 'laguito' # 'saladito_muy_corto'
    nombre = str(fn) + '.shp'
    print(f"{nombre = }")

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




    print('Creando lista de hojas...')
    L = crear_lista_de_hojas(F)
    pprint(L)
    distancias = extraer_distancias(L)
    print(f'{distancias = }')
    print("Etiquetando distancias...")
    etiquetar_distancias(distancias)

    print('Guardando shapefile de hojas...')
    nombre_hojas = str(fn) + '_hojas.shp'
    helpers.shapefile_from_data(L, crs='EPSG:32721', fn=nombre_hojas)

    # hojas = [hoja['geometry'] for hoja in L]
    # pprint(hojas)
    # helpers.plot_polygon(MultiPolygon(hojas))

    print('Obteniendo diferencias...')
    diferencias = obtener_diferencias(R, L)
    etiquetadas = etiquetar_cuellos(diferencias, [80, 800])

    print('Guardando shapefile de diferencias etiquetadas...')
    nombre_difs = str(fn) + '_difs.shp'
    helpers.shapefile_from_data(etiquetadas, crs='EPSG:32721', fn=nombre_difs)

    # Extraer los cuellos con todas sus etiquetas
    print("Extrayendo cuellos...")
    cuellos = [e for e in etiquetadas if e['es_cuello']]

    # pprint(cuellos)
    # cuellos_geoms = [c['geometry'] for c in cuellos]
    # helpers.plot_polygon(MultiPolygon(cuellos_geoms))
    print('Guardando shapefile de cuellos...')
    nombre_cuellos = str(fn) + '_cuellos.shp'
    helpers.shapefile_from_data(etiquetadas, crs='EPSG:32721', fn=nombre_cuellos)

    # Extraer las lineas que son intersección en cada cuello
    print('Extrayendo lineas...')
    lineas = extraer_lineas(cuellos, hojas=L, umbrales=[80, 800])

    print('Guardando shapefile de linestrings...')
    nombre_lineas = str(fn) + '_lineas.shp'
    helpers.shapefile_from_data(lineas, crs='EPSG:32721', fn=nombre_lineas)

    # Rectificar las líneas
    print('Rectificando lineas...')
    rectificadas = rectificar_lineas(R, lineas)

    print('Guardando shapefile de rectificadas...')
    nombre_rectif = str(fn) + '_rectificadas.shp'
    helpers.shapefile_from_data(rectificadas, crs='EPSG:32721', fn=nombre_rectif)

    # Dividir el polígono por las líneas rectificadas
    print('Dividiendo poligonos...')
    divididos = dividir_poligono(R, rectificadas)

    print('Guardando shapefile de divididos...')
    nombre_divididos = str(fn) + '_divididos.shp'
    helpers.shapefile_from_data(divididos, crs='EPSG:32721', fn=nombre_divididos)


    # Una vez que se obtienen los cuellos, las partes significativas son
    #  la intersección y la diferencia entre R y cuellos.
    # print('Obteniendo partes significativas...')
    # partes = obtener_partes_significativas(R, cuellos)

    # print('Guardando shapefile de descomposición...')
    # nombre_desc = str(fn) + '_desc.shp'
    # helpers.shapefile_from_geom(partes, crs='EPSG:32721', fn=nombre_desc)
    # helpers.plot_polygon(MultiPolygon(partes))


if __name__ == '__main__':
    main()
