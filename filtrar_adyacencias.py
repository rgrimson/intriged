#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analizar adyacencias en una filtración."""

# %% Librerías
from pathlib import Path

from rtree import index
import numpy as np
from shapely.geometry.linestring import LineString
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.ops import polygonize
from shapely import (
    unary_union,
    line_merge,
    union_all,
)

from aux import helpers
from aux import exceptions


# %% Calcular Filtración Recursiva
def calcular_filtracion_recursiva(P, cod='0', r=0, r_step=1, verb=0, eps=0.001):
    """Calcular la filtración recursiva de P.

    Returns: Un diccionario {'cod': cod, 'd': d, 'P': P, 'F': F}, donde F puede
    ser None si P es indivisible y se absorbe en d, o una lista de dos o más
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

# %% Calcular Filtración Recursiva a partir de lista
def calcular_filtracion_recursiva_a_partir_de_lista(P, cod='0', radios=None, verb=0, eps=0.001):
    """Calcular la filtración recursiva de P.

    Returns: Un diccionario {'cod': cod, 'd': d, 'P': P, 'F': F}, donde F puede
    ser None si P es indivisible y se absorbe en d, o una lista de dos o más
    elementos, donde cada elemento es un diccionario correspondiente a cada
    parte de la descomposición de P en d.
    """

    if not isinstance(radios, list) or radios == []:
        raise exceptions.FiltrationError('Radios vacios.')

    #defino r y r_step
    r = radios[0]
    radios_pendientes = radios.copy()

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
    d = radios_pendientes[0]

    # Q es la intersección entre P y su buffereado
    #  (prácticamente el mismo buffereado, que debería estar contenido en P
    #  excepto porque tiene vértices que P no tiene).
    buffered = helpers.reinflado(P, d, eps)
    Q = (P.intersection(buffered)).normalize()
    topo_Q = helpers.topo(Q)
    radios_pendientes.pop(0)

    # Mientras que P y Q tengan misma cantidad de partes:
    while len(topo_P) == len(topo_Q) and radios_pendientes:
        d = radios_pendientes[0]
        radios_pendientes.pop(0)
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
            frQ = calcular_filtracion_recursiva_a_partir_de_lista(Q, cod+str(i), radios_pendientes, verb,
                                                 eps)
            rta['F'].append(frQ)

        return rta


    # Si no, LQ == [Polygon()] o se acabaron los radios pendientes y LQ == [P] y sale de la recursión devolviendo el
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

        # Calcular el ratio len_inter / perimetro.
        D['ratio'] = len_inter / D['geometry'].exterior.length

        # Agregar el índice de Miller.
        D['miller'] = helpers.get_miller(D['geometry'])

    return diff


# %% Etiquetar cuellos
def etiquetar_cuellos(diff, umbrales, max_ratio, max_miller):
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

            # Comprobar si todas las distancias se encuentran dentro de un
            #  mismo rango, en cuyo caso no es cuello.
            for umbral in umbrales:
                if not (all(distancias < umbral) or
                        all(distancias >= umbral)):
                    d['es_cuello'] = True
                    break

    return diff


# %% Extraer lineas de cuellos
def extraer_lineas_de_cuellos(cuellos, umbrales, eps=0.001, min_length=1.0):
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
            #  para mantener las líneas de adyacencias con distancia menor
            #  al umbral máximo.
            try:
                umbral_max = max(umb
                             for umb in umbrales
                             # Menor o igual porque puede haber un umbral
                             #  justo en la distancia máxima.
                             if umb <= cuello['dists'][max_index])
            except Exception as e:
                print(f'{cuello = }')
                raise e

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
    lineas = [{'geometry': helpers.extender_linea(line, eps),
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
        extendida = helpers.extender_linea({'geometry': rectificada}, eps)
        intersecciones.append({'id': linea['id'],
                               'cuello_id': linea['cuello_id'],
                               'cod': linea['cod'],
                               'dist': linea['dist'],
                               'geometry': extendida})

    return intersecciones


# %% Dividir polígono
def dividir_poligono(P, lineas):
    """Dividir un polígono por una lista de lineas."""

    # Implementación modificada de split.
    # Extraer las geometrías de las líneas.
    geoms = [linea['geometry'] for linea in lineas]
    # Agregarle los anillos del polígono.
    anillos = helpers.extraer_anillos(P)
    geoms.extend(anillos)
    # Unir en una multilinestring
    unioned = union_all(geoms)
    # Poligonizar y quedarse con los polígonos cuyo punto representativo
    #  esté dentro de P.
    subpoligonos = [{'geometry': poly}
                    for poly in polygonize(unioned)
                    if poly.representative_point().intersects(P)]
    return subpoligonos


# %% Etiquetar divididos
def etiquetar_divididos(subpolis, radios, max_ratio, eps=0.001):
    """Etiquetar los subpolígonos divididos.

    Returns: Una lista de diccionarios [{`id`: id, `d`: d, `geometry`: P}, ...]
    donde P es un subpolígono de la división y se absorbe en d.
    """
    # Iterar sobre la lista de supolígonos divididos.
    for i, subpol in enumerate(subpolis):
        # Agregar id
        subpol['id'] = i
        # Definir P
        P = subpol['geometry']
        # P tiene que ser un polígono singlepart, si no eleva un error.
        if not isinstance(P, Polygon):
            textos = ['P no es un poligono singlepart.']
            msg = '\n'.join(textos)
            raise exceptions.NotAPolygonError(msg)

        # P tiene que ser un polígono válido, si no eleva un error.
        if not P.is_valid:
            textos = ['P es invalido.']
            msg = '\n'.join(textos)
            raise exceptions.InvalidGeometryError(msg)

        # P no puede ser un polígono vacío, si no eleva un error.
        if P.is_empty:
            textos = ['P es vacío.']
            msg = '\n'.join(textos)
            raise exceptions.EmptyGeometryError(msg)

        # Distancia de filtración.
        radios_pendientes = radios.copy()
        d = radios_pendientes[0]
        radios_pendientes.pop(0)

        # Q es el buffer negativo de P.
        Q = helpers.buffer_negativo(P, d, eps)

        # Mientras que Q no sea un polígono vacío:
        while radios_pendientes and not Q.is_empty:
            d = radios_pendientes[0]
            radios_pendientes.pop(0)
            Q = helpers.buffer_negativo(P, d, eps)

        # Agregar d a subpol.
        subpol['d'] = d

    # Analizar adyacencias.
    idx = index.Index()
    for subpol in subpolis:
        subpol['n'] = 0
        subpol['ids_adyac'] = []
        subpol['dists_ady'] = []
        subpol['lines'] = []
        idx.insert(subpol['id'], subpol['geometry'].bounds)

    # Recorrer cada subpolígono S y analizar cuáles son sus adyacentes
    for S in subpolis:
        # len_inter lleva el largo de las intersecciones.
        len_inter = 0

        # j es el índice de los supolígonos que intersecan a este mismo.
        for j in idx.intersection(S['geometry'].bounds):  # pylint: disable=not-an-iterable
            if (j != S['id'] and
                S['geometry'].intersects(subpolis[j]['geometry'])):

                # Sumar una adyacencia.
                S['n'] += 1
                # Agregar el id del S adyacente.
                S['ids_adyac'].append(subpolis[j]['id'])
                # Agregar la distancia de filtración de la hoja a la lista de
                #  distancias.
                S['dists_ady'].append(subpolis[j]['d'])

                # Calcular el largo de la intersección y sumarlo al total.
                line = subpolis[j]['geometry'].intersection(S['geometry'])
                S['lines'].append(line)
                len_inter += helpers.get_length(line)

        # Si n == 0, hubo un problema con la filtración. Analizarlo.
        if S['n'] == 0:
            textos = ['No hay adyacencias en subpoligono dividido',
                      f'id = {S["id"]}.']
            msg = '\n'.join(textos)
            raise exceptions.FiltrationError(msg)


        # Calcular el ratio len_inter / perimetro.
        S['ratio'] = len_inter / S['geometry'].exterior.length

        # Agregar el índice de Miller.
        S['miller'] = helpers.get_miller(S['geometry'])

    # Decidir si un subpolígono debe ser unido a otro.
    # Para que un polígono se una a otro debe tener un ratio alto.
    for S in subpolis:
        if S['ratio'] > max_ratio:
            S['unir'] = True
        else:
            S['unir'] = False

    return subpolis


# %% Unir subpolígonos divididos
def unir_divididos(divididos):
    """Unir subpolígonos divididos."""
    # Crear dos listas, una de subpolígonos que quedan y otra de los que
    #  se unen.
    subpolis_quedan = []
    subpolis_unen = []
    for S in divididos:
        if S['unir']:
            subpolis_unen.append(S)
        else:
            # Redefinir el polígono que va a quedar.
            queda = {}
            queda['id'] = S['id']
            queda['geometry'] = S['geometry']
            queda['d'] = S['d']
            queda['unidos'] = []
            subpolis_quedan.append(queda)

    # Unir subpolígonos.
    for se_une in subpolis_unen:
        # Unirse al adyacente de mayor distancia.
        max_index = se_une['dists_ady'].index(max(se_une['dists_ady']))
        id_adyac = se_une['ids_adyac'][max_index]

        # Para redefinir el polígono que va a quedar, extraigo su índice de la
        #  lista de subpolígonos que quedan.
        i_lista = [i for i, s
                    in enumerate(subpolis_quedan)
                    if s['id'] == id_adyac]

        if len(i_lista) != 1:
            texts = ['No se encontro el subpoligono a unir.',
                     f'id_adyac = {id}.']
            msg = '\n'.join(texts)
            raise exceptions.FiltrationError(msg)
        else:
            i = i_lista[0]

        # Unir a la geometría del polígono que queda y anotar la union.
        subpolis_quedan[i]['geometry'] = (subpolis_quedan[i]['geometry']
                                            .union(se_une['geometry'])
        )
        subpolis_quedan[i]['unidos'].append(se_une['id'])

    return subpolis_quedan


# %% Main
def main(nombre, radios, umbrales, max_ratio, max_miller):
    """Leer un shapefile, filtrarlo y verificar los radios."""
    home_dir = Path.home()
    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'
    fn = wdir / nombre
    nombre_shp = str(fn) + '.shp'
    print(f"{nombre = }")

    if nombre == 'pol_single_hole' or nombre == 'pol_rafa':
        R = helpers.gen_poly(tipo='sintetico', nombre='pol_single_hole')
    else:
        R = helpers.gen_poly(tipo='fn', nombre=nombre_shp)

    # helpers.plot_polygon(R)

    print('Calculando filtración recursiva...')
    # Filtración usando r y r-step
    # F = calcular_filtracion_recursiva(R, r_step=1)

    # Filtración usando lista de umbrales
    F = calcular_filtracion_recursiva_a_partir_de_lista(R, radios=radios)

    print('Guardando shapefile de filtración...')
    D = crear_lista_de_diccionarios(F)
    nombre_filtr = str(fn) + '_filtr.shp'
    helpers.shapefile_from_data(D, crs='EPSG:32721', fn=nombre_filtr)

    print('Creando lista de hojas...')
    L = crear_lista_de_hojas(F)
    # pprint(L)

    print('Guardando shapefile de hojas...')
    nombre_hojas = str(fn) + '_hojas.shp'
    helpers.shapefile_from_data(L, crs='EPSG:32721', fn=nombre_hojas)

    # hojas = [hoja['geometry'] for hoja in L]
    # pprint(hojas)
    # helpers.plot_polygon(MultiPolygon(hojas))

    print('Obteniendo diferencias...')
    diferencias = obtener_diferencias(R, L)
    etiquetadas = etiquetar_cuellos(diferencias, umbrales, max_ratio, max_miller)

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
    lineas = extraer_lineas_de_cuellos(cuellos, umbrales=umbrales)

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
    print('Etiquetando poligonos divididos...')
    divididos_etiquetados = etiquetar_divididos(divididos,
                                                radios=radios,
                                                max_ratio=max_ratio)

    print('Guardando shapefile de divididos...')
    nombre_divididos = str(fn) + '_divididos.shp'
    helpers.shapefile_from_data(divididos_etiquetados,
                                crs='EPSG:32721',
                                fn=nombre_divididos)

    # Unir divididos
    print('Uniendo subpolígonos divididos...')
    unidos = unir_divididos(divididos_etiquetados)

    print('Guardando shapefile de unidos...')
    nombre_unidos = str(fn) + '_unidos.shp'
    helpers.shapefile_from_data(unidos,
                                crs='EPSG:32721',
                                fn=nombre_unidos)

    return None

# %% Constantes
if __name__ == '__main__':
    # Nombre del polígono a filtrar, puede ser:
    # 'pol_single_hole', 'pol_rafa', 'saladito_muy_corto', 'laguito', 'saladito_muy_corto'
    NOMBRE = 'laguito'
    # Radios de filtración
    RADIOS = list(range(20, 1000, 20))
    # Umbrales de agrupación
    UMBRALES = [100, 200, 400, 600]
    # Umbrales de decisión
    MAX_RATIO = 0.3  # Ratio (largo de intersección) / (perímetro)
    MAX_MILLER = 0.6  # Coeficiente de Miller
    main(NOMBRE, RADIOS, UMBRALES, MAX_RATIO, MAX_MILLER)
