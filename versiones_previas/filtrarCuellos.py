#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aplicar descomposición de partes a partir de detección de cuellos."""

import os

#import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
from tqdm import tqdm

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
#from shapely import intersects
#from shapely import union
#import itertools
from sklearn.cluster import DBSCAN


guardar_intermedios = True

#%%
def lpolys(PMP):
    """Generar una lista de polígonos singlepart."""
    #devuelve una lista con los polígonos simples que componen el P o MP.
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
def GC2MP(P):
    """Generar un MultiPolyogn a partir de una GeometryCollection."""
    if type(P) == GeometryCollection:
        return MultiPolygon([p for p in list(P.geoms) if type(p)==Polygon])
    else:
        return P


#%%
def topo(P):
    """Generar la topología de un polígono o multipolígono."""
    #devuelve una descripcion de la topología del polígono o del multip
    #es recursiva para multipolígonos
    if P:
        if type(P)==Polygon:
            return [len(P.interiors)]
        elif type(P)==MultiPolygon:
            return [topoP(Q) for Q in list(P.geoms)]
        elif type(P) == GeometryCollection:
            return [p for p in list(P.geoms) if type(p)==Polygon]
        else:
            print("Tipo no reconocido")


    else:
        return []

#%%
def topoP(P):
    """Generar la topologia de un polígono."""
    if P and type(P)==Polygon:
        return len(P.interiors)
    else:
        return []


#%%

# Funciones requeridas para la descomposición a partir de cuellos

def get_area(polygon):
    """Obtener el área de un polígono."""
    return polygon.area


def es_cuello(part_v, dif_element):
    """Verificar si dif_element es un cuello de part_v."""
    diferencia = (part_v.buffer(0)).difference(dif_element.buffer(0))
    return topo(diferencia) != topo(part_v)


def identificar_cuello(part_v, dif_list):
    """Identificar el elemento de dif_list que es el cuello de part_v."""
    # Ordenar las partes por área de mayor a menor
    dif_list = sorted(dif_list, key=get_area, reverse=True)

    # Iterar sobre cada parte, de mayor a menor área, hasta encontrar el cuello.
    for dif_element in dif_list:
        if es_cuello(part_v, dif_element):
            return dif_element.buffer(0)

    # Si no se encontró el cuello entre todos los elementos de dif_list:
    return None


#%% calcular radios donde cambia la topología

def calcular_saltos_filtracion(P, r_min=0, r_step=1):
    """Calcular los radios que producen cambio de topología."""
    # Lista de radios a devolver
    D = []

    # Topología del polígono original.
    ta = topo(P)

    # Primer distancia es el radio mínimo
    d = r_min

    # Mientras que la topología del buffereado no resulte vacía
    t = ta
    while t != []:

        # Buffer de cero para limpiar la topología, mencionado en:
        #  https://shapely.readthedocs.io/en/stable/manual.html#object.buffer
        Q = (P.buffer(-d).buffer(d)).buffer(0)
        t = topo(Q)

        # Analizar únicamente si cambia la cantidad de partes
        if len(t) != len(ta):
            D.append(d)

        d += r_step

    return D


#%%

# Algoritmo de agrupamiento provisto por inteligencia artificial.
# No estoy seguro que esto sea realmente lo que necesitamos.
def agrupar_radios(radios, eps, min_samples):
    """Agrupar radios y extraer radios máximos por cada grupo."""
    # Convertir la lista de radios en una matriz 2D para DBSCAN
    radios = np.array(radios).reshape(-1, 1)

    # Aplicar DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(radios)

    # Obtener las etiquetas de los grupos
    etiquetas = db.labels_

    # Agrupar los radios según las etiquetas
    grupos = {}
    for etiqueta, radio in zip(etiquetas, radios):
        if etiqueta not in grupos:
            grupos[etiqueta] = []
        grupos[etiqueta].append(radio[0])

    print(f'{grupos = }')

    # Excluir el ruido (etiquetas con valor -1)
    # y tomar el valor de radio máximo de cada grupo
    radios_definitivos = [max(grupo)
                          for etiqueta, grupo
                          in grupos.items() if etiqueta != -1]

    return radios_definitivos


#%% descomponer a partir de extraer los cuellos

# La hipótesis es que para un radio que hace cambiar la topología, entre las
#  geometrías que resultan de la diferencia entre un polígono y su buffereado,
#  la geometría que tiene mayor área es una parte significativa que separa a
#  otras dos.
# En una figura de cabeza, cuello y cuerpo, el buffereado de radio mayor que
#  el ancho mínimo del cuello / 2, contendría principalmente la cabeza y el
#  cuerpo. La diferencia entre la figura original y el buffereado contendría
#  muchos ángulos redondeados y el cuello, coleccionados en un multipolígono.
# El cuello, que sería la parte de mayor área, es una parte significativa que
#  separa a otras dos: en este caso la cabeza y el resto del cuerpo.

# Editado después de pruebas:
# La hipótesis resulta falsa en algunos casos, por lo que se analizan las
#  partes de la diferencia, de mayor a menor área, para saber cuál es el cuello
#  (el que hace cambiar la topología al restarse de la figura original).

# La función almacena el cuello como parte definitiva y la cabeza y el cuerpo
#  como geometrías modificadas, y vuelve a evaluar las geometrías modificadas
#  por si encuentra que vuelven a cambiar su topología para este radio.

def descomponer(P, D):
    """Descomponer P en partes significativas para los radios en D."""
    # Lista de listas de radio y polígono singlepart definitivos
    # [[radio, Polygon], ...]
    definitivos = []

    # Evaluar una lista de los polígonos singlepart que componen a P.
    partes_a_evaluar = [parte.buffer(0) for parte in lpolys(P)]

    # Iterar sobre los radios.
    for d in D:
        print(f'Analizando radio {d}.')

        # Lista de polígonos singlepart para evaluar en el próximo radio.
        proximo_radio = []

        # Recursividad mientras haya partes a evaluar.
        while partes_a_evaluar:
            print(f'Cantidad de partes a evaluar: {len(partes_a_evaluar)}.')

            # Iniciar una lista de geometrías modificadas (cabeza y cuerpo),
            #  para ser evaluadas en la próxima recursión.
            modificadas = []

            # Iterar sobre cada parte a evaluar.
            for part in partes_a_evaluar:

                # Obtener la topología de esta parte.
                topo_part = topo(part)

                # Hacer el buffer-in/buffer-out.
                buffered = (part.buffer(-d).buffer(d)).buffer(0)

                # Obtener la topología de buffered.
                topo_buffered = topo(buffered)

                # Si la topología de buffered es un conjunto vacío
                #  (es decir que buffered es una geometría vacía):
                if topo_buffered == []:
                    # Agregar esta parte a la lista de polígonos definitivos
                    #  porque con este radio, esta parte se pierde.
                    definitivos.append([d, part])

                # Si la topología no es vacía, pero no cambió respecto de la
                #  parte evaluada:
                elif topo_part == topo_buffered:
                    # Agregar esta parte a la lista a evaluar en el próximo radio.
                    proximo_radio.append(part)

                # Si no, entonces cambió la topología. Extraer el cuello a
                #  geometrías definitivas y las modificadas al próximo while
                #  para ser reevaluadas.
                else:
                    # Incluir vértices de intersección a las geometrías originales.
                    p_n_b = (part.intersection(buffered)).buffer(0)
                    p_b = (part.difference(buffered)).buffer(0)
                    b_p = (buffered.difference(part)).buffer(0)
                    part_v = p_n_b.union(p_b).buffer(0)
                    buffered_v = p_n_b.union(b_p).buffer(0)

                    # Calcular la diferencia.
                    dif = (part_v.difference(buffered_v)).buffer(0)

                    # Convertir a lista de polígonos singlepart.
                    dif_list = lpolys(dif)

                    # Encontrar el cuello entre la dif_list.
                    cuello = identificar_cuello(part_v, dif_list)
                    # Si no se encuentra cuello es porque todas las partes que
                    #  componen dif_list son protuberancias de part_v.
                    # No debería suceder (ya se analizó que este buffered
                    #  tiene distinta topología que part), pero puede suceder...
                    if not cuello:
                        print("No se encontró un cuello!")
                        proximo_radio.append(part)
                        continue

                    # Agregar el cuello a la lista de polígonos definitivos
                    #  para este radio (si se reevaluara, igual debería
                    #  perderse para este radio).
                    definitivos.append([d, cuello])

                    # Y agregar la diferencia entre la parte original que se
                    #  está evaluando y el cuello, a la lista de modificadas
                    #  (una vez que se terminen de evaluar todas las partes
                    #  originales, se reevalúan las partes modificadas).
                    p_c = (part_v.difference(cuello)).buffer(0)
                    modificadas.extend(lpolys(p_c))

            # Una vez que termine el for loop sobre las partes originales,
            #  evaluar las modificadas en el while, que termina cuando no hayan
            #  modificadas.
            partes_a_evaluar = modificadas

        # Una vez que se analizó todo el radio, las partes a evaluar son las
        #  que se agregaron a proximo_radio (no cambiaron topología).
        partes_a_evaluar = proximo_radio

    return definitivos


#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%
home_dir = os.path.expanduser("~")
wdir = os.path.join(home_dir, 'Projects/2024 - Filtracion/salado/')
fn = wdir + 'laguito'#'saladito_muy_corto'

gdf = gpd.read_file(fn + '.shp')
#miro solo el primer polígono
R=gdf.iloc[0].geometry #union(gdf.iloc[0]['geometry']...)

coords = ((0., 0.), (0., 3.), (3., 3.), (3., 2.), (10.,2.), (10.,4), (15.,4.),(15.,-1.),(10.,-1.), (10.,1.),  (3., 1.), (3., 0.), (0., 0.))
S=Polygon(coords)

#%%
import geopandas as gpd
#gdf = gpd.geodataframe.GeoDataFrame(pd.DataFrame(D.items(),columns=['cod','geometry']))

P = R

r_min = 0
r_step = 100
D = calcular_saltos_filtracion(P, r_min, r_step)
print(f'{D = }')

# Agrupar radios
eps = r_step * 1.5
min_samples = 1
radios = agrupar_radios(D, eps, min_samples)

# Resulta un problema si el polígono es muy irregular, casi todos los radios
#  son consecutivos
print(f'{radios = }')

definitivos = descomponer(P, D)

gdfo = gpd.geodataframe.GeoDataFrame(
    definitivos,
    columns=['radio', 'geometry'],
    crs=gdf.crs)

print('gdfo = ')
print(gdfo)
gdfo.to_file(wdir + 'agrupada.shp')
