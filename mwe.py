#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genera un minimal working example para el que falle el buffer in."""

# %% Librerías
from pathlib import Path
from pprint import pprint

import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection

import shapely


# %% Imprimir Polígono
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


# %% Topologia de Polygon, MultiPolygon o GeometryCollection
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


# %% Topologia de Polygon
def topoP(P):
    """Generar la topologia de un polígono."""
    if P and type(P)==Polygon:
        return len(P.interiors)
    else:
        return []


# %% Crear shapefile a partir de una lista de geometrías
def shapefile_from_geom(geoms, crs, fn):
    """Guardar el shapefile de una lista de geometrías."""

    gdfo = gpd.geodataframe.GeoDataFrame(geoms, columns=['geometry'])
    gdfo = gdfo.set_crs(crs)
    gdfo.to_file(str(fn) + '_invalid.shp')


# %% Sacar i-esimo vertice
def sacar_i(R,i=1):
    """Sacar el i-esimo vertice de un polígono."""
    coords = list(R.exterior.coords)
    new_coords = coords[0:i]
    new_coords.extend(coords[i+1:])
    new_R = shapely.geometry.Polygon(new_coords)
    return new_R

# %% Simplificar
def simplificar(R):
    """Simplificar R hasta encontrar el mwe."""
    i=1
    coords = list(R.exterior.coords)
    l = len(coords)
    B=R.buffer(-437, resolution=5, cap_style=1, join_style=1, mitre_limit=2.0, single_sided=False)
    topoOrig=topo(B)

    while i<l-1:
        new_R=sacar_i(R,i)
        B=new_R.buffer(
            -437,
            resolution=5,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            single_sided=False
        )

        if topo(B)==topoOrig:
            R=new_R
            print('s',end='')
            coords = list(R.exterior.coords)
            l = len(coords)
        else:
            print('x',end='')
            i+=1
    print()

    return R


# %% Analizar segmentos
def analizar_segmentos(simplificado, min_seg, max_seg):
    """Calcular la cantidad de segmentos que hace que el mwe no falle."""
    results = {}
    for resolution in range(min_seg, max_seg):
        buffered = simplificado.buffer(-437, resolution=resolution)
        topo_buffered = topo(buffered)
        results[resolution] = topo_buffered

    return results

# %% Main
def main():
    """Leer un shapefile, simplificarlo y crear un nuevo shapefile."""
    #print(f'{shapely.__version__ = }')
    home_dir = Path.home()

    wdir = home_dir / 'Projects/2024 - Filtracion/salado/'

    #Gaby, con el MWE1 no me salen dos polígonos, no se porque. Vuelvo al falladito_invalid_P
    fn = wdir / 'mwe2' #'falladito_invalid_P' #'mwe1'

    gdf2 = gpd.read_file(str(fn) + '.shp')

    R = gdf2.iloc[0].geometry # union(gdf.iloc[0]['geometry']...)
    print(f'{R = }')
    #R = shapely.set_precision(R, 1/1024)

    #plot_polygon(R) #print(f'{topo(R) = }')
    B = R.buffer(-437)
    #plot_polygon(B)
    #print(f'{type(B) = }')
    topoOrig = topo(B)

    print(f'{topoOrig = }')

    r_simpl = simplificar(R)

    print(f'{r_simpl = }')

    shapefile_from_geom([r_simpl], gdf2.crs, fn)

    # Análisis de segmentos del buffer
    analisis_seg = analizar_segmentos(r_simpl, 5, 20)
    pprint(analisis_seg)


if __name__ == '__main__':
    main()
