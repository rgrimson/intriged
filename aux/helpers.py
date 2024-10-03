#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Proveer funciones auxiliares para las filtraciones."""

# %% Librerías
import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection


# %% Listar Polígonos
def lpolys(PMP):
    """Generar una lista de polígonos singlepart.

    Devuelve una lista con los polígonos simples que componen el P o MP.
    """
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


# %% Convertir GeometryCollection a MultiPolygon
def GC2MP(P):
    """Generar un MultiPolyogn a partir de una GeometryCollection."""
    if type(P) == GeometryCollection:
        return MultiPolygon([p for p in list(P.geoms) if type(p)==Polygon])
    else:
        return P


# %% Calcular topología
def topo(P):
    """Calcular la topología de un polígono.

    El tipo de geometría puede ser Polygon, MultiPolygon, o GeometryCollection.
    """
    def topoP(P):
        if P and type(P)==Polygon:
            return len(P.interiors)
        else:
            return []

    # Devuelve una descripcion de la topología del polígono, multipolígono
    #  o componentes poligonales de una geometry collection.
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


# %% Crear shapefile a partir de una lista de geometrías
def shapefile_from_geom(geoms, crs, fn):
    """Guardar el shapefile de una lista de geometrías."""

    gdfo = gpd.geodataframe.GeoDataFrame(geoms, columns=['geometry'])
    gdfo = gdfo.set_crs(crs)
    gdfo.to_file(str(fn) + '_invalid.shp')

    return 0


# %% Imprimir Polígono
def plot_polygon(polygon):
    """Imprimir un polígono o multipolígono."""
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

    return 0
