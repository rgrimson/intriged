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
    gdfo.to_file(str(fn))

    return None


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

    return None


# %% Generar polígono
def gen_poly(tipo='sintetico', nombre='pol_single_hole'):
    """Generar un polígono.

    Tipos y nombres implementados:
        tipo == 'sintetico', nombre == 'pol_single_hole' | 'pol_rafa'.
        tipo == 'fn', nombre == path-like object.
    """
    if tipo == 'sintetico':
        if nombre == 'pol_single_hole':
            ext_coords = ((0., 0.), (0., 15.), (15., 15.), (15., 9.), (18., 9.),
                    (21., 15.), (36., 15.), (36., 0.), (27., 0.), (27., 9.),
                    (33., 9.), (33., 12.), (24., 12.), (24., 0.), (18., 0.),
                    (18., 3.), (15., 6.), (15., 0.), (0., 0.))
            pol_single = Polygon(ext_coords)
            hole_coords = ((6., 3.), (6., 9.), (12., 9.), (9., 3.), (6., 3.))
            hole = Polygon(hole_coords)
            poly = pol_single.difference(hole)

        elif nombre == 'pol_rafa':
            coords = ((0., 0.), (0., 3.), (3., 3.), (3., 2.), (10.,2.), (10.,4.),
                    (15.,4.), (15., -1.), (10., -1.), (10., 1.), (3., 1.),
                    (3., 0.), (0., 0.))
            poly = Polygon(coords)

        else:
            msg = f'No está implementado el poligono sintético solicitado: {nombre}.'
            raise NotImplementedError(msg)

    elif tipo == 'fn':
        gdf = gpd.read_file(str(nombre))
        poly = gdf.iloc[0].geometry

    else:
        msg = f'No está implementado el tipo solicitado: {tipo}.'
        raise NotImplementedError(msg)

    return poly
