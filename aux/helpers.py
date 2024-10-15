#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Proveer funciones auxiliares para las filtraciones."""

# %% Librerías
import geopandas as gpd
import matplotlib.pyplot as plt

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
from shapely import unary_union

#import exceptions
from aux import exceptions


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


# %% Reinflar
def reinflado(P, d, eps=0.001, quad_segs=16):
    """Reinflar un polígono en una distancia `d`."""
    buffered = (P.buffer(-d, quad_segs=quad_segs)
                .buffer(d+eps, quad_segs=quad_segs))
    if not buffered.is_valid:
        textos = [
            'Buffered es invalido.',
            f'{P = }',
            f'{d = }',
            f'{buffered = }',
        ]
        msg = '\n'.join(textos)
        raise exceptions.InvalidGeometryError(msg)

    return buffered


# %% Agregar Vertices necesarios
def agregar_vertices(A,B,C=None):
    """Agregar los vértices necesarios para operaciones de solapamiento."""
    A_n_B = (A & B).buffer(0)
    if C:
        A_dif_B = (A - B - C).buffer(0)
    else:
        A_dif_B = (A - B).buffer(0)
    B_dif_A = (B - A).buffer(0)

    A_v = (A_n_B | A_dif_B).buffer(0)
    B_v = (A_n_B | B_dif_A).buffer(0)
    return A_v, B_v


# %% Obtener intersección y diferencia
def get_inter_diff(P, hojas, eps, quad_segs=16):
    """Obtener la intersección y la diferencia.

    Dados un polígono `P` y una lista de diccionarios de hojas `hojas`,
    obtener P & list y P - list.
    `hojas` tiene la estructura: [{'cod': cod, 'd': d, 'geometry': H}, ...].

    Devolver una lista `inter` de diccionarios de la forma:
    [{'cod': cod, 'd': d, 'geometry': inter_geom}, ...],
    y una lista `diff` de diccionarios de la forma:
    [{'geometry': diff_geom}, ...] para ser luego llenada con los atributos
    de las adyacencias.
    """
    # Calcular la intersección una por una para no perder la estructura de
    #  los diccionarios.
    def _intersecar(P, hoja_geom, eps, quad_segs):
        inter_geom = (hoja_geom
                      .buffer(eps, quad_segs=quad_segs)
                      .intersection(P)
        )
        if not inter_geom.is_valid:
            msg = "`intersecar` no es valido."
            raise exceptions.InvalidGeometryError(msg)

        return inter_geom

    # Agregar además un identificador numérico a cada hoja
    inter = [{
        'cod': hoja['cod'],
        'd': hoja['d'],
        'geometry': _intersecar(P, hoja['geometry'], eps, quad_segs)
    } for hoja in hojas]

    # Para las diferencias, extraer las geometrías de las hojas y crear un
    #  multipolígono por unaray_union, hacerles un buffer de `eps` y calcular
    #  las diferencias.
    # Devolver una estructura similar (lista de diccionarios) para después
    #  completar con datos de adyacencias.
    lista = [hoja['geometry'] for hoja in hojas]
    unidos = unary_union(lista)
    buffered = unidos.buffer(eps, quad_segs=quad_segs)

    if not buffered.is_valid:
        textos = [
            'Geometría no válida al hacer buffer de la lista unida.',
            f'{lista = }',
            f'{buffered = }',
        ]
        msg = '\n'.join(textos)
        raise exceptions.InvalidGeometryError(msg)

    diferencias = lpolys(P.difference(buffered))

    diff = [{
        'geometry': dif,
        'n': 0,  # Cuenta de adyacencias.
        'cods': [],  # Códigos de las hojas adyacentes.
        'dists': [],  # Distancias de filtracion de las hojas adyacentes.
        'Miller': 0
        } for dif in diferencias]

    return inter, diff


# %% Obtener una hoja por su id
def get_hoja(hoja_id, hojas):
    """Obtener una hoja por su código."""
    # Generar un iterador con las hojas que satisfacen que
    #  hoja['hoja_id'] == hoja_id
    #  y devolver el primero de ese iterador o None si está vacío.
    return next((hoja for hoja in hojas if hoja['hoja_id'] == hoja_id), None)


# %% Save polygon
def save_poly(P, fn='poly.shp'):
    """Guardar un polígono a shapefile."""
    gpd.GeoDataFrame(geometry=P).to_file(fn)


# %% Save polygons list
def save_plist(L, fn='poly.shp'):
    """Guardar una lista de polígonos a shapefile."""
    gpd.GeoDataFrame(geometry=L).to_file(fn)


# %% Crear shapefile a partir de lista de diccionarios
def shapefile_from_data(data, crs, fn):
    """Guardar el shapefile de una lista de geometrías."""
    gdf = gpd.geodataframe.GeoDataFrame(data, crs=crs)
    gdf.to_file(str(fn))


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
