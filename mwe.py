

from pathlib import Path
from pprint import pprint


from packaging import version
#import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
#import numpy as np
#import seaborn as sns
#from tqdm import tqdm

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
#from shapely import wkt
#from shapely import intersects
#from shapely import union

import shapely


#%%
# Imprimir Polígono
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


#%%
# Topologia de Polygon, MultiPolygon o GeometryCollection
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


#%%
# Topologia de Polygon
def topoP(P):
    """Generar la topologia de un polígono."""
    if P and type(P)==Polygon:
        return len(P.interiors)
    else:
        return []


#%%
# Crear shapefile a partir de geometría
def shapefile_from_geom(geoms):
    """Guardar el shapefile de un polígono."""

    gdfo = gpd.geodataframe.GeoDataFrame(geoms, columns=['geometry'])
    gdfo = gdfo.set_crs(gdf.crs)
    gdfo.to_file(str(fn) + '_invalid.shp')


#%%
def sacar_i(R,i=1):
    """Sacar el i-esimo vertice de un polígono."""
    coords = list(R.exterior.coords)
    new_coords = coords[0:i]
    new_coords.extend(coords[i+1:])
    new_R = shapely.geometry.Polygon(new_coords)
    return new_R

##%
def simplificar(R):
    """Simplificar R hasta encontrar el mwe."""
    i=1
    coords = list(R.exterior.coords)
    l = len(coords)
    while i<len(coords)-1:
        new_R=sacar_i(R,i)
        B=new_R.buffer(-437)
        if topo(B)==topoOrig:
            R=new_R
            print('s',end='')
            coords = list(R.exterior.coords)
            l = len(coords)
        else:
            print('x',end='')
            i+=1


##%

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
    R = shapely.set_precision(R, 1/1024)

    #plot_polygon(R) #print(f'{topo(R) = }')
    B=R.buffer(-437)
    #plot_polygon(B)
    #print(f'{type(B) = }')
    topoOrig=topo(B)

    print(f'{topoOrig = }')



if __name__ == '__main__':
    main()
