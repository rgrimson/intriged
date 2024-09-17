#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 10:36:32 2024

@author: rgrimson
"""

#import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
from tqdm import tqdm

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.collection import GeometryCollection
from shapely import intersects
#from shapely import union
#import itertools


guardar_intermedios = True

#%%
def lpolys(PMP):
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
    if type(P) == GeometryCollection:
       return MultiPolygon([p for p in list(P.geoms) if type(p)==Polygon])
    else:
        return P

    
#%%
def topo(P):
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
    if P and type(P)==Polygon:
        return len(P.interiors)
    else:
        return []

#%% para sacar errores de tanto partir y unir
def clean_interiors_P(P,eps=0.0001):
        list_interiors = []
        
        for interior in P.interiors:
            p = Polygon(interior)

            if p.area > eps:
                list_interiors.append(interior)

        return Polygon(P.exterior.coords, holes=list_interiors)
        
def clean_interiors_MP(P,eps=0.0001):
    list_parts = []

    for polygon in P.geoms:
        list_parts.append(clean_interiors_P(polygon,eps))
    return MultiPolygon(list_parts)

def clean_interiors(PMP,eps=0.0001):
    if type(PMP)==Polygon:
        return clean_interiors_P(PMP,eps)
    elif type(PMP)==MultiPolygon:
        return clean_interiors_MP(PMP,eps)
    
    
#%% calcular filtración de manera recursiva
def calcular_filtracion_recursiva(P,r,cod,r_step=1):
    if type(P)!=Polygon:
        a = 0
        b = a/a
        return b
    ta = topo(P)
    if ta==[]:
        return Polygon()
    print(cod,r,ta)

    d = r + r_step    
    Q=(P.intersection(P.buffer(-d).buffer(d))).normalize()
    t = topo(Q)
    while len(t)==len(ta): #busco el primer cambio en la topología
        d+=r_step
        Q=(P.intersection(P.buffer(-d).buffer(d))).normalize()
        t = topo(Q)


    LQ = lpolys(Q)
    if len(LQ)>1:
        rta = []
        for i,Q in enumerate(LQ):
            frQ = calcular_filtracion_recursiva(Q,d,cod+str(i),r_step)
            if frQ:
                rta.append((Q,frQ))
            else:
                rta.append(Q)
        return rta
    elif LQ[0]:
        return calcular_filtracion_recursiva(LQ[0],d,cod,r_step)
    else:
        return Polygon()
        
#%%
def calcular_filtracion(P):
    return (P,calcular_filtracion_recursiva(P,0,'',r_step=1))

#%%
def dicc_filtracion(FP,cod=''):
    if type(FP)==Polygon:
        d = {cod:FP}
    else:
        d = {cod:FP[0]}
        for i,p in enumerate(FP[1]):
            d.update(dicc_filtracion(p,cod+str(i)))
    return d
#%%
def dicc_patriarcas(FP,cod=''):
    if type(FP)==Polygon:
        d = {cod:FP}
    else:
        #d = {cod:FP[0]}
        d={}
        if len(FP[1])>1:
            for i,p in enumerate(FP[1]):
                d.update(dicc_patriarcas(p,cod+str(i)))
    return d

#%%
import shapely


def agrupar_filtracion(F):
    if type(F)==Polygon:
        return [F]
    else:
        PG = F[0] #poligono grande
        DescPG = F[1] #descomposición del PG
        DescPG = [agrupar_filtracion(P) for P in DescPG]
        #lista de poligonos chicos
        LPC = []
        for LP in DescPG:
            LPC.extend(LP)
        D = PG-shapely.union_all(LPC) #componentes del grande que no estan en el chico
        LD = lpolys(D) #como lista
        for i,p in enumerate(LD): #las miro una a una
            J=[]
            print(i,end=': ')
            for j,q in enumerate(LPC): #me fijo que PChicos tocas
                if intersects(p, q):
                    #print(i,j)
                    J.append(j)
            if len(J)>1:
                print(J)
            elif len(J)==1:
                j=J[0]
                q=LPC[j]
                LPC[j]=q.union(p)
            else:
                print("PROBLEMA")
        return LPC
                    
#%%                    

            
#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%
    
wdir = '/home/rgrimson/Projects/2024 - Filtracion/salado/'
fn = wdir + 'laguito'#'saladito_muy_corto'

gdf = gpd.read_file(fn + '.shp')
#miro solo el primer polígono
R=gdf.iloc[0].geometry #union(gdf.iloc[0]['geometry']...)

coords = ((0., 0.), (0., 3.), (3., 3.), (3., 2.), (10.,2.), (10.,4), (15.,4.),(15.,-1.),(10.,-1.), (10.,1.),  (3., 1.), (3., 0.), (0., 0.))
S=Polygon(coords)

#%%
import geopandas as gpd
#gdf = gpd.geodataframe.GeoDataFrame(pd.DataFrame(D.items(),columns=['cod','geometry']))

P=R
FP = calcular_filtracion(P)
#%%
#guardar filtración
D=dicc_filtracion(FP,cod='')
gdfo = gpd.geodataframe.GeoDataFrame(D.items(),columns=['cod','geometry'])
gdfo=gdfo.set_crs(gdf.crs)
gdfo.to_file(wdir + '_filt.shp')    
#%%
# guardar patriarcas

DP=dicc_patriarcas(FP,cod='')
gdfo = gpd.geodataframe.GeoDataFrame(DP.items(),columns=['cod','geometry'])
gdfo=gdfo.set_crs(gdf.crs)
gdfo.to_file(fn + '_leaves.shp')    

#%%

FA=agrupar_filtracion(FP)
gdfo = gpd.geodataframe.GeoDataFrame(FA,columns=['geometry'])
gdfo=gdfo.set_crs(gdf.crs)
gdfo.to_file(wdir + 'agrupada.shp')    

#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
