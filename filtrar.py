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

#%%
def topo(P):
    #devuelve una descripcion de la topología del polígono o del multip
    #es recursiva para multipolígonos
    if P:
      if type(P)==Polygon:
        return [len(P.interiors)]
      elif type(P)==MultiPolygon:
        return [topoP(Q) for Q in list(P.geoms)]
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
    
    
#%% calcular radios donde cambia la topología
def calcular_saltos_filtracion(P,r_min = 0, r_max=1000, r_step=1):
    ta=[] # inicializo la topología para que guarde el primer radio
    D = []# radios guardados
    T=[]  # topologías
    Qo=P
    for d in tqdm(range(r_min, r_max, r_step)):
        Q=(Qo.buffer(-d).buffer(d)).normalize()
        t = topo(Q)
        if t!=ta:
            D.append(d)
            T.append(t)
            ta=t
            Qo=Q
            #print(d,t,len(t))
    
    if t!=[]: #si no llegúe al conjunto vación, agrego un radio grande.
        point1, point2 = P.bounds[:2], P.bounds[2:]
        dist = round(np.linalg.norm(np.array(point1) - np.array(point2))/2+1)
        D.append(dist) #agrego un radio grande
        T.append([]) #y la topología del conjunto vacío
        
    return D,T
#%% computar y guardar la filtracion (los polígonos)
#from shapely import union_all 

def calcular_filtracion(P, D):
    #dado un polígono y un conjunto de distancias
    #calcula la filtración de P para esas D
    F=[] # inicializo la filtración vacía
    K=[] # y los radios como vacios
    
    Qo=P #inicializo el polígono anterior, es el primero
    for d in tqdm(D):
        #print(d)
        Q=Qo.buffer(-d).buffer(d).normalize()
        dif=(Qo-Q).normalize() #diferencias entre etapas (radios)
        Qo=Q
        if d>0:
            L=lpolys(dif) #separo la diferencia en polígonos simples 
            F+=L          #los agrego a la filtración 
            K+=[d]*len(L) #indicando el radio en que aparecen
    return F,K

#%% def compute intersection matrix
from rtree import index

def compute_intersections(F):
    #calcula las intersecciones de elementos de F
    #devuelve la lista de intersecciones con elementos de mayor indice
    Int=[[] for f in F]
    idx = index.Index()

    # Completar el indice R-tree usando el bounds de cada poly
    print("Armando indice")
    i = 0
    for p in tqdm(F):
        idx.insert(i, p.bounds)
        i+=1
    
    #armo la lista de intersecciones con elementos de mayor indice
    print("Calculando intersecciones")
    Int=[[] for f in F] #inicializo vacía
    i = 0
    for p in tqdm(F):
        for j in idx.intersection(p.bounds):
            pi = p
            if j>i:
                q = F[j]
                if intersects(p, q):
                    Int[i].append(j) #si hay interseccion, agregarla
        if Int[i]==[]: #si no le encontré intersecciones, reviso
          for j in idx.intersection(p.bounds):
            pi = p.buffer(0.001) #lo inflo un milimetro
            if j>i:
                q = F[j]
                if intersects(pi, q):
                    Int[i].append(j)
        i+=1
    return Int
#%% funcion calcular_relaciones

def calcular_relaciones(F,K,Int):
    #armo una lista de padres de cada polígono
    #un polígono es padre de otro si se tocan y es mas viejo (tiene mayor indice)
    

    Padres = [[]]*len(F) #inicializo lista de padres como vacía
    RadioPatriarcas = {} #dicc de radios para los patriarcas
    K = np.array(K)
    R=sorted(np.unique(K),reverse=True)
    N=np.zeros(len(F),dtype=int) # cuantos padres tiene
    for d in tqdm(R):
        Kd=np.where(K==d)[0]
        for j in Kd: #para los j de radio d...
            #print("j: ",j)
            I = Int[j] 
            for i in I: #agrego todos los padres que tenga
                if N[j]==0:
                    Padres[j]=Padres[i]
                    N[j]=N[i]
                if N[j]>=1 and Padres[j]!=Padres[i]:
                    Padres[j] = list(np.unique(Padres[j]+Padres[i]))
                    N[j]=len(Padres[j])
                        
            if N[j]==0: #si no tiene padres, es patriarca
                    Padres[j]=[j] #es su propio padre
                    N[j] = 1
                    RadioPatriarcas[j]=d  #lo anoto en el dicc, con su radio
    return Padres, RadioPatriarcas
            
#%%
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#%%
    
wdir = '/home/rgrimson/Downloads/estructura/salado/'
fn = wdir + 'poly1_utm21s'#'saladito_muy_corto'

gdf = gpd.read_file(fn + '.shp')
#miro solo el primer polígono
P=gdf.iloc[0].geometry #union(gdf.iloc[0]['geometry']...)

#%% calcular y mostrar puntos para la filtración
print("calcular saltos para filtracion")
D,T = calcular_saltos_filtracion(P)
for d,t in zip(D,T):
    print(d,t)
#%% plotear la cantidad de componentes conexas para cada d, salvo el ultimo
plt.scatter(D[:-1],[len(t) for t in T[:-1]])
plt.xlabel("radio en m del buffer")
plt.ylabel("cantidad de componentes")


#%%
print("calcular filtracion")
F,K = calcular_filtracion(P, D)
data =  {'d': K, 'geometry': F}

#%%
if guardar_intermedios:
    print("guardar filtracion")
    gdf = gpd.GeoDataFrame(data,crs='epsg:32721')
    gdf.to_file( fn + '_filt.shp')

#%%
print("Calculando adyacencias entre partes")
Int = compute_intersections(F)

#%%
print("Defino Patriarcas y Padres de cada poligono")
Padres, RadioPatriarcas = calcular_relaciones(F,K,Int)

#%% guardar la relacion de paternidad   
if guardar_intermedios:
    print("")
    D =  {'d': K, 'geometry': F, 'P':Padres}
    gdf = gpd.GeoDataFrame(D,crs='epsg:32721')
    gdf.to_file( fn + '_filt_Paternidad.shp')


#ACA TERMINA LA FILTRACIÓN. AHORA EMPIEZO A UNIR


#%% como prueba,me quedo con los que tienen un solo padre
# print("un solo padre")
# # un dissolve y lo guardo
# N=np.array([len(Padres[i]) for i in range(len(F))]) #cuantos padres
# Padre=[l[0] for l in np.array(Padres,dtype=object)[N==1]] #los que tienen uno
# dg = gdf[N==1].copy() #me quedo con los bien definidos (un solo padre)
# dg['C'] = Padre #lo anoto
# dg['radio']=[RadioPatriarcas[p] for p in dg['C']]
# gdp=dg.dissolve('C')
# gdp.to_file( fn + '_filt_P2_Diss.shp')

# plt.hist(RadioPatriarcas.values(),bins=100,log=True)
#%% ahora uso umbrales para armar clases
# la idea es juntar poligonos adyacentes por debajo de cada umbral

print('defino clases a partir de los umbrales')
umbrales = [230,600]
en más tráfico que los ca
PU=[]
ClasePatriarcas = {}
uo = 0 #umbral anterior, comienza en 0
clase = 0

for u in umbrales:
    P=[p for p in RadioPatriarcas if (RadioPatriarcas[p]<=u) and (RadioPatriarcas[p]>uo)]
    PU.append(P)
    for p in P:
        ClasePatriarcas[p]=clase 
    clase += 1
    
    uo=u

P=[p for p in RadioPatriarcas if (RadioPatriarcas[p]>uo)]
PU.append(P)
for p in P:
    ClasePatriarcas[p]=clase 
clase += 1
    
Clase = [max([ClasePatriarcas[p] for p in L]) for L in Padres]
gdf['C'] = Clase #lo anoto en el gdf

#%%
# #guardo intermedio
# gdq=gdf.dissolve('C')
# gdq.geometry = [clean_interiors(P) for P in gdq.geometry]
# gdq.to_file( fn + '_filt_P2_Diss_U_clean.shp')


#%%
from itertools import combinations

AdyPatriarcas = {} #a cada patriarca le asocio sus adyacentes

for L in Padres:
    if len(L) == 1:
        p = L[0]
        if not p in AdyPatriarcas:
          AdyPatriarcas[p] = []
    else:
        pairs = combinations(L,2)
        for p,q in pairs:
            if not p in AdyPatriarcas:
                AdyPatriarcas[p]=[q]
            elif not q in AdyPatriarcas[p]:
                AdyPatriarcas[p].append(q)
            if not q in AdyPatriarcas:
                AdyPatriarcas[q]=[p]
            elif not p in AdyPatriarcas[q]:
                AdyPatriarcas[q].append(p)
                
#%% defino una etiqueta para cada patriarca
# uso la misma para patriarcas adyacentes de la misma clase

clasePatriarca = [[ClasePatriarcas[q],q] for q in ClasePatriarcas]
clasePatriarca.sort() 
EtiqPatriarcas={}
etiqueta=0
for c,p in clasePatriarca:
    print('---------')
    print(c,p)
    if p not in EtiqPatriarcas:
        etiqueta+=1
        EtiqPatriarcas[p] = etiqueta
        print(f'nueva etiq: {EtiqPatriarcas[p]}')
        
        Ady = AdyPatriarcas[p]
        print(Ady)
        nuevo = True
        
        while nuevo:
            nuevo = False
            for q in Ady:
                if (ClasePatriarcas[q] == c) and (q not in EtiqPatriarcas):
                    EtiqPatriarcas[q] = etiqueta
                    print('  ', q,'-->',etiqueta)
                    Ady+=AdyPatriarcas.get(q)
                    nuevo = True
            #if nuevo:
            #    Ady = list(itertools.chain.from_iterable([AdyPatriarcas.get(q) for q in AdyPatriarcas[p] if AdyPatriarcas.get(q)]))
            print(Ady)
    else:
        print(f'ya esta {EtiqPatriarcas[p]}')
#%% a cada polígono le asigno una etiqueta. 
Etiq = []
for L in Padres:
    if len(L) == 1:
        p = L[0]
        Etiq.append(EtiqPatriarcas[p])
    else:
        CP = [[ClasePatriarcas[q],q] for q in L]
        clasePatriarca.sort() 
        p = CP[0][1]
        Etiq.append(EtiqPatriarcas[p])
        
#%% disuelvo, limpio  y guardo. Plotear categorizado por "E"
print('disuelvo')
gdf['E'] = Etiq #lo anoto
gde=gdf.dissolve('E')

print('limpio')
gde.geometry = [clean_interiors(P) for P in gde.geometry]

print('guardo')
gde.to_file( fn + '_filt_etiq.shp')


