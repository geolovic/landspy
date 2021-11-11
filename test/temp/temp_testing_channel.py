# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Network, Channel, Flow, BNetwork
from osgeo import ogr
import numpy as np

#def rivers_to_channels(path, net, idfield=""):
    
path = "../data/in/jebja_channels.shp"
idfield=""
net = Network("../data/in/jebja30_net.dat")

# Open river shapefile
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(path)
layer = dataset.GetLayer()
geom_type = layer.GetGeomType()
lydef = layer.GetLayerDefn()
id_fld = lydef.GetFieldIndex(idfield)

points = []

for feat in layer:
    geom = feat.GetGeometryRef()
    if geom.GetGeometryCount() > 1:
        continue
    
    head = geom.GetPoint(0)
    mouth = geom.GetPoint(geom.GetPointCount()- 1)
    points.append([head, mouth])


canales = []
for canal in points:
    head = canal[0]
    mouth = canal[1]
    
    canales.append(net.get_channel(head, mouth))
    
import matplotlib.pyplot as plt
canal = canales[0]
fig = plt.figure()
ax = fig.add_subplot(111)
dir(canal)
canal.get_xy()
for canal in canales:
    xy = canal.get_xy()
    ax.plot(xy[:,0], xy[:,1])