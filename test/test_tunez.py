# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network, BNetwork
import numpy as np
import ogr
import matplotlib.pyplot as plt

## Get network
fld = Flow("C:/Users/Usuario/Desktop/tunez/srtm30_fd.tif")
net = Network("C:/Users/Usuario/Desktop/tunez/network_4k.net")
# Open outlets
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open("C:/Users/Usuario/Desktop/tunez/outlets.shp")
layer = dataset.GetLayer()
outlets = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    outlets.append((geom.GetX(), geom.GetY(), idx))
dataset = None
layer = None
#
## Open heads and sort them
#driver = ogr.GetDriverByName("ESRI Shapefile")
#dataset = driver.Open("data/in/heads.shp")
#layer = dataset.GetLayer()
#heads = []
#for feat in layer:
#    geom = feat.GetGeometryRef()
#    idx = feat.GetField("id")
#    heads.append((geom.GetX(), geom.GetY(), idx))
#heads = np.array(heads)
#heads = heads[np.argsort(heads[:,2])]
#dataset = None
#layer = None
#
outlets = net.snap_points(np.array(outlets))
outlets = np.array(outlets)

outlets = outlets[:2, :]
basins = fld.get_drainage_basins(outlets)
basins.save("C:/Users/Usuario/Desktop/tunez/cuencas.tif")
#
#bnet = BNetwork(net, basins, heads, 5)


#dem = DEM("C:/Users/Usuario/Desktop/tunez/srtm30_dem.tif")
#fld = Flow(dem, True)
#fld.save("C:/Users/Usuario/Desktop/tunez/srtm30_fd.tif")
#fac = fld.get_flow_accumulation()
#fac.save("C:/Users/Usuario/Desktop/tunez/srtm30_fac.tif")
#net = Network(fld, 4000, 0.45, 25)
#net.save("C:/Users/Usuario/Desktop/tunez/network_4k.net")
