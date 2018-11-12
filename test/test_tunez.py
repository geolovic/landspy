# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network, BNetwork
import numpy as np
import ogr
import matplotlib.pyplot as plt

# Get network
dem = DEM("/Users/vicen/Desktop/Tunez/gisdata/srtm30_dem.tif")
fld = Flow (dem, verbose=True)
#fld = Flow("/Users/vicen/Desktop/Tunez/gisdata/srtm30_fd.tif")
net = Network("/Users/vicen/Desktop/Tunez/gisdata/network_4k.net")
#
# Open outlets
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open("/Users/vicen/Desktop/Tunez/gisdata/outlets.shp")
layer = dataset.GetLayer()
outlets = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    outlets.append((geom.GetX(), geom.GetY(), idx))
dataset = None
layer = None

## Snap outlets and get drainage basins
outlets = net.snap_points(np.array(outlets))
outlets = np.array(outlets)
basins = fld.get_drainage_basins(outlets)
basins.save("/Users/vicen/Desktop/Tunez/gisdata/cuencas.tif")


#bnet = BNetwork(net, basins, heads, 5)


#dem = DEM("C:/Users/Usuario/Desktop/tunez/srtm30_dem.tif")
#fld = Flow(dem, True)
#fld.save("C:/Users/Usuario/Desktop/tunez/srtm30_fd.tif")
#fac = fld.get_flow_accumulation()
#fac.save("C:/Users/Usuario/Desktop/tunez/srtm30_fac.tif")
#net = Network(fld, 4000, 0.45, 25)
#net.save("C:/Users/Usuario/Desktop/tunez/network_4k.net")
