# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Flow, DEM, Network, BNetwork, Grid
import numpy as np
import ogr
import matplotlib.pyplot as plt

# Get network and basins
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"
dem = DEM(basedir + "/srtm30_dem.tif")
fld = Flow(dem, verbose=True)
#net = Network(basedir + "/net_4k.net")

# Open outlets and snap them
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(basedir + "/medjerda_outlet.shp")
layer = dataset.GetLayer()
outlets = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    outlets.append((geom.GetX(), geom.GetY(), idx))
dataset = None
layer = None

outlets = net.snap_points(np.array(outlets))
basins = fld.get_drainage_basins(outlets)
basins.save(basedir + "/medjerda_basin.tif")
