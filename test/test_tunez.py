# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network, BNetwork
import numpy as np
import ogr
import matplotlib.pyplot as plt

basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"

# Get DEM, Flow and Network and save them
dem = DEM(basedir + "/srtm30_dem.tif")
fld = Flow (dem, verbose=True)
fld.save(basedir + "/srtm30_fd.tif")
net = Network(fld, 4000)
net.save(basedir + "/net_4k.net")
fac = fld.get_flow_accumulation()
fac.save(basedir + "/srtm30_fac.tif")
net.export_to_shp(basedir + "/network.shp", con=True)
streams = net.get_stream_order()
streams.save(basedir + "/red_strahler.tif")

## Load Flow and Network
#fld = Flow(basedir + "/srtm30_fd.tif")
#net = Network(basedir + "/net_4k.net")
#
#
#
#
#
## Open outlets and snap them
#driver = ogr.GetDriverByName("ESRI Shapefile")
#dataset = driver.Open(basedir + "/outlets.shp")
#layer = dataset.GetLayer()
#outlets = []
#for feat in layer:
#    geom = feat.GetGeometryRef()
#    idx = feat.GetField("id")
#    outlets.append((geom.GetX(), geom.GetY(), idx))
#dataset = None
#layer = None
#outlets = net.snap_points(np.array(outlets))
#basins = fld.get_drainage_basins(outlets)
#basins.save(basedir + "/basins.tif")
