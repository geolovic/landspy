#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:15:55 2020

@author: vicen
"""

from topopy import Grid, DEM, Flow, Network
import ogr
import numpy as np

basin_path = "../data/in/jebja30_basins.tif"
myheads_path = "../data/in/jebja30_myheads.shp"

# Cargamos xy de cabeceras
dataset = ogr.Open(myheads_path)
layer = dataset.GetLayer()
xy = []
for feat in layer:
    geom = feat.geometry()
    xy.append([geom.GetX(), geom.GetY(), feat["id"]])
    
xy = np.array(xy)
if xy.shape[1] > 2:
    pos = np.argsort(xy[:, 2])
    xy = xy[pos]

# Objeto network
net = Network("../data/in/jebja30_network.net")
# Objeto basin
basin = Grid(basin_path)
# Cabeceras
heads = xy
# Id
bid = 1

# Get heads inside the basin

# Sort heads
# If "id" field is present

# If there is not "id" field, head by elevation

# Snap heads to network heads





