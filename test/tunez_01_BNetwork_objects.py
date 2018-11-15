# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Grid, Network, BNetwork
import numpy as np
import ogr
import os

basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"

# Load network, basins and heads
net = Network(basedir + "/net_4k.net")
main_basins = Grid(basedir + "/basins.tif")
med_basin = Grid(basedir + "/medjerda_basin.tif")

# Get main heads
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(basedir + "/heads.shp")
layer = dataset.GetLayer()
main_heads = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    main_heads.append((geom.GetX(), geom.GetY(), idx))
main_heads = np.array(main_heads)

# Get medjerda heads
dataset = driver.Open(basedir + "/medjerda_head.shp")
layer = dataset.GetLayer()
med_heads = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    med_heads.append((geom.GetX(), geom.GetY(), idx))
med_heads = np.array(med_heads)

# Create main BNetwork objects and save them
if not os.path.exists(basedir + "/BNetworks"):
    os.mkdir(basedir + "/BNetworks")
    
bids = np.unique(main_basins.read_array())[1:] # To skip zero
for bid in bids:
    pos = np.where(main_heads[:,2] == bid)
    head = main_heads[pos]
    bnet = BNetwork(net, main_basins, head, bid)
    bnet.save(basedir + "/BNetworks/bnet_{0:02}.net".format(bid))
    
# Create Medjerda BNetwork object and save it
med_heads = med_heads[np.argsort(med_heads[:, 2])]
bnet = BNetwork(net, med_basin, med_heads, 1)
bnet.save(basedir + "/BNetworks/bnet_medjerda.net")