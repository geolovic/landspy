# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Network, BNetwork, Grid
import numpy as np
import ogr
import matplotlib.pyplot as plt

# INPUT PARAMETERS
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"

# Parameters from "tunez_02_thetaref_analysis.py
r2 = [0.923, 0.941, 0.912]
mn = [0.2, 0.3, 0.4]

#########################################################
#########################################################

# Get network and basins
net = Network(basedir + "/net_4k.net")
basins = Grid(basedir + "/medjerda_basin.tif")

# Get basin heads and sort them by id
driver = ogr.GetDriverByName("ESRI Shapefile")
dataset = driver.Open(basedir + "/medjerda_head.shp")
layer = dataset.GetLayer()
heads = []
for feat in layer:
    geom = feat.GetGeometryRef()
    idx = feat.GetField("id")
    heads.append((geom.GetX(), geom.GetY(), idx))
dataset = None
layer = None

heads = np.array(heads)
heads = heads[np.argsort(heads[:, 2])]

bnet = BNetwork(net, basins, heads)
bnet.save(basedir + "/chi_analysis/bnet_medjerda.net")
fig = plt.figure()
n = 1

for thetaref in mn:
    bnet.calculate_chi(thetaref)
    main_ch = bnet.get_main_channel()
    ax = fig.add_subplot(1, 3, n)
#    ax.plot(bnet._chi, bnet._zx, color="0.75", ls="None", marker=".", ms=1)
    ax.plot(main_ch[:, 5], main_ch[:, 2], ls="-", c="0.3", lw=1)
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0)
    if n ==1:
        ax.set_xlabel("Chi [m]")
        ax.set_ylabel("Elevation [m]")
    else:
        ax.set_yticklabels([])
        ax.set_xlabel("Chi [m]")
    ax.set_title("$m/n = {0} (R^2 = {1})$".format(thetaref, r2[n-1]))
    n += 1
