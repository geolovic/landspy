# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Network, BNetwork, Grid
import numpy as np
import ogr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# INPUT PARAMETERS
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"

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

bnet = BNetwork(net, basins, heads, 1)

desv = []

for thetaref in np.arange(0, 1, 0.01):
    bnet.calculate_chi(thetaref)
    chi = bnet._chi
    zi = bnet._zx
    poli, SCR = np.polyfit(chi, zi, deg = 1, full = True)[:2]
    r2 = float(1 - SCR/(zi.size * zi.var()))
    desv.append((thetaref, r2))
    
desv = np.array(desv)

fig, ax = plt.subplots()


ax.set_xlabel("$m/n$")
ax.set_ylabel("$R^2$")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

# Get maximum R2
pos = np.argmax(desv[:,1])

ax.plot([0, 0.3, 0.3], [0.9407316, 0.9407316, 0], linestyle="-", c="orange", lw=0.8)
ax.grid(True, which='major', axis="both", linestyle ="--", c="0.7", lw=0.5)
ax.plot(desv[:,0], desv[:, 1])