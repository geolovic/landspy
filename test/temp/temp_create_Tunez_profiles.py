# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Network, Channel
import ogr
import numpy as np
import matplotlib.pyplot as plt

# Test longitudinal profile
# Load profiles (TProfile)
net = Network("../data/in/tunez_net.dat")
net.calculate_gradients(20)
shp_path = "../data/in/tunez_channels.shp"

dataset = ogr.Open(shp_path)
layer = dataset.GetLayer(0)
canales =[] 

for n in range(layer.GetFeatureCount()):
    feat = layer.GetFeature(n)
    geom = feat.GetGeometryRef()
    head = geom.GetPoint(0)
    mouth = geom.GetPoint(geom.GetPointCount() - 1)
    canal = net.get_channel(head, mouth)
    canales.append(canal)
    
fig, ax = plt.subplots()
for canal in canales:
    xy = canal.get_xy()
    ax.plot(xy[:,0], xy[:,1], c="k")
    
np.save("../data/in/canales_tunez.npy", np.array(canales))