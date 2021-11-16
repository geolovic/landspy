# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import Network, Channel, Flow, BNetwork
import ogr
import numpy as np
import matplotlib.pyplot as plt


fd = Flow("../data/in/jebja30_fd.tif")
net = Network("../data/in/jebja30_net.dat")
cuencas = fd.get_drainage_basins()
bnet1 = BNetwork(net, cuencas, None, 3)
bnet2 = BNetwork(net, cuencas, None, 10)

shp_path = "../data/out/canales_basin.shp"
bnet1.export_to_shp(shp_path, True)

dataset = ogr.Open(shp_path)
layer = dataset.GetLayer(0)

canales = []
for n in range(layer.GetFeatureCount()):
    feat = layer.GetFeature(n)
    geom = feat.GetGeometryRef()
    head = geom.GetPoint(0)
    mouth = geom.GetPoint(geom.GetPointCount() - 1)
    canal = net.get_channel(head, mouth)
    canales.append(canal)
    
canales = np.array(canales)
np.save("../data/out/canales-b3_jebja30.npy", canales)
