# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network, BNetwork
import numpy as np
import matplotlib.pyplot as plt


dem = DEM("data/in/ttunez.tif")
flow = Flow(dem)
net = Network(flow, 4000)

streams = net.get_streams()
streams.save("data/in/ttunez_str.tif")

outlet = [548253, 4106795]
head = np.array([531292, 4096471]).reshape(1, 2)

basin = flow.get_drainage_basins([548253, 4106795], asgrid=False)
bnet = BNetwork(net, basin, head)
bnet.calculate_distances()
#fig, ax = plt.subplots()
#
canal = bnet.get_main_channel()

#row, col = bnet.ind_2_cell(canal)
#xi, yi = bnet.cell_2_xy(row, col)
#plt.plot(xi, yi)