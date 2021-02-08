#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 march, 2020 (COVID19 quarantine)
Testing suite for BNetwork class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 26 march, 2020 (COVID19 quarantine)
"""

import unittest
import numpy as np
import os
from topopy import Flow, Basin, Network, BNetwork, DEM, extract_points
infolder = "data/in"
outfolder = "data/out"

file = "jebja30"
basin = Basin("{}/{}_sbasin.tif".format(infolder, file))
net = Network("{}/{}_net.dat".format(infolder, file))

# Generamos 50 puntos aleatorios dentro de la extensión del objeto Network
# Estos 50 puntos se usaran como cabeceras
# xmin, xmax, ymin, ymax = net.get_extent()
# xi = np.random.randint(xmin, xmax, 50)
# yi = np.random.randint(ymin, ymax, 50)
# heads = np.array((xi, yi)).T

# np.savetxt("data/out/cabeceras_total.txt", np.array((xi, yi)).T, delimiter=";", header="X;Y\n",
#            comments="")

heads = extract_points("data/out/test.shp", "id")

bnet = BNetwork(net, basin, heads)

row, col = bnet.ind_2_cell(bnet._heads)
x, y = bnet.cell_2_xy(row, col)

# snap_heads = net.snap_points(heads, "heads")
# x, y = snap_heads[:,0], snap_heads[:,1]
# inside = basin.is_inside(x, y)

# x = x[inside]
# y = y[inside]

np.savetxt("data/out/cabeceras_test.txt", np.array((x, y)).T, delimiter=";", header="X;Y\n",
            comments="")