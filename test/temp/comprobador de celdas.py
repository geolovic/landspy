#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:24:47 2019

@author: vicen
"""

# Comprobador de celdas
from topopy import DEM, Flow, Network
import numpy as np

threshold = 250
celdas = [248581, 178025, 168315, 247404, 98686, 44155, 108239]

dem = DEM("../data/in/morocco.tif")
fld = Flow(dem)
net = Network(fld, threshold)
canales = net.get_streams()
canales.save("../data/out/canales_celdas.tif")

confs = net.get_stream_poi(kind="confluences", coords="XY")

row, col = dem.ind_2_cell(celdas)
xi, yi = dem.cell_2_xy(row, col)
celdas = np.array((xi, yi)).T

np.savetxt("../data/out/celdas.txt", np.array(celdas), delimiter=";", comments="", header="x;y")
np.savetxt("../data/out/confs.txt", confs, delimiter=";", comments="", header="x;y")




