# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 12:14:38 2018

@author: VicenPC
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from topopy import DEM, Flow, Network


fd = Flow("../data/in/morocco_fd.tif")
dem = DEM("../data/in/morocco.tif")

threshold = 1500

fac = fd.get_flow_accumulation(nodata=False, asgrid=False)
w = fac > threshold
w = w.ravel()
I   = w[fd._ix]
ix  = fd._ix[I]
ixc = fd._ixc[I]


di = np.zeros(fd._ncells)
for n in np.arange(ix.size)[::-1]:
    grow, gcol = fd.ind_2_cell(ix[n])
    rrow, rcol = fd.ind_2_cell(ixc[n])
    gx, gy = fd.cell_2_xy(grow, gcol)
    rx, ry = fd.cell_2_xy(rrow, rcol)
    d_gr = np.sqrt((gx - rx)**2 + (gy - ry)**2)
    di[ix[n]] = di[ixc[n]] + d_gr

dii = di.reshape(fd.get_dims())
