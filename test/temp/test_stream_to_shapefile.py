#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:53:42 2018

@author: vicen
"""

import sys
sys.path.append("../../")
from topopy import Flow, Grid
import numpy as np
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt

fd = Flow()
fd.load_gtiff("../data/in/fd_tunez.tif")
threshold = 1000

ixcix = np.zeros(fd._ncells, np.int)
ixcix[fd._ix] = np.arange(len(fd._ix))

x = 507194.338
y = 4060151.087
row, col = fd.xy_2_cell(x, y)

channel_ix = fd.cell_2_ind(row, col)
channel_points = [channel_ix]

new_ind = channel_ix
while ixcix[new_ind] != 0:
    new_ind = fd._ixc[ixcix[new_ind]]
    channel_points.append(new_ind)

#while ixcix[channel_points[-1]] != 0:
#    channel_points.append(fd._ixc[ixcix[channel_points[-1]]])

marr = np.zeros(fd._ncells, np.int)
marr[channel_points] = 1
marr = marr.reshape(fd._dims)
plt.imshow(marr)


#channel_ix = fd.cell_2_ind(row, col)
#channel_points = []
#
#add_ind = channel_ix
#channel_points.append(add_ind)
#
#while ixcix[add_ind] != 0:
#    add_ind = ixcix[add_ind]
#    channel_points.append(fd._ixc[add_ind])