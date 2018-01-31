#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 21:52:56 2018

@author: vicen
"""


import sys
sys.path.append("../../")
from topopy import Flow, Grid, DEM
import numpy as np
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt


arr = np.array([[49, 50, 50, 50],
                [49, 47, 47, 46],
                [49, 48, 45, 45],
                [49, 47, 46, 43],
                [46, 47, 47, 47]], dtype=np.int16)

dem = DEM()
dem.set_array(arr)
flow = Flow(dem)

ixcix = np.zeros(flow._ncells, np.int)
ixcix[flow._ix] = np.arange(len(flow._ix))

new_ind = 0
channel_points = [new_ind]
while ixcix[new_ind] != 0:
    new_ind = flow._ixc[ixcix[new_ind]]
    channel_points.append(new_ind)

print(flow._ix)
print(flow._ixc)
print(ixcix)
print(channel_points)