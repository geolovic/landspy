#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 10:06:47 2018

@author: vicen
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid

MY_GRID = "data/small25.tif"
ids = np.load("data/MY_GRID_100rnd_id.npy")
rows = np.load("data/MY_GRID_100rnd_row.npy")
cols = np.load("data/MY_GRID_100rnd_col.npy")
xi = np.load("data/MY_GRID_100rnd_X.npy")
yi = np.load("data/MY_GRID_100rnd_Y.npy")
zi = np.load("data/MY_GRID_100rnd_Z.npy")

dem = Grid(MY_GRID)

row = rows[10:20].tolist()
col = cols[10:20].tolist()
idx = ids[10:20].tolist()
print(idx)
print(dem.cell_2_ind(row, col))
