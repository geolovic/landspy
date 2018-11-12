# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Network, BNetwork, Grid
import numpy as np
import ogr
import matplotlib.pyplot as plt

# Get network and basins
net = Network("/Users/vicen/Desktop/Tunez/gisdata/network_4k.net")
basins = Grid("/Users/vicen/Desktop/Tunez/gisdata/cuencas.tif").read_array()
idx = np.unique(basins)[1:]

bid = 1
basin = basins
heads = None

net2 = Network("/Users/vicen/Desktop/Tunez/gisdata/network_4k.net")
self = net2

# CREATE A BASIN OBJECT


# Set basin cells to 1
basin = np.where(basin==bid, 1, 0)
                
# Get basin extent
c1 = basin.max(axis=0).argmax() - 1
r1 = basin.max(axis=1).argmax() - 1
c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax() + 1 
r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax() + 1

# Check boundaries conditions
if c1 < 0:
    c1 = 0
if r1 < 0:
    r1 = 0
if c2 >= basin.shape[1]:
    c2 = basin.shape[1] - 1
if r2 >= basin.shape[0]:
    r2 = basin.shape[0] - 1
    
basin_cl = basin[r1:r2, c1:c2] # clipped basin

# Create the new grid
self._size = (basin_cl.shape[1], basin_cl.shape[0])
self._dims = (basin_cl.shape[0], basin_cl.shape[1])
geot = net._geot
ULx = geot[0] + geot[1] * c1
ULy = geot[3] + geot[5] * r1
self._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
self._cellsize = (geot[1], geot[5])
self._proj = net._proj
self._ncells = basin_cl.size
self._threshold = net._threshold
self._thetaref = net._thetaref
self._slp_npoints = net._slp_npoints
self._ksn_npoints = net._ksn_npoints

# Get only points inside the basin
basin = basin.astype(np.bool).ravel()
I = basin[net._ix]
self._ax = net._ax[I]
self._zx = net._zx[I]
self._dd = net._dd[I]
self._dx = net._dx[I]
self._chi = net._chi[I]
self._slp = net._slp[I]
self._ksn = net._ksn[I]
self._r2slp = net._r2slp[I]
self._r2ksn = net._r2ksn[I]

# Get new indices for the new grid
ix = net._ix[I]
ixc = net._ixc[I]
rowix, colix = net.ind_2_cell(ix)
rowixc, colixc = net.ind_2_cell(ixc)
nrowix = rowix - r1
ncolix = colix - c1
nrowixc = rowixc - r1
ncolixc = colixc - c1
self._ix = self.cell_2_ind(nrowix, ncolix)
self._ixc = self.cell_2_ind(nrowixc, ncolixc)

# Get basin heads
bheads = self.get_stream_poi(coords="IND")
aux_z = np.zeros(self._ncells)
aux_z[self._ix] = self._zx
belev = aux_z[bheads]
bheads = bheads[np.argsort(belev)]

# If heads are provided, append them first
if heads is not None:
    heads = self.snap_points(heads)
    bheads = bheads.tolist()



        