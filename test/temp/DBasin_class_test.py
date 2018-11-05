#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 October, 2018
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 23 October, 2018
"""

from topopy import Network, DEM, Flow, Grid
import matplotlib.pyplot as plt
import numpy as np
infolder = "../data/in"
outfolder = "../data/out"

dem = DEM(infolder + "/morocco.tif")
fd = Flow(dem)
net = Network(dem, fd, 1500)

streams = net.get_streams()
streams.save(outfolder + "/canales_orig.tif")

outlet = np.array([579213, 504282]).reshape(1, 2)

basin = fd.get_drainage_basins(net.snap_points(outlet), asgrid=False)
c1 = basin.max(axis=0).argmax()
r1 = basin.max(axis=1).argmax()
c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax()
r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax()
basin_cl = basin[r1:r2, c1:c2]

nrow = basin_cl.shape[0]
ncol = basin_cl.shape[1]
outgrid = Grid()

outgrid._size = (ncol, nrow)
outgrid._dims = (nrow, ncol)
geot = net._geot
ULx = geot[0] + geot[1] * c1
ULy = geot[3] + geot[5] * r1
outgrid._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
outgrid._cellsize = geot[1]
outgrid._proj = net._proj
outgrid._ncells = basin_cl.size
outgrid._nodata = 0
outgrid._array = basin_cl
outgrid._tipo = str(basin.dtype)

fd.get_drainage_basins(net.snap_points(outlet)).save(outfolder + "/basin_orig.tif")
outgrid.save(outfolder + "/basin_test.tif")


# GET ONLY POINTS INSIDE BASIN
basin = basin.astype(np.bool).ravel()
I = basin[net._ix]

ix = net._ix[I]
ixc = net._ixc[I]

# Get grid channel cells
w = np.zeros(net._ncells, dtype=np.int8)
w[ix] = 1
w[ixc] = 1
w = w.reshape(net._dims)

plt.imshow(w)

#I = basin.astype(np.bool).ravel()[net._ix]
#
#ix = net._ix[I]
#ixc = net._ixc[I] 
#
rowix, colix = net.ind_2_cell(ix)
nrowix = rowix - r1
ncolix = colix - c1
newix = outgrid.cell_2_ind(nrowix, ncolix)

#rowixc, colixc = net.ind_2_cell(ixc)
#rowixc -= r1
#colixc -= c1
#newixc = outgrid.cell_2_ind(rowixc, colixc)
#
w = np.zeros(basin_cl.size, dtype=np.int8)
w[newix] = 1
#w[newixc] = 1
w = w.reshape(outgrid._dims)

streams = Grid()
streams.copy_layout(outgrid)
streams._nodata = 0
streams._array = w
streams._tipo = str(w.dtype)

streams.save(outfolder + "/canales_prueba.tif")


