#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 march, 2020 (COVID19 quarantine)
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 28 march, 2020 (COVID19 quarantine)
"""

from topopy import Network, Grid, Basin, Flow, DEM
import numpy as np

# def __init__(self, net, basingrid=None, heads=None, bid=1):
file = "jebja30"

# Cargamos Flow y Network
fd = Flow("../data/in/{}_fd.tif".format(file))
net = Network("../data/in/{}_network.net".format(file))
dem = DEM("../data/in/{}.tif".format(file))
cuencas = fd.get_drainage_basins(min_area = 0.0025)   

# Parameteres to fake a BNetwork
class Foo:
    pass

bnet = Foo()
self = bnet

bids, counts = np.unique(cuencas.read_array(), return_counts=True)
if 0 in bids:
    bids = bids[1:]
    counts = counts[1:]
idmax = bids[np.argmax(counts)]
idmin = bids[np.argmin(counts)]
idrdn = np.random.choice(bids)
bids = [idmin, idmax, idrdn]

basingrid = cuencas
bid= bids[0]

xmin, xmax, ymin, ymax = net.get_extent()
xi = np.random.randint(xmin, xmax, 50)
yi = np.random.randint(ymin, ymax, 50)

heads = np.array((xi, yi)).T

# Get a basin grid with net dimensions (basin cells will set to 1)
if isinstance(basingrid, Basin):
    basin = np.zeros(net.get_dims(), dtype=np.int8)
    x = basingrid.get_geotransform()[0] + (basingrid.get_geotransform()[1]/2)
    y = basingrid.get_geotransform()[3] + (basingrid.get_geotransform()[5]/2)
    row, col = net.xy_2_cell(x, y)
    arr = np.where(basingrid.read_array()==basingrid.get_nodata(), 0, 1)
    basin[row:row+basingrid.get_dims()[0], col:col+basingrid.get_dims()[1]] = arr        
elif isinstance(basingrid, Grid):
    basin = np.where(basingrid.read_array()==bid, 1, 0)

# # Get limits for the input basin
# c1 = basin.max(axis=0).argmax()
# r1 = basin.max(axis=1).argmax()
# c2 = basin.shape[1] - np.fliplr(basin).max(axis=0).argmax()
# r2 = basin.shape[0] - np.flipud(basin).max(axis=1).argmax()

# # Cut basin by those limits
# basin_cl = basin[r1:r2, c1:c2]

# # Create Grid
# self._size = (basin_cl.shape[1], basin_cl.shape[0])
# self._dims = (basin_cl.shape[0], basin_cl.shape[1])
# geot = basingrid._geot
# ULx = geot[0] + geot[1] * c1
# ULy = geot[3] + geot[5] * r1
# self._geot = (ULx, geot[1], 0.0, ULy, 0.0, geot[5])
# self._cellsize = (geot[1], geot[5])
# self._proj = basingrid._proj
# self._ncells = basin_cl.size
# self._ncells = basin_cl.size
# self._threshold = net._threshold
# self._thetaref = net._thetaref
# self._slp_npoints = net._slp_npoints
# self._ksn_npoints = net._ksn_npoints

# # Get only points inside the basin
# basin_bool = basin.astype(np.bool).ravel()
# I = basin_bool[net._ix]
# self._ax = net._ax[I]
# self._zx = net._zx[I]
# self._dd = net._dd[I]
# self._dx = net._dx[I]
# self._chi = net._chi[I]
# self._slp = net._slp[I]
# self._ksn = net._ksn[I]
# self._r2slp = net._r2slp[I]
# self._r2ksn = net._r2ksn[I]

# # Get new indices for the new grid
# ix = net._ix[I]
# ixc = net._ixc[I]
# # Givers (ix)
# row, col = net.ind_2_cell(ix)
# x, y = net.cell_2_xy(row, col)
# newrow, newcol = self.xy_2_cell(x, y)
# self._ix = self.cell_2_ind(newrow, newcol)
# # Receivers (ixc)
# row, col = net.ind_2_cell(ixc)
# x, y = net.cell_2_xy(row, col)
# newrow, newcol = self.xy_2_cell(x, y)
# self._ixc = self.cell_2_ind(newrow, newcol)
# self._heads = np.array([])

# if heads is not None:
#     # Sort heads if "id" field is present
#     if heads.shape[1] > 2:
#         pos = np.argsort(heads[:, 2])
#         heads = heads[pos]
    
#     # Snap heads to network heads
#     heads = net.snap_points(heads, "heads")                   
    
#     # Get heads inside the basin (taking into account nodata)
#     heads = heads[basingrid.is_inside(heads[:,0], heads[:,1])]
#     row, col = basingrid.xy_2_cell(heads[:,0], heads[:,1])
#     pos = np.where(basin[row, col] > 0)
#     heads = heads[pos]
#     row, col = self.xy_2_cell(heads[:,0], heads[:,1])
#     self._heads = self.cell_2_ind(row, col)
#     self._heads = np.unique(self._heads)
    
# elif self._heads.size == 0:
#     heads = self.get_stream_poi("heads", "IND")
#     ixcix = np.zeros(self._ncells, np.int)
#     ixcix[self._ix] = np.arange(self._ix.size)
#     pos = np.argsort(-self._zx[ixcix[heads]])
#     self._heads = heads[pos][0:1] #Take only the first (highest) head