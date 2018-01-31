#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:32:58 2018

@author: vicen
"""
import sys
import numpy as np
from scipy import ndimage
import time
# Add to the path code folder and data folder
sys.path.append("../../")
from topopy import DEM
from skimage import graph

start = time.time()

# Get DEM and fill sinks
dem = DEM("../data/tunez2.tif")
dims = dem._array.shape
ncells = dem._array.size

time_lap = time.time()
# 01 Fill sinks
fill = dem.fill_sinks()
topodiff = fill.read_array() -  dem.read_array()
topodiff = topodiff.astype(np.float32)
dem = fill
print("Sinks filled -- {0:.3f} seconds".format(time.time() - time_lap))
    
time_lap = time.time()  
# 02 Get flats and sills
flats, sills = dem.identify_flats(False)
print("Flats identified -- {0:.3f} seconds".format(time.time() - time_lap))

time_lap = time.time()  
# 03 Create auxiliar topography for 'carve' option
tweight = 2
carvemin = 0.1                       
struct = np.ones((3, 3), dtype="int8")
lbl_arr, nlbl = ndimage.label(flats.read_array(), structure=struct)
lbls = np.arange(1, nlbl + 1)

for lbl in lbls:
    topodiff[lbl_arr == lbl] = (topodiff[lbl_arr==lbl].max() - topodiff[lbl_arr == lbl])**tweight + carvemin
    
del lbl_arr
print("Auxiliar topography created -- {0:.3f} seconds".format(time.time() - time_lap))

time_lap = time.time()  
# 04 Identify presill pixels, i.e. pixels immediately upstream to sill pixels *
zvals = dem.read_array()
flats_arr = flats.read_array().astype("bool")
dims = zvals.shape
row, col = sills.find()

rowadd = np.array([-1, -1, 0, 1, 1,  1,  0, -1])
coladd = np.array([ 0,  1, 1, 1, 0, -1, -1, -1])
ps_rows = np.array([], dtype=np.int32)
ps_cols = np.array([], dtype=np.int32)

for n in range(8):
    rowp = row + rowadd[n]
    colp = col + coladd[n]
    # Avoid neighbors outside array (remove cells and their neighbors)
    valid_rc = (rowp >= 0) & (colp >= 0) & (rowp < dims[0]) & (colp < dims[1])
    rowp = rowp[valid_rc]
    colp = colp[valid_rc]
    # Discard cells (row-col pairs) that do not fullfill both conditions
    cond01 = zvals[row[valid_rc], col[valid_rc]] == zvals[rowp, colp]
    cond02 = flats_arr[rowp, colp]
    valid_pix = np.logical_and(cond01, cond02)
    ps_rows = np.append(ps_rows, rowp[valid_pix])
    ps_cols = np.append(ps_cols, colp[valid_pix])

ps_pos = zip(ps_rows, ps_cols)
print("Presills identified -- {0:.3f} seconds".format(time.time() - time_lap))

time_lap = time.time()  
# 05 Calculate weights for the cost-distance analysis
flats_arr = np.invert(flats_arr)
topodiff[flats_arr] = 99999
lg = graph.MCP_Geometric(topodiff)
topodiff = lg.find_costs(starts=ps_pos)[0] + 1
topodiff[flats_arr] = -99999
print("Weights calculated -- {0:.3f} seconds".format(time.time() - time_lap))

time_lap = time.time() 
# 06 Sort pixels 
# Sort the flat areas
rdem = zvals.ravel()
topodiff = topodiff.ravel()
ix_flats = np.argsort(-topodiff, kind='mergesort')

# Sort the rest of the pixels from the DEM
ndx = np.arange(ncells, dtype=np.int)
ndx = ndx[ix_flats]
ix = ndx[np.argsort(-rdem[ndx], kind='mergesort')]
print("Pixels sorted -- {0:.3f} seconds".format(time.time() - time_lap))

time_lap = time.time()
# 07 Look for receivers cells
pp = np.zeros(dims, dtype=np.int32)
IX = np.arange(ncells, dtype=np.int32)
pp = pp.ravel()
pp[ix] = IX
pp = pp.reshape(dims)

demarr = dem.fill_sinks().read_array()
f_dem = demarr.ravel()

# Get cardinal neighbors
footprint= np.array([[0, 1, 0],
                     [1, 1, 1],
                     [0, 1, 0]], dtype=np.int)
IXC1 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
xxx1 = np.copy(IXC1)
IX = IXC1.ravel()[ix]
IXC1 = ix[IX]
G1   = (rdem[ix]-rdem[IXC1])/(dem.get_cellsize())

# Get diagonal neighbors
footprint= np.array([[1, 0, 1],
                     [0, 1, 0],
                     [1, 0, 1]], dtype=np.int)
IXC2 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
xxx2 = np.copy(IXC2)
IX = IXC2.ravel()[ix]
IXC2 = ix[IX]
G2   = (rdem[ix]-rdem[IXC2])/(dem.get_cellsize() * np.sqrt(2))

del IX

# Get the steepest one
I  = (G1<=G2) & (xxx2.ravel()[ix]>xxx1.ravel()[ix])
ixc = IXC1
ixc[I] = IXC2[I]
I = ixc == ix
I = np.invert(I)
ix = ix[I]
ixc = ixc[I]
print("Receivers pixels identified -- {0:.3f} seconds".format(time.time() - time_lap))


# 08 Calculate flow accumulation
time_lap = time.time() 
nix = len(ix)

A = np.ones(ncells)
for n in range(nix):
    A[ixc[n]] = A[ix[n]] + A[ixc[n]]

A = A.reshape(dims)

acc = DEM()
acc.copy_layout(dem)
acc.set_array(A)
print("Flow accumulation calculated -- {0:.3f} seconds".format(time.time() - time_lap))

acc.save("flow_Arr3.tif")

