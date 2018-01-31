#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:32:58 2018

@author: vicen
"""
import warnings
warnings.filterwarnings('ignore')

import sys
import numpy as np
from scipy import ndimage
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import DEM
from topopy.ext.distance import cost
from mcompare import compare_indexes, compare_row_cols, compare_arrays, load_array
from skimage import graph
import scipy.io as sio
dem = DEM("../data/tunez.tif")

# USED TO DEBUG
dims = dem._array.shape
ncells = dem._array.size
   

def _get_aux_topography(dem):
    """
    This function calculate the auxiliary topography by using the topography in 
    depressions to derive the most realistic flow paths. It uses as weight sweights the 
    differences between the filled and the raw DEM.
    
    References:
    -----------
    This algoritm is adapted to Python from TopoToolbox matlab codes by Wolfgang Schwanghart 
    (version 2017-09-02). It is equivalent to use the preprocessing option
    'carve' in that code (default preprocessing option).
    
    Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
    for topographic analysis. Environ. Model. Softw. 25, 770–781. 
    https://doi.org/10.1016/j.envsoft.2009.12.002
    
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
    """
    
    # Fill the DEM *
    fill = dem.fill_sinks()
    diff = fill.read_array() - dem.read_array() 
    dem = fill
        
    # Identify flats and sills *
    flats, sills = dem.identify_flats(nodata=False)
    
    # Derive the cost of routing through sills *
    carvemin = 0.1                       
    struct = np.ones((3, 3), dtype="int8")
    
    lbl_arr, nlbl = ndimage.label(flats.read_array(), structure=struct)
    lbls = np.arange(1, nlbl + 1)
    
    tweight = 2
    
    for lbl in lbls:
        diff[lbl_arr == lbl] = (diff[lbl_arr==lbl].max() - diff[lbl_arr == lbl])**tweight + carvemin
        
    del lbl_arr
    
    
    # Get presill pixels i.e. pixels immediately upstream to sill pixels *
    zvals = dem.read_array()
    flats_arr = flats.read_array().astype("bool")
    dims = zvals.shape
    row, col = sills.find()
    
    rowadd = np.array([-1, -1, 0, 1, 1,  1,  0, -1])
    coladd = np.array([ 0,  1, 1, 1, 0, -1, -1, -1])
    presill_rows = np.array([], dtype="int")
    presill_cols = np.array([], dtype="int")
    
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
        presill_rows = np.append(presill_rows, rowp[valid_pix])
        presill_cols = np.append(presill_cols, colp[valid_pix])
    
    starts = [xx for xx in zip(presill_rows, presill_cols)]
    
    
    # Calulate auxiliary topography, ie. the cost surface seeded at presill pixels *
    flats_arr = np.invert(flats_arr)
    diff[flats_arr] = 99999
    lg = graph.MCP_Geometric(diff)
    diff = lg.find_costs(starts=starts)[0] + 1
    diff[flats_arr] = -99999
    
    return diff

def sort_pixels(dem, aux_topo):
    dem = dem.fill_sinks().read_array()
    rdem = dem.ravel()
#    if aux_topo:
    diff = aux_topo.ravel()
    ix_sorted_flats = np.argsort(-diff, kind='mergesort')
    
#    ### TEST
#    ix_sorted_flats = sio.loadmat("IXSortedFlats.mat")["IXSortedFlats"] - 1
#    row, col = np.unravel_index(ix_sorted_flats, dims, "F")
#    ix_sorted_flats = np.ravel_multi_index((row, col), dims)
#    ix_sorted_flats = ix_sorted_flats.reshape((ncells,))
#    ### TEST
    
#    ix_sorted_flats_z = np.append(rdem[ix_sorted_flats].reshape(ncells, 1),
#                              ix_sorted_flats.reshape(ncells, 1), axis=1).reshape((ncells, 2))
#    
#    ind = np.argsort(ix_sorted_flats_z[:,0])[::-1]
#    ix = np.array(ix_sorted_flats_z[:,1][ind], dtype=np.int32)
    ndx = np.arange(ncells)
    ndx = ndx[ix_sorted_flats]
    
    # This F#@#@ line took me three days
    # Numpy has not 'descent', but you CANNOT invert the index array! 
    # Since you'll mess up all the previous sorting!!!
    #ix = np.argsort(rdem[ndx])[::-1] # My F@#@¢@ mistake!!
    ix = np.argsort(-rdem[ndx], kind='mergesort')
    ix = ndx[ix]
  
    return ix

aux_topo = _get_aux_topography(dem)
ix = sort_pixels(dem, aux_topo)
#del aux_topo

#ix = sio.loadmat('ix.mat')['ix'] - 1
#ix = ix.reshape((ncells, ))
#row, col = np.unravel_index(ix, dims, "F")
#ix = np.ravel_multi_index((row, col), dims)

pp = np.zeros(dims, dtype=np.int32)
IX = np.arange(ncells, dtype=np.int32)

pp = pp.ravel()
pp[ix] = IX
pp = pp.reshape(dims)

demarr = dem.fill_sinks().read_array()
f_dem = demarr.ravel()

# CARDINAL NEIGHBORS
footprint= np.array([[0, 1, 0],
                     [1, 1, 1],
                     [0, 1, 0]], dtype=np.int)
IXC1 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
xxx1 = np.copy(IXC1)
IX = IXC1.ravel()[ix]
IXC1 = ix[IX]
G1   = (f_dem[ix]-f_dem[IXC1])/(dem.get_cellsize())

"""
% cardinal neighbors
IXC1 = imdilate(pp,[0 1 0; 1 1 1; 0 1 0]>0);
xxx1 = IXC1;
IX   = IXC1(FD.ix);
IXC1 = FD.ix(IX);
G1   = (dem(FD.ix)-dem(IXC1))/(FD.cellsize);
G1(FD.ix == IXC1) = -inf;
"""

# DIAGONAL NEIGHBORS
footprint= np.array([[1, 0, 1],
                     [0, 1, 0],
                     [1, 0, 1]], dtype=np.int)
IXC2 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
xxx2 = np.copy(IXC2)
IX = IXC2.ravel()[ix]
IXC2 = ix[IX]
G2   = (f_dem[ix]-f_dem[IXC2])/(dem.get_cellsize() * np.sqrt(2))

"""
IXC2 = imdilate(pp,[1 0 1; 0 1 0; 1 0 1]>0);
xxx2 = IXC2;
IX   = IXC2(FD.ix);
IXC2 = FD.ix(IX);
G2   = (dem(FD.ix)-dem(IXC2))/(norm([FD.cellsize,FD.cellsize]));
"""

I  = (G1<=G2) & (xxx2.ravel()[ix]>xxx1.ravel()[ix])
ixc = IXC1
ixc[I] = IXC2[I]
I = ixc == ix
I = np.invert(I)
ix = ix[I]
ixc =ixc[I]


"""
% choose the steeper one

I    = G1<=G2 & xxx2(FD.ix)>xxx1(FD.ix);
FD.ixc = IXC1;
FD.ixc(I) = IXC2(I);

I = FD.ixc == FD.ix;
FD.ix(I) = [];
FD.ixc(I) = [];
"""
#ix = sio.loadmat('ix.mat')['ix'] - 1
#ix = ix.reshape((304980, ))
#row, col = np.unravel_index(ix, dims, "F")
#ix = np.ravel_multi_index((row, col), dims)

#ixc = sio.loadmat('ixc.mat')['ixc'] - 1
#ixc = ixc.reshape((304980, ))
#row, col = np.unravel_index(ixc, dims, "F")
#ixc = np.ravel_multi_index((row, col), dims)

nix = len(ix)

A = np.ones(ncells)
for n in range(nix):
    A[ixc[n]] = A[ix[n]] + A[ixc[n]]

A = A.reshape(dims)

acc = DEM()
acc.copy_layout(dem)
acc.set_array(A)

acc.save("flow_Arr.tif")


##del aux_topo
#
#aux_topo = _get_aux_topography(dem)
#
#
#dem = dem.fill_sinks().read_array()
#f_dem = dem.ravel()
#
#IXSortedFlats = np.argsort(aux_topo.ravel())[::-1]
##del aux_topo
#
#ndx = np.arange(ncells, dtype=np.int32)
#ndx = ndx[IXSortedFlats]
##del IXSortedFlats
#
#ix = np.argsort(f_dem[ndx])[::-1]
#ix = ndx[ix]
#
#resta = ix2 - ix


#ix = np.argsort(rdem[ix_sorted_flats])[::-1]
#
#sdem = dem.ravel()[ix]


#fdem = np.ravel(mdem, "F")
#ix = np.argsort(fdem[ix_sorted_flats])[::-1]

## Sort pixels (givers)
#aux_topo = diff.ravel()
#dem = dem.fill_sinks()
#
#f_dem = dem.read_array().ravel("F")

#ix = np.argsort(f_dem[ix_sorted_flats])[::-1]

#ndx = np.arange(ncells, dtype = np.int32)[ix_sorted_flats]

#ix = ndx[np.argsort(f_dem[ndx])] ## FD.IX (Givers)


#
#r, c = np.unravel_index(ix, dims)
#ind = np.ravel_multi_index((r, c), dims, order = "F")
#
##
##ncells = dims[0] * dims[1]
##
#pp = np.zeros(dims, dtype=np.int32)
#iixx = np.arange(dims[0] * dims[1])
#
#pp[np.unravel_index(ix, dims)] = iixx
##pp = pp.reshape(dims)






                    


