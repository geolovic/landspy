#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:55:46 2018

@author: vicen
"""
import numpy as np
from scipy import ndimage
cimport numpy
cimport scipy

def sort_pixels(dem, cellsize):
    # 01 Fill sinks
    fill = dem.fill_sinks()
    topodiff = fill.read_array() -  dem.read_array()
    topodiff = topodiff.astype(np.float32)
    dem = fill
    
    # 02 Get flats and sills
    flats, sills = dem.identify_flats(False)
    
    # 03 Get presills (i.e. pixels immediately upstream to sill pixels)
    presill_pos = self._get_presills(flats, sills, dem)
    
    del sills

    # 04 Get the auxiliar topography for the flats areas
    topodiff = self._get_topodiff(topodiff, flats)

    # 05 Get the weights inside the flat areas (for the cost-distance analysis)
    weights = self._get_weights(flats, topodiff, presill_pos)
    
    del flats, topodiff, presill_pos
    self.weights = weights
    # 06 Sort pixels (givers)
    ix = self._sort_pixels(dem, weights)
    
    # 07 Get receivers
    ixc = self._get_receivers(ix, dem)

    I = ixc == ix
    I = np.invert(I)
    self._ix = ix[I]
    self._ixc = ixc[I]

    
def _get_presills(self, flats, sills, dem):
    """
    This functions extracts the presill pixel locations (i.e. pixelimmediately 
    upstream to sill pixels)- Adapted from TopoToolbox matlab codes.
    """
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
    
    return list(ps_pos)


def _get_topodiff(self, topodiff, flats):
    """
    This function calculate an auxiliar topography to sort the flats areas
    """
    tweight = 2
    carvemin = 0.1                       
    struct = np.ones((3, 3), dtype="int8")
    lbl_arr, nlbl = ndimage.label(flats.read_array(), structure=struct)
    lbls = np.arange(1, nlbl + 1)
    
    for lbl in lbls:
        topodiff[lbl_arr == lbl] = (topodiff[lbl_arr==lbl].max() - topodiff[lbl_arr == lbl])**tweight + carvemin
                
    return topodiff

def _get_weights(self, flats, topodiff, ps_pos):
    """
    This function calculate weights in the flats areas by doing a cost-distance analysis.
    It uses presill positions as seed locations, and an auxiliar topography as friction
    surface.
    """
    flats_arr = flats.read_array().astype("bool")
    flats_arr = np.invert(flats_arr)
    topodiff[flats_arr] = 99999
    lg = graph.MCP_Geometric(topodiff)
    topodiff = lg.find_costs(starts=ps_pos)[0] + 1
    topodiff[flats_arr] = -99999
    
    return topodiff

def _sort_pixels(self, dem, weights):
    """
    Sort the cells of a DEM in descending order. It uses weights for the flats areas
    """
    # Sort the flat areas
    rdem = dem.read_array().ravel()
    rweights = weights.ravel()
    ix_flats = np.argsort(-rweights, kind='mergesort')
    
    # Sort the rest of the pixels from the DEM
    ndx = np.arange(self._ncells, dtype=np.int)
    ndx = ndx[ix_flats]
    ix = ndx[np.argsort(-rdem[ndx], kind='mergesort')]
    
    return ix

def _get_receivers(self, ix, dem):
    
    dem = dem.read_array()
    rdem = dem.ravel()
    
    pp = np.zeros(self._dims, dtype=np.int32)
    IX = np.arange(self._ncells, dtype=np.int32)
    pp = pp.ravel()
    pp[ix] = IX
    pp = pp.reshape(self._dims)
            
    # Get cardinal neighbors
    footprint= np.array([[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]], dtype=np.int)
    IXC1 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
    xxx1 = np.copy(IXC1)
    IX = IXC1.ravel()[ix]
    IXC1 = ix[IX]
    G1   = (rdem[ix]-rdem[IXC1])/(self._cellsize)
    
    # Get diagonal neighbors
    footprint= np.array([[1, 0, 1],
                         [0, 1, 0],
                         [1, 0, 1]], dtype=np.int)
    IXC2 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
    xxx2 = np.copy(IXC2)
    IX = IXC2.ravel()[ix]
    IXC2 = ix[IX]
    G2   = (rdem[ix]-rdem[IXC2])/(self._cellsize * np.sqrt(2))
    
    del IX
    
    # Get the steepest one
    I  = (G1<=G2) & (xxx2.ravel()[ix]>xxx1.ravel()[ix])
    ixc = IXC1
    ixc[I] = IXC2[I]
    
    return ixc

def fill_sinks(dem):
    pass
