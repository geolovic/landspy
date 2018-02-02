#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 17:14:39 2018
Temp testing script for Network
@author: vicen
"""
import sys
sys.path.append("../../")
from topopy import Flow, DEM
from scipy.sparse import csc_matrix
import numpy as np
import matplotlib.pyplot as plt


def get_stream_poi(fd, threshold, kind="heads"):
    """
    This function find point of interest of the drainage network. These points of interest
    can be 'heads', 'confluences' or 'outlets'.
    
    Parameters:
    ===========
    threshold : *int* 
      Area threshold to initiate a channels (in cells)
    kind : *str* {'heads', 'confluences', 'outlets'}
      Kind of point of interest to return. 
      
    Returns:
    ==========
    (row, col) : *tuple*
      Tuple of numpy nd arrays with the location of the points of interest
      
    References:
    -----------
    The algoritms to extract the point of interest have been adapted to Python 
    from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
    August, 2017). These smart algoritms use sparse arrays with giver-receiver indexes, to 
    derive point of interest in a really efficient way. Cite:
            
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
    """
    # Check input parameters
    if kind not in ['heads', 'confluences', 'outlets']:
        return
    
    # Get the Flow Accumulation and select cells that meet the threholds
    fac = fd.get_flow_accumulation(nodata = False, asgrid=False)
    w = fac > threshold
    del fac
    
    # Build sparse array
    w = w.ravel()
    I   = w[fd._ix]
    ix  = fd._ix[I]
    ixc = fd._ixc[I]
    aux_vals = np.ones(ix.shape, dtype=np.int8)

    sp_arr = csc_matrix((aux_vals, (ix, ixc)), shape=(fd._ncells, fd._ncells))
    del I, ix, ixc, aux_vals # Clean up (Don't want to leave stuff in memory)
    
    if kind == 'heads':
        sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
        out_pos = (sum_arr == 0) & w
    elif kind == 'confluences':
        sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
        out_pos = sum_arr > 1
    elif kind == 'outlets':
        sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
        out_pos = (sum_arr == 0) & w  
        
    out_pos = out_pos.reshape(fd._dims)
    row, col = np.where(out_pos)
    
    return row, col

dem = DEM("../data/tunez.tif")
fd = Flow()
fd.load_gtiff("../data/tunez_fd.tif")
fig, ax = plt.subplots()
dem.plot(ax)

streams = fd.get_network(1000)
streams.plot(ax)

# Heads
row, col = get_stream_poi(fd, 1000, 'heads')
ax.plot(col, row, "ro")
#
# Confluences
row, col = get_stream_poi(fd, 1000, 'confluences')
ax.plot(col, row, "bs")

