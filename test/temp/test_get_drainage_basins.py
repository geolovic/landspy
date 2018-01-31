#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:42:34 2018

@author: vicen
"""
import sys
sys.path.append("../../")
from topopy import Flow, Grid
import numpy as np

fd = Flow()
fd.load_gtiff("../data/tunez_fd.tif")

def get_drainage_basins(fd, outlets=[], asgrid=True):
    """
    This function extracts the drainage basins for the Flow object and returns a Grid object that can
    be saved into the disk.
    
    Parameters:
    ===========
    outlets : *list* or *tuple*
      List or tuple with (xi, yi) coordinate for outlets. xi and xi can be numbers, lists, or numpy.ndarrays
      If outlets is an empty list (default) the basins will be extracted for all the outlets in the Flow object.
    asgrid : *bool*
      Indicates if the network is returned as topopy.Grid (True) or as a numpy.array

    Return:
    =======
    basins : *topopy.Grid* object with the different drainage basins.

    Usage:
    =====
    basins = fd.drainage_basins() # Extract all the basins in the Flow objects
    basins = fd.drainage_basins([520359.7, 4054132.2]) # Creates the basin for the outlet
    xi = [520359.7, 519853.5]
    yi = [4054132.2, 4054863.5]
    basins = fd.drainage_basins((xi, yi)) # Create two basins according xi and yi coordinates

    References:
    -----------
    The algoritms to extract the drainage basins have been adapted to Python 
    from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
    August, 2017). Cite:
            
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
    """
    temp_ix = fd._ix
    temp_ixc = fd._ixc
    
    # If outlets are not specified, basins for all the outlets will be extracted
    if outlets == []:
        nbasins = 0
        D = np.zeros(fd._ncells, np.int)
        nix = len(temp_ix)
        for n in range(nix-1,-1,-1):
            if D[temp_ixc[n]] == 0:
                nbasins += 1
                D[temp_ixc[n]] = nbasins
                outlets.append(temp_ixc[n])
            D[temp_ix[n]] = D[temp_ixc[n]]
                  
    # Outlets coordinates are provided
    else:
        x, y = outlets
        row, col = fd.xy_2_cell(x, y)
        inds = fd.cell_2_ind(row, col)
        D = np.zeros(fd._ncells, np.int)
        for n, inds in enumerate(inds):
            D[inds] = n+1
        
        nix = len(temp_ix)
        for n in range(nix-1,-1,-1):
            if (D[temp_ixc[n]] != 0) & (D[temp_ix[n]] == 0):
                D[temp_ix[n]] = D[temp_ixc[n]]
    
    D = D.reshape(fd._dims)
    if asgrid:
        return fd._create_output_grid(D, 0)
    else:
        return D

x = [520359.726010026, 519853.453259635]
y = [4054132.16717364, 4054863.45003531]

basins = get_drainage_basins(fd, (x, y))
basins.plot()

