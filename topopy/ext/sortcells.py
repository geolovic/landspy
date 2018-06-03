#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:48:28 2018

@author: vicen
"""
import numpy as np
from scipy import ndimage
from skimage import graph

def sort_pixels(dem, auxtopo=False, filled=False, verbose=False, order="C"):
    
    # Get DEM properties
    cellsize = dem.get_cellsize()
    nodata_val = dem.get_nodata()
    if nodata_val is None:
        nodata_val = -9999
    if filled:
        auxtopo=False
        
    # 01 Fill sinks
    # If filled is True, input DEM was previously pit-filled
    if filled:
        fill = dem
        dem_arr = fill.read_array()
        topodiff = np.zeros(dem_arr.shape, dem_arr.dtype)
        del(dem)
    else:
        fill = dem.fill_sinks()
        dem_arr = dem.read_array()
        fill_arr = fill.read_array()
        topodiff = fill_arr - dem_arr
        dem_arr = fill_arr
        del(dem)     
    if verbose:
        print("1/7 - DEM filled")
    
    # 02 Get flats and sills
    flats, sills = fill.identify_flats(as_array=True)
    del(fill)
    if verbose:
        print("2/7 - Flats and sills identified")
    
    # 03 Get presills (i.e. pixels immediately upstream to sill pixels)
    presills_pos = get_presills(dem_arr, flats, sills)
    if verbose:
        print("3/7 - Presills identified")
    
    # 04 Get the auxiliar topography for the flats areas
    if auxtopo:
        topodiff = get_aux_topography(topodiff.astype(np.float32), flats.astype(np.int8))
    else:
        topodiff = np.zeros(dem_arr.shape, dtype=np.int8)
        topodiff[flats] = 1
        
    if verbose:
        print("4/7 - Auxiliar topography generated")
    
    # 05 Get the weights inside the flat areas (for the cost-distance analysis)
    weights = get_weights(flats, topodiff, presills_pos)
    if verbose:
        print("5/7 - Weights calculated")
    
    del flats, sills, presills_pos, topodiff

    # 06 Sort pixels (givers)
    ix = sort_dem(dem_arr, weights)
    if verbose:
        print("6/7 - Pixels sorted")
    
    # 07 Get receivers
    ixc = get_receivers(ix, dem_arr, cellsize)
    if verbose:
        print("7/7 - Receivers calculated")

    # 08 Remove givers==receivers
    ind = np.invert(ixc == ix) # givers == receivers
    ix = ix[ind]
    ixc = ixc[ind]
    
    # 09 Remove receivers marked as nodatas
    w = dem_arr != nodata_val
    w = w.ravel()
    I   = w[ixc]
    ix  = ix[I]
    ixc = ixc[I]
    
    return ix, ixc


def get_presills(filldem, flats, sills, as_positions=True):
    """
    This functions extracts the presill pixel locations (i.e. pixel immediately 
    upstream to sill pixels)- Adapted from TopoToolbox matlab codes.
    
    Parameters:
    -----------
    filldem : *np.ndarray*
      Array of values representing a filled DEM
    flats : *np.ndarray*
      Numpy logical array with location of the flats (cells without downward 
      neighboring cells)
    sills: *np.ndarray*
      Numpy logical array with location of the sills (cells where flat regions 
      spill over into lower terrain)
    
    Return:
    -------
    presills_pos : list
      List of tuples (row, col) with the location of the presill pixels
      
    References:
    -----------
    This algoritm is adapted from TopoToolbox matlab codes by Wolfgang Schwanghart 
    
    Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
    for topographic analysis. Environ. Model. Softw. 25, 770–781. 
    https://doi.org/10.1016/j.envsoft.2009.12.002
    
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014  
    """
    dims = filldem.shape
    row, col = np.where(sills)
    
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
        cond01 = filldem[row[valid_rc], col[valid_rc]] == filldem[rowp, colp]
        cond02 = flats[rowp, colp]
        valid_pix = np.logical_and(cond01, cond02)
        ps_rows = np.append(ps_rows, rowp[valid_pix])
        ps_cols = np.append(ps_cols, colp[valid_pix])
    
    if as_positions:
        ps_pos = list(zip(ps_rows, ps_cols))
        return ps_pos
    else:
        presills = np.zeros(dims, dtype=np.bool)
        presills[ps_rows, ps_cols] = True
        return presills

    
def get_aux_topography(topodiff, flats):
    """
    This function calculate an auxiliar topography to sort the flats areas

    Parameters:
    -----------
    flats : *numpy.array* [dtype=int8]
      Numpy array [np.int8 type] with the location of the flats
    topodiff : *numpy.array* [dtype=float32]
      Numpy array [np.float32 type] with the auxiliar topography (diference between filldem and dem)
    
    Return:
    -------
    aux_topography : *np.array*
      Auxiliar topography to sort flats areas
    """              
    struct = np.ones((3, 3), dtype=np.int8)
    lbl_arr, nlbl = ndimage.label(flats, structure=struct)
    lbls = np.arange(1, nlbl + 1)
    #aux_topo = np.copy(topodiff) # Unmark to preserve topodiff (uses more memory)
    aux_topo = topodiff
    
    for lbl in lbls:
        aux_topo[lbl_arr == lbl] = (aux_topo[lbl_arr==lbl].max() - aux_topo[lbl_arr == lbl])**2 + 0.1
                
    return aux_topo


def get_weights(flats, aux_topo, presills_pos):
    """
    This function calculate weights in the flats areas by doing a cost-distance analysis.
    It uses presill positions as seed locations, and an auxiliar topography as friction
    surface.

    Parameters:
    -----------
    flats : *numpy.array* [dtype = np.bool]
      Numpy array with the location of the flats surfaces
    aux_topo: *numpy.array* [dtype = np.float32]
      Numpy array with the auxiliar topography
    presill_pos *list*
      List of tuples (row, col) with the location of the presills

    Returns:
    --------
    weigths : *numpy.array*
      Numpy array with the cost of routing throught the flat areas
    
    References:
    -----------
    This algoritm is adapted from the TopoToolbox matlab codes by Wolfgang Schwanghart.
    
    Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
    for topographic analysis. Environ. Model. Softw. 25, 770–781. 
    https://doi.org/10.1016/j.envsoft.2009.12.002
    
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
    """
    flats = np.invert(flats)
    aux_topo[flats] = 99999
    if len(presills_pos) > 0:
        lg = graph.MCP_Geometric(aux_topo)
        aux_topo = lg.find_costs(starts=presills_pos)[0] + 1
    aux_topo[flats] = -99999
    
    return aux_topo


def sort_dem(dem_arr, weights, order="C"):
    """
    Sort the cells of a DEM in descending order. It uses a weights array to
    sort the flats areas
    
    Parameters:
    -----------
    dem_arr : *numpy.ndarray* 
      Numpy array representing a filled DEM
    weights :  *numpy.ndarray* 
      Numpy array with the weights to sort the flats areas
    order : *str*
      Order of the returned indexes ("C" row-major (C-style) or "F", column-major 
      (Fortran-style) order.
    """
    ncells = dem_arr.shape[0] * dem_arr.shape[1]
    # Sort the flat areas
    rdem = dem_arr.ravel(order=order)
    rweights = weights.ravel(order=order)
    ix_flats = np.argsort(-rweights, kind='mergesort')
    
    # Sort the rest of the pixels from the DEM
    ndx = np.arange(ncells, dtype=np.int)
    ndx = ndx[ix_flats]
    ix = ndx[np.argsort(-rdem[ndx], kind='mergesort')]
    
    return ix


def get_receivers(ix, dem_arr, cellsize, order="C"):
    """
    This function obtain the receiver cells for an array of "givers" cells 
    represented by linear indexes.
    
    Parameters:
    -----------
    ix : *numpy.array*
      Linear indexes for the givers cells
    dem_arr : *numpy.array*
      Numpy array representing the DEM to obtain receivers
    cellsize : *float* / *int*
      Cellsize of the DEM
    order : *str*
      Order of the returned indexes ("C" row-major (C-style) or "F", column-major 
      (Fortran-style) order
      
    Returns:
    --------
    ixc : *numpy.array*
      Linear indexes for the receivers cells
      
    References:
    -----------
    This algoritm is adapted from the TopoToolbox matlab codes by Wolfgang Schwanghart.
    
    Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
    for topographic analysis. Environ. Model. Softw. 25, 770–781. 
    https://doi.org/10.1016/j.envsoft.2009.12.002
    
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014    
    """
    
    ncells = dem_arr.shape[0] * dem_arr.shape[1]
    dims = dem_arr.shape
    rdem = dem_arr.ravel(order=order)
    
    pp = np.zeros(dims, dtype=np.int32)
    IX = np.arange(ncells, dtype=np.int32)
    pp = pp.ravel(order=order)
    pp[ix] = IX
    pp = pp.reshape(dims, order=order)
            
    # Get cardinal neighbors
    footprint= np.array([[0, 1, 0],
                         [1, 1, 1],
                         [0, 1, 0]], dtype=np.int)
    IXC1 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
    xxx1 = np.copy(IXC1)
    IX = IXC1.ravel(order=order)[ix]
    IXC1 = ix[IX]
    G1   = (rdem[ix]-rdem[IXC1])/(cellsize)
    
    # Get diagonal neighbors
    footprint= np.array([[1, 0, 1],
                         [0, 1, 0],
                         [1, 0, 1]], dtype=np.int)
    IXC2 = ndimage.morphology.grey_dilation(pp, footprint=footprint)
    xxx2 = np.copy(IXC2)
    IX = IXC2.ravel(order=order)[ix]
    IXC2 = ix[IX]
    G2   = (rdem[ix]-rdem[IXC2])/(cellsize * np.sqrt(2))
    
    # Get the steepest one
    I  = (G1<=G2) & (xxx2.ravel(order=order)[ix]>xxx1.ravel(order=order)[ix])
    ixc = IXC1
    ixc[I] = IXC2[I]
    
    return ixc