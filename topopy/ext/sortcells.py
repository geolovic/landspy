#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 16:48:28 2018

@author: vicen
"""
import numpy as np
from scipy import ndimage
from skimage import graph

def sort_pixels(dem, order="C"):

    cellsize = dem.get_cellsize()
    nodata_val = dem.get_nodata()
    if nodata_val is None:
        nodata_val = -9999

    # 01 Fill sinks
    dem = dem.read_array()
    fill = fill_sinks(dem, nodata_val)
    topodiff = fill - dem
    dem = fill
    
    # 02 Get flats and sills
    flats, sills = identify_flats(dem, nodata_val)
    
    # 03 Get presills (i.e. pixels immediately upstream to sill pixels)
    presills_pos = get_presills(dem, flats, sills)
    
    # 04 Get the auxiliar topography for the flats areas
    topodiff = get_aux_topography(topodiff.astype(np.float32), flats.astype(np.int8))
    
    # 05 Get the weights inside the flat areas (for the cost-distance analysis)
    weights = get_weights(flats, topodiff, presills_pos)
    
    del flats, sills, presills_pos, topodiff

    # 06 Sort pixels (givers)
    ix = sort_dem(dem, weights)
    
    # 07 Get receivers
    ixc = get_receivers(ix, dem, cellsize)

    # 08 Remove nodata cells
    ind = ixc == ix
    ind = np.invert(ind)
    ix = ix[ind]
    ixc = ixc[ind]
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


def fill_sinks(input_array, nodata_val):
    """
    Fill sinks method adapted from  fill depressions/sinks in floating point array
    
    Parameters:
    ----------
    input_array : *numpy.ndarray* 
      Input array representing the DEM to be filled
    nodata_val : *int* or *float*
      Value for nodata (cells with values of nodata will be maintained)
    
    Returns:
    ----------
    fill_dem : *numpy.ndarray* 
      Array with depressions/sinks filled

    This algorithm has been adapted (with minor modifications) from the 
    Charles Morton slow fill algorithm (with ndimage and python 3 was not slow
    at all). 
    
    References
    ----------
    Soile, P., Vogt, J., and Colombo, R., 2003. Carving and Adaptive
    Drainage Enforcement of Grid Digital Elevation Models.
    Water Resources Research, 39(12), 1366
    
    Soille, P., 1999. Morphological Image Analysis: Principles and
    Applications, Springer-Verlag, pp. 173-174

    """    
    # Change nodata values for very low values
    copyarr = np.copy(input_array)
    nodata_pos = np.where(input_array==nodata_val)
    copyarr[nodata_pos] = -9999.
    
    # Set h_max to a value larger than the array maximum to ensure
    #   that the while loop will terminate
    h_max = copyarr.max() + 100

    # Build mask of cells with data not on the edge of the image
    # Use 3x3 square Structuring element
    inside_mask = ndimage.morphology.binary_erosion(
        np.isfinite(copyarr),
        structure=np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool))

    # Initialize output array as max value test_array except edges
    output_array = np.copy(copyarr)
    output_array[inside_mask] = h_max

    # Array for storing previous iteration
    output_old_array = np.copy(copyarr)
    output_old_array[:] = 0

    # Cross structuring element
    el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)

    # Iterate until marker array doesn't change
    while not np.array_equal(output_old_array, output_array):
        output_old_array = np.copy(output_array)
        output_array = np.maximum(
            copyarr,
            ndimage.grey_erosion(output_array, size=(3, 3), footprint=el))

    # Put back nodata values and return
    output_array[nodata_pos] = nodata_val
    
    return output_array


def identify_flats(fill_array, nodata_val):
    """
    This functions returns two binary Grid objects (values 0, 1) with flats and sills. 
    Flats are defined  as cells without downward neighboring cells. Sills are cells where 
    flat regions spill over into lower terrain. It the DEM has nodata, they will be maintained
    in output grids.
    
    Parameters:
    -----------
    fill_arr : *numpy.ndarray*
      Array of values representing a filled DEM
    
    Returns:
    ----------
    [flats, sills] : *numpy.ndarray*
      Logical arrays with flats and sill location as True
    
    References:
    -----------
    This algoritm is adapted from identifyflats.m by Wolfgang Schwanghart 
    (version of 17. August, 2017) included in TopoToolbox matlab codes.
    
    Schwanghart, W., Kuhn, N.J., 2010. TopoToolbox: A set of Matlab functions 
    for topographic analysis. Environ. Model. Softw. 25, 770–781. 
    https://doi.org/10.1016/j.envsoft.2009.12.002
    
    Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
    MATLAB-based software for topographic analysis and modeling in Earth 
    surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014

    """
    
    # Change nodata to very low values
    nodata_ids = np.where(fill_array==nodata_val)
    fill_array[nodata_ids] = -9999
    
    footprint = np.ones((3, 3), dtype=np.int8)
    
    # Identify flats throught a image binary erosion
    # Flats will be True where cells don't have lower neighbors
    flats = ndimage.morphology.grey_erosion(fill_array, footprint=footprint) == fill_array
    
    # Remove flats from the borders
    flats[:,0] = False
    flats[:,-1] = False
    flats[0,:] = False
    flats[-1,:] = False
    
    # Remove flats for nodata values and cells bordering them
    flats[nodata_ids] = False
    auxmat = np.zeros(flats.shape, dtype="bool")
    auxmat[nodata_ids] = True
    nodata_bord = ndimage.morphology.grey_dilation(auxmat, footprint=footprint)
    flats[nodata_bord] = False
    
    # Identify sills
    sills = np.empty(fill_array.shape)
    sills.fill(-9999.)
    sills[flats] = fill_array[flats]
    aux_dil = ndimage.morphology.grey_dilation(sills, footprint=footprint)
    sills = np.logical_and(aux_dil == fill_array, np.logical_not(flats))
    sills[nodata_ids] = False
    
    # Return
    return flats, sills
    
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