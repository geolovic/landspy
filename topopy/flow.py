# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# December 26, 2017
#
# Last modified 31 May, 2018

import numpy as np
import gdal
from scipy import ndimage
from skimage import graph
from . import Grid, PRaster

class Flow(PRaster):
    
    def __init__(self, dem="", auxtopo=False, filled=False, verbose=False):
        """
        Class that define a network object (topologically sorted giver-receiver cells)
        
        Parameters:
        ===========
        dem : *DEM object* or *str*
          topopy.DEM instance with the input Digital Elevation Model, or path to a previously saved Flow object. If the 
          parameter is an empty string, it will create an empty Flow instance.
        auxtopo : boolean
          Boolean to determine if a auxiliar topography is used (much slower). The auxiliar
          topography is calculated with elevation differences between filled and un-filled dem.
        filled : boolean
          Boolean to check if input DEM was already pit-filled. The fill algoritm implemented in the DEM object, 
          althoug fast, consumes a lot of memory. In some cases could be necessary fill the DEM with alternative GIS tools.
        verbose : boolean
          Boolean to show processing messages in console to known the progress. Usefull with large DEMs to se the evolution.
        
        References:
        -----------
        The algoritm to created the topologically sorted network has been adapted to Python from FLOWobj.m 
        by Wolfgang Schwanghart (version of 17. August, 2017) included in TopoToolbox matlab codes (really
        smart algoritms there!). If use, please cite:
                
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        """
        
        if dem == "":
            # Creates an empty Flow object
            self._size = (1, 1)
            self._dims = (1, 1)
            self._geot = (0., 1., 0., 0., 0., -1.)
            self._cellsize = self._geot[1]
            self._proj = ""
            self._ncells = 1
            self._nodata_pos = np.array([], dtype=np.int32)
            self._ix = np.array([], dtype=np.int32)
            self._ixc = np.array([], dtype=np.int32)
        
        elif type(dem) == str:
            # Loads the Flow object in GeoTiff format
            try:
                self.load_gtiff(dem)
            except:
                raise FlowError("Error opening the Geotiff")
        else:
            #try:
            # Set Network properties
            self._size = dem.get_size()
            self._dims = dem.get_dims()
            self._geot = dem.get_geotransform()
            self._cellsize = dem.get_cellsize()
            self._proj = dem.get_projection()
            self._ncells = dem.get_ncells()
            self._nodata_pos = np.ravel_multi_index(dem.get_nodata_pos(), self._dims)            
            # Get topologically sorted nodes (ix - givers, ixc - receivers)
            self._ix, self._ixc = sort_pixels(dem, auxtopo=auxtopo, filled=filled, verbose=verbose)
#            except:
#                raise FlowError("Unexpected Error creating the Flow object")
    
    def save_gtiff(self, path):
        """
        Saves the flow object as a geotiff. The geotiff file it wont have any
        sense if its open with GIS software.
        
        Parameters:
        ===========
        path : *str* 
          Path to store the geotiff file with the Flow data

        The organization of this geotiff is as follow::
            
        * Band 1 --> Givers pixels reshaped to self._dims
        * Band 2 --> Receiver pixels reshaped to self._dims
        * Band 3 --> Nodata band (pixels with nodata == 1, pixels with data == 0)
        """
        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._dims[1], self._dims[0], 3, 4)
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)

        no_cells = self._ncells - len(self._ix)

        miss_cells = np.zeros(no_cells, np.uint32)
        ix = np.append(self._ix, miss_cells)
        ixc = np.append(self._ixc, miss_cells)
        nodata_arr = np.zeros(self._dims, np.uint32)
        nodata_pos = np.unravel_index(self._nodata_pos, self._dims)
        nodata_arr[nodata_pos] = 1

        ix = ix.reshape(self._dims)
        ixc = ixc.reshape(self._dims)

        raster.GetRasterBand(1).WriteArray(ix)
        raster.GetRasterBand(2).WriteArray(ixc)
        raster.GetRasterBand(3).WriteArray(nodata_arr)
        raster.GetRasterBand(3).SetNoDataValue(no_cells)

    def load_gtiff(self, path):
        """
        Load a geotiff file with flow direction information. This geotiff must
        have been saved with the save_gtiff() function.
        
        Parameters:
        ===========
        path : *str* 
          Path for the Flow geotiff.
        """
        raster = gdal.Open(path)

        # Set Network properties        
        self._size = (raster.RasterXSize, raster.RasterYSize)
        self._dims = (raster.RasterYSize, raster.RasterXSize)
        self._geot = raster.GetGeoTransform()
        self._cellsize = self._geot[1]
        self._proj = raster.GetProjection()
        self._ncells = raster.RasterYSize * raster.RasterXSize

        # Load nodata values
        banda = raster.GetRasterBand(3)
        no_cells = banda.GetNoDataValue()
        arr = banda.ReadAsArray()
        self._nodata_pos = np.ravel_multi_index(np.where(arr==1), self._dims)

        # Load ix, ixc
        banda = raster.GetRasterBand(1)
        arr = banda.ReadAsArray()
        self._ix = arr.ravel()[0:int(self._ncells - no_cells)]
        banda = raster.GetRasterBand(2)
        arr = banda.ReadAsArray()
        self._ixc = arr.ravel()[0:int(self._ncells - no_cells)]
    
    def get_flow_accumulation(self, weights=None, nodata=True, asgrid=True):
        """
        Calculates the flow accumulation from the topologically sorted pixels of the
        Flow object. As pixels of the Flow objects are sorted topologically, the flow
        accumulation can be obtained very fast with a computational time that is linearly
        dependent on the number of cell of the DEM.
        
        Parameters:
        ===========  
        weights : *topopy.Grid*
          Grid with weights for the flow accumulation (p.e. precipitation values)
        nodata : *bool*
          Boolean flag that indicates if the output flow accumulation grid will maintain NoData values. 
          If nodata=False, nodata values will be filled with 0 and NoDataValue will set to None. 
        asgrid : *bool*
          Indicates if the network is returned as topopy.Grid (True) or as a numpy.array
        
        Usage:
        ======
        flowacc = fd.get_flow_accumulation() # Create a flow accumulation Grid object
        flowacc.save("C:/Temp/flow_acc.tif") # Save the flow accumulation in the disk
        
        Reference:
        ----------
        Braun, J., Willett, S.D., 2013. A very efficient O(n), implicit and parallel 
        method to solve the stream power equation governing fluvial incision and landscape 
        evolution. Geomorphology 180–181, 170–179. 
        """
        if weights:
            facc = weights.read_array()
        else:
            facc = np.ones(self._ncells, np.uint32)
        
        nix = len(self._ix)
        for n in range(nix):
            facc[self._ixc[n]] += facc[self._ix[n]]
        
        facc = facc.reshape(self._dims)
        if nodata:
            nodata_val = np.iinfo(np.uint32).max
        else:
            nodata_val = 0
        
        row, col = np.unravel_index(self._nodata_pos, self._dims)
        facc[row, col] = nodata_val
        
        # Get the output in form of a Grid object
        if not nodata:
            nodata_val = None
        if asgrid:
            return self._create_output_grid(facc, nodata_val)
        else:
            return facc

    def get_drainage_basins(self, outlets=[], asgrid=True):
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
        basins = fd.drainage_basins() # Extract all the basins in the Flow object
        basins = fd.drainage_basins([520359.7, 4054132.2]) # Creates the basin for the specified outlet
        xi = [520359.7, 519853.5]
        yi = [4054132.2, 4054863.5]
        basins = fd.drainage_basins((xi, yi)) # Create two basins according xi and yi coordinates

        References:
        -----------
        The algoritms to extract the drainage basins have been adapted to Python 
        from Topotoolbox matlab codes developed by Wolfgang Schwanghart (version of 17. 
        August, 2017).
                
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        """

        
        # If outlets are not specified, all the basins will be extracted
        if outlets == []:
            temp_ix = self._ix
            temp_ixc = self._ixc
            nbasins = 0
            basin_arr = np.zeros(self._ncells, np.int)
            nix = len(temp_ix)
            for n in range(nix-1,-1,-1):
                # If receiver is zero, add a new basin
                if basin_arr[temp_ixc[n]] == 0:
                    nbasins += 1
                    basin_arr[temp_ixc[n]] = nbasins
                
                # Mark basin giver with the id of the basin receiver
                basin_arr[temp_ix[n]] = basin_arr[temp_ixc[n]]
                      
        # Outlets coordinates are provided
        else:
            temp_ix = self._ix
            temp_ixc = self._ixc
            x, y = outlets
            row, col = self.xy_2_cell(x, y)
            inds = self.cell_2_ind(row, col)
            basin_arr = np.zeros(self._ncells, np.int)

            # Change basin array outlets by the basin id (starting to 1)
            if inds.size == 1:
                basin_arr[inds] = 1
            else:
                for n, inds in enumerate(inds):
                    basin_arr[inds] = n + 1
                
            nix = len(temp_ix)
            # Loop by all the sorted cells
            for n in range(nix-1,-1,-1):
                # If the basin receiver is not Zero and basin giver is zero
                if (basin_arr[temp_ixc[n]] != 0) & (basin_arr[temp_ix[n]] == 0):
                    # Mark giver with the basin id of the receiver
                    basin_arr[temp_ix[n]] = basin_arr[temp_ixc[n]]
        
        # Reshape and return
        basin_arr = basin_arr.reshape(self._dims)  
        
        if asgrid:
            return self._create_output_grid(basin_arr, 0)
        else:
            return basin_arr
            
    def _create_output_grid(self, array, nodata_value=None):
        """
        Convenience function that creates a Grid object from an input array. The array
        must have the same shape that self._dims and will maintain the Flow object 
        properties as dimensions, geotransform, reference system, etc.
        
        Parameters:
        ===========
        array : *numpy.ndarray*
          Array to convert to a Grid object
        nodata_value _ *int* / *float*
          Value for NoData values
          
        Returns:
        ========
        Grid object with the same properties that Flow
        """
        grid = Grid()
        grid.copy_layout(self)
        grid._nodata = nodata_value
        grid._array = array
        grid._tipo = str(array.dtype)
        return grid  

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


class FlowError(Exception):
    pass