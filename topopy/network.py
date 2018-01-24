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
# Last modified January 24, 2017

import numpy as np
import gdal
from scipy import ndimage
from skimage import graph
from scipy.sparse import csc_matrix
from . import Grid

class Flow():
    
    def __init__(self, dem=""):
        """
        Class that define a network object (topologically sorted giver-receiver cells)
        
        Parameters:
        ===========
        dem : *DEM object* 
          DEM object with the input Digital Elevation Model
          
        preprocess : *str* {'carve', 'fill'}
          Type of preprocessing anaysis. These preprocessing algoritms sort pixels topologically 
          in flats areas to ensure flow routing through them.        
        
        References:
        -----------
        The algoritm to created the topologically sorted network has been adapted to Python from FLOWobj.m 
        by Wolfgang Schwanghart (version of 17. August, 2017) included in TopoToolbox matlab codes (really
        smart algoritms there!). If use, please cite:
                
        Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 - 
        MATLAB-based software for topographic analysis and modeling in Earth 
        surface sciences. Earth Surf. Dyn. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
        """
        
        # Creates an empty Network
        if dem == "":
            return
        
        # Set Network properties
        self._dims = (dem._array.shape)
        self._ncells = dem.get_ncells()
        self._cellsize = dem.get_cellsize()
        self._geot = dem.get_geotransform()
        self._proj = dem.get_projection()
        self._cellsize = dem.get_cellsize()
        self._nodata_pos = np.ravel_multi_index(dem.get_nodata_pos(), self._dims)

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

    def xy_2_cell(self, x, y):
        """
        Get row col indexes coordinates from xy coordinates
        
        Parameters:
        ===========
        x : *number*, *list of values*, *array of values*
            Value (or list of values) with x coordinates
        y : *number*, *list of values*, *array of values*
            Value (or list of values) with y coordinates
            
        Return:
        =======
        (row, col) : Tuple with row and column indexes (or list/array of values depending on the input type)
        """
        is_number = False
        is_list = False
        
        if type(x) == int or type(x) == float:
            is_number = True
        
        if type(x) == list:
            is_list = True
            x = np.array(x)
            y = np.array(y)
            
        row = (self._geot[3] - y) / self._geot[1]
        col = (x - self._geot[0]) / self._geot[1]
        
        if is_number:
            return int(row), int(col)
        elif is_list:
            return row.astype("int").tolist(), col.astype("int").tolist()
        else:
            return row.astype("int"), col.astype("int")

    def cell_2_xy(self, row, col):
        """
        Get xy coordinates from cells (row and col values)
        
        Parameters:
        ===========
        row : *number*, *list of values*, *array of values*
            Value (or list of values) with row indexes
        col : *number*, *list of values*, *array of values*
            Value (or list of values) with column indexes
            
        Return:
        =======
        (x, y) : Tuple of coordinates (or list/array of coordinates depending on the input type)
        """
        is_list = False                
        if type(row) == list:
            is_list = True
            row = np.array(row)
            col = np.array(col)
        
        x = self._geot[0] + self._geot[1] * col + self._geot[1] / 2
        y = self._geot[3] - self._geot[1] * row - self._geot[1] / 2
        
        if is_list:
            return x.tolist(), y.tolist()
        else:
            return x, y
        
    def ind_2_cell(self, ind):
        """
        Get row col indexes from cells ids (cells are indexed from left-to-right)
        
        Parameters:
        ===========
        ind : *number*, *list of values*, *array of values*
            Cell index
            
        Return:
        =======
        (row, col) : Tuple with row and column indexes (or list/array of values depending on the input type)
        """
        is_list = False                
        if type(ind) == list:
            is_list = True
        
        row, col = np.unravel_index(ind, self._dims) 
        if is_list:
            return (row.tolist(), col.tolist())
        else:
            return (row, col)
    
    def cell_2_ind(self, row, col):
        """
        Get cell ids from row and column indexes (cells are indexed from left-to-right)
        
        Parameters:
        ===========
        row : *number*, *list of values*, *array of values*
            Value (or list of values) with row indexes
        col : *number*, *list of values*, *array of values*
            Value (or list of values) with column indexes
            
        Return:
        =======
        ind : Index of the cell at (row, col) (or list/array of ids depending on the input type)
        """
        is_list = False                
        if type(row) == list:
            is_list = True
        
        ind = np.ravel_multi_index((row, col), self._dims)
        if is_list:
            return ind.tolist()
        else:
            return ind
    
    def get_network(self, threshold):
        """
        This function extract a drainage network by using a determined area threshold

        Parameters:
        ===========
        threshold : *int* 
          Area threshold to initiate a channels (in cells)
        """
        # Get the Flow Accumulation and select cells that meet the threholds
        fac = self.flow_accumulation(nodata = False).read_array()
        w = fac > threshold
        w = w.astype(np.int8)
        grid = self._create_output_grid(w, 0)
        return grid

    def flow_accumulation(self, nodata=True):
        """
        Calculates the flow accumulation from the topologically sorted pixels of the
        Flow object. As pixels of the Flow objects are sorted topologically, the flow
        accumulation can be obtained very fast with a computational time that is linearly
        dependent on the number of cell of the DEM.
        
        Parameters:
        ===========
        nodata : *bool*
          Boolean flag that indicates if the output flow accumulation grid will maintain NoData values. 
          If nodata=False, nodata values will be filled with 0 and NoDataValue will set to None. 

        Reference:
        ----------
        Braun, J., Willett, S.D., 2013. A very efficient O(n), implicit and parallel 
        method to solve the stream power equation governing fluvial incision and landscape 
        evolution. Geomorphology 180–181, 170–179. 
        """
        # TODO -- Test function
        facc = np.ones(self._ncells, np.uint32)
        nix = len(self._ix)
        for n in range(nix):
            facc[self._ixc[n]] = facc[self._ix[n]] + facc[self._ixc[n]]
        
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
        acc = _create_output_grid(facc, nodata_val)
        return acc

    def get_stream_poi(self, threshold, kind="heads"):
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
        # TODO -- Test function
        # Check input parameters
        if kind not in ['heads', 'confluences', 'outlets']:
            return np.array([]), np.array([])
        
        # Get the Flow Accumulation and select cells that meet the threholds
        fac = self.flow_accumulation(nodata = False).read_array()
        w = fac > threshold
        del fac
        
        # Build a sparse array with giver-receivers cells
        w = w.ravel()
        I   = w[self._ix]
        ix  = self._ix[I]
        ixc = self._ixc[I]
        aux_vals = np.ones(ix.shape, dtype=np.int8)
    
        sp_arr = csc_matrix((aux_vals, (ix, ixc)), shape=(self._ncells, self._ncells))
        del I, ix, ixc, aux_vals # Clean up (Don't want to leave stuff in memory)
        
        if kind == 'heads':
            # Heads will be channel cells marked only as givers (ix) but not as receivers (ixc) 
            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
            out_pos = (sum_arr == 0) & w
        elif kind == 'confluences':
            # Confluences will be channel cells with two or givers
            sum_arr = np.asarray(np.sum(sp_arr, 0)).ravel()
            out_pos = sum_arr > 1
        elif kind == 'outlets':
            # Outlets will be channel cells marked only as receivers (ix) but not as givers (ixc) 
            sum_arr = np.asarray(np.sum(sp_arr, 1)).ravel()
            out_pos = (sum_arr == 0) & w  
            
        out_pos = out_pos.reshape(self._dims)
        row, col = np.where(out_pos)
        
        return row, col

    def save_gtiff(self, path):
        """
        Saves the flow object as a geotiff. The geotiff file it wont have any
        sense if its open with GIS software.
        The organization of this geotiff is as follow:
        Band 1 --> Givers pixels reshaped to self._dims
        Band 2 --> Receiver pixels reshaped to self._dims
        Band 3 --> Nodata band (pixels with nodata == 1, pixels with data == 0)
        """
        # TODO -- Test function
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
        # TODO -- Test function
        raster = gdal.Open(path)

        # Set Network properties
        self._dims = (raster.RasterYSize, raster.RasterXSize)
        self._ncells = raster.RasterYSize * raster.RasterXSize
        self._geot = raster.GetGeoTransform()
        self._proj = raster.GetProjection()
        self._cellsize = self._geot[1]

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
      
    def drainage_basins(self, outlets = []):
        """
        This function extracts the drainage basins for the Flow object and returns a Grid object that can
        be saved into the disk.
        
        Parameters:
        ===========
        outlets : *list* or *tuple*
          List or tuple with (xi, yi) coordinate for outlets. xi and xi can be numbers, lists, or numpy.ndarrays
          If outlets is an empty list (default) the basins will be extracted for all the outlets in the Flow object.
        
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
        temp_ix = self._ix
        temp_ixc = self._ixc
        
        # If outlets are not specified, basins for all the outlets will be extracted
        if outlets == []:
            nbasins = 0
            D = np.zeros(self._ncells, np.int)
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
            if not len(x) == len(y):
                return
            row, col = self.xy_2_cell(x, y)
            inds = self.cell_2_ind(row, col)
            D = np.zeros(self._ncells, np.int)
            for n, inds in enumerate(inds):
                D[inds] = n+1
            
            nix = len(temp_ix)
            for n in range(nix-1,-1,-1):
                if (D[temp_ixc[n]] != 0) & (D[temp_ix[n]] == 0):
                    D[temp_ix[n]] = D[temp_ixc[n]]
        
        D = D.reshape(self._dims)  
        basins = self._create_output_grid(self, D, 0)
        return basins  
            
    def _create_output_grid(self, array, nodata_value=None):
        """
        This convenience function creates a Grid object from an input array. The array
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
        topopy.Grid object with the same properties that Flow
        """
        # TODO -- Test Function
        grid = Grid()
        grid._size = (self._dims[1], self._dims[0])
        grid._geot = self._geot
        grid._cellsize = self._geot[1]
        grid._proj = self._proj
        grid._nodata = nodata_value
        
        grid._array = array
        grid._tipo = str(array.dtype)
        
        return grid     
            
            
            