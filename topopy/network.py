# -*- coding: utf-8 -*-

# network.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# Version: 1.0
# December 26, 2017
#
# Last modified December 26, 2017

import numpy as np
from scipy import ndimage
from skimage import graph
from . import Grid

class Network():
    
    def __init__(self, dem="", preprocess='carve'):
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
        This algoritm has been adapted to Python from FLOWobj.m by Wolfgang Schwanghart 
        (version of 17. August, 2017) included in TopoToolbox matlab codes. If use, please
        cite:
                
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
            
            
    def flow_accumulation(self):
        """
        TODO -- Escribir documentacion
        """
        facc = np.ones(self._ncells, np.uint32)
        nix = len(self._ix)
        for n in range(nix):
            facc[self._ixc[n]] = facc[self._ix[n]] + facc[self._ixc[n]]
        
        facc = facc.reshape(self._dims)
        maxval = np.iinfo(np.uint32).max
        row, col = np.unravel_index(self._nodata_pos, self._dims)
        facc[row, col] = maxval
        
        acc = Grid()
        acc._size = (self._dims[1], self._dims[0])
        acc._geot = self._geot
        acc._cellsize = self._geot[1]
        acc._proj = self._proj
        acc._nodata = maxval
        acc._array = facc
        acc._tipo = str(facc.dtype)
        
        return acc


        
            
            
            
            
            