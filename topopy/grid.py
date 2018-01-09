# -*- coding: utf-8 -*-

# Grid.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# Version: 1.0
# December 26, 2017
#
# Last modified December 26, 2017


import gdal
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}

class Grid():
        
    def __init__(self, path="", band=1):
        """
        Class to manipulate rasters
        
        Parameters:
        ================
        path : *str* 
          Path to the raster
        band : *str*
          Band to be open (usually don't need to be modified)
        """
        raster = gdal.Open(path)
        
        if raster:
            banda = raster.GetRasterBand(band)
            self._size = (banda.XSize, banda.YSize)
            self._geot = raster.GetGeoTransform()
            self._cellsize = self._geot[1]
            self._proj = raster.GetProjection()
            self._nodata = banda.GetNoDataValue()
            self._array = banda.ReadAsArray()
            self._tipo = str(self._array.dtype)
            self._put_nans()
        else:
            self._geot = (0., 1., 0., 0., 0., -1.)
            self._cellsize = self._geot[1]
            self._proj = ""
            self._nodata = None
            self._array = np.array([[0]], dtype="float")
            self._tipo = str(self._array.dtype)
            self._size = self._array.shape
            
    def _put_nans(self):
        """
        Changes nodata values by nans
        """
        self._array = self._array.astype("float")
        if self._nodata:
            inds = np.where(self._array == self._nodata)
            self._array[inds] = np.nan
    
    def copy_layout(self, grid):
        """
        Copy all the parameters from another Grid instance except grid data (internal array)
        
        Parameters:
        ================
        grid : *Grid* 
          Grid instance from which parameters will be copied
        """
        self._size = grid.get_size()
        self._geot = grid.get_geotransform()
        self._cellsize = self._geot[1]
        self._proj = grid.get_projection()
        self._nodata = grid.get_nodata()
        
    def set_array(self, array):
        """
        Set data array for a Grid Object. 
        If Grid is an empty Grid [get_size() = (1, 1)], data are set and its size recalculated
        If Grid is not an empty Grid, the input array should match with Grid Size
        
        Parameters:
        ================
        array : *numpy array* 
          Numpy array with the data
        """
        # Si hemos creado un grid vacio, le podemos especificar el array
        if self._size == (1, 1):          
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
            self._size = (array.shape[1], array.shape[0])
            self._put_nans()       
        # Solo se podra cambiar el array interno si las dimensiones coinciden    
        elif array.shape == (self._size[1], self._size[0]):    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
            self._put_nans()
        else:
            return 0
        
    def read_array(self, nans=True):
        """
        Return the internal array of the Grid
        
        Parameters:
        ==========
        nans : *bool*
          If True, the array is read with nan values. If False, the array is read in its original 
          type and with nodata values instead of nans
        """
        if nans:
            return self._array
        else:
            copyarr = np.copy(self._array)
            if self._nodata:
                copyarr[np.where(np.isnan(copyarr))] = self._nodata
            return copyarr.astype(self._tipo)  
    
    def find(self):
        """
        Find the non-zero elements in the array. Return a tuple of arrays with
        row and col positios.
        """
        return np.where(self._array > 0)
    
    def max(self):
        """
        Return the maximun value of the Grid
        """
        return np.nanmax(self._array)
    
    def min(self):
        """
        Return the minimun value of the Grid
        """
        return np.nanmin(self._array)
    
    def mean(self):
        """
        Return the mean value of the Grid
        """
        return np.nanmean(self._array)
    
    def set_value(self, row, col, value):
        """
        Set the value for a cell of the grid at (row, col)
        Input value will be converted to array datatype
        
        Parameters:
        ================
        row, col : *int* 
          Row and column indexes
        value : *number*
          Value for the cell (row, col)
        """
        self._array[row, col] = value
    
    def get_value(self, row, col):
        """
        Get the value for a cell of the grid at (row, col)
        
        Parameters:
        ================
        row, col : *int* 
          Row and column indexes
        """
        if type(row) == list:
            return self._array[row, col].tolist()
        else:    
            return self._array[row, col]
    
    def get_size(self):
        """
        Return a tuple with the size of the grid (Xsize, Ysize)
        """
        return self._size
    
    def get_projection(self):
        """
        Return a string with the projection of the grid in WKT
        """
        return self._proj
    
    def get_nodata(self):
        """
        Return the value for NoData in the grid. This value could be None if NoData
        are not defined in the grid.
        """
        return self._nodata
    
    def get_nodata_pos(self):
        """
        Return the position of the NoData values as a tuple of two arrays (rows, columns)
        """
        return np.where(np.isnan(self._array))

    def get_cellsize(self):
        """
        Return the grid cellsize
        """
        return self._cellsize
    
    def get_geotransform(self):
        """
        Return the GeoTranstorm matrix of the grid. This matrix has the form:
        *(ULx, Cx, Tx, ULy, Ty, Cy)*
        
        * ULx = Upper-Left X coordinate (upper-left corner of the pixel)
        * ULy = Upper-Left Y coordinate (upper-left corner of the pixel)
        * Cx = X Cellsize
        * Cy = Y Cellsize (negative value)
        * Tx = Rotation in X axis
        * Ty = Rotation in Y axis
        """
        return self._geot            

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
        
        row, col = np.unravel_index(ind, self._array.shape) 
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
        
        ind = np.ravel_multi_index((row, col), self._array.shape)
        if is_list:
            return ind.tolist()
        else:
            return ind
    
    def set_nodata(self, value):
        """
        Sets the nodata value for the Grid
        """
        self._nodata = float(value)
        self._put_nans()
        
    def values_2_nodata(self, value):
        """
        Changes values to nodata,  if Grid nodata is defined (whether is not None). 
        
        Parameters:
        ===========
        value : *number*, *list of values*, *array of values*
          Value or values that will change to NoData
        """
        if not self._nodata:
            return
        
        is_number = False
        if type(value) == int or type(value) == float:
            is_number = True
        
        if not is_number:
            inds = np.array([])
            for num in value:
                idx = np.where(self._array == num)
                self._array[idx] = self._nodata        
        else:
            inds = np.where(self._array == value)
            self._array[inds] = self._nodata
        
        self._put_nans()
        
    def plot(self, ax=None):
        """
        Plots the grid in a new Axes or in a existing one
        
        Parameters:
        ===========
        ax : *matplotlib.Axe*
          If is not defined, the function will use plt.imshow()
        """
        arr = self._array.astype("float")
        arr = np.copy(arr)
        if self._nodata:
            nodata = float(self._nodata)
            ids = np.where(arr==nodata)
            arr[ids] = np.nan
        if ax:
            ax.imshow(arr)
        else:
            plt.imshow(arr)
    
    def save(self, path):
        """
        Saves the grid in the disk
        
        Parameters:
        ================
        path : *str* 
          Path where new raster will be saved
        """
        # Check if the type of the internal array is compatible with gdal
        if str(self._tipo) not in NTYPES.keys():
            return 0
        else:
            tipo = NTYPES[str(self._tipo)]

        # Put back nodata values and change type
        copyarr = np.copy(self._array)
        if self._nodata:
            copyarr[np.where(np.isnan(copyarr))] = self._nodata        
        copyarr = copyarr.astype(self._tipo)
        
        # Prepare driver to write raster
        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._size[0], self._size[1], 1, tipo)
        if not raster:
            return
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)
        if self._nodata:
            raster.GetRasterBand(1).SetNoDataValue(self._nodata)
        
        raster.GetRasterBand(1).WriteArray(copyarr)

class DEM(Grid):
    
    def identify_flats(self, nodata=True):
        """
        This functions returns two binary Grid objects (values 0, 1) with flats and sills. 
        Flats are defined  as cells without downward neighboring cells. Sills are cells where 
        flat regions spill over into lower terrain. It the DEM has nodata, they will be maintained
        in output grids.
        
        Returns:
        ----------
        [flats, dills] : Grid objects
        
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

        z_arr = self._array
        
        # Handle nodata (nan values)
        nans = False
        if np.isnan(z_arr.max()):
            nans = True
            nans_ids = np.where(np.isnan(z_arr))
            z_arr[nans_ids] = -9999.
        
        footprint = np.ones((3, 3), dtype="int")
        
        # Identify flats throught a image binary erosion
        # Flats will be True where cells don't have lower neighbors
        flats = ndimage.morphology.grey_erosion(z_arr, footprint=footprint) == z_arr
        
        # Remove flats from the borders
        flats[:,0] = False
        flats[:,-1] = False
        flats[0,:] = False
        flats[-1,:] = False
        
        # Remove flats for nodata values and cells bordering them
        if nans:
            flats[nans_ids] = False
            auxmat = np.zeros(flats.shape, dtype="bool")
            auxmat[nans_ids] = True
            nodata_bord = ndimage.morphology.grey_dilation(auxmat, footprint=footprint)
            flats[nodata_bord] = False
        
        # Identify sills
        sills = np.empty(z_arr.shape)
        sills.fill(-9999.)
        sills[flats] = z_arr[flats]
        aux_dil = ndimage.morphology.grey_dilation(sills, footprint=footprint)
        sills = np.logical_and(aux_dil == z_arr, np.logical_not(flats))
        if nans:
            sills[nans_ids] = False
        
        # Prepare outputs (as Grid objects)
        res = []
        for arr in [flats, sills]:
            grid = Grid()
            grid.copy_layout(self)
            grid.set_array(arr)
            if nodata and nans:
                grid._array[nans_ids] = np.nan
            res.append(grid)
        
        return res

    
    def fill_sinks(self, four_way=False):
        """
        Fill sinks method adapted from  fill depressions/sinks in floating point array
        
        Parameters:
        ----------
        input_array : [ndarray] Input array to be filled
        four_way : [bool] Searchs the 4 (True) or 8 (False) adjacent cells
        
        Returns:
        ----------
        [ndarray] Filled array
    
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
        # Change nan values to a high value (will not have)
        copyarr = np.copy(self._array)
        nan_inds = np.where(np.isnan(copyarr))
        copyarr[nan_inds] = -9999.
        
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
        if four_way:
            el = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]).astype(np.bool)
        else:
            el = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]).astype(np.bool)
    
        # Iterate until marker array doesn't change
        while not np.array_equal(output_old_array, output_array):
            output_old_array = np.copy(output_array)
            output_array = np.maximum(
                copyarr,
                ndimage.grey_erosion(output_array, size=(3, 3), footprint=el))

        # Put back nodata values and change type
        output_array[nan_inds] = self._nodata
        output_array = output_array.astype(self._tipo)
        # Create output filled DEM
        filled_dem = DEM()
        filled_dem.copy_layout(self)
        filled_dem.set_array(output_array)
        return filled_dem

                   