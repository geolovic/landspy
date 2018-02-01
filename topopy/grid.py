# -*- coding: utf-8 -*-

# grid.py
# Jose Vicente Perez Pena
# Dpto. Geodinamica-Universidad de Granada
# 18071 Granada, Spain
# vperez@ugr.es // geolovic@gmail.com
#
# MIT License (see LICENSE file)
# Version: 1.0
# December 26, 2017
#
# Last modified January 24, 2018


import gdal
import numpy as np
from scipy import ndimage

# This import statement avoid issues with matplotlib in Mac when using Python not as a Framework
# If matplotlib is not imported, Grid.plot() will not work.
try: 
    import matplotlib.pyplot as plt
    PLT = True
except Exception as e:
    PLT = False


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
        
        if path:
            raster = gdal.Open(path)
            if not raster:
                raise FileNotFoundError
            
            banda = raster.GetRasterBand(band)
            self._size = (banda.XSize, banda.YSize)
            self._geot = raster.GetGeoTransform()
            self._cellsize = self._geot[1]
            self._proj = raster.GetProjection()
            self._nodata = banda.GetNoDataValue()
            self._array = banda.ReadAsArray()
            self._tipo = str(self._array.dtype)
        
        else:
            self._geot = (0., 1., 0., 0., 0., -1.)
            self._cellsize = self._geot[1]
            self._proj = ""
            self._nodata = None
            self._array = np.array([[0]], dtype=np.float)
            self._tipo = str(self._array.dtype)
            self._size = (self._array.shape[1], self._array.shape[0])
    
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
        # Si el grid no está vacio, solo se podra cambiar el array interno por otro de
        # dimensiones equivalentes    
        elif array.shape == (self._size[1], self._size[0]):    
            self._array = np.copy(array)
            self._tipo = str(self._array.dtype)
        else:
            return 0
        
    def read_array(self, ascopy=False):
        """
        Return the internal array of the Grid
        
        Parameters:
        ==========
        ascopy : *bool*
          If True, the returned array is a memory view of the Grid original array.
        """
        if ascopy:
            return np.copy(self._array)
        else:
            return self._array
    
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
        datapos = np.where(self._array != self._nodata)
        return np.max(self._array[datapos])
    
    def min(self):
        """
        Return the minimun value of the Grid
        """
        datapos = np.where(self._array != self._nodata)
        return np.min(self._array[datapos])
    
    def mean(self):
        """
        Return the mean value of the Grid
        """
        datapos = np.where(self._array != self._nodata)
        return np.mean(self._array[datapos])
    
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
        Return a tuple with the size of the grid (XSize, YSize)
        """
        return self._size
    
    def get_dims(self):
        """
        Return a tuple with the size of the internal array (nrow, ncol)
        """
        return self._array.shape
    
    def get_ncells(self):
        """
        Return the total number of cells of the Grid
        """
        return self._array.size
    
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
        if self._nodata is None:
            return (np.array([], dtype=np.int), np.array([], dtype=np.int))
        else:
            return np.where(self._array == self._nodata)

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
        Get row col indexes coordinates from XY coordinates
        
        Parameters:
        ===========
        x : X coordinates (number, list, or numpy.ndarray)
        y : Y coordinates (number, list, or numpy.ndarray)
            
        Return:
        =======
        (row, col) : Tuple of ndarray with row and column indexes
        """
        x = np.array(x)
        y = np.array(y)       
        row = (self._geot[3] - y) / self._geot[1]
        col = (x - self._geot[0]) / self._geot[1]
        return row.astype(np.int32), col.astype(np.int32)

    def cell_2_xy(self, row, col):
        """
        Get XY coordinates from row and column cell indexes
        
        Parameters:
        ===========
        row : row indexes (number, list, or numpy.ndarray)
        col : column indexes (number, list, or numpy.ndarray)
            
        Return:
        =======
        (x, y) : Tuple of ndarray with X and Y coordinates
        """
        row = np.array(row)
        col = np.array(col)
        x = self._geot[0] + self._geot[1] * col + self._geot[1] / 2
        y = self._geot[3] - self._geot[1] * row - self._geot[1] / 2
        return x, y
        
    def ind_2_cell(self, ind):
        """
        Get row col indexes from cells linear indexes (row-major, C-style)
        
        Parameters:
        ===========
        ind : linear indexes (number, list, or numpy.ndarray)
        
        Return:
        =======
        (row, col) : Tuple of ndarray with row and column indexes
        """
        row, col = np.unravel_index(ind, self._array.shape) 
        return (row, col)
    
    def cell_2_ind(self, row, col):
        """
        Get cell linear indexes (row-major, C-style) from row and column indexes
        
        Parameters:
        ===========
        row : row indexes (number, list, or numpy.ndarray)
        col : column indexes (number, list, or numpy.ndarray)
            
        Return:
        =======
        ind : Linear indexes (row-major, C-style)
        """
        ind = np.ravel_multi_index((row, col), self._array.shape)
        return ind
    
    def set_nodata(self, value):
        """
        Sets the nodata value for the Grid
        """
        self._nodata = value
        
    def values_2_nodata(self, value):
        """
        Change specific values to NoData (if Grid nodata is defined). 
        
        Parameters:
        ===========
        value : Value or values that will be changed to NoData
        """
        if self._nodata is None:
            return
        value = np.array(value)
        for val in value:
            idx = np.where(self._array == val)
            self._array[idx] = self._nodata        
    
    def plot(self, ax=None):
        """
        Plots the grid in a new Axes or in a existing one
        
        Parameters:
        ===========
        ax : *matplotlib.Axe*
          If is not defined, the function will use plt.imshow()
        """
        if not PLT:
            return
        
        arr = self.read_array()
        if self._nodata:
            mask = self._array == self._nodata
            arr = np.ma.array (self._array, mask = mask)
        
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
        
        # Prepare driver to write raster
        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._size[0], self._size[1], 1, tipo)
        if not raster:
            return
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)
        
        if not self._nodata is None:
            raster.GetRasterBand(1).SetNoDataValue(self._nodata)
        
        raster.GetRasterBand(1).WriteArray(self._array)

class DEM(Grid):
    
    def identify_flats(self, nodata=True):
        """
        This functions returns two binary Grid objects (values 0, 1) with flats and sills. 
        Flats are defined  as cells without downward neighboring cells. Sills are cells where 
        flat regions spill over into lower terrain. It the DEM has nodata, they will be maintained
        in output grids.
        
        Returns:
        ----------
        [flats, sills] : Grid objects
        
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
        
        # Change nodata to very low values
        nodata_ids = self.get_nodata_pos()
        z_arr[nodata_ids] = -9999
        
        footprint = np.ones((3, 3), dtype=np.int8)
        
        # Identify flats throught a image binary erosion
        # Flats will be True where cells don't have lower neighbors
        flats = ndimage.morphology.grey_erosion(z_arr, footprint=footprint) == z_arr
        
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
        sills = np.empty(z_arr.shape)
        sills.fill(-9999.)
        sills[flats] = z_arr[flats]
        aux_dil = ndimage.morphology.grey_dilation(sills, footprint=footprint)
        sills = np.logical_and(aux_dil == z_arr, np.logical_not(flats))
        sills[nodata_ids] = False
        
        # Prepare outputs (as Grid objects)
        res = []
        for arr in [flats, sills]:
            grid = Grid()
            grid.copy_layout(self)
            grid.set_nodata(-1)
            grid._tipo = 'int8'
            grid.set_array(arr.astype(np.int8))
            
            if nodata:
                grid._array[nodata_ids] =  -1
            else:
                grid._array[nodata_ids] =  0
                grid.set_nodata(None)
            
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
        # Change nan values to a very low value
        copyarr = np.copy(self._array)
        nodata_pos = self.get_nodata_pos()
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
        if self._nodata:
            output_array[nodata_pos] = self._nodata
        # Create output filled DEM
        filled_dem = DEM()
        filled_dem.copy_layout(self)
        filled_dem.set_array(output_array)
        return filled_dem

                   