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
import matplotlib.pyplot as plt


NTYPES = {'int8': 3, 'int16': 3, 'int32': 5, 'int64': 5, 'uint8': 1, 'uint16': 2,
          'uint32': 4, 'uint64': 4, 'float16': 6, 'float32': 6, 'float64': 7}
GTYPES = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32', 6: 'float32', 
          7: 'float64'}

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
            self._tipo = NTYPES[str(self._array.dtype)]
            if self._tipo <= 5 and self._nodata:
                self._nodata = int(self._nodata)
        else:
            self._size = (0, 0)
            self._geot = (0., 1., 0., 0., 0., -1.)
            self._cellsize = self._geot[1]
            self._proj = ""
            self._nodata = None
            self._array = np.empty((0,0))
            self._tipo = NTYPES[str(self._array.dtype)]
    
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
        
    def set_data(self, array):
        """
        Set data array for a Grid Object. 
        If Grid is an empty Grid [get_size() = (0, 0)], data are set and its size recalculated
        If Grid is not an empty Grid, the input array should match with Grid Size
        
        Parameters:
        ================
        array : *numpy array* 
          Numpy array with the data
        """
        if self._size == (0, 0):
            # Si hemos creado un grid vacio, le podemos especificar el array
            self._array = np.copy(array)
            self._tipo = NTYPES[str(self._array.dtype)]
            self._size = (array.shape[1], array.shape[0])
            
        elif array.shape == (self._size[1], self._size[0]):
            # Solo se podra cambiar el array interno si las dimensiones coinciden
            self._array = np.copy(array)
            self._tipo = NTYPES[str(self._array.dtype)]
            if self._tipo <= 5 and self._nodata:
                self._nodata = int(self._nodata)
            elif self._tipo > 5 and self._nodata:
                self._nodata = float(self._nodata)
        else:
            return
        
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
    
    def read_array(self):
        """
        Return the internal array of the Grid
        """
        return self._array

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
        TODO Escribir documentacion
        TODO Crear test
        TODO Añadir funcionalidad para listas y numeros (no arrays)
        """
        return np.unravel_index(ind, self._array.shape)
    
    def cell_2_ind(self, row, col):
        """
        TODO Escribir documentacion
        TODO Crear test
        TODO Añadir funcionalidad para listas y numeros (no arrays)
        """
        return np.ravel_multi_index((row, col), self._array.shape, mode='clip')
    
    def fill(self, output=""):
        #TODO Crear funcion
        pass
    
    def save(self, path):
        """
        Saves the grid in the disk
        
        Parameters:
        ================
        path : *str* 
          Path where new raster will be saved
        """
        if str(self._array.dtype) not in NTYPES.keys():
            return
        else:
            tipo = NTYPES[str(self._array.dtype)]

        driver = gdal.GetDriverByName("GTiff")
        raster = driver.Create(path, self._size[0], self._size[1], 1, tipo)
        raster.SetGeoTransform(self._geot)
        raster.SetProjection(self._proj)
        if self._nodata:
            raster.GetRasterBand(1).SetNoDataValue(self._nodata)
        
        raster.GetRasterBand(1).WriteArray(self._array)
        