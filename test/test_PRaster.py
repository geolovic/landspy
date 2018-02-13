# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""

import unittest
import sys
import numpy as np
import gdal
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import PRaster

class TestPRaster00(unittest.TestCase):
    
    def test_properties(self):
        # Test PRaster getters
        files = ["small25.tif", "tunez.tif", "tunez2.tif"]
        for file in files:
            raster = PRaster("data/{0}".format(file))
            size = raster.get_size()
            dims = raster.get_dims()
            ncells = raster.get_ncells()
            cellsize = raster.get_cellsize()
            geot = raster.get_geotransform()
            computed = (size, dims, ncells, cellsize, geot)
            graster = gdal.Open("data/{0}".format(file))
            band = graster.GetRasterBand(1)
            arr = band.ReadAsArray()
            ggeot = graster.GetGeoTransform()
            expected = ((band.XSize, band.YSize), arr.shape, arr.size, ggeot[1], ggeot)
            self.assertEqual(computed, expected)
            
    def test_projections(self):
        # Test PRaster getters
        files = ["small25.tif", "tunez.tif", "tunez2.tif"]
        
        for file in files:
            # Test PRaster getters for a raster
            raster = PRaster("data/{0}".format(file))
            graster = gdal.Open("data/{0}".format(file))
            self.assertEqual(raster.get_projection(), graster.GetProjection())
            
    def test_empty(self):
        # Test PRaster getters in an empty PRaster
        raster = PRaster()
        size = raster.get_size()
        dims = raster.get_dims()
        ncells = raster.get_ncells()
        cellsize = raster.get_cellsize()
        geot = raster.get_geotransform()
        proj = raster.get_projection()
        computed = (size, dims, ncells, cellsize, geot, proj)
        expected = ((1, 1), (1, 1), 1, 1., (0., 1., 0., 0., 0., -1.), "")
        self.assertEqual(computed, expected)
        
    def test_copy_layout(self):
        # Test copy_layout for some rasters
        files = ["small25.tif", "tunez.tif", "tunez2.tif"]
        for file in files:
            b_raster = PRaster()
            c_raster = PRaster("data/{0}".format(file))
            b_raster.copy_layout(c_raster)
            
            size = c_raster.get_size()
            dims = c_raster.get_dims()
            ncells = c_raster.get_ncells()
            cellsize = c_raster.get_cellsize()
            geot = c_raster.get_geotransform()
            proj = c_raster.get_projection()
            computed = (size, dims, ncells, cellsize, geot, proj)
    
            size = b_raster.get_size()
            dims = b_raster.get_dims()
            ncells = b_raster.get_ncells()
            cellsize = b_raster.get_cellsize()
            geot = b_raster.get_geotransform()
            proj = b_raster.get_projection()
            expected = (size, dims, ncells, cellsize, geot, proj)
            
            self.assertEqual(computed, expected)

#class TestPRaster01(unittest.TestCase):
#    
#    def setUp(self):        
#        # Load test data
#        self.ids = np.load("data/small25_100rnd_id.npy")
#        self.rows = np.load("data/small25_100rnd_row.npy")
#        self.cols = np.load("data/small25_100rnd_col.npy")
#        self.xi = np.load("data/small25_100rnd_X.npy")
#        self.yi = np.load("data/small25_100rnd_Y.npy")
#        
#    def test_xy_2_cell_01(self):        
#        raster = PRaster("data/small25.tif")
#        xi = self.xi
#        yi = self.yi
#        rows = self.rows
#        cols = self.cols
#        c_rows, c_cols = raster.xy_2_cell(xi, yi)
#        res = (np.array_equal(rows, c_rows), np.array_equal(cols, c_cols))
#        self.assertEqual(res, (True, True))
#        
#    def test_xy_2_cell_02(self):
#        raster = PRaster("data/small25.tif")
#        x = 471927
#        y = 4116048
#        row, col = raster.xy_2_cell(x, y)
#        self.assertEqual((43, 71), (row, col))
#        
        
if __name__ == "__main__":
    unittest.main()