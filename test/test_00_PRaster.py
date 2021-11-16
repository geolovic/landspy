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
from osgeo import gdal
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import PRaster
infolder = "data/in"

class TestPRaster00(unittest.TestCase):
    
    def test_properties(self):
        # Test PRaster getters
        files = ["small25.tif", "tunez.tif", "jebja30.tif"]
        for file in files:
            raster = PRaster(infolder + "/{0}".format(file))
            size = raster.get_size()
            dims = raster.get_dims()
            ncells = raster.get_ncells()
            cellsize = raster.get_cellsize()
            geot = raster.get_geotransform()
            computed = (size, dims, ncells, cellsize, geot)
            graster = gdal.Open(infolder + "/{0}".format(file))
            band = graster.GetRasterBand(1)
            arr = band.ReadAsArray()
            ggeot = graster.GetGeoTransform()
            expected = ((band.XSize, band.YSize), arr.shape, arr.size, (ggeot[1], ggeot[5]), ggeot)
            self.assertEqual(computed, expected)
            
    def test_projections(self):
        # Test PRaster getters
        files = ["small25.tif", "tunez.tif", "jebja30.tif"]
        
        for file in files:
            # Test PRaster getters for a raster
            raster = PRaster(infolder + "/{0}".format(file))
            graster = gdal.Open(infolder + "/{0}".format(file))
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
        expected = ((1, 1), (1, 1), 1, (1.0, -1.0), (0., 1., 0., 1., 0., -1.), "")
        self.assertEqual(computed, expected)
        
    def test_copy_layout(self):
        # Test copy_layout for some rasters
        files = ["small25.tif", "tunez.tif", "jebja30.tif"]
        for file in files:
            b_raster = PRaster()
            c_raster = PRaster(infolder + "/{0}".format(file))
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

class TestPRaster01(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load(infolder + "/np_files/small25_100rnd_id.npy")
        self.rows = np.load(infolder + "/np_files/small25_100rnd_row.npy")
        self.cols = np.load(infolder + "/np_files/small25_100rnd_col.npy")
        self.xi = np.load(infolder + "/np_files/small25_100rnd_X.npy")
        self.yi = np.load(infolder + "/np_files/small25_100rnd_Y.npy")
        
    def test_xy_2_cell_01(self):   
        raster = PRaster(infolder + "/small25.tif")
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        c_rows, c_cols = raster.xy_2_cell(xi, yi)
        res = (np.array_equal(rows, c_rows), np.array_equal(cols, c_cols))
        self.assertEqual(res, (True, True))
        
    def test_xy_2_cell_02(self):
        raster = PRaster(infolder + "/small25.tif")
        x = 471927
        y = 4116048
        row, col = raster.xy_2_cell(x, y)
        self.assertEqual((43, 71), (row, col))
        
    def test_xy_2_cell_03(self):
        raster = PRaster(infolder + "/small25.tif")
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows
        cols = self.cols
        c_rows, c_cols = raster.xy_2_cell(xi, yi)
        res = (np.array_equal(rows, c_rows), np.array_equal(cols, c_cols))
        self.assertEqual(res, (True, True))
        
        
    def test_cell_2_xy_01(self):   
        raster = PRaster(infolder + "/small25.tif")
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        cxi, cyi = raster.cell_2_xy(rows, cols)
        res = (np.array_equal(cxi, xi), np.array_equal(cyi, yi))
        self.assertEqual(res, (True, True))
        
    def test_cell_2_xy_02(self):
        raster = PRaster(infolder + "/small25.tif")
        x = 471927.1
        y = 4116048.5
        row =  43
        col = 71
        cx, cy = raster.cell_2_xy(row, col)
        self.assertEqual((x, y), (cx, cy))
        
    def test_cell_2_xy_03(self):
        raster = PRaster(infolder + "/small25.tif")
        xi = self.xi
        yi = self.yi
        rows = self.rows.tolist()
        cols = self.cols.tolist()
        cxi, cyi = raster.cell_2_xy(rows, cols)
        res = (np.array_equal(xi, cxi), np.array_equal(yi, cyi))
        self.assertEqual(res, (True, True))      
        
    def test_cell_2_ind_01(self):   
        raster = PRaster(infolder + "/small25.tif")
        rows = self.rows
        cols = self.cols
        ids = self.ids
        cids = raster.cell_2_ind(rows, cols)
        self.assertEqual(np.array_equal(ids, cids), True)
        
    def test_cell_2_ind_02(self):
        raster = PRaster(infolder + "/small25.tif")
        row =  25
        col = 11
        ids = 4986
        cids = raster.cell_2_ind(row, col)
        self.assertEqual(ids, cids)
    
    def test_cell_2_ind_03(self):   
        raster = PRaster(infolder + "/small25.tif")
        rows = self.rows.tolist()
        cols = self.cols.tolist()
        ids = self.ids
        cids = raster.cell_2_ind(rows, cols)
        self.assertEqual(np.array_equal(ids, cids), True)
        
    def test_is_inside(self):
        # points 1, 2, 4, 8 --> -Fuera de ráster
        # points 3, 5 --> Dentro de ráster, pero en NoData
        # points 6, 7, 9 --> Dentro de ráster
        
        puntos = np.array([[476154., 4115084.],
                          [472289., 4112838.],
                          [471317., 4114050.],
                          [472874., 4117717.],
                          [472205., 4114091.],
                          [470795., 4116411.],
                          [472257., 4115565.],
                          [469572., 4115376.],
                          [473877., 4114844.]])
        x = puntos[:,0]
        y = puntos[:,1]
        raster = PRaster(infolder + "/small25.tif")
        computed = raster.is_inside(x, y)
        expected = np.array([False, False, True, False, True, True, True ,False, True])
        self.assertEqual(np.array_equal(computed, expected), True)
        
if __name__ == "__main__":
    unittest.main()