
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
from topopy import Grid
infolder = "data/in"
outfolder = "data/out"

class TestGrid01(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load(infolder + "/np_files/small25_100rnd_id.npy")
        self.rows = np.load(infolder + "/np_files/small25_100rnd_row.npy")
        self.cols = np.load(infolder + "/np_files/small25_100rnd_col.npy")
        self.xi = np.load(infolder + "/np_files/small25_100rnd_X.npy")
        self.yi = np.load(infolder + "/np_files/small25_100rnd_Y.npy")
        self.zi = np.load(infolder + "/np_files/small25_100rnd_Z.npy")
        
    def test_set_read_array_00(self):
        # Set data to an empty array and to an array with differnt dimensions
        arr = np.arange(25).reshape(5,5).astype(np.int8)
        dem = Grid()
        res01 = dem.set_array(arr)
        res02 = dem.set_array(np.arange(9).reshape(3,3).astype(np.int8))
        computed = (res01, res02)
        self.assertEqual(computed, (None, 0))
    
    def test_set_read_array_01(self):
        # Read array (NOT Sharing memory)
        arr = np.arange(25).reshape(5,5).astype(np.int8)
        dem = Grid()
        dem.set_array(arr)
        out_arr = dem.read_array(True)
        out_arr[0, 0] = 99
        computed = dem.get_value(0, 0)
        expected = 0
        self.assertEqual(computed, expected)    
        
    def test_set_read_array_02(self):
        # Read array (Sharing memory)
        arr = np.arange(25).reshape(5,5).astype(np.int8)
        dem = Grid()
        dem.set_array(arr)
        out_arr = dem.read_array(False)
        out_arr[0, 0] = 99
        computed = dem.get_value(0, 0)
        expected = 99
        self.assertEqual(computed, expected)    
            
    def test_set_data_00(self):
        arr = np.arange(9).reshape((3, 3))
        arr[[1, 2], [2, 1]] = 8
        dem = Grid()
        dem.set_nodata(8)
        dem.set_array(arr)
        row, col = dem.get_nodata_pos()
        computed = dem.get_value(row, col)
        expected = np.array([8, 8, 8])
        self.assertEqual(np.array_equal(computed, expected), True)
          
    def test_save(self):
        dem = Grid(infolder + "/small25.tif")
        dem.save(outfolder + "/a_dummy_dem.tif")
        dem2 = Grid(outfolder + "/a_dummy_dem.tif")
        expected = dem.get_value([20, 30, 40, 50], [20, 30, 40, 50])
        computed = dem2.get_value([20, 30, 40, 50], [20, 30, 40, 50])
        self.assertEqual(np.array_equal(computed, expected), True)

    def test_save_02(self):
        # Testing nodata with value of Zero
        dem = Grid()
        np.random.seed(1)
        arr = np.random.randint(0, 100, (25, 25))
        arr[np.where(arr%5==0)] = 0
        dem.set_array(arr)
        dem.set_nodata(0)
        dem.save(outfolder + "/a_dummy_dem2.tif")
        # Open with gdal
        raster = gdal.Open(outfolder + "/a_dummy_dem2.tif")
        banda = raster.GetRasterBand(1)
        nodata = banda.GetNoDataValue()
        self.assertEqual(nodata, 0)
           
    def test_values_2_nodata(self):
        dem = Grid()
        dem.set_array(np.arange(25).reshape((5, 5)))
        dem.set_nodata(-99.)
        dem.values_2_nodata([10, 11, 12, 13, 14])
        row, col = dem.ind_2_cell([10, 11, 12, 13, 14])
        computed = dem.get_value(row, col)
        expected = np.array([-99, -99, -99, -99, -99])
        res = np.array_equal(computed, expected)
        self.assertEqual(res, True)

    def test_get_value_01(self):
        dem = Grid(infolder + "/small25.tif")
        # Taking row, col in a nan position (88)
        ind = 88
        row, col = self.rows[ind], self.cols[ind]
        computed = dem.get_value(row, col)
        self.assertEqual(computed, -9999)

    def test_get_value_02(self):
        dem = Grid(infolder + "/small25.tif")
        # Taking row, col in other position (with value)
        ind = 25
        row, col = self.rows[ind], self.cols[ind]
        expected = self.zi[ind]
        computed = dem.get_value(row, col)
        self.assertEqual(computed, expected)
        
    def test_get_value_03(self):
        dem = Grid(infolder + "/small25.tif")
        # Taking row, col outside array
        row, col = 199, 133
        self.assertRaises(IndexError, dem.get_value, row, col)
    
    def test_get_value_04(self):
        dem = Grid(infolder + "/small25.tif")
        # Taking row, col as numpy arrays
        expected = self.zi
        computed = dem.get_value(self.rows, self.cols)
        self.assertEqual(np.nansum(computed),np.nansum(expected))
        
    def test_get_value_05(self):
        dem = Grid(infolder + "/small25.tif")
        # Taking row, col as lists
        expected = np.nansum(self.zi.tolist())
        res =  dem.get_value(self.rows.tolist(), self.cols.tolist())
        computed = np.nansum(res)
        self.assertEqual(computed, expected)
        
    def test_get_nodata_pos(self):
        dem = Grid(infolder + "/small25.tif")
        arr = dem._array        
        row, col = dem.get_nodata_pos()
        crow, ccol= np.where(arr == -9999)
        self.assertEqual((np.array_equal(row, crow), np.array_equal(col, ccol)), (True, True))
        
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
        raster = Grid(infolder + "/small25.tif")
        computed = raster.is_inside(x, y)
        expected = np.array([False, False, False, False, False, True, True ,False, True])
        self.assertEqual(np.array_equal(computed, expected), True)
        
        
    
if __name__ == "__main__":
    unittest.main()