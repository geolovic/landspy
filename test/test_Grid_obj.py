#!/usr/bin/env python3
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
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid

MY_GRID = "data/small25.tif"

class GridPropertyTests(unittest.TestCase):
    
    def setUp(self):        
        # Create a Grid object
        self.dem = Grid(MY_GRID)
        
    def test_get_size(self):
        expected = (199, 133)
        computed = self.dem.get_size()
        self.assertEqual(computed, expected)
   
    def test_get_projection(self):
        expected = 'PROJCS["ETRS89 / UTM zone 30N",GEOGCS["ETRS89",DATUM["European_Terrestrial_Reference_System_1989",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6258"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4258"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","25830"]]'
        computed = self.dem.get_projection()
        self.assertEqual(computed, expected) 

    def test_get_nodata(self):
        expected = (-9999, int)
        computed = (self.dem.get_nodata(), type(self.dem.get_nodata()))
        self.assertEqual(computed, expected)
        
    def test_get_nodata2(self):
        dem2 = Grid()
        dem2.copy_layout(self.dem)
        ncol, nrow = self.dem.get_size()
        fix_arr = np.empty((nrow, ncol), dtype="float16")
        fix_arr.fill(25.)
        dem2.set_data(fix_arr)
        expected = (-9999.0, float)
        computed = (dem2.get_nodata(), type(dem2.get_nodata()))
        self.assertEqual(computed, expected)
        
    def test_get_cellsize(self):
        expected = 25.
        computed = self.dem.get_cellsize()
        self.assertEqual(computed, expected)
        
    def test_get_geot(self):
        expected = (470139.6, 25.0, 0.0, 4117136.0, 0.0, -25.0)
        computed = self.dem.get_geotransform()
        self.assertEqual(computed, expected)

class GridCopySaveTests(unittest.TestCase):
    
    def setUp(self):        
        # Create a Grid object
        self.dem = Grid(MY_GRID)
        
    def test_copy_layout(self):
        dem2 = Grid()
        dem2.copy_layout(self.dem)
        expected = (self.dem._size, self.dem._geot, self.dem._cellsize, 
                    self.dem._proj, self.dem._nodata)
        
        computed = (dem2._size, dem2._geot, dem2._cellsize, 
            dem2._proj, dem2._nodata)
        
        self.assertEqual(computed, expected)
        
    def test_set_data_00(self):
        dem = Grid()
        data = np.arange(9).reshape((3, 3))
        dem.set_data(data)
        
        expected = ((0., 1., 0., 0., 0., -1.), 4)
        computed = (dem.get_geotransform(), dem.get_value(1, 1))
        self.assertEqual(computed, expected)
        
    def test_set_data_01(self):
        dem = Grid()
        dem.copy_layout(self.dem)
        # try with an array with different size
        dem.set_data(np.arange(9).reshape((3, 3)))
        
        expected = np.empty((0,0))
        computed = dem._array
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)
        
    def test_set_data_02(self):
        dem = Grid()
        dem.copy_layout(self.dem)
        # try with an array with different size
        ncol, nrow = self.dem.get_size()
        fix_arr = np.empty((nrow, ncol), dtype="float")
        fix_arr.fill(25.)
        dem.set_data(fix_arr)
        computed = dem.get_value(100, 100)
        expected = 25.
        self.assertEqual(computed, expected)
        
    def test_set_value_read_array(self):
        # First create a copy of the grid
        dem = Grid()
        dem.copy_layout(self.dem)
        dem.set_data(self.dem.read_array())
        
        row_ids = np.random.randint(0, dem.get_size()[1], 25)
        col_ids = np.random.randint(0, dem.get_size()[0], 25)
        
        computed = []
        expected = []
        
        for n in range(25):
            value = np.random.randint(888, 999)
            dem.set_value(row_ids[n], col_ids[n], value) 
            expected.append(value)
            computed.append(dem.get_value(row_ids[n], col_ids[n]))
            
        self.assertEqual(computed, expected)
        
    def test_save(self):
        # First create a copy of the grid
        dem = Grid()
        dem.copy_layout(self.dem)
        dem.set_data(self.dem.read_array())
        
        # Save it
        dem.save("data/dummy_grid.tif")
        
        # Open it
        dem = Grid("data/dummy_grid.tif")
        
        # Take a value
        expected = self.dem.get_value(25, 25)
        computed = dem.get_value(25, 25)
        
        self.assertEqual(computed, expected)
        

class GridValueTests(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load("data/MY_GRID_100rnd_id.npy")
        self.rows = np.load("data/MY_GRID_100rnd_row.npy")
        self.cols = np.load("data/MY_GRID_100rnd_col.npy")
        self.xi = np.load("data/MY_GRID_100rnd_X.npy")
        self.yi = np.load("data/MY_GRID_100rnd_Y.npy")
        self.zi = np.load("data/MY_GRID_100rnd_Z.npy")
        
        # Create a Grid object
        self.dem = Grid(MY_GRID)  

    def test_get_value01(self):
        # Taking row, col as integer
        ind = np.random.randint(0, 100)
        row, col = self.rows[ind], self.cols[ind]
        expected = self.zi[ind]
        computed = self.dem.get_value(row, col)
        self.assertEqual(computed, expected)
    
    def test_get_value02(self):
         # Taking row, col as numpy arrays
        expected = self.zi
        computed = self.dem.get_value(self.rows, self.cols)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)
        
    def test_get_value03(self):
        # Taking row, col as lists
        expected = self.zi.tolist()
        computed = self.dem.get_value(self.rows.tolist(), self.cols.tolist())
        self.assertEqual(computed, expected)

    def test_xy2cell01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_xy2cell02(self):
        ind = np.random.randint(0, 100)
        x = self.xi[ind]
        y = self.yi[ind]
        row = self.rows[ind]
        col = self.cols[ind]
        expected = (row, col)
        computed = self.dem.xy_2_cell(x, y)
        self.assertEqual(computed, expected)
        
    def test_xy2cell03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows.tolist()
        cols = self.cols.tolist()
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        self.assertEqual(computed, expected)   

    
suite1 = unittest.TestLoader().loadTestsFromTestCase(GridPropertyTests)
suite2 = unittest.TestLoader().loadTestsFromTestCase(GridCopySaveTests)
suite3 = unittest.TestLoader().loadTestsFromTestCase(GridValueTests)

suite = unittest.TestSuite([suite1, suite2, suite3])
unittest.TextTestRunner(verbosity=2).run(suite)  