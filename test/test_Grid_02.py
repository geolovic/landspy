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
        
if __name__ == "__main__":
    unittest.main()