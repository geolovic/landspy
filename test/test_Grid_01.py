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

    def test_set_nodata(self):
        dem = Grid()
        dem.set_data(np.arange(25).reshape((5, 5)))
        dem.set_nodata_value(-99.)
        expected = (-99, int)
        computed = (dem.get_nodata(), type(dem.get_nodata()))
        self.assertEqual(computed, expected)
        
    def test_set_nodata2(self):
        dem = Grid()
        dem.set_data(np.arange(25).reshape((5, 5)))
        dem.set_nodata_value(-99.)
        dem.set_nodata(10)
        expected = -99
        row, col = dem.ind_2_cell(10)
        computed = dem.get_value(row, col)
        self.assertEqual(computed, expected)
        
    def test_set_nodata3(self):
        dem = Grid()
        dem.set_data(np.arange(25).reshape((5, 5)))
        dem.set_nodata_value(-99.)
        dem.set_nodata([10, 11, 12, 13, 14])
        row, col = dem.ind_2_cell([10, 11, 12, 13, 14])
        computed = dem.get_value(row, col)
        expected = [-99, -99, -99, -99, -99]
        self.assertEqual(computed, expected)
        
if __name__ == "__main__":
    unittest.main()