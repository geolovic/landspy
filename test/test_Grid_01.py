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
    
    def test_empty_grid(self):
        # Test que crea un objeto grid vacio
        dem = Grid()
        computed = (dem.get_size(), dem.get_geotransform(), dem.get_cellsize(),
                    dem.get_projection(), dem.get_nodata(), dem._array, str(dem._tipo))
        
        expected = ((1, 1), (0., 1., 0., 0., 0., -1.), 1., "", None, np.array([[0]], dtype="float"), "float64")
        self.assertEqual(computed, expected)
    
    def test_create_grid(self):
        dem = Grid(MY_GRID)
        computed = (dem.get_size(), dem.get_geotransform(), dem.get_cellsize(),
                    dem.get_projection(), dem.get_nodata(), str(dem._tipo))
        proj_wkt = 'PROJCS["ETRS89 / UTM zone 30N",GEOGCS["ETRS89",DATUM["European_Terrestrial_Reference_System_1989",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6258"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4258"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","25830"]]'
        expected = ((199, 133), (470139.6, 25.0, 0.0, 4117136.0, 0.0, -25.0), 25., proj_wkt, -9999., "int16")
        self.assertEqual(computed, expected)

    def test_set_data_00(self):
        arr = np.arange(9).reshape((3, 3))
        arr[[1, 2], [2, 1]] = 8
        dem = Grid()
        dem.set_nodata(8)
        dem.set_array(arr)
        computed = (np.isnan(dem._array[1, 2]), np.isnan(dem._array[2, 1]))
        expected = (True, True)
        self.assertEqual(computed, expected)
    
    def test_read_array_00(self):
        # read_array(nans = True)
        arr = np.arange(9).reshape((3, 3))
        arr[[1, 2], [2, 1]] = 8
        dem = Grid()
        dem.set_array(arr)
        dem.set_nodata(8)
        dem_arr = dem.read_array()
        computed = (np.isnan(dem_arr[1, 2]), np.isnan(dem_arr[2, 1]))
        expected = (True, True)
        self.assertEqual(computed, expected)
        
    def test_read_array_01(self):
        # read_array(nans = False)
        arr = np.arange(9).reshape((3, 3))
        arr[[1, 2], [2, 1]] = 8
        dem = Grid()
        dem.set_array(arr)
        dem.set_nodata(8)
        dem_arr = dem.read_array(False)
        computed = (dem_arr[1, 2], dem_arr[2, 1], dem_arr[2, 2])
        expected = (8, 8, 8)
        self.assertEqual(computed, expected)
        
    def test_save(self):
        dem = Grid(MY_GRID)
        dem.save("data/dummy_dem.tif")
        dem2 = Grid("data/dummy_dem.tif")
        
        expected = dem._array[[20, 30, 40, 50], [20, 30, 40, 50]]
        computed = dem2._array[[20, 30, 40, 50], [20, 30, 40, 50]]
        self.assertEqual(computed.sum(), expected.sum())
             
    def test_values_2_nodata(self):
        dem = Grid()
        dem.set_array(np.arange(25).reshape((5, 5)))
        dem.set_nodata(-99.)
        dem.values_2_nodata([10, 11, 12, 13, 14])
        computed = np.isnan(dem._array[2]).tolist()
        expected = [True, True, True, True, True]
        self.assertEqual(computed, expected)

 
if __name__ == "__main__":
    unittest.main()