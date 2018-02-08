#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thrusday 08 Feb 2018
Testing suite for topopy Flow class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import DEM, Flow

class CreateFlow(unittest.TestCase):
    
    def test_create_load(self):
        # Create Flow object
        dem_files = ['tunez', 'tunez2', 'small25']        
        for filename in dem_files:
            dem = DEM('data/' + filename + ".tif")
            fd = Flow(dem)
            fd2 = Flow('data/fd_' + filename + ".tif")
            computed = np.array_equal(fd2._ix, fd._ix)
            self.assertEqual(computed, True)
           
    def test_flow_properties_01(self):
        dem_files = ['tunez', 'tunez2', 'small25'] 
        for filename in dem_files:
            dem = DEM('data/' + filename + ".tif")
            fd = Flow('data/fd_' + filename + ".tif")
            computed = (fd.get_size(), fd.get_dims(), fd.get_ncells(), fd.get_projection(), fd.get_cellsize(), fd.get_geotransform())
            expected = (dem.get_size(), dem.get_dims(), dem.get_ncells(), dem.get_projection(), dem.get_cellsize(), dem.get_geotransform())
            self.assertEqual(computed, expected)
    
    def test_flow_properties_02(self):
        # Testing an empty Flow object
        dem = DEM()
        fd = Flow()
        computed = (fd.get_size(), fd.get_dims(), fd.get_ncells(), fd.get_projection(), fd.get_cellsize(), fd.get_geotransform())
        expected = (dem.get_size(), dem.get_dims(), dem.get_ncells(), dem.get_projection(), dem.get_cellsize(), dem.get_geotransform())
        self.assertEqual(computed, expected)
        
   
if __name__ == "__main__":
    unittest.main()