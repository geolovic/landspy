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
from topopy import DEM

MY_GRID = "data/small25.tif"

class DEMTests(unittest.TestCase):
    
    def setUp(self):               
        # Create a DEM object
        self.dem = DEM(MY_GRID)  

    def test_fill_01(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 17, 20],
                        [19, 18, 31, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.set_array(arr)
        fill = dem.fill_sinks()
        computed = fill.read_array().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 31, 39],
                        [19, 29, 42, 51]])
        expected = arr.tolist()
        
        self.assertEqual(computed, expected)
    

        
    def test_fill_02(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 17, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.set_array(arr)
        fill = dem.fill_sinks()
        computed = fill.read_array().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, 51]])
        
        expected = arr.tolist()
    
        self.assertEqual(computed, expected)
    


if __name__ == "__main__":
    unittest.main()