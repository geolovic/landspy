#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 25, 2018
Testing suite for topopy.Flow.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: September 25, 2018
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Flow, DEM, Network
infolder = "data/in"

class StreamPoiTest(unittest.TestCase):
    
#    def test_stream_poi_01(self):
#        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]
#        for file in dem_files:
#            dem = DEM(infolder + "/" + file)
#            fd = Flow(dem)
#            thr = int(fd.get_ncells() * 0.01)
#            net = Network(fd, dem, thr)
#            
#            out01 = fd.get_stream_poi(thr, "heads", "CELL")
#            out02 = net.get_stream_poi("heads", "CELL")
#            
#            computed = np.array_equal(out01, out02)
#            self.assertEqual(computed, True)
            
    def test_stream_poi_02(self):
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]
        for file in dem_files:
            dem = DEM(infolder + "/" + file)
            fd = Flow(dem)
            thr = int(fd.get_ncells() * 0.01)
            net = Network(fd, dem, thr)
            
            out01 = fd.get_stream_poi(thr, "confluences", "CELL")
            out02 = net.get_stream_poi("confluences", "CELL")
            
            computed = np.array_equal(out01, out02)
            print(file)
            self.assertEqual(computed, True)

#            
#    def test_stream_poi_03(self):
#        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]
#        for file in dem_files:
#            dem = DEM(infolder + "/" + file)
#            fd = Flow(dem)
#            thr = int(fd.get_ncells() * 0.01)
#            net = Network(fd, dem, thr)
#            
#            out01 = fd.get_stream_poi(thr, "outlets", "CELL")
#            out02 = net.get_stream_poi("outlets", "CELL")
#            
#            computed = np.array_equal(out01, out02)
#            self.assertEqual(computed, True)


if __name__ == "__main__":
    unittest.main()