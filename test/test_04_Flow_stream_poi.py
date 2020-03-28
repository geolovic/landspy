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
from topopy import Flow
infolder = "data/in"
outfolder = "data/out"

class StreamPoiTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.get_stream_poi(thr, kind, "XY")
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)

    def test_stream_poi_02(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.get_stream_poi(thr, kind, "CELL")
                x, y = fd.cell_2_xy(poi[:, 0], poi[:, 1])
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
                
    def test_stream_poi_03(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.get_stream_poi(thr, kind, "IND")
                row, col = fd.ind_2_cell(poi)
                x, y = fd.cell_2_xy(row, col)
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True) 

if __name__ == "__main__":
    unittest.main()