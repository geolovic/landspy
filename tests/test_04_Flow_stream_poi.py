#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 25, 2018
Testing suite for landspy.Flow.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: September 25, 2018
@last_modified: 19 september, 2022
"""

import unittest
import numpy as np
from landspy import Flow

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

class StreamPoiTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.getNCells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.streamPoi(thr, kind, "XY")
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)

    def test_stream_poi_02(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.getNCells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.streamPoi(thr, kind, "CELL")
                x, y = fd.cellToXY(poi[:, 0], poi[:, 1])
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
                
    def test_stream_poi_03(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.getNCells() * 0.0025)
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.streamPoi(thr, kind, "IND")
                row, col = fd.indToCell(poi)
                x, y = fd.cellToXY(row, col)
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True) 

if __name__ == "__main__":
    unittest.main()