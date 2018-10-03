#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 October, 2018
Testing suite for topopy.Network.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 02 October, 2018
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Flow, DEM, Network
infolder = "data/in"
outfolder = "data/out"

class StreamPoiTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            dem_path = infolder +  "/{0}.tif".format(file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            net = Network(dem, fd, thr)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.get_stream_poi(kind, "XY")
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)

    def test_stream_poi_02(self):    
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            dem_path = infolder +  "/{0}.tif".format(file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            net = Network(dem, fd, thr)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.get_stream_poi(kind, "CELL")
                x, y = fd.cell_2_xy(poi[:, 0], poi[:, 1])
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
    
    def test_stream_poi_03(self):                 
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            dem_path = infolder +  "/{0}.tif".format(file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            net = Network(dem, fd, thr)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.get_stream_poi(kind, "IND")
                row, col = fd.ind_2_cell(poi)
                x, y = fd.cell_2_xy(row, col)
                poi = np.array((x, y)).T                 
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
#
if __name__ == "__main__":
    unittest.main()