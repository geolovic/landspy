#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 October, 2018
Testing suite for landspy.Network.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 02 October, 2018
@last_modified: 19 september, 2022
"""

import unittest
import numpy as np
import sys
# Add to the path code folder and data folder
sys.path.append("../src/")
from landspy import Network
infolder = "data/in"
outfolder = "data/out"


class StreamPoiTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.streamPoi(kind, "XY")
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)

    def test_stream_poi_02(self):    
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.streamPoi(kind, "CELL")
                x, y = net.cellToXY(poi[:, 0], poi[:, 1])
                poi = np.array((x, y)).T
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
    
    def test_stream_poi_03(self):                 
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.streamPoi(kind, "IND")
                row, col = net.indToCell(poi)
                x, y = net.cellToXY(row, col)
                poi = np.array((x, y)).T                 
                spoi = np.loadtxt(infolder +  "/{0}_{1}.txt".format(file, kind), delimiter=";", skiprows=1)
                compare = np.array_equal(poi, spoi)
                self.assertEqual(compare, True)
#
if __name__ == "__main__":
    unittest.main()