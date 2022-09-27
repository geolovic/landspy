#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 3, 2018
Testing suite for landspy.Network.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: October 3, 2018
@last_modified: 19 september, 2022
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../src/")
from landspy import Network
infolder = "data/in"
outfolder = "data/out"

class SnapPoiTest(unittest.TestCase):
    
    def test_snap_poi_01(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = net.getExtent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Hacemos snap a los stream poi
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.streamPoi(kind, "XY")
                snap_pp = net.snapPoints(puntos, kind)
                # Comprobamos que punto esta entre los POI
                for row in snap_pp:
                    res = row in poi
                    self.assertEqual(res, True)

    def test_snap_poi_02(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = net.getExtent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Hacemos snap a celdas de canal
            snap_pp = net.snapPoints(puntos, kind="channel")
            row, col = net.xyToCell(snap_pp[:, 0], snap_pp[:, 1])
            inds = net.cellToInd(row, col)
            # Comprobamos que punto esta entre los POI
            for ind in inds:
                res = ind in net._ix
                self.assertEqual(res, True)


if __name__ == "__main__":
    unittest.main()