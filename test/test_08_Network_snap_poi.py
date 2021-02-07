#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 3, 2018
Testing suite for topopy.Network.get_stream_poi() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: October 3, 2018
@modified:  07 february, 2021
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Network
infolder = "data/in"
outfolder = "data/out"

class SnapPoiTest(unittest.TestCase):
    
    def test_snap_poi_01(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = net.get_extent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Hacemos snap a los stream poi
            for kind in ["heads", "confluences", "outlets"]:
                poi = net.get_stream_poi(kind, "XY")
                snap_pp = net.snap_points(puntos, kind)
                # Comprobamos que punto esta entre los POI
                for row in snap_pp:
                    res = row in poi
                    self.assertEqual(res, True)

    def test_snap_poi_02(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = net.get_extent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Hacemos snap a celdas de canal
            snap_pp = net.snap_points(puntos, kind="channel")
            row, col = net.xy_2_cell(snap_pp[:, 0], snap_pp[:, 1])
            inds = net.cell_2_ind(row, col)
            # Comprobamos que punto esta entre los POI
            for ind in inds:
                res = ind in net._ix
                self.assertEqual(res, True)


if __name__ == "__main__":
    unittest.main()