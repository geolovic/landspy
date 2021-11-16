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

class SnapPoiTest(unittest.TestCase):
    
    def test_snap_poi_01(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.0025)
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = fd.get_extent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Hacemos snap a los stream poi
            for kind in ["heads", "confluences", "outlets"]:
                poi = fd.get_stream_poi(thr, kind, "XY")
                snap_pp = fd.snap_points(puntos, thr, kind)
                # Comprobamos que punto esta entre los POI
                for row in snap_pp:
                    res = row in poi
                    self.assertEqual(res, True)


    def test_snap_poi_02(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.0025)
            # Obtenemos 20 puntos aleatorios
            x1, x2, y1, y2 = fd.get_extent()
            xi = (x2 - x1) * np.random.random(25) + x1
            yi = (y2 - y1) * np.random.random(25) + y1
            puntos = np.array((xi, yi)).T
            
            # Obtenemos lista con Flow Accumulations cells
            fac = fd.get_flow_accumulation(nodata=False, asgrid=False)
            ch_cells = np.where(fac.ravel() >= thr)[0]
            
            # Hacemos snap a celdas de canal
            snap_pp = fd.snap_points(puntos, thr, kind="channel")
            row, col = fd.xy_2_cell(snap_pp[:, 0], snap_pp[:, 1])
            inds = fd.cell_2_ind(row, col)
            # Comprobamos que punto esta entre los POI
            for ind in inds:
                res = ind in ch_cells
                self.assertEqual(res, True)


if __name__ == "__main__":
    unittest.main()