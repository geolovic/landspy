#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 October, 2018
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 04 October, 2018
"""

import unittest
from topopy import Flow, DEM, Network
infolder = "data/in"
outfolder = "data/out"

class NetworkClassTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            dem_path = infolder +  "/{0}.tif".format(file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            # Simplemente probamos que no hay fallos al crear el objeto net
            net = Network(dem, fd, thr)
            net.export_2_points(outfolder + "/{0}_chandata.txt".format(file))



if __name__ == "__main__":
    unittest.main()