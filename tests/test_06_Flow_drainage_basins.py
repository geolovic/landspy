#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 25, 2018
Testing suite for landspy.Flow.get_drainage_basin() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: September 25, 2018
@last_modified: 19 september, 2022
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../src/")
from landspy import Flow
infolder = "data/in"
outfolder = "data/out"

class DrainageBasinTest(unittest.TestCase):
    
    def test_drainage_basins_00(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            
            # Creamos 10 cuencas aleatorias 
            ri = np.random.randint(0, fd.getDims()[0], 10)
            ci = np.random.randint(0, fd.getDims()[1], 10)
            rdn_ids = np.random.randint(100, 700, 10)
            xi, yi = fd.cellToXY(ri, ci)
            # Extract basins
            outlets = np.array((xi, yi, rdn_ids)).T
            threshold = int(fd.getNCells() * 0.05)
            snap_outlets = fd.snapPoints(outlets, threshold, kind="channel")
            snap = np.append(snap_outlets, rdn_ids.reshape(rdn_ids.size, 1), 1)
            basins = fd.drainageBasins(snap)
            
            basins.save(outfolder + "/rnd_snap_basins_{0}.tif".format(file))
    
    def test_drainage_basins_01(self):
        # Test 10 random basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            
            # Creamos 10 cuencas aleatorias con 10 ids aleatorios
            ri = np.random.randint(0, fd.getDims()[0], 10)
            ci = np.random.randint(0, fd.getDims()[1], 10)
            rdn_ids = np.random.randint(100, 700, 10)
            xi, yi = fd.cellToXY(ri, ci)
            # Extract basins
            outlets = np.array((xi, yi, rdn_ids)).T
            basins = fd.drainageBasins(outlets)
            basins.save(outfolder + "/rnd_basins_{0}.tif".format(file))

    def test_drainage_basins_02(self):
        # Test extracting the biggest basin and input as a list [x, y]
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            # Get coords of the highest flow accumulation
            fac = fd.flowAccumulation()
            maxval = fac.max()
            rowpos, colpos = np.where(fac.readArray()==maxval)
            xi, yi = fd.cellToXY(rowpos, colpos)
            # Extract basin
            xi, yi = xi[0], yi[0]
            outlets = [xi, yi]
            basins = fd.drainageBasins(outlets)
            basins.save(outfolder + "/max_basin_{0}.tif".format(file))       

    def test_drainage_basins_03(self):
        # Test extracting the biggest basin and input as a list [x, y]
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            # Get coords of the highest flow accumulation
            fac = fd.flowAccumulation()
            maxval = fac.max()
            rowpos, colpos = np.where(fac.readArray()==maxval)
            xi, yi = fd.cellToXY(rowpos, colpos)
            # Extract basin
            outlets = np.array((xi, yi)).T
            basins = fd.drainageBasins(outlets)
            basins.save(outfolder + "/max_basin2_{0}.tif".format(file))
            
    def test_drainage_basins_04(self):
        # Test extracting all basins (min_area = 10%)
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            # Extract all basin with min_area
            basins = fd.drainageBasins()
            basins.save(outfolder + "/all_basins-minarea_{0}.tif".format(file))

    def test_drainage_basins_05(self):
        # Test extracting all basins
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            fd = Flow(flw_path)
            # Extract all basin without min_area
            basins = fd.drainageBasins(min_area = 0)
            basins.save(outfolder + "/all_basins{0}.tif".format(file))
   
if __name__ == "__main__":
    unittest.main()