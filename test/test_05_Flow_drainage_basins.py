#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 25, 2018
Testing suite for topopy.Flow.get_drainage_basin() function
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

class DrainageBasinTest(unittest.TestCase):
    
    def test_create_Flow_01(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            
            # Creamos 10 cuencas aleatorias 
            ri = np.random.randint(0, fd._dims[0], 10)
            ci = np.random.randint(0, fd._dims[1], 10)
            xi, yi = fd.cell_2_xy(ri, ci)
            # Extract basins
            outlets = np.array((xi, yi)).T
            basins = fd.get_drainage_basins(outlets)
            basins.save(outfolder + "/rnd_basins_{0}".format(file))

    def test_create_Flow_02(self):
        # Test extracting the biggest basin and input as a list [x, y]
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            # Get coords of the highest flow accumulation
            fac = fd.get_flow_accumulation()
            maxval = fac.max()
            rowpos, colpos = np.where(fac.read_array()==maxval)
            xi, yi = fd.cell_2_xy(rowpos, colpos)
            # Extract basin
            xi, yi = xi[0], yi[0]
            outlets = [xi, yi]
            basins = fd.get_drainage_basins(outlets)
            basins.save(outfolder + "/max_basin_{0}".format(file))       

    def test_create_Flow_03(self):
        # Test extracting the biggest basin and input as a list [x, y]
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            # Get coords of the highest flow accumulation
            fac = fd.get_flow_accumulation()
            maxval = fac.max()
            rowpos, colpos = np.where(fac.read_array()==maxval)
            xi, yi = fd.cell_2_xy(rowpos, colpos)
            # Extract basin
            outlets = np.array((xi, yi)).T
            basins = fd.get_drainage_basins(outlets)
            basins.save(outfolder + "/max_basin2_{0}".format(file))
            
    def test_create_Flow_04(self):
        # Test extracting all basins (min_area = 10%)
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            # Extract all basin with min_area
            basins = fd.get_drainage_basins()
            basins.save(outfolder + "/all_basins-minarea_{0}".format(file))

    def test_create_Flow_05(self):
        # Test extracting all basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            # Extract all basin without min_area
            basins = fd.get_drainage_basins(min_area = 0)
            basins.save(outfolder + "/all_basins{0}".format(file))
   
if __name__ == "__main__":
    unittest.main()