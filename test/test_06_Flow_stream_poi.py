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
import sys, os
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Flow
infolder = "data/in"
outfolder = "data/out"

class StreamPoiTest(unittest.TestCase):
    
    def test_stream_poi_01(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "heads", "CELL")
            row, col = outlets[:, 0], outlets[:, 1]
            xi, yi = fd.cell_2_xy(row, col)
            out_arr = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/heads_CELL_{0}.txt".format(file), out_arr, delimiter=";", header=cabecera, comments="")

    def test_stream_poi_02(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "heads", "XY")
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/heads_XY_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="")            

    def test_stream_poi_03(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "heads", "IND")
            row, col = fd.ind_2_cell(outlets)
            xi, yi = fd.cell_2_xy(row, col)
            outlets = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/heads_IND_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="") 

    def test_stream_poi_04(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "confluences", "CELL")
            row, col = outlets[:, 0], outlets[:, 1]
            xi, yi = fd.cell_2_xy(row, col)
            out_arr = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/confluences_CELL_{0}.txt".format(file), out_arr, delimiter=";", header=cabecera, comments="")

    def test_stream_poi_05(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "confluences", "XY")
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/confluences_XY_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="")            

    def test_stream_poi_06(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "confluences", "IND")
            row, col = fd.ind_2_cell(outlets)
            xi, yi = fd.cell_2_xy(row, col)
            outlets = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/confluences_IND_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="") 

    def test_stream_poi_07(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "outlets", "CELL")
            row, col = outlets[:, 0], outlets[:, 1]
            xi, yi = fd.cell_2_xy(row, col)
            out_arr = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/outlets_CELL_{0}.txt".format(file), out_arr, delimiter=";", header=cabecera, comments="")

    def test_stream_poi_08(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "outlets", "XY")
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/outlets_XY_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="")            

    def test_stream_poi_09(self):
        # Test 10 random basins
        dem_files = ['tunez.tif', 'small25.tif',  "jebja30.tif"]        
        for file in dem_files:
            flw_path = infolder +  "/fd_{0}".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            outlets = fd.get_stream_poi(thr, "outlets", "IND")
            row, col = fd.ind_2_cell(outlets)
            xi, yi = fd.cell_2_xy(row, col)
            outlets = np.array((xi, yi)).T
            cabecera = "x;y"
            file = os.path.splitext(file)[0]
            np.savetxt(outfolder + "/outlets_IND_{0}.txt".format(file), outlets, delimiter=";", header=cabecera, comments="") 
   
   
if __name__ == "__main__":
    unittest.main()