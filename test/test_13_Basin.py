#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 February, 2020
Testing suite for Network.get_chi_shapefile() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 11 February, 2020
"""

import unittest
from topopy import Grid, Basin
infolder = "data/in"
outfolder = "data/out"

class BasinClassTest(unittest.TestCase):
    
    def test_BasinClass01(self):
        dem_path = infolder + "/jebja30.tif"
        basin_path = infolder + "/jebja30_basins.tif"
        
        basin = Basin(dem_path, basin_path)
        basin.save(outfolder + "/jebja30_basin01.tif")
        
    def test_BasinClass02(self):
        dem_path = infolder + "/jebja30.tif"
        basin_path = infolder + "/jebja30_basins.tif"
        
        basingrid = Grid(basin_path)   
        basin = Basin(dem_path, basingrid)
        basin.save(outfolder + "/jebja30_basin02.tif")
       
    def test_BasinClass03(self):
        dem_path = infolder + "/jebja30_basin.tif"  
        basin = Basin(dem_path)
        basin.save(outfolder + "/jebja30_basin03.tif")
        
    def test_BasinClass04(self):
        dem_path = infolder + "/tunez.tif"
        basin_path = infolder + "/tunez_basins.tif"
        basin = Basin(dem_path, basin_path, idx=15)
        basin.save(outfolder + "/tunez_basin01.tif")
            
if __name__ == "__main__":
    unittest.main()