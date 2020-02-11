#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21 January, 2020
Testing suite for Network.get_chi_shapefile() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 21 January, 2020
"""

import unittest
from topopy import DEM, Grid, Basin
infolder = "data/in"
outfolder = "data/out"

class GridBasinClassTest(unittest.TestCase):
    
    def test_BasinClass(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        basin_idx = [(3, 5, 8, 6), (4, 10, 12, 21), (3, 12, 13, 14), (3, 4, 6, 12)]
        
        for n, file in enumerate(files):
            # Cargamos cuencas y DEM
            basin_path = outfolder +  "/all_basins-minarea_{0}.tif".format(file)
            basins = Grid(basin_path)
            dem_path = infolder +  "/{0}.tif".format(file)
            basin_ids = basin_idx[n]
            
            for n, idbasin in enumerate(basin_ids):
                # Extract basins
                basin = Basin(dem_path, basins, idbasin)
                out_path = outfolder +  "/{0}_basin{1}.tif".format(file, n)
                basin.save(out_path)
            
            
if __name__ == "__main__":
    unittest.main()