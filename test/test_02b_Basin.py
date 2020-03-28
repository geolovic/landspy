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
import numpy as np
from topopy import Flow, Basin, DEM
infolder = "data/in"
outfolder = "data/out"

class BasinClassTest(unittest.TestCase):
    
    def test_BasinClass01(self):
        files = ["small25", "jebja30", "tunez"]
        for file in files:
            # Cargamos flow
            flw_path = "{}/{}_fd.tif".format(infolder, file)
            dem_path = "{}/{}.tif".format(infolder, file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            # Obtenemos cuencas
            cuencas = fd.get_drainage_basins(min_area = 0.0025)
            # Mayor, menor, y random
            bids, counts = np.unique(cuencas.read_array(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]
                
            idmax = bids[np.argmax(counts)]
            idmin = bids[np.argmin(counts)]
            idrdn = np.random.choice(bids)
            
            idxs = [idmax, idmin, idrdn]
            lbls = ["max", "min", "rdn"]
            
            for n in range(3):
                basin = Basin(dem, cuencas, idxs[n])
                basin_path = "{}/{}_basin{}.tif".format(outfolder, file, lbls[n])
                basin.save(basin_path)
                basin2 = Basin(basin_path)
                computed = np.array_equal(basin.read_array(), basin2.read_array())
                self.assertEqual(computed, True)  

if __name__ == "__main__":
    unittest.main()