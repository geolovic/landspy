#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20 January, 2020
Testing suite for Network.get_chi_shapefile() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 20 January, 2020
@modified:  07 february, 2021
"""

import unittest, os
from topopy import Network
infolder = "data/in"
outfolder = "data/out"

class NetworkGetChiShapefile(unittest.TestCase):
    
    def test_get_chi_shp(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Exportamos canales a shapefile
            out_shp = outfolder +  "/{0}_chi.shp".format(file)
            net.get_chi_shapefile(out_shp, 250)
            
            computed = os.path.exists(out_shp)
            self.assertEqual(computed, True)
            
            
if __name__ == "__main__":
    unittest.main()