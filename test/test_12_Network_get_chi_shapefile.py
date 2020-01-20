#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20 January, 2020
Testing suite for Network.get_chi_shapefile() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 20 January, 2020
"""

import unittest
from topopy import Network
infolder = "data/in"
outfolder = "data/out"

class NetworkGetChiShapefile(unittest.TestCase):
    
    def test_get_chi_shp(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_network.net".format(file)
            net = Network(net_path)
            
            # Exportamos canales a shapefile
            out_shp = outfolder +  "/{0}_chi.shp".format(file)
            net.get_chi_shapefile(out_shp, 250)
            
if __name__ == "__main__":
    unittest.main()