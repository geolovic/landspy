#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 October, 2018
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 23 October, 2018
"""

import unittest
from topopy import Network
infolder = "data/in"
outfolder = "data/out"

class NetworkExportToShp(unittest.TestCase):
    
    def test_get_streams(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_network.net".format(file)
            net = Network(net_path)
            
            # Exportamos canales a shapefile
            out_shp = outfolder +  "/{0}_str.shp".format(file)
            out_con_shp = outfolder +  "/{0}_strc.shp".format(file)
            net.export_to_shp(out_shp)
            net.export_to_shp(out_con_shp, True)
            

if __name__ == "__main__":
    unittest.main()