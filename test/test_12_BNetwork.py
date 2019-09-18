#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18 September, 2019
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 18 September, 2019
"""

import unittest
from topopy import Network, BNetwork, extract_points, Flow
import numpy as np
infolder = "data/in"
outfolder = "data/out"

class BNetwork_test(unittest.TestCase):
    
    def test_create_BNetwork(self):
        # Cargamos objeto network guardado previamente
        net_path = infolder +  "/tunez_network.net"
        fld_path = infolder + "/tunez_fd.tif"
        net = Network(net_path)
        fld = Flow(fld_path)  
        
        outlets = extract_points(infolder + "/tunez_outlets.shp", "id")
        heads = extract_points(infolder + "/tunez_heads.shp", "id")
        basins = fld.get_drainage_basins(outlets)
        
        for n in np.unique(basins._array):
            if n == 0: continue
            bnet = BNetwork(net, basins, heads, n)
            canales = bnet.get_streams()
            canales.save(outfolder + "/{}_basin.tif".format(str(n).zfill(2)))
         

if __name__ == "__main__":
    unittest.main()