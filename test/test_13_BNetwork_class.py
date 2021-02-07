#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 march, 2020 (COVID19 quarantine)
Testing suite for BNetwork class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 26 march, 2020 (COVID19 quarantine)
"""

import unittest
import numpy as np
from topopy import Flow, Basin, Network, BNetwork, DEM
infolder = "data/in"
outfolder = "data/out"

class BNetworkClassTest(unittest.TestCase):
    
    def test_BNetworkClass(self):
        files = ["small25", "jebja30", "tunez"]
        print("Executing test for BNetwork class...")
        for file in files:
            # Cargamos DEM, Flow, Network
            fd = Flow("{}/{}_fd.tif".format(infolder, file))
            net = Network("{}/{}_network.net".format(infolder, file))
            dem = DEM("{}/{}.tif".format(infolder, file))
            
            # Generamos todas las cuencas
            cuencas = fd.get_drainage_basins(min_area = 0.0025)
            
            # Generamos 20 puntos aleatorios dentro de la extensión del objeto Network
            # Estos 20 puntos se usara como cabeceras
            xmin, xmax, ymin, ymax = net.get_extent()
            xi = np.random.randint(xmin, xmax, 50)
            yi = np.random.randint(ymin, ymax, 50)
            heads = np.array((xi, yi)).T
            
            # Cogemos 5 cuencas aleatorias
            bids = np.random.choice(np.unique(cuencas.read_array())[1:], 5)
            print("Test with DEM: ", file)
            for bid in bids:
                if np.random.randint(101) < 70:
                    print("Creating BNetwork from grid, bid={}".format(bid))
                    bnet = BNetwork(net, cuencas, heads, bid)
                else:
                    print("Creating BNetwork from Basin, bid={}".format(bid))
                    basin = Basin(dem, cuencas, bid)
                    bnet = BNetwork(net, basin, heads, bid)
                
                bnet_path = "{}/{}_{}_bnetwork.net".format(outfolder, file, bid)
                bnet.save(bnet_path)
                bnet2 = BNetwork(bnet_path)
                
                computed = np.array_equal(bnet._ix, bnet2._ix)
                self.assertEqual(computed, True)  
        

if __name__ == "__main__":
    unittest.main()