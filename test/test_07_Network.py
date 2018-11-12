#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 October, 2018
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 10 October, 2018
"""

import unittest
import numpy as np
from topopy import Flow, Network
infolder = "data/in"
outfolder = "data/out"

class NetworkClassTest(unittest.TestCase):
    
    def test_create_network(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            net_path = infolder + "/{0}_network.net".format(file)
            fd = Flow(flw_path)
            thr = int(fd.get_ncells() * 0.01)
            
            # Creamos objeto network
            net = Network(fd, thr)

            # Cargamos objeto network guardado
            net2 = Network(net_path)
            
            # Comparamos propiedades,
            prop_net = [net._size, net._dims, net._geot, net._cellsize, net._ncells, net._proj, net._threshold]
            prop_net2 = [net2._size, net2._dims, net2._geot, net2._cellsize, net2._ncells, net2._proj, net._threshold]
            self.assertEqual(prop_net, prop_net2)
            
            # Comparamos los datos
            arr1 = np.array((net._ix, net._ixc, net._ax, net._dx, net._zx,
                             net._chi, net._slp, net._ksn, net._r2slp, 
                             net._r2ksn, net._dd))
            arr2 = np.array((net2._ix, net2._ixc, net2._ax, net2._dx, net2._zx,
                             net2._chi, net2._slp, net2._ksn, net2._r2slp, 
                             net2._r2ksn, net2._dd))            
            res = np.array_equal(arr1, arr2)
            self.assertEqual(res, True)


if __name__ == "__main__":
    unittest.main()