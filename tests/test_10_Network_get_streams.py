#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 October, 2018
Testing suite for Network get_streams functions
@author: J. Vicente Perez
@email: geolovic@gmail.com
@date: 08 October, 2018
@last_modified: 15 october, 2025
"""

import unittest
from landspy import Network

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

class NetworkGetStreams(unittest.TestCase):
    
    def test_get_streams(self):
        files = ["small25", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = outfolder +  "/{0}_net.dat".format(file)
            net = Network(net_path)
            
            # Probamos que no haya fallos en las funciones de get_streams
            streams = net.getStreams()
            streams.save(outfolder + "/{0}_streams.tif".format(file))
            computed = [os.path.exists(outfolder + "/{0}_streams.tif".format(file))]
            segments = net.getStreamSegments()
            segments.save(outfolder + "/{0}_segments.tif".format(file))
            computed.append(os.path.exists(outfolder + "/{0}_segments.tif".format(file)))
            orders = net.getStreamOrders()
            orders.save(outfolder + "/{0}_ord.tif".format(file))
            computed.append(os.path.exists(outfolder + "/{0}_ord.tif".format(file)))
            
            self.assertEqual(computed, [True, True, True])

if __name__ == "__main__":
    unittest.main()