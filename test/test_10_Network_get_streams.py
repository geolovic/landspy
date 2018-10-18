#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 October, 2018
Testing suite for Network get_streams functions
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 08 October, 2018
"""

import unittest
from topopy import Network
infolder = "data/in"
outfolder = "data/out"

class NetworkGetStreams(unittest.TestCase):
    
    def test_get_streams(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            # Cargamos objeto network guardado previamente
            net_path = infolder +  "/{0}_network.net".format(file)
            net = Network(net_path)
            
            # Probamos que no haya fallos en las funciones de get_streams
            streams = net.get_streams()
            streams.save(outfolder + "/{0}_streams.tif".format(file))
            segments = net.get_stream_segments()
            segments.save(outfolder + "/{0}_segments.tif".format(file))
            orders = net.get_stream_order()
            orders.save(outfolder + "/{0}_ord.tif".format(file))

if __name__ == "__main__":
    unittest.main()