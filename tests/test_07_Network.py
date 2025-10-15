#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 August, 2018
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 11 August, 2018
@last_modified: 19 september, 2022
"""

import unittest
import numpy as np
from landspy import Flow, Network

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

class NetworkClassTest(unittest.TestCase):
    
    def test_empty_network(self):
        net = Network()
        
        # Test PRaster properties
        self.assertEqual(net.getCellSize(), (1.0, -1.0))
        self.assertEqual(net.getDims(), (1, 1))
        self.assertEqual(net.getSize(), (1, 1))
        self.assertEqual(net.getExtent(), (0.0, 1.0, 0.0, 1.0))
        self.assertEqual(net.getGeot(), (0.0, 1.0, 0.0, 1.0, 0.0, -1.0))
        self.assertEqual(net.getNCells(), 1)
        self.assertEqual(net.getCRS(), "")
        
        # Test PRaster functions
        self.assertEqual(net.cellToInd(0, 0), .0)
        self.assertEqual(net.cellToXY(0, 0), (0.5, 0.5))
        self.assertEqual(net.xyToCell(1, 1), (0, 1))
        self.assertEqual(net.indToCell(0), (0, 0))
        
        # Test saving functions
        path = outfolder + "/net_delete.dat"
        net.save(path)
        self.assertEqual(os.path.exists(path), True)
        
        path = outfolder + "/points_delete.txt"
        net.exportPoints(path)
        self.assertEqual(os.path.exists(path), True)

        # Test other functions
        self.assertEqual(np.array_equal(net.getStreams(False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.getStreamSegments(False), np.array([0]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.getStreamOrders("strahler", False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.getStreamOrders('shreeve', False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.getStreamOrders('dinosaur', False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.streamPoi("heads", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("heads", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("heads", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.streamPoi("confluences", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("confluences", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("confluences", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.streamPoi("outlets", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("outlets", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.streamPoi("outlets", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.streamPoi("dinosaur", "dinosaur"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.snapPoints(np.array((5, 5)).reshape(1, 2)), np.array([[0.5, 0.5]])), True)
        self.assertEqual(net.isInside(0, 0), False)
        net.calculateGradients(0)
    
    def test_create_network(self):
        files = ["small25", "tunez", "jebja30"]
    
        for file in files:
            flw_path = outfolder +  "/{0}_fd.tif".format(file)
            net_path = outfolder + "/{0}_net.dat".format(file)
            fd = Flow(flw_path)
            
            # Creamos objeto network
            net = Network(fd, gradients=True)
            # Guardamos objeto network y cargamos
            net.save(net_path)
            net2 = Network(net_path)
            
            # Comparamos propiedades,
            prop_net = [net._size, net._geot, net._proj]
            prop_net2 = [net2._size, net2._geot, net2._proj]
            self.assertEqual(prop_net, prop_net2)
            
            # Comparamos los datos
            arr1 = np.array((net._ix, net._ixc, net._ax, net._dx, net._zx,
                              net._chi, net._slp, net._ksn, net._dd))
            arr2 = np.array((net2._ix, net2._ixc, net2._ax, net2._dx, net2._zx,
                              net2._chi, net2._slp, net2._ksn, net2._dd))            
            res = np.array_equal(arr1, arr2)
            self.assertEqual(res, True)



if __name__ == "__main__":
    unittest.main()