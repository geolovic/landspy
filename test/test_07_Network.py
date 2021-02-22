#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 August, 2018
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 11 August, 2018
@modified:  16 february, 2021
"""

import unittest
import numpy as np
from topopy import Flow, Network
import os
infolder = "data/in"
outfolder = "data/out"

class NetworkClassTest(unittest.TestCase):
    
    def test_empty_network(self):
        net = Network()
        
        # Test PRaster properties
        self.assertEqual(net.get_cellsize(), (1.0, -1.0))
        self.assertEqual(net.get_dims(), (1, 1))
        self.assertEqual(net.get_size(), (1, 1))
        self.assertEqual(net.get_extent(), (0.0, 1.0, 0.0, 1.0))
        self.assertEqual(net.get_geotransform(), (0.0, 1.0, 0.0, 1.0, 0.0, -1.0))
        self.assertEqual(net.get_ncells(), 1)
        self.assertEqual(net.get_projection(), "")
        
        # Test PRaster functions
        self.assertEqual(net.cell_2_ind(0, 0), .0)
        self.assertEqual(net.cell_2_xy(0, 0), (0.5, 0.5))
        self.assertEqual(net.xy_2_cell(1, 1), (0, 1))
        self.assertEqual(net.ind_2_cell(0), (0, 0))
        
        # Test saving functions
        path = outfolder + "/net_delete.dat"
        net.save(path)
        self.assertEqual(os.path.exists(path), True)
        if os.path.exists(path):
            os.remove(path)
        
        path = outfolder + "/points_delete.txt"
        net.export_to_points(path)
        self.assertEqual(os.path.exists(path), True)
        if os.path.exists(path):
            os.remove(path)
        
        path = outfolder +"/shp_delete.shp"
        net.export_to_shp(path)
        self.assertEqual(os.path.exists(path), True)
        if os.path.exists(path):
            os.remove(path)
        
        path = outfolder +"/chi_delete.shp"
        net.get_chi_shapefile(path, 0)
        self.assertEqual(os.path.exists(path), True)
        if os.path.exists(path):
            os.remove(path)

        # Test other functions
        self.assertEqual(np.array_equal(net.get_streams(False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.get_stream_segments(False), np.array([0]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.get_stream_orders("strahler", False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.get_stream_orders('shreeve', False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.get_stream_orders('dinosaur', False), np.array([1]).reshape(1, 1)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("heads", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("heads", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("heads", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("confluences", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("confluences", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("confluences", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("outlets", "XY"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("outlets", "CELL"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("outlets", "IND"), np.array([])), True)
        self.assertEqual(np.array_equal(net.get_stream_poi("dinosaur", "dinosaur"), np.array([]).reshape(0, 2)), True)
        self.assertEqual(np.array_equal(net.snap_points(np.array((5, 5)).reshape(1, 2)), np.array([[0.5, 0.5]])), True)
        self.assertEqual(net.is_inside(0, 0), False)
        net.calculate_gradients(0)
    
    def test_create_network(self):
        files = ["small25", "tunez", "jebja30"]
    
        for file in files:
            flw_path = infolder +  "/{0}_fd.tif".format(file)
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