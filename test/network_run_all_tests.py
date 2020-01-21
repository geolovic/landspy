#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05 November, 2018
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 21 January, 2020
"""
import unittest

from test_07_Network import NetworkClassTest
from test_08_Network_snap_poi import SnapPoiTest
from test_09_Network_stream_poi import StreamPoiTest
from test_10_Network_get_streams import NetworkGetStreams
from test_11_Network_export_to_shp import NetworkExportToShp
from test_12_Network_get_chi_shapefile import NetworkGetChiShapefile
from test_13_Basin import GridBasinClassTest

suite1 = unittest.TestLoader().loadTestsFromTestCase(NetworkClassTest)
suite2 = unittest.TestLoader().loadTestsFromTestCase(SnapPoiTest)
suite3 = unittest.TestLoader().loadTestsFromTestCase(StreamPoiTest)
suite4 = unittest.TestLoader().loadTestsFromTestCase(NetworkGetStreams)
suite5 = unittest.TestLoader().loadTestsFromTestCase(NetworkExportToShp)
suite6 = unittest.TestLoader().loadTestsFromTestCase(NetworkGetChiShapefile)
suite7 = unittest.TestLoader().loadTestsFromTestCase(GridBasinClassTest)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4, 
                            suite5, suite6, suite7])
unittest.TextTestRunner(verbosity=2).run(suite)