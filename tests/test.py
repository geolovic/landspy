#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 19 september, 2022
"""
import unittest

from test_00_PRaster import TestPRaster00, TestPRaster01
from test_01_Grid import TestGrid01
from test_02_DEM import DEM_class, DEMFlatTest
from test_03_Flow import FlowValueTest, FlowPropertiesTest, FlowCreateTest
from test_03b_Basin import BasinClassTest
from test_04_Flow_stream_poi import StreamPoiTest
from test_05_Flow_snap_poi import SnapPoiTest
from test_06_Flow_drainage_basins import DrainageBasinTest
from test_07_Network import NetworkClassTest
from test_08_Network_snap_poi import SnapPoiTest as netSnapPoiTest
from test_09_Network_stream_poi import StreamPoiTest as netStreamPoiTest
from test_10_Network_get_streams import NetworkGetStreams
from test_11_Network_export_to_shp import NetworkExportToShp
from test_12_Network_get_chi_shapefile import NetworkGetChiShapefile
from test_13_BNetwork_class import BNetworkClassTest
from test_13_Network_get_channel import NetworkGetChannelTest
from test_14_BNetwork_main_channel import BNetworkGetMainChannelTest
from test_14_Channel_class import ChannelClassTest
from test_15_Channels_save_load import ChannelSaveTest

# Create output directory if does not exits
import os
if not os.path.exists("data/out"):
    os.mkdir("data/out")

# Grid DEM tests
suite1 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster00)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster01)
suite3 = unittest.TestLoader().loadTestsFromTestCase(TestGrid01)
suite4 = unittest.TestLoader().loadTestsFromTestCase(DEM_class)
suite5 = unittest.TestLoader().loadTestsFromTestCase(DEMFlatTest)

# Flow Tests
suite6 = unittest.TestLoader().loadTestsFromTestCase(FlowValueTest)
suite7 = unittest.TestLoader().loadTestsFromTestCase(FlowPropertiesTest)
suite8 = unittest.TestLoader().loadTestsFromTestCase(FlowCreateTest)
suite9 = unittest.TestLoader().loadTestsFromTestCase(StreamPoiTest)
suite10 = unittest.TestLoader().loadTestsFromTestCase(SnapPoiTest)
suite11 = unittest.TestLoader().loadTestsFromTestCase(DrainageBasinTest)
suite12 = unittest.TestLoader().loadTestsFromTestCase(BasinClassTest)

# Network Tests
suite13 = unittest.TestLoader().loadTestsFromTestCase(NetworkClassTest)
suite14 = unittest.TestLoader().loadTestsFromTestCase(netSnapPoiTest)
suite15 = unittest.TestLoader().loadTestsFromTestCase(netStreamPoiTest)
suite16 = unittest.TestLoader().loadTestsFromTestCase(NetworkGetStreams)
suite17 = unittest.TestLoader().loadTestsFromTestCase(NetworkExportToShp)
suite18 = unittest.TestLoader().loadTestsFromTestCase(NetworkGetChiShapefile)

# BNetwork and Channel Tests
suite19 = unittest.TestLoader().loadTestsFromTestCase(BNetworkClassTest)
suite20 = unittest.TestLoader().loadTestsFromTestCase(NetworkGetChannelTest)
suite21 = unittest.TestLoader().loadTestsFromTestCase(BNetworkGetMainChannelTest)
suite22 = unittest.TestLoader().loadTestsFromTestCase(ChannelClassTest)
suite23 = unittest.TestLoader().loadTestsFromTestCase(ChannelSaveTest)

# Running tests
suite = unittest.TestSuite([suite1, suite2, suite3, suite4, suite5, suite6, suite7,
                            suite8, suite9, suite10, suite11, suite12, suite13, suite14,
                            suite15, suite16, suite17, suite18, suite19, suite20, suite21, 
                            suite22, suite23])

unittest.TextTestRunner(verbosity=2).run(suite)

# Clean output directory
for ff in os.listdir("data/out"):
    os.remove("data/out/{}".format(ff))
    