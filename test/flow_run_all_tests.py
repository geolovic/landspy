#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 01 October, 2018
"""
import unittest

from test_03_Flow import FlowValueTest, FlowPropertiesTest, FlowCreateTest
from test_04_Flow_stream_poi import StreamPoiTest
from test_05_Flow_snap_poi import SnapPoiTest
from test_06_Flow_drainage_basins import DrainageBasinTest

suite1 = unittest.TestLoader().loadTestsFromTestCase(FlowValueTest)
suite2 = unittest.TestLoader().loadTestsFromTestCase(FlowPropertiesTest)
suite3 = unittest.TestLoader().loadTestsFromTestCase(FlowCreateTest)
suite4 = unittest.TestLoader().loadTestsFromTestCase(StreamPoiTest)
suite5 = unittest.TestLoader().loadTestsFromTestCase(SnapPoiTest)
suite6 = unittest.TestLoader().loadTestsFromTestCase(DrainageBasinTest)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4, suite5, suite6])
unittest.TextTestRunner(verbosity=2).run(suite)