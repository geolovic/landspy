#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""
import unittest

from test_Grid import TestGrid01
from test_PRaster import TestPRaster00, TestPRaster01
from test_DEM import DEMFillTest, DEMFlatTest

suite1 = unittest.TestLoader().loadTestsFromTestCase(TestGrid01)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster00)
suite3 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster01)
suite4 = unittest.TestLoader().loadTestsFromTestCase(DEMFillTest)
suite5 = unittest.TestLoader().loadTestsFromTestCase(DEMFlatTest)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4, suite5])
unittest.TextTestRunner(verbosity=2).run(suite)