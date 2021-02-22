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

from test_00_PRaster import TestPRaster00, TestPRaster01
from test_01_Grid import TestGrid01
from test_02_DEM import DEM_class, DEMFlatTest
from test_02b_Basin import BasinClassTest

suite1 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster00)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestPRaster01)
suite3 = unittest.TestLoader().loadTestsFromTestCase(TestGrid01)
suite4 = unittest.TestLoader().loadTestsFromTestCase(DEM_class)
suite5 = unittest.TestLoader().loadTestsFromTestCase(DEMFlatTest)
suite6 = unittest.TestLoader().loadTestsFromTestCase(BasinClassTest)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4, suite5, suite6])
unittest.TextTestRunner(verbosity=2).run(suite)