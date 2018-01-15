#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""
import unittest

from test_Grid_01 import GridPropertyTests
from test_Grid_02 import GridValueTests
from test_DEM_fill import DEMFillTest
from test_DEM_flats import DEMFlatTest
from test_distance_cython import DistanceTest, CostTest

suite1 = unittest.TestLoader().loadTestsFromTestCase(GridPropertyTests)
suite2 = unittest.TestLoader().loadTestsFromTestCase(GridValueTests)
suite3 = unittest.TestLoader().loadTestsFromTestCase(DEMFillTest)
suite4 = unittest.TestLoader().loadTestsFromTestCase(DEMFlatTest)
suite5 = unittest.TestLoader().loadTestsFromTestCase(DistanceTest)
suite6 = unittest.TestLoader().loadTestsFromTestCase(CostTest)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4, suite5, suite6])
unittest.TextTestRunner(verbosity=2).run(suite)