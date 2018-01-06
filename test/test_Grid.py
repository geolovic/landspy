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
from test_Grid_02 import GridCopySaveTests
from test_Grid_03 import GridValueTests
from test_DEM_fill import DEMTests

suite1 = unittest.TestLoader().loadTestsFromTestCase(GridPropertyTests)
suite2 = unittest.TestLoader().loadTestsFromTestCase(GridCopySaveTests)
suite3 = unittest.TestLoader().loadTestsFromTestCase(GridValueTests)
suite4 = unittest.TestLoader().loadTestsFromTestCase(DEMTests)

suite = unittest.TestSuite([suite1, suite2, suite3, suite4])
unittest.TextTestRunner(verbosity=2).run(suite)  