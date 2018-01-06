#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""
import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid
from test_Grid_01 import GridPropertyTests
from test_Grid_02 import GridCopySaveTests
from test_Grid_03 import GridValueTests


