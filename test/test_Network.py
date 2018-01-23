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
import scipy.io as sio
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import DEM, Network


class PresillTest(unittest.TestCase):
    
    def test_presills_01(self):
        # Open DEM and get flats, sills and presills
        dem = DEM('data/tunez.tif')
        flats, sills = dem.identify_flats(False)
        fd = Network()
        presills = fd._get_presills(flats, sills, dem)
        # Get MatLab Presills
        mlab_file = 'data/mlab_files/presill_tunez.mat'
        mpresills = sio.loadmat(mlab_file)['PreSillPixel']
        mpresills = mpresills.reshape((mpresills.shape[0],)) - 1
        row, col = np.unravel_index(mpresills, dem.get_dims(), "F")
        mpresills = np.zeros(dem.get_dims(), np.int8)
        mpresills[row, col] = 1
        computed = np.array_equal(presills, mpresills)
        self.assertEqual(computed, True)
    
    
if __name__ == "__main__":
    unittest.main()