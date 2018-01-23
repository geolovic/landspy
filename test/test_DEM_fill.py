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
from topopy import DEM

class DEMFillTest(unittest.TestCase):

    def test_fill_01(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 17, 20],
                        [19, 18, 31, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.set_array(arr)
        fill = dem.fill_sinks()
        computed = fill.read_array().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 31, 39],
                        [19, 29, 42, 51]])
        expected = arr.tolist()
        
        self.assertEqual(computed, expected)
    
  
    def test_fill_02(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 17, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.set_array(arr)
        fill = dem.fill_sinks()
        computed = fill.read_array().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, 51]])
        
        expected = arr.tolist()
    
        self.assertEqual(computed, expected)
      
    def test_fill_03(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 19, 39],
                        [19, 29, 42, -99]])
        dem = DEM()
        dem.set_nodata(-99)
        dem.set_array(arr)
        fill = dem.fill_sinks()
        computed = fill.read_array(False).tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, -99]])
        
        expected = arr.tolist()
    
        self.assertEqual(computed, expected)
          
    def test_fill_04(self):
        dem = DEM("data/small25.tif")
        fill = dem.fill_sinks().read_array().astype("int16")
        
        mfill = sio.loadmat('data/mlab_files/fill_small25.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.get_nodata_pos()
        mfill[nodatapos] = dem.get_nodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        
    def test_fill_05(self):
        dem = DEM("data/tunez.tif")
        fill = dem.fill_sinks().read_array().astype("int16")
        
        mfill = sio.loadmat('data/mlab_files/fill_tunez.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.get_nodata_pos()
        mfill[nodatapos] = dem.get_nodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        
    def test_fill_06(self):
        dem = DEM("data/tunez2.tif")
        fill = dem.fill_sinks().read_array().astype("int16")
        
        mfill = sio.loadmat('data/mlab_files/fill_tunez2.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.get_nodata_pos()
        mfill[nodatapos] = dem.get_nodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        

if __name__ == "__main__":
    unittest.main()