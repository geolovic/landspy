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
infolder = "data/in"
outfolder = "data/out"


class DEM_class(unittest.TestCase):
      
    def test_empty_dem(self):
        # test an empty DEM
        dem = DEM()
        computed = (dem._size, dem._geot, dem._proj, dem._nodata)
        expected = ((1, 1), (0.0, 1.0, 0.0, 1.0, 0.0, -1.0), "", -9999.0)
        self.assertEqual(computed, expected)
        
    def test_empty_dem02(self):
        dem = DEM()
        computed = np.array([[0]], dtype=np.float)
        expected = dem._array       
        self.assertEqual(np.array_equal(computed, expected), True)
        
    def test_copyDEM(self):
        files = ["small25", "jebja30", "tunez"]
        for file in  files:
            dem = DEM("{}/{}.tif".format(infolder, file))
            dem2 = dem.copy()
            computed = np.array_equal(dem._array, dem2._array)
            self.assertEqual(computed, True)
        
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
        dem = DEM(infolder + "/small25.tif")
        fill = dem.fill_sinks().read_array().astype("int16")
        
        mfill = sio.loadmat(infolder + '/mlab_files/fill_small25.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.get_nodata_pos()
        mfill[nodatapos] = dem.get_nodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        
    def test_fill_05(self):
        dem = DEM(infolder + "/tunez.tif")
        fill = dem.fill_sinks().read_array().astype("int16")
        
        mfill = sio.loadmat(infolder + '/mlab_files/fill_tunez.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.get_nodata_pos()
        mfill[nodatapos] = dem.get_nodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        

class DEMFlatTest(unittest.TestCase):
    
    def test_identify_flats_00(self):               
        # Create a DEM object and make fill
        dem = DEM(infolder + "/tunez.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=False)
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat(infolder + "/mlab_files/flats_tunez.mat")['flats']
        m_sills = sio.loadmat(infolder + "/mlab_files/sills_tunez.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      
    def test_identify_flats_02(self):               
        # Create a DEM object and make fill
        dem = DEM(infolder + "/small25.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=False)
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat(infolder + "/mlab_files/flats_small25.mat")['flats']
        m_sills = sio.loadmat(infolder + "/mlab_files/sills_small25.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      

if __name__ == "__main__":
    unittest.main()