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
 

class DEMFlatTest(unittest.TestCase):
    
    def test_identify_flats_00(self):               
        # Create a DEM object and make fill
        dem = DEM("data/tunez.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=False)
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat("data/mlab_files/flats_tunez.mat")['flats']
        m_sills = sio.loadmat("data/mlab_files/sills_tunez.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
     
    def test_identify_flats_01(self):               
        # Create a DEM object and make fill
        dem = DEM("data/tunez2.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=False)
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat("data/mlab_files/flats_tunez2.mat")['flats']
        m_sills = sio.loadmat("data/mlab_files/sills_tunez2.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      
    def test_identify_flats_02(self):               
        # Create a DEM object and make fill
        dem = DEM("data/small25.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=False)
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat("data/mlab_files/flats_small25.mat")['flats']
        m_sills = sio.loadmat("data/mlab_files/sills_small25.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
     
    def test_identify_flats_03(self):               
        # Create a DEM object and make fill
        dem = DEM("data/tunez2.tif")
        fill = dem.fill_sinks()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identify_flats(nodata=True)
        flats_nodata = flats.get_nodata_pos()
        sills_nodata = sills.get_nodata_pos()
        flats = flats.read_array()
        sills = sills.read_array()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat("data/mlab_files/flats_tunez2.mat")['flats']
        m_sills = sio.loadmat("data/mlab_files/sills_tunez2.mat")['sills']
        
        m_flats = m_flats.astype(np.int8)
        m_sills = m_sills.astype(np.int8)
        
        # Put nodata inplace
        m_flats[flats_nodata] = -1
        m_sills[sills_nodata] = -1
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      

if __name__ == "__main__":
    unittest.main()