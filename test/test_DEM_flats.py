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
from topopy import DEM
 

class DEMFlatTest(unittest.TestCase):
    
    def test_identify_flats_00(self):               
        # Create a DEM object
        dem = DEM("data/tunez.tif")
        
        # Load matlab-derivated flats and sills
        m_flats = np.load("data/tunez_mlab_flats.npy")
        m_sills = np.load("data/tunez_mlab_sills.npy")
        
        # Compute dem flats and sills
        flats, sills = dem.identify_flats()
        flats = flats.read_array().astype("int8")
        sills = sills.read_array().astype("int8")
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
        
    def test_identify_flats_01(self):               
        # Create a DEM object (DEM Tunez with NoData)
        dem = DEM("data/tunez2.tif")
        
        # Load matlab-derivated flats and sills
        m_flats = np.load("data/tunez2_mlab_flats.npy")
        m_sills = np.load("data/tunez2_mlab_sills.npy")
        
        # Compute dem flats and sills
        # Flats and sills from Grid class will have nodata values
        # (Matlab do not have)
        flats, sills = dem.identify_flats()
        flat_arr = flats.read_array()
        sill_arr = sills.read_array()
        flat_arr[np.where(np.isnan(flat_arr))] = 0
        sill_arr[np.where(np.isnan(sill_arr))] = 0
        
        flats = flat_arr.astype("int8")
        sills = sill_arr.astype("int8")
        
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))


if __name__ == "__main__":
    unittest.main()