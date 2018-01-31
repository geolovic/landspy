#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for topopy Flow class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""

import unittest
import sys
import numpy as np
import scipy.io as sio
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import DEM, Flow


class CreateFlow(unittest.TestCase):
    
    def test_presill(self):
        # Create Flow object
        in_dems = ['tunez', 'tunez2', 'small25']
        
        for path_dem in in_dems:
            dem = DEM('data/' + path_dem + ".tif")
            dem = dem.fill_sinks()
            flats, sills = dem.identify_flats(False)
            fd = Flow()
            fd._dims = dem.get_dims()
            rowcol = fd._get_presills(flats, sills, dem)
            rowcol = np.array(rowcol)
            row = rowcol[:,0]
            col = rowcol[:,1]
            presills = np.zeros(fd._dims, np.int)
            presills[row, col]  = 1

            # Get matlab-derived data
            mlab_file = 'data/mlab_files/presills_{0}.mat'.format(path_dem)
            mpresills = sio.loadmat(mlab_file)['PreSillPixel']
            row, col = np.unravel_index(mpresills-1, fd._dims, "F")
            mpresills = np.zeros(fd._dims, np.int)
            mpresills[row, col]  = 1
            self.assertEqual(np.array_equal(presills, mpresills), True)
    
    def test_auxtopo(self):
        # Create Flow object
        in_dems = ['tunez', 'tunez2', 'small25']
        
        for path_dem in in_dems:
            dem = DEM('data/' + path_dem + ".tif")
            fill = dem.fill_sinks()
            topodiff = fill.read_array() - dem.read_array()
            topodiff = topodiff.astype(np.float32)
            flats, sills = fill.identify_flats(False)
            fd = Flow()
            fd._dims = dem.get_dims()
            auxtopo = fd._get_topodiff(topodiff, flats)
            
            # Get matlab-derived data
            mlab_file = 'data/mlab_files/auxtopo_{0}.mat'.format(path_dem)
            mauxtopo = sio.loadmat(mlab_file)['D'].astype(np.int16)
          
            self.assertEqual(np.array_equal(auxtopo, mauxtopo), True)
    
    def test_weights(self):
        # Create Flow object
        in_dems = ['tunez', 'tunez2', 'small25']
        
        for path_dem in in_dems:
            dem = DEM('data/' + path_dem + ".tif")
            fill = dem.fill_sinks()
            topodiff = fill.read_array() - dem.read_array()
            topodiff = topodiff.astype(np.float32)
            dem = fill
            flats, sills = fill.identify_flats(False)
            fd = Flow()
            fd._dims = dem.get_dims()
            auxtopo = fd._get_topodiff(topodiff, flats)
            presill_pos = fd._get_presills(flats, sills, dem)
            weights = fd._get_weights(flats, auxtopo, presill_pos)
            
            # Get matlab-derived data
            mlab_file = 'data/mlab_files/weights_{0}.mat'.format(path_dem)
            mweights = sio.loadmat(mlab_file)['D']
            
            resta = np.abs(weights - mweights)
          
            self.assertEqual(np.all(resta < 0.01), True)
            
if __name__ == "__main__":
    unittest.main()