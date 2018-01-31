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
import gdal
# Add to the path code folder and data folder
sys.path.append("../")
from sortcells import identify_flats
 

class FlatsSillsTest(unittest.TestCase):
    
    def load_matlab_array(self, path, key, nptype, nodata_val):
        marray = sio.loadmat(path)[key]        
        if nodata_val:
            nodatamask = np.isnan(marray)
            marray[nodatamask] = nodata_val
        marray = marray.astype(nptype)
        return marray

    def load_raster(self, path):
        raster = gdal.Open(path)
        banda = raster.GetRasterBand(1)
        arr = banda.ReadAsArray()
        nodata = banda.GetNoDataValue()
        return arr, nodata
    
    def test_identify_flats_00(self):

        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]
        
        for idx, file in enumerate(files):
            
            # Locate data
            fill_path = "data/fill_{0}.npy".format(file)
            flats_mlab_path = "data/flats_{0}.mat".format(file)
            sills_mlab_path = "data/sills_{0}.mat".format(file)
            
            # Load numpy data
            fill_arr = np.load(fill_path)
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999
            
            # Get Flats and Sills
            flats, sills = identify_flats(fill_arr, nodata)
            
            # Get MatLab results
            mflats = self.load_matlab_array(flats_mlab_path, 'flats', np.bool, nodata)
            msills = self.load_matlab_array(sills_mlab_path, 'sills', np.bool, nodata)
            
            # Compare
            computed = [np.array_equal(flats, mflats), np.array_equal(sills, msills)]
            self.assertEqual(computed, [True, True])             
    

if __name__ == "__main__":
    unittest.main()