#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for get_presills() function
@author: J. Vicente Perez
@email: geolovic@hotmail.com

Last modified: May 23, 2018
"""

import unittest
import sys
import numpy as np
import scipy.io as sio
import gdal
# Add to the path code folder and data folder
sys.path.append("../")
from sortcells import get_presills
infolder = "data"

class PreSillTest(unittest.TestCase):
    
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
    
    def test_get_presills(self):
        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]
        
        for idx, file in enumerate(files):
            # Locate data
            fill_path = infolder + "/fill_{0}.npy".format(file)
            flats_path = infolder + "/flats_{0}.npy".format(file)
            sills_path = infolder + "/sills_{0}.npy".format(file)
            presills_mlab_path = infolder + "/mlab_files/presills_{0}.mat".format(file)
            
            # Load numpy data
            fill_arr = np.load(fill_path)
            flats_arr = np.load(flats_path)
            sills_arr = np.load(sills_path)
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999
            
            # Get Presills
            presills = get_presills(fill_arr, flats_arr, sills_arr, False)
            
            # Get MatLab results
            mpresill_inds = self.load_matlab_array(presills_mlab_path, 'PreSillPixel', np.int32, nodata)
            mpresill_inds = mpresill_inds.ravel() - 1
            row, col = np.unravel_index(mpresill_inds, fill_arr.shape, "F")
            mpresills = np.zeros(fill_arr.shape, dtype=np.bool)
            mpresills[row, col] = True
            
            # Compare
            computed = np.array_equal(presills, mpresills)
            self.assertEqual(computed, True)            
        

if __name__ == "__main__":
    unittest.main()