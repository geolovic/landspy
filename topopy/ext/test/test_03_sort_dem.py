#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for sort_dem() function
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
from sortcells import sort_dem
infolder = "data"

class SortDEMTest(unittest.TestCase):
    
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

    def test_auxtopo(self):
        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]

        for idx, file in enumerate(files):
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999
            # Load numpy data
            pfill = np.load(infolder + "/fill_{0}.npy".format(file))
            # Change nodata values to large value (to mimic Matlab sorting)
            pfill[np.where(pfill==nodata)] = np.iinfo(pfill.dtype).max
            pweights = np.load(infolder + "/weights_{0}.npy".format(file))
            
            # Load matlab data
            mix = sio.loadmat(infolder + "/mlab_files/ix0_{0}.mat".format(file))['ix0']
            mix = mix.ravel() - 1
            
            # Sort DEM pixels
            ix = sort_dem(pfill, pweights, order="F")

            # Compare
            res = np.array_equal(ix, mix)
            self.assertEqual(res, True) 


if __name__ == "__main__":
    unittest.main()