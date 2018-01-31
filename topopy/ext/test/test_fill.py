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
from sortcells import fill_sinks

#(input_array, nodata_val=None):

class FillTest(unittest.TestCase):

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

    def test_fill_01(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 17, 20],
                        [19, 18, 31, 39],
                        [19, 29, 42, 51]])

        fill = fill_sinks(arr, -9999)
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 31, 39],
                        [19, 29, 42, 51]])
        computed = np.array_equal(fill, arr)
        expected = True       
        self.assertEqual(computed, expected)
    
    def test_fill_02(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 17, 39],
                        [19, 29, 42, 51]])

        fill = fill_sinks(arr, -9999)
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, 51]])
        
        computed = np.array_equal(fill, arr)
        expected = True
        self.assertEqual(computed, expected)
    
    def test_fill_03(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 19, 39],
                        [19, 29, 42, -99]])
        
        fill = fill_sinks(arr, -99)
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, -99]])
        
        computed = np.array_equal(fill, arr)
        expected = True
        self.assertEqual(computed, expected)

    def test_fill_04(self):
        # Testing real dems
        files = ['tunez', 'small25', 'tunez2']
        for file in files:
            dempath = "data/{0}.tif".format(file)
            mlab_path = "data/fill_{0}.mat".format(file)
        
            dem_arr, nodata = self.load_raster(dempath)
            if not nodata:
                nodata = -9999
            fill = fill_sinks(dem_arr, nodata)
            mfill = self.load_matlab_array(mlab_path, 'fill', dem_arr.dtype, nodata)
            
            computed = np.array_equal(fill, mfill)
            self.assertEqual(computed, True)


if __name__ == "__main__":
    unittest.main()