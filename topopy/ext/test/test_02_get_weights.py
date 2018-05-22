#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for get_weights() function
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
from sortcells import get_weights
infolder = "data"

class WeightsTest(unittest.TestCase):
    
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
    
    def test_get_weights(self):

        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]
        
        for idx, file in enumerate(files):
            
            # Locate data
            auxtopo_path = infolder + "/auxtopo_{0}.npy".format(file)
            flats_path = infolder + "/flats_{0}.npy".format(file)
            presills_path = infolder + "/presills_{0}.npy".format(file)
            weights_mlab_path = infolder + "/mlab_files/weights_{0}.mat".format(file)
            
            # Load numpy data
            auxtopo_arr = np.load(auxtopo_path)
            flats_arr = np.load(flats_path)
            presills_arr = np.load(presills_path)
            presills_pos = [(n[0], n[1]) for n in presills_arr]
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999

            # Get weights
            weights = get_weights(flats_arr, auxtopo_arr, presills_pos)

            # Load MatLab data
            mweights = self.load_matlab_array(weights_mlab_path, "D", np.float32, nodata)

            # Compare
            resta = np.abs(weights - mweights)
            res = np.all(resta<0.001)
            self.assertEqual(res, True)

if __name__ == "__main__":
    unittest.main()