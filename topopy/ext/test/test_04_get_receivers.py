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
# Add to the path code folder and data folder
sys.path.append("../")
from sortcells import sort_dem, get_receivers
infolder = "data"

class GetReceiversTest(unittest.TestCase):
    
    def test_get_receivers(self):
        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]
        cellsizes = [28.126263910642397, 25.0, 28.126263910642397]

        for idx, file in enumerate(files):
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999
            cellsize = cellsizes[idx]
            
            # Load numpy data
            pfill = np.load(infolder + "/fill_{0}.npy".format(file))
            # Change nodata values to large value (to mimic Matlab sorting)
            pfill[np.where(pfill==nodata)] = np.iinfo(pfill.dtype).max
            pweights = np.load(infolder + "/weights_{0}.npy".format(file))
            
            # Load matlab data
            mixc = sio.loadmat(infolder + "/mlab_files/ixc0_{0}.mat".format(file))['ixc']
            mixc = mixc.ravel() - 1
            
            # Sort DEM pixels
            ix = sort_dem(pfill, pweights, order="F")
            # Get receivers
            # Put again nodata in place
            nodataval = np.iinfo(pfill.dtype).max
            pfill = pfill.astype(np.float32)
            pfill[np.where(pfill==nodataval)] = np.nan
            ixc = get_receivers(ix, pfill, cellsize, order="F")
            
            # Compare
            res = np.array_equal(ixc, mixc)
            self.assertEqual(res, True) 


if __name__ == "__main__":
    unittest.main()