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
from sortcells import get_aux_topography
 

class AuxTopoTest(unittest.TestCase):
    
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
    
    def test_aux_topography(self):

        # Data for testing
        files = ['tunez', 'small25', 'tunez2']
        nodatas = [None, -9999.0, -9999.0]
        
        for idx, file in enumerate(files):
            
            # Locate data
            topodiff_path = "data/topodiff_{0}.npy".format(file)
            flats_path = "data/flats_{0}.npy".format(file)
            auxtopo_mlab_path = "data/auxtopo_{0}.mat".format(file)
            
            # Load numpy data
            topodiff_arr = np.load(topodiff_path)
            flats_arr = np.load(flats_path)
            nodata = nodatas[idx]
            if not nodata:
                nodata = -9999

            # Get aux topography (change first input array types)
            topodiff_arr = topodiff_arr.astype(np.float32)
            flats_arr = flats_arr.astype(np.int8)
            auxtopo = get_aux_topography(topodiff_arr, flats_arr)
            
            # Get matlab-derived data
            mlab_file = 'data/auxtopo_{0}.mat'.format(file)
            mauxtopo = sio.loadmat(mlab_file)['D']
            # Reemplazamos nan por ceros
            nanmask = np.isnan(mauxtopo)
            mauxtopo[nanmask] = 0.0

            # Comparamos
            res = np.array_equal(auxtopo, mauxtopo)
            self.assertEqual(res, True)
            
            # # Get Auxiliar topography [Changing types first !!!]
            # topodiff_arr = topodiff_arr.astype(np.float32)
            # flats_arr = flats_arr.astype(np.int8)
            # auxtopo = get_aux_topography(topodiff_arr, flats_arr)
            
            # # Get MatLab results
            # mlab_file = 'data/auxtopo_{0}.mat'.format(file)
            # mauxtopo = sio.loadmat(mlab_file)['D']
            # # Reemplazamos nan por ceros
            # nanmask = np.isnan(mauxtopo)
            # mauxtopo[nanmask] = 0.0
            
            # # Compare
            # computed = np.array_equal(auxtopo, mauxtopo)
            #self.assertEqual(computed, True) 
    # def test_auxtopo(self):
    #     # Create Flow object
    #     in_dems = ['tunez', 'tunez2', 'small25']
        
    #     for path_dem in in_dems:
    #         # Get DEM
    #         dempath = "data/{0}.tif".format(path_dem)        
    #         dem_arr, nodata = self.load_raster(dempath)
    #         if not nodata:
    #             nodata = -9999
    #         # Fill DEM and get flats and sills
    #         fill_arr = fill_sinks(dem_arr, nodata)
    #         flats, sills = identify_flats(fill_arr, nodata)
    #         # Get auxiliar topography (fill - dem)
    #         topodiff = fill_arr - dem_arr
    #         # Change types and get auxiliar topography
    #         topodiff = topodiff.astype(np.float32)
    #         flats = flats.astype(np.int8)
    #         auxtopo = get_aux_topography(topodiff, flats)
            
    #         # Get matlab-derived data
    #         mlab_file = 'data/auxtopo_{0}.mat'.format(path_dem)
    #         mauxtopo = sio.loadmat(mlab_file)['D']
    #         # Reemplazamos nan por ceros
    #         nanmask = np.isnan(mauxtopo)
    #         mauxtopo[nanmask] = 0.0
          
    #         self.assertEqual(np.array_equal(auxtopo, mauxtopo), True) 
    

if __name__ == "__main__":
    unittest.main()