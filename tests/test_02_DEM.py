#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:58:02 2017
Testing suite for landspy Grid class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 19 september, 2022
"""

import unittest
import numpy as np
import scipy.io as sio
from skimage.morphology import reconstruction
from landspy import DEM

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"


class DEM_class(unittest.TestCase):
      
    def test_empty_dem(self):
        # test an empty DEM
        dem = DEM()
        computed = (dem._size, dem._geot, dem._proj, dem._nodata)
        expected = ((1, 1), (0.0, 1.0, 0.0, 1.0, 0.0, -1.0), "", -9999.0)
        self.assertEqual(computed, expected)
        
    def test_empty_dem02(self):
        dem = DEM()
        computed = np.array([[0]], dtype=np.float64)
        expected = dem._array       
        self.assertEqual(np.array_equal(computed, expected), True)
        
    def test_copyDEM(self):
        files = ["small25", "jebja30", "tunez"]
        for file in  files:
            dem = DEM("{}/{}.tif".format(infolder, file))
            dem2 = dem.copy()
            computed = np.array_equal(dem._array, dem2._array)
            self.assertEqual(computed, True)
        
    def test_fill_01(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 17, 20],
                        [19, 18, 31, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.setArray(arr)
        fill = dem.fill()
        computed = fill.readArray().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 31, 39],
                        [19, 29, 42, 51]])
        expected = arr.tolist()
        
        self.assertEqual(computed, expected)

  
    def test_fill_02(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 17, 39],
                        [19, 29, 42, 51]])
        dem = DEM()
        dem.setArray(arr)
        fill = dem.fill()
        computed = fill.readArray().tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, 51]])
        
        expected = arr.tolist()
    
        self.assertEqual(computed, expected)
      
    def test_fill_03(self):
        arr = np.array([[49, 36, 29, 29],
                        [32, 17, 12, 20],
                        [19, 17, 19, 39],
                        [19, 29, 42, -99]])
        dem = DEM()
        dem.setNodata(-99)
        dem.setArray(arr)
        fill = dem.fill()
        computed = fill.readArray(False).tolist()
        arr = np.array([[49, 36, 29, 29],
                        [32, 19, 19, 20],
                        [19, 19, 19, 39],
                        [19, 29, 42, -99]])
        
        expected = arr.tolist()
    
        self.assertEqual(computed, expected)
          
    def test_fill_04(self):
        dem = DEM(infolder + "/small25.tif")
        fill = dem.fill().readArray().astype("int16")
        
        mfill = sio.loadmat(infolder + '/mlab_files/fill_small25.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.getNodataPos()
        mfill[nodatapos] = dem.getNodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        
    def test_fill_05(self):
        dem = DEM(infolder + "/tunez.tif")
        fill = dem.fill().readArray().astype("int16")
        
        mfill = sio.loadmat(infolder + '/mlab_files/fill_tunez.mat')['fill']
        mfill = mfill.astype("int16")
        # Matlab files contain "nan" in the nodata positions
        nodatapos = dem.getNodataPos()
        mfill[nodatapos] = dem.getNodata()
        
        computed = np.array_equal(fill, mfill)
        self.assertEqual(computed, True)
        

class DEMFillMemoryTest(unittest.TestCase):

    def make_dem(self, dtype):
        dem = DEM()
        dem.setArray(np.array([[9, 9, 9], [9, 1, 9], [9, 9, 9]], dtype=dtype))
        dem._geot = (100, 30, 0, 200, 0, -30)
        dem._proj = 'test projection'
        return dem

    def test_fill_preserves_dtype_layout_and_original(self):
        for dtype in ('int16', 'int32', 'float32', 'float64'):
            with self.subTest(dtype=dtype):
                dem = self.make_dem(dtype)
                original = dem.readArray().copy()
                result = dem.fill()
                np.testing.assert_array_equal(result.readArray(), np.full((3, 3), 9))
                np.testing.assert_array_equal(dem.readArray(), original)
                self.assertFalse(np.shares_memory(result.readArray(), dem.readArray()))
                self.assertEqual(result.readArray().dtype, original.dtype)
                self.assertEqual(result.getSize(), dem.getSize())
                self.assertEqual(result.getGeot(), dem.getGeot())
                self.assertEqual(result.getCRS(), dem.getCRS())
                self.assertEqual(result.getNodata(), dem.getNodata())

    def test_fill_inplace(self):
        dem = self.make_dem('float32')
        original_array = dem.readArray()
        result = dem.fill(inplace=True)
        self.assertIs(result, dem)
        np.testing.assert_array_equal(dem.readArray(), np.full((3, 3), 9))
        self.assertEqual(original_array[1, 1], 1)
        self.assertEqual(dem.readArray().dtype, np.dtype('float32'))
        self.assertEqual(dem.getGeot(), (100, 30, 0, 200, 0, -30))

    def test_fill_array_return_modes(self):
        for inplace in (False, True):
            with self.subTest(inplace=inplace):
                dem = self.make_dem('float64')
                result = dem.fill(as_array=True, inplace=inplace)
                np.testing.assert_array_equal(result, np.full((3, 3), 9))
                if inplace:
                    self.assertIs(result, dem.readArray())
                else:
                    self.assertEqual(dem.readArray()[1, 1], 1)
                    self.assertFalse(np.shares_memory(result, dem.readArray()))

    def test_fill_nodata_outlet(self):
        for inplace in (False, True):
            with self.subTest(inplace=inplace):
                dem = DEM()
                dem.setNodata(-99)
                dem.setArray(np.array([[9, 9, 9], [9, 1, -99], [9, 9, 9]], dtype='int16'))
                expected = dem.readArray().copy()
                np.testing.assert_array_equal(dem.fill(inplace=inplace).readArray(), expected)

    def test_fill_without_nodata(self):
        dem = self.make_dem('float32')
        dem.setNodata(None)
        result = dem.fill(inplace=True)
        self.assertIsNone(result.getNodata())
        np.testing.assert_array_equal(result.readArray(), np.full((3, 3), 9))

    def test_fill_all_nodata(self):
        dem = DEM()
        dem.setArray(np.full((5, 5), dem.getNodata(), dtype='float32'))
        expected = dem.readArray().copy()
        np.testing.assert_array_equal(dem.fill(inplace=True).readArray(), expected)


class DEMPriorityFloodTest(unittest.TestCase):

    def test_matches_reconstruction_on_random_terrain(self):
        rng = np.random.default_rng(71)
        for dtype in ('int16', 'int64', 'uint16', 'float16', 'float32', 'float64'):
            for shape in ((1, 1), (1, 19), (19, 1), (2, 13), (17, 23)):
                with self.subTest(dtype=dtype, shape=shape):
                    dem = DEM()
                    arr = rng.integers(0, 100, shape).astype(dtype)
                    if np.dtype(dtype).kind == 'f':
                        arr /= 8
                    dem.setArray(arr)
                    dem.setNodata(None)
                    seed = dem.readArray().copy()
                    seed[1:-1, 1:-1] = seed.max()
                    expected = reconstruction(seed, dem.readArray(), 'erosion').astype(dem.readArray().dtype)
                    actual = dem.fill(as_array=True)
                    np.testing.assert_array_equal(actual, expected)
                    np.testing.assert_array_equal(dem.readArray(), arr)

    def test_diagonal_outlet(self):
        dem = DEM()
        dem.setArray(np.array([[1, 9, 9], [9, 1, 9], [9, 9, 9]], dtype='float32'))
        np.testing.assert_array_equal(dem.fill(as_array=True), dem.readArray())

    def test_negative_sentinel_is_processed_as_elevation(self):
        dem = DEM()
        arr = np.full((5, 5), 9, dtype='float32')
        arr[1:4, 1:4] = 1
        arr[2, 2] = dem.getNodata()
        dem.setArray(arr)
        expected = np.full((5, 5), 9, dtype='float32')
        np.testing.assert_array_equal(dem.fill(as_array=True), expected)

    def test_large_integer_elevations_remain_exact(self):
        dem = DEM()
        arr = np.full((3, 3), 2**54 + 3, dtype='int64')
        arr[1, 1] = 2**54 + 1
        dem.setArray(arr)
        np.testing.assert_array_equal(dem.fill(as_array=True), np.full_like(arr, 2**54 + 3))

    def test_nan_leaves_input_unchanged(self):
        dem = DEM()
        arr = np.array([[1, 1, 1], [1, np.nan, 1], [1, 1, 1]])
        dem.setArray(arr)
        with self.assertRaises(ValueError):
            dem.fill(inplace=True)
        np.testing.assert_array_equal(dem.readArray(), arr)


class DEMFlatTest(unittest.TestCase):
    
    def test_identify_flats_00(self):               
        # Create a DEM object and make fill
        dem = DEM(infolder + "/tunez.tif")
        fill = dem.fill()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identifyFlats(nodata=False)
        flats = flats.readArray()
        sills = sills.readArray()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat(infolder + "/mlab_files/flats_tunez.mat")['flats']
        m_sills = sio.loadmat(infolder + "/mlab_files/sills_tunez.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      
    def test_identify_flats_02(self):               
        # Create a DEM object and make fill
        dem = DEM(infolder + "/small25.tif")
        fill = dem.fill()
        
        # Identify flats and sills and load arrays
        flats, sills = fill.identifyFlats(nodata=False)
        flats = flats.readArray()
        sills = sills.readArray()
        
        # Load matlab flats and sills
        m_flats = sio.loadmat(infolder + "/mlab_files/flats_small25.mat")['flats']
        m_sills = sio.loadmat(infolder + "/mlab_files/sills_small25.mat")['sills']
        
        # Compare
        computed = (np.array_equal(m_flats, flats),
                    np.array_equal(m_sills, sills))
        self.assertEqual(computed, (True, True))
      

if __name__ == "__main__":
    unittest.main()
