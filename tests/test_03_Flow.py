#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thrusday 08 Feb 2018
Testing suite for landspy Flow class
@author: J. Vicente Perez
@email: geolovic@gmail.com
@last_modified: 19 september, 2022
"""

import unittest
from unittest.mock import patch
import numpy as np
from landspy import Flow, DEM
from landspy.flow import get_weights, sort_dem

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

class FlowFillReuseTest(unittest.TestCase):
    def test_unit_weights_without_cost_surface(self):
        flats = np.array([[True, True], [False, True]], dtype=bool)
        weights = get_weights(flats, None, [(0, 0)])
        self.assertEqual(weights.dtype, np.dtype('float32'))
        np.testing.assert_array_equal(weights[0], [1, 2])
        self.assertEqual(weights[1, 0], -99999)
        np.testing.assert_allclose(weights[1, 1], np.float32(1 + np.sqrt(2)), rtol=1e-7)
        fallback = get_weights(flats, None, [])
        np.testing.assert_array_equal(fallback, [[1, 1], [-99999, 1]])


    def test_weight_barriers_and_unreachable_flats(self):
        flats = np.array([[True, False, True],
                          [True, False, True]], dtype=bool)
        costs = np.ones(flats.shape, dtype=np.float64, order='F')
        weights = get_weights(flats, costs, [(0, 0)])
        self.assertEqual(weights.dtype, np.dtype('float64'))
        np.testing.assert_array_equal(weights[:, 0], [1, 2])
        np.testing.assert_array_equal(weights[:, 1], [-99999, -99999])
        self.assertTrue(np.isinf(weights[:, 2]).all())

    def test_weight_diagonal_distance(self):
        flats = np.ones((2, 2), dtype=bool)
        costs = np.ones(flats.shape, dtype=np.float64, order='F')
        weights = get_weights(flats, costs, [(0, 0)])
        np.testing.assert_allclose(weights, [[1, 2], [2, 1 + np.sqrt(2)]],
                                   rtol=1e-15)

    def test_weights_without_presills(self):
        flats = np.array([[True, False]], dtype=bool)
        weights = get_weights(flats, np.ones((1, 2), order='F'), [])
        np.testing.assert_array_equal(weights, [[1, -99999]])
        self.assertEqual(weights.dtype, np.dtype('float64'))

    def test_stable_sort_order(self):
        rng = np.random.default_rng(19)
        elevations = rng.integers(0, 5, (17, 23)).astype('float32')
        weights = rng.integers(0, 4, elevations.shape).astype('float64')
        for order in ('C', 'F'):
            with self.subTest(order=order):
                expected = np.lexsort((np.arange(elevations.size),
                                       -weights.ravel(order=order),
                                       -elevations.ravel(order=order)))
                np.testing.assert_array_equal(sort_dem(elevations, weights, order), expected)

    def test_unconnected_cells(self):
        flow = Flow()
        flow._size = (3, 3)
        flow._ix = np.array([0, 2, 7], dtype='uint32')
        flow._ixc = np.array([1, 3, 8], dtype='uint32')
        np.testing.assert_array_equal(flow._get_nodata_pos(), [4, 5, 6])

    def test_fill_count_and_elevations(self):
        original = np.array([[9, 8, 7, 6, 5],
                             [8, 7, 6, 5, 4],
                             [7, 6, 0, 4, 3],
                             [6, 5, 4, 3, 2],
                             [5, 4, 3, 2, 1]], dtype='float32')
        for filled in (False, True):
            for raw_z in (False, True):
                for auxtopo in (False, True):
                    with self.subTest(filled=filled, raw_z=raw_z, auxtopo=auxtopo):
                        dem = DEM()
                        dem.setArray(original)
                        filled_array = dem.fill(as_array=True)
                        self.assertGreater(filled_array[2, 2], original[2, 2])
                        if filled:
                            dem.setArray(filled_array)
                        before = dem.readArray().copy()
                        with patch.object(DEM, 'fill', autospec=True,
                                          side_effect=DEM.fill) as fill_call:
                            flow = Flow(dem, filled=filled, raw_z=raw_z, auxtopo=auxtopo)
                        self.assertEqual(fill_call.call_count, 0 if filled else 1)
                        expected = before if raw_z else filled_array
                        np.testing.assert_array_equal(flow._zx, expected.ravel()[flow._ix])
                        np.testing.assert_array_equal(dem.readArray(), before)


class FlowValueTest(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load(infolder + "/np_files/small25_100rnd_id.npy")
        self.rows = np.load(infolder + "/np_files/small25_100rnd_row.npy")
        self.cols = np.load(infolder + "/np_files/small25_100rnd_col.npy")
        self.xi = np.load(infolder + "/np_files/small25_100rnd_X.npy")
        self.yi = np.load(infolder + "/np_files/small25_100rnd_Y.npy")
        self.zi = np.load(infolder + "/np_files/small25_100rnd_Z.npy")
        # Load Flow object
        self.fd = Flow(infolder + "/small25_fd.tif")
    
    def test_xy_2_cell_01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (rows, cols)
        computed = self.fd.xyToCell(xi, yi)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_xy_2_cell_02(self):
        ind = np.random.randint(0, 100)
        x = self.xi[ind]
        y = self.yi[ind]
        row = self.rows[ind]
        col = self.cols[ind]
        expected = (row, col)
        computed = self.fd.xyToCell(x, y)
        self.assertEqual(computed, expected)
        
    def test_xy_2_cell_03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows
        cols = self.cols
        crows, ccols = self.fd.xyToCell(xi, yi)
        computed = (np.array_equal(rows, crows), np.array_equal(cols, ccols))
        self.assertEqual(computed, (True, True))
        
    def test_ind_2_cell_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = (row, col)
        computed = self.fd.indToCell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind_2_cell_02(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        crow, ccol = self.fd.indToCell(ind)
        computed = (np.array_equal(row, crow), np.array_equal(col, ccol))
        self.assertEqual(computed, (True, True))
        
    def test_ind_2_cell_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = (row, col)
        computed = self.fd.indToCell(ind)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_cell_2_ind_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = ind
        computed = self.fd.cellToInd(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell_2_ind_02(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.fd.cellToInd(row, col)
        res = np.array_equal(expected, computed)
        self.assertEqual(res, True)
        
    def test_cell_2_ind_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.fd.cellToInd(row, col)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)
        
    def test_cell_2_xy_01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (xi, yi)
        computed = self.fd.cellToXY(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))
        
    def test_cell_2_xy_02(self):
        ind = np.random.randint(0, 100)
        xi = self.xi[ind]
        yi = self.yi[ind]
        rows = self.rows[ind]
        cols = self.cols[ind]
        expected = (xi, yi)
        computed = self.fd.cellToXY(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))
        
    def test_cell_2_xy_03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows
        cols = self.cols
        expected = (xi, yi)
        computed = self.fd.cellToXY(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))


class FlowPropertiesTest(unittest.TestCase):
    
    def test_flow_properties_01(self):
        files = ["small", "small25", "tunez", "jebja30"]
        for file in files:
            dem = DEM(infolder + "/{0}.tif".format(file))
            fd = Flow(dem)
            computed = (fd.getSize(), fd.getDims(), fd.getNCells(), fd.getCRS(), fd.getCellSize(), fd.getGeot())
            expected = (dem.getSize(), dem.getDims(), dem.getNCells(), dem.getCRS(), dem.getCellSize(), dem.getGeot())
            self.assertEqual(computed, expected)
    
    def test_flow_properties_02(self):
        # Testing an empty Flow object
        dem = DEM()
        fd = Flow()
        computed = (fd.getSize(), fd.getDims(), fd.getNCells(), fd.getCRS(), fd.getCellSize(), fd.getGeot())
        expected = (dem.getSize(), dem.getDims(), dem.getNCells(), dem.getCRS(), dem.getCellSize(), dem.getGeot())
        self.assertEqual(computed, expected)    


class FlowCreateTest(unittest.TestCase):
    
    def test_create_Flow_01(self):
        # Create Flow object with and without auxiliar topography
        files = ["small", "small25", "tunez", "jebja30"]
        for file in files:
            dem_path = infolder + "/{0}.tif".format(file)
            dem = DEM(dem_path)
            # Flow obj without auxiliar topography
            fd = Flow(dem, auxtopo=False, filled=False, verbose=False)
            # Save the flow object for later tests
            fd.save(outfolder + "/{0}_fd.tif".format(file))
            computed = (fd.getSize(), fd.getDims(), fd.getNCells(), fd.getCRS(), fd.getCellSize(), fd.getGeot())
            expected = (dem.getSize(), dem.getDims(), dem.getNCells(), dem.getCRS(), dem.getCellSize(), dem.getGeot())
            self.assertEqual(computed, expected)    
    
    def test_create_Flow_02(self):
        # Create Flow object with and without auxiliar topography
        files = ["small", "small25", "tunez", "jebja30"]
        for file in files:
            dem_path = infolder + "/{0}.tif".format(file)
            dem = DEM(dem_path)
            # Flow obj without auxiliar topography
            fd = Flow(dem, auxtopo=True, filled=False, verbose=False)
            computed = (fd.getSize(), fd.getDims(), fd.getNCells(), fd.getCRS(), fd.getCellSize(), fd.getGeot())
            expected = (dem.getSize(), dem.getDims(), dem.getNCells(), dem.getCRS(), dem.getCellSize(), dem.getGeot())
            self.assertEqual(computed, expected)    
            

    def test_create_Flow_03(self):
        # Create Flow object with previously filled DEM (auxtopo is ignored)
        files = ["small", "small25", "tunez", "jebja30"]
        for file in files:
            dem_path = infolder + "/{0}.tif".format(file)
            dem = DEM(dem_path)
            dem = dem.fill()
            # Flow obj with filled dem and without auxiliar topography
            fd = Flow(dem, auxtopo=False, filled=True, verbose=False)
            computed = (fd.getSize(), fd.getDims(), fd.getNCells(), fd.getCRS(), fd.getCellSize(), fd.getGeot())
            expected = (dem.getSize(), dem.getDims(), dem.getNCells(), dem.getCRS(), dem.getCellSize(), dem.getGeot())
            self.assertEqual(computed, expected) 

if __name__ == "__main__":
    unittest.main()
