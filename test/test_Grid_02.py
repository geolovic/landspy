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
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid

MY_GRID = "data/small25.tif"

class GridValueTests(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load("data/small25_100rnd_id.npy")
        self.rows = np.load("data/small25_100rnd_row.npy")
        self.cols = np.load("data/small25_100rnd_col.npy")
        self.xi = np.load("data/small25_100rnd_X.npy")
        self.yi = np.load("data/small25_100rnd_Y.npy")
        self.zi = np.load("data/small25_100rnd_Z.npy")
        
        # Create a DEM
        self.dem = Grid(MY_GRID)

    def test_get_value_01(self):
        # Taking row, col in a nan position (88)
        ind = 88
        row, col = self.rows[ind], self.cols[ind]
        computed = self.dem.get_value(row, col)
        self.assertEqual(computed, -9999)

    def test_get_value_02(self):
        # Taking row, col in other position (with value)
        ind = 25
        row, col = self.rows[ind], self.cols[ind]
        expected = self.zi[ind]
        computed = self.dem.get_value(row, col)
        self.assertEqual(computed, expected)
        
    def test_get_value_03(self):
        # Taking row, col outside array
        row, col = 199, 133
        self.assertRaises(IndexError, self.dem.get_value, row, col)
    
    def test_get_value_04(self):
        # Taking row, col as numpy arrays
        expected = self.zi
        computed = self.dem.get_value(self.rows, self.cols)
        self.assertEqual(np.nansum(computed),np.nansum(expected))
        
    def test_get_value_05(self):
        # Taking row, col as lists
        expected = (np.nansum(self.zi.tolist()), list)
        res =  self.dem.get_value(self.rows.tolist(), self.cols.tolist())
        computed = (np.nansum(res), type(res))
        self.assertEqual(computed, expected)

    def test_xy_2_cell_01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_xy_2_cell_02(self):
        ind = np.random.randint(0, 100)
        x = self.xi[ind]
        y = self.yi[ind]
        row = self.rows[ind]
        col = self.cols[ind]
        expected = (row, col)
        computed = self.dem.xy_2_cell(x, y)
        self.assertEqual(computed, expected)
        
    def test_xy_2_cell_03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows.tolist()
        cols = self.cols.tolist()
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        self.assertEqual(computed, expected)
        
    def test_ind_2_cell_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind_2_cell_02(self):
        ind = self.ids.tolist()
        row = self.rows.tolist()
        col = self.cols.tolist()
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind_2_cell_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_cell_2_ind_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell_2_ind_02(self):
        ind = self.ids.tolist()
        row = self.rows.tolist()
        col = self.cols.tolist()
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell_2_ind_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)

if __name__ == "__main__":
    unittest.main()