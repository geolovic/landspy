#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thrusday 08 Feb 2018
Testing suite for topopy Flow class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
"""

import unittest
import sys
import numpy as np
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Flow

class FlowValueTest(unittest.TestCase):
    
    def setUp(self):        
        # Load test data
        self.ids = np.load("data/small25_100rnd_id.npy")
        self.rows = np.load("data/small25_100rnd_row.npy")
        self.cols = np.load("data/small25_100rnd_col.npy")
        self.xi = np.load("data/small25_100rnd_X.npy")
        self.yi = np.load("data/small25_100rnd_Y.npy")
        self.zi = np.load("data/small25_100rnd_Z.npy")
        
        # Load Flow object
        self.fd = Flow("data/fd_small25.tif")
    
    def test_xy_2_cell_01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (rows, cols)
        computed = self.fd.xy_2_cell(xi, yi)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_xy_2_cell_02(self):
        ind = np.random.randint(0, 100)
        x = self.xi[ind]
        y = self.yi[ind]
        row = self.rows[ind]
        col = self.cols[ind]
        expected = (row, col)
        computed = self.fd.xy_2_cell(x, y)
        self.assertEqual(computed, expected)
        
    def test_xy_2_cell_03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows
        cols = self.cols
        crows, ccols = self.fd.xy_2_cell(xi, yi)
        computed = (np.array_equal(rows, crows), np.array_equal(cols, ccols))
        self.assertEqual(computed, (True, True))
        
    def test_ind_2_cell_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = (row, col)
        computed = self.fd.ind_2_cell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind_2_cell_02(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        crow, ccol = self.fd.ind_2_cell(ind)
        computed = (np.array_equal(row, crow), np.array_equal(col, ccol))
        self.assertEqual(computed, (True, True))
        
    def test_ind_2_cell_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = (row, col)
        computed = self.fd.ind_2_cell(ind)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_cell_2_ind_01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = ind
        computed = self.fd.cell_2_ind(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell_2_ind_02(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.fd.cell_2_ind(row, col)
        res = np.array_equal(expected, computed)
        self.assertEqual(res, True)
        
    def test_cell_2_ind_03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.fd.cell_2_ind(row, col)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)
        
    def test_cell_2_xy_01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (xi, yi)
        computed = self.fd.cell_2_xy(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))
        
    def test_cell_2_xy_02(self):
        ind = np.random.randint(0, 100)
        xi = self.xi[ind]
        yi = self.yi[ind]
        rows = self.rows[ind]
        cols = self.cols[ind]
        expected = (xi, yi)
        computed = self.fd.cell_2_xy(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))
        
    def test_cell_2_xy_03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows
        cols = self.cols
        expected = (xi, yi)
        computed = self.fd.cell_2_xy(rows, cols)
        res = (np.array_equal(expected[0], computed[0]), np.array_equal(expected[1], computed[1]))
        self.assertEqual(res, (True, True))

        
if __name__ == "__main__":
    unittest.main()