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
        self.ids = np.load("data/MY_GRID_100rnd_id.npy")
        self.rows = np.load("data/MY_GRID_100rnd_row.npy")
        self.cols = np.load("data/MY_GRID_100rnd_col.npy")
        self.xi = np.load("data/MY_GRID_100rnd_X.npy")
        self.yi = np.load("data/MY_GRID_100rnd_Y.npy")
        self.zi = np.load("data/MY_GRID_100rnd_Z.npy")
        
        # Create a Grid object
        self.dem = Grid(MY_GRID)  

    def test_get_value01(self):
        # Taking row, col as integer
        ind = np.random.randint(0, 100)
        row, col = self.rows[ind], self.cols[ind]
        expected = self.zi[ind]
        computed = self.dem.get_value(row, col)
        self.assertEqual(computed, expected)
    
    def test_get_value02(self):
         # Taking row, col as numpy arrays
        expected = self.zi
        computed = self.dem.get_value(self.rows, self.cols)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)
        
    def test_get_value03(self):
        # Taking row, col as lists
        expected = self.zi.tolist()
        computed = self.dem.get_value(self.rows.tolist(), self.cols.tolist())
        self.assertEqual(computed, expected)

    def test_xy2cell01(self):
        xi = self.xi
        yi = self.yi
        rows = self.rows
        cols = self.cols
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_xy2cell02(self):
        ind = np.random.randint(0, 100)
        x = self.xi[ind]
        y = self.yi[ind]
        row = self.rows[ind]
        col = self.cols[ind]
        expected = (row, col)
        computed = self.dem.xy_2_cell(x, y)
        self.assertEqual(computed, expected)
        
    def test_xy2cell03(self):
        xi = self.xi.tolist()
        yi = self.yi.tolist()
        rows = self.rows.tolist()
        cols = self.cols.tolist()
        expected = (rows, cols)
        computed = self.dem.xy_2_cell(xi, yi)
        self.assertEqual(computed, expected)
        
    def test_ind2cell01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind2cell02(self):
        ind = self.ids.tolist()
        row = self.rows.tolist()
        col = self.cols.tolist()
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        self.assertEqual(computed, expected)
        
    def test_ind2cell03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = (row, col)
        computed = self.dem.ind_2_cell(ind)
        comparison = (computed[0] == expected[0]).all() and (computed[1] == expected[1]).all()
        self.assertEqual(comparison, True)
        
    def test_cell2ind01(self):
        idx = np.random.randint(0, 100)
        ind = self.ids[idx]
        row = self.rows[idx]
        col = self.cols[idx]
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell2ind02(self):
        ind = self.ids.tolist()
        row = self.rows.tolist()
        col = self.cols.tolist()
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        self.assertEqual(computed, expected)
        
    def test_cell2ind03(self):
        ind = self.ids
        row = self.rows
        col = self.cols
        expected = ind
        computed = self.dem.cell_2_ind(row, col)
        comparison = (computed == expected).all()
        self.assertEqual(comparison, True)

if __name__ == "__main__":
    unittest.main()