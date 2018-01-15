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
from topopy.c_distance import distance, cost


class DistanceTest(unittest.TestCase):
    
    def test_distance_00(self):
        in_arr = np.array([[1, 0, 0, 0], 
                           [0, 0, 0, 0],
                           [0, 0, 0, 0], 
                           [0, 0, 0, 0]], dtype=np.int)
        
        rs_arr = np.array([[0.0000, 1.0000, 2.0000, 3.0000], 
                           [1.0000, 1.4142, 2.4142, 3.4142],
                           [2.0000, 2.4142, 2.8284, 3.8284], 
                           [3.0000, 3.4142, 3.8284, 4.2426]], dtype=np.float32) 
        
        dist = distance(in_arr)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001))

    def test_distance_01(self):
        in_arr = np.array([[0, 0, 0, 0, 0], 
                           [0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0], 
                           [0, 0, 0, 0, 0]], dtype=np.int)
        
        rs_arr = np.array([[2.4142, 1.4142, 1.0000, 1.4142, 2.4142], 
                           [2.0000, 1.0000, 0.0000, 1.0000, 2.0000],
                           [2.0000, 1.0000, 0.0000, 1.0000, 2.0000], 
                           [2.4142, 1.4142, 1.0000, 1.4142, 2.4142]], dtype=np.float32) 
        
        dist = distance(in_arr)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001)) 
        
    def test_distance_02(self):
        in_arr = np.array([[1, 0, 0, 0, 0], 
                           [0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0], 
                           [0, 0, 0, 0, 1]], dtype=np.int)
        
        rs_arr = np.array([[0.0000, 1.0000, 2.0000, 3.0000, 3.0000], 
                           [1.0000, 1.4142, 2.4142, 2.4142, 2.0000],
                           [2.0000, 2.4142, 2.4142, 1.4142, 1.0000], 
                           [3.0000, 3.0000, 2.0000, 1.0000, 0.0000]], dtype=np.float32) 
        
        dist = distance(in_arr)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001))        
    
    def test_distance_03(self):
        in_arr = np.zeros((500, 500), dtype=np.int)
        in_arr[[0, 1, 0, 250, 251, 499, 499, 498], [0, 0, 1, 250, 250, 499, 498, 499]] = 1
        
        rs_arr = np.load("data/res_disttest03.npy")
        dist = distance(in_arr)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001)) 


class CostTest(unittest.TestCase):
    
    def test_cost_00(self):
        in_arr = np.array([[1, 0, 0, 0], 
                           [0, 0, 0, 0],
                           [0, 0, 0, 0], 
                           [0, 0, 0, 0]], dtype=np.int)
    
        cost_surf = np.array([[1, 1, 2, 2], 
                              [1, 1, 2, 2],
                              [1, 1, 2, 2], 
                              [3, 3, 3, 3]], dtype=np.float)
        
        rs_arr = np.array([[0.0000, 1.0000, 2.5000, 4.5000], 
                           [1.0000, 1.4142, 2.9142, 4.9142],
                           [2.0000, 2.4142, 3.5355, 5.5355], 
                           [4.0000, 4.4142, 5.2426, 7.0711]], dtype=np.float) 
        
        dist = cost(in_arr, cost_surf)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001))

    def test_cost_01(self):
        in_arr = np.array([[0, 0, 0, 0, 0], 
                           [0, 0, 1, 0, 0],
                           [0, 0, 1, 0, 0], 
                           [0, 0, 0, 0, 0]], dtype=np.int)
    
        cost_surf = np.array([[1, 1, 2, 2, 2], 
                              [1, 1, 1, 2, 2],
                              [1, 1, 2, 2, 3], 
                              [3, 3, 3, 3, 4]], dtype=np.float)
        
        rs_arr = np.array([[2.4142, 1.4142, 1.5000, 2.1213, 4.1213], 
                           [2.0000, 1.0000, 0.0000, 1.5000, 3.5000],
                           [2.4142, 1.4142, 0.0000, 2.0000, 4.5000], 
                           [4.2426, 3.4142, 2.5000, 3.5355, 6.2426]], dtype=np.float32) 
        
        dist = cost(in_arr, cost_surf)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 0.001)) 
    
    def test_cost_02(self):
        in_arr = np.zeros((500, 500), dtype=np.int)
        in_arr[[0, 1, 498, 499], [0, 0, 499, 499]] = 1
        np.random.seed(1)
        cost_surf = np.random.randint(1, 5, (500, 500)) 
        rs_arr = np.load("data/res_costtest02.npy")
        dist = cost(in_arr, cost_surf)    
        diff = dist - rs_arr
        self.assertTrue(np.all(diff < 20)) 
        
if __name__ == "__main__":
    unittest.main()