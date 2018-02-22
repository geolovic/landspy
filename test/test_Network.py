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
import matplotlib.pyplot as plt
# Add to the path code folder and data folder
sys.path.append("../")
from topopy import Grid, Flow, Network

class CreateNetwork(unittest.TestCase):
    
    def test_create_load(self):
        dem_files = ['tunez', 'tunez2', 'small25']        
        for filename in dem_files: 
           fd = Flow("data/fd_{0}.tif".format(filename))
           st = Network(fd, 1000)
           computed = [st.get_dims(), st.get_size(), st.get_ncells(), 
                       st.get_cellsize(), st.get_geotransform(), st.get_projection()]
           expected = [fd.get_dims(), fd.get_size(), fd.get_ncells(), 
                       fd.get_cellsize(), fd.get_geotransform(), fd.get_projection()] 
           self.assertTrue(computed, expected)
       
    def test_streams(self):
        dem_files = ['tunez', 'tunez2', 'small25']        
        for filename in dem_files:
           fd = Flow("data/fd_{0}.tif".format(filename))
           st = Network(fd, 1000)
           streams = st.get_streams()
           st01 = streams.read_array()
           st02 = Grid("data/str_{0}.tif".format(filename)).read_array()
           self.assertTrue(np.array_equal(st01, st02), True)
           
    def test_streampoi(self):
        dem_files = ['tunez', 'tunez2', 'small25']        
        for filename in dem_files:
            fd = Flow("data/fd_{0}.tif".format(filename))
            st = Network(fd, 1000)
            kinds = ['heads', 'confluences', 'outlets']
            for kind in kinds:
                poi = st.get_stream_poi(kind)
                
                
   
           
if __name__ == "__main__":
    unittest.main()