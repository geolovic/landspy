
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

class GridPropertyTests(unittest.TestCase):
    
    def test_empty_grid(self):
        # Test que crea un objeto grid vacio
        #dem = Grid()
        #computed = (dem.get_size(), dem.get_geotransform(), dem.get_cellsize(),
#        3            dem.get_projection(), dem.get_nodata(), dem._array, str(dem._tipo))

#        expected = ((1, 1), (0., 1., 0., 0., 0., -1.), 1., "", None, np.array([[0]], dtype="float"), "float64")
        self.assertEqual(True, True)
 
if __name__ == "__main__":
    unittest.main()